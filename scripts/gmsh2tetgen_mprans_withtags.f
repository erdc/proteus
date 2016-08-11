	program gmsh2tetgen_withtags
!=============================================================
! program to read in the original input(*.mesh) file and 
! translates to tetgen format
! also output vtk files for visualization of generated mesh 
! this version includes element material types
!=============================================================
	implicit none
      
	integer :: i, j, ier, meshf, itmp
	integer :: NSHL, NSHLb, NSD, NNode, NElem, NEdge, EType
	real(kind=8), allocatable :: xg(:,:)
	real(kind=8) theta, xtmp(3), st,ct
	integer :: rNFace, NFACE, hNFace, FACE_ID
	integer, allocatable ::  Face_IEN(:,:), FaceID(:)
	integer, allocatable :: hFace_IEN(:,:)
	integer, allocatable :: IEN(:,:)
! hold element material types assuming these come from the physical group in gmsh
	integer, allocatable :: IEMAT(:)

	character(len=128) :: fname
	character(len=30) :: length_scale_string
	character(len=30) :: perm_x_string
	character(len=30) :: perm_y_string
	character(len=30) :: perm_z_string
  	character(len=128) :: ctmp
	real(kind=8) length_scale
	integer :: perm(3)

	EType = 2
	NSHL  = 4
	NSHLb = 3

!==================================
! Read mesh in .mesh format
!==================================
	  meshf = 11      
	  length_scale=1.0
	  perm(1) = 1
	  perm(2) = 2
	  perm(3) = 3      
	  i = iargc()
	  if (i.eq.1) then  
	    call getarg(1,fname)
	 else if (i.eq.2) then
	    call getarg(1,fname)
	    call getarg(2,length_scale_string)
	    read (length_scale_string,*) length_scale  
	 else if (i.eq.5) then
	    call getarg(1,fname)
	    call getarg(2,length_scale_string)
	    read (length_scale_string,*) length_scale  
	    call getarg(3,perm_x_string)           
	    read (perm_x_string,*) perm(1)  
	    call getarg(4,perm_y_string)           
	    read (perm_y_string,*) perm(2)  
	    call getarg(5,perm_z_string)           
	    read (perm_z_string,*) perm(3)  
	  else
	    write(*,*) 'Wrong input'
	    stop
	  endif    
    
 	 open(meshf, file = fname, status = 'old')
 	 write(*,*) "Reading ", fname   
  
 	 read(meshf,*) ctmp, itmp  
 	 read(meshf,*) ctmp
 	 read(meshf,*) NSD
  	write(*,*) "Dimension", NSD
  
! Vertices
  	read(meshf,*) ctmp  
 	 read(meshf,*) NNode
 	 write(*,*) ctmp, " ", NNode
 	 allocate(xg(NNode,NSD))
	  do i = 1, NNode
	    read(meshf,*) ( xg(i,j), j = 1, NSD), itmp    
	  end do      
	  read(meshf,*) ctmp

! Edges
	  if (ctmp.eq.'Edges') then  
	    read(meshf,*) NEdge
	    write(*,*) ctmp, " ", Nedge
	    do i = 1, NEdge
	      read(meshf,*) itmp,itmp,itmp
	    end do    
	    read(meshf,*) ctmp
	  endif

! Triangles  (in temporary variables)
	  if (ctmp.eq.'Triangles') then
	    read(meshf,*) NFace
	    write(*,*) ctmp, " ", NFace
	    allocate(Face_IEN(NFACE,NSHLb))
	    allocate(FaceID(NFace))  
	    do i = 1, NFace
	      read(meshf,*) (Face_IEN(i,j), j = 1, NSHLb),FaceID(i)  
	    end do 
	    read(meshf,*) ctmp
	  endif
    
! Tetrahedra  
	  read(meshf,*) NElem
	  write(*,*) ctmp, " ", NElem  
	  allocate(IEN(NElem,NSHL)) 
	  allocate(IEMAT(NElem))
	  do i = 1,NElem
	    IEMAT(i) = 0
	  enddo
	  do i = 1,NElem
   	     read(meshf,*) IEN(i,1),IEN(i,2),IEN(i,4),IEN(i,3), itmp
   	     if (itmp.ge.0) then
      	       IEMAT(i)=itmp
  	     endif
 	  end do 
  
  
!==================================
! output elements/nodes/faces
! tetgen format 
!==================================
  	meshf = 11      
  	fname = 'mesh.ele' 
  	open(meshf, file = fname, status = 'unknown')
  	write(*,*) "Writing ", fname  
  	write(meshf,*) NElem, NSHL, 1
  	do i = 1, NElem
	    write(meshf,'(6I11)') i, (IEN(i,j), j = 1, NSHL), IEMAT(i)
	end do      
	close(meshf)

  	fname = 'mesh.node'
  	write(*,*) "Writing ", fname  
  	open(meshf, file = fname, status = 'unknown')
  	write(meshf,*) NNode, NSD, 0,0 
  	do i = 1, NNode
    	  write(meshf,'(I11,x, 3E17.9)') i, (  xg(i,perm(j))*length_scale, j = 1, NSD)
  	end do  
  	close(meshf)

  	fname = 'mesh.face'
  	write(*,*) "Writing ", fname   
  	open(meshf, file = fname, status = 'unknown')
  	write(meshf,*) NFace,1  
  	do i = 1, NFace
 	   write(meshf,'(5I11)') i,(Face_IEN(i,j), j=1,NSHLb),FaceID(i)
  	end do       
  	close(meshf)

!==============================
! output mesh in vtk
!==============================
  	write(*,*) "Generate mesh.vtk file ..."
	meshf = 99
	fname = trim('mesh.vtk')
	open(meshf, file=fname, status='unknown', form='formatted')
      
	write(meshf,'(a)') '# vtk DataFile Version 3.0'
	write(meshf,'(a)') 'vtk output'
  	write(meshf,'(a)') 'ASCII'
	write(meshf,'(a)') 'DATASET UNSTRUCTURED_GRID'
  
 	write(meshf,'(a,x,I11,x,a)') 'POINTS ',NNode, 'float'

  	do i = 1, NNode
  	  write(meshf,'(3E17.8)') (real(xg(i,j),4), j = 1, NSD)    
  	end do

  	write(meshf,'(a,x,I11,x,I11)') 'CELLS ',NElem, NElem*5
  	do i = 1, NElem
  	  write(meshf,'(5I11)') 4, (IEN(i,j)-1, j = 1, 4)
  	end do

  	write(meshf,'(a,x,I11)') 'CELL_TYPES',NElem
  	do i = 1, NElem
  	  write(meshf,'(I11)') 10
  	enddo

  	close(meshf)
  
!==============================
! output boundary in vtk
!==============================
  	allocate(hFace_IEN(NFACE,NSHLb))
   
  	do FACE_ID = 1, maxval(faceID)
    
  	   hNFace  = 0
   	   do i = 1, NFace
      	      if (FaceID(i).eq.FACE_ID) then
       	         hNFace =  hNFace + 1   
          	 Face_IEN(hNFace,:) = Face_IEN(i,:)
      	      endif  
    	   enddo
    
    	   if (hNFace.ge.1) then
    	      write(*,*) "Generate hull.vtk file ..."
      	      meshf = 99
    	      write(fname,'(I11)') FACE_ID
      	      fname = trim('face.'// trim(adjustl(fname))//'.vtk')
      	      open(meshf, file=fname, status='unknown', form='formatted')
      
     	      write(meshf,'(a)') '# vtk DataFile Version 3.0'
      	      write(meshf,'(a)') 'vtk output'
      	      write(meshf,'(a)') 'ASCII'
       	      write(meshf,'(a)') 'DATASET UNSTRUCTURED_GRID'
  
      	      write(meshf,'(a,x,I11,x,a)') 'POINTS ',NNode, 'float'

     	      do i = 1, NNode
     	         write(meshf,'(3E17.8)') (real(xg(i,j)), j = 1, NSD)    
     	      end do

     	      write(meshf,'(a,x,I11,x,I11)') 'CELLS ',hNFace, hNFace*4
     	      do i = 1, hNFace
      	        write(meshf,'(4I11)') 3, (hFace_IEN(i,j)-1, j = 1, NSHLb)
     	      end do

      	      write(meshf,'(a,x,I11)') 'CELL_TYPES',hNFace
      	      do i = 1, hNFace
       	         write(meshf,'(I11)') 5
      	      enddo

      	      close(meshf)
    
   	  endif
  	enddo
  
!==============================
      
  	deallocate(IEN,Face_IEN,FaceID,xg,IEMAT)
   
	end program gmsh2tetgen_withtags
