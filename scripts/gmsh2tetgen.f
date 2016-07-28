program gmsh2tetgen
!=============================================================
! program to read in the original input(*.mesh) file and 
! translates to tetgen format
! also output vtk files for visualization of generated mesh 
!=============================================================
  implicit none
      
  integer :: i, j, ier, meshf, itmp
  integer :: NSHL, NSHLb, NSD, NNode, NElem, NEdge, EType
  real(kind=8), allocatable :: xg(:,:)
  real(kind=8) theta, xtmp(3), st,ct
  integer :: rNFace, NFACE, hNFace, FACE_ID
  integer, allocatable ::  Face_IEN(:,:), FaceID(:)
  integer, allocatable :: hFace_IEN(:,:)
  integer, allocatable :: IEN(:,:), o2n(:), n2o(:)
  logical, allocatable :: check(:)
  
  character(len=30) :: fname
  character(len=40) :: ctmp

  EType = 2
  NSHL  = 4
  NSHLb = 3

  !==================================
  ! Read mesh in .mesh format
  !==================================
  meshf = 11      
  i = iargc()
  if (i.eq.1) then  
    call getarg(1,fname)
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
  do i = 1,NElem
    read(meshf,*) IEN(i,1),IEN(i,2),IEN(i,4),IEN(i,3), itmp
  end do 
  
  
  !==================================
  ! Remove useless nodes
  !==================================  
  allocate(check(NNode),o2n(NNode),n2o(NNode))
  check=.false.
  do i = 1,NElem
    check(IEN(i,1)) = .true. 
    check(IEN(i,2)) = .true. 
    check(IEN(i,4)) = .true. 
    check(IEN(i,3)) = .true. 
  end do   
  
  itmp =0
  do i = 1, NNode
    if (.not.check(i)) then
       write(*,*) "Node: ", i
    else
      itmp = itmp  + 1
      n2o(itmp) = i
      o2n(i)    = itmp
    endif 
    
  end do   
  write(*,*) "Gross nodes: ",NNode       
  NNode = itmp      	
  write(*,*) "Net nodes: ",NNode       
	
  !==================================
  ! output elements/nodes/faces
  ! tetgen format 
  !==================================
  meshf = 11      
  fname = 'mesh.ele' 
  open(meshf, file = fname, status = 'unknown')
  write(*,*) "Writing ", fname  
  write(meshf,*) NElem, NSHL, 0         
  do i = 1, NElem
    write(meshf,'(5I8)') i, (o2n(IEN(i,j)), j = 1, NSHL)
  end do      
  close(meshf)

  fname = 'mesh.node'
  write(*,*) "Writing ", fname  
  open(meshf, file = fname, status = 'unknown')
  write(meshf,*) NNode, NSD, 0,0 
  do i = 1, NNode
    write(meshf,'(I8,x, 3E17.9)') i, (  xg(n2o(i),j)/1000.0, j = 1, NSD)
  end do  
  close(meshf)

  fname = 'mesh.face'
  write(*,*) "Writing ", fname   
  open(meshf, file = fname, status = 'unknown')
  write(meshf,*) NFace,1  
  do i = 1, NFace
    write(meshf,'(5I8)') i, (o2n(Face_IEN(i,j)), j = 1, NSHLb), FaceID(i)
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
  
  write(meshf,'(a,x,I10,x,a)') 'POINTS ',NNode, 'float'

  do i = 1, NNode
    write(meshf,'(3E17.8)') (real(xg(n2o(i),j),4), j = 1, NSD)    
  end do

  write(meshf,'(a,x,I10,x,I10)') 'CELLS ',NElem, NElem*5
  do i = 1, NElem
    write(meshf,'(5I10)') 4, (o2n(IEN(i,j))-1, j = 1, 4)
  end do

  write(meshf,'(a,x,I10)') 'CELL_TYPES',NElem
  do i = 1, NElem
    write(meshf,'(I10)') 10
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
        hFace_IEN(hNFace,:) = Face_IEN(i,:)
      endif  
    enddo
    
    if (hNFace.ge.1) then
      write(*,*) "Generate hull.vtk file ..."
      meshf = 99
      write(fname,'(I8)') FACE_ID
      fname = trim('face.'// trim(adjustl(fname))//'.vtk')
      open(meshf, file=fname, status='unknown', form='formatted')
      
      write(meshf,'(a)') '# vtk DataFile Version 3.0'
      write(meshf,'(a)') 'vtk output'
      write(meshf,'(a)') 'ASCII'
      write(meshf,'(a)') 'DATASET UNSTRUCTURED_GRID'
  
      write(meshf,'(a,x,I10,x,a)') 'POINTS ',NNode, 'float'

      do i = 1, NNode
        write(meshf,'(3E17.8)') (real(xg(n2o(i),j)), j = 1, NSD)    
      end do

      write(meshf,'(a,x,I10,x,I10)') 'CELLS ',hNFace, hNFace*4
      do i = 1, hNFace
        write(meshf,'(4I10)') 3, (o2n(hFace_IEN(i,j))-1, j = 1, NSHLb)
      end do

      write(meshf,'(a,x,I10)') 'CELL_TYPES',hNFace
      do i = 1, hNFace
        write(meshf,'(I10)') 5
      enddo

      close(meshf)
    
    endif
  enddo
  
  !==============================
      
  deallocate(IEN,Face_IEN,FaceID,xg)
   
end program gmsh2tetgen
