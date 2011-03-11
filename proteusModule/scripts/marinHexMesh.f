!------------------------------------------------------------------------
!                                                                        
!        Main routine to call all the subroutines                        
!                                                                        
!------------------------------------------------------------------------
      program marinMesh 

      implicit none

      integer xsn,ysn,zsn
      real*8  xs(4),ys(4),zs(3)  
      real*8  dx(3),dy(3),dz(2)  
      integer xn(3),yn(3),zn(2) 
      integer xnn,ynn,znn 
      
      integer meshf
      logical ex

      real*8, allocatable :: x(:),y(:),z(:)

      integer NNODZ, NEL
      real*8, allocatable :: xg(:,:)
      integer, allocatable :: IEN(:,:)

      integer NNODZ2, NEL2
      real*8, allocatable :: xg2(:,:)
      integer, allocatable :: IEN2(:,:)

      xsn=4
      xs(1)=0.0d0
      xs(2)=2.3955d0  
      xs(3)=2.5565d0
      xs(4)=3.22d0

      ysn=4
      ys(1)=0.0d0
      ys(2)=0.2985d0
      ys(3)=0.7015d0
      ys(4)=1.0d0

      zsn=3
      zs(1)=0.0d0
      zs(2)=0.161d0
      zs(3)=1.0d0


      meshf=11
      inquire(file='size.mesh', exist=ex) 

      
      if (ex) then
        open(meshf, file='size.mesh', status='unknown')  
	read(meshf,*) dx
	read(meshf,*) dy
	read(meshf,*) dz	
	close(meshf)
      else
        dx(:) = 0.15d0
        dy(:) = 0.15d0
        dz(:) = 0.15d0
      endif


      write(*,*) "Get Refinement"
      call getrefine(xs,dx,xsn,xn)
      call getrefine(ys,dy,ysn,yn)
      call getrefine(zs,dz,zsn,zn)
      xnn = sum(xn)+1
      ynn = sum(yn)+1
      znn = sum(zn)+1
      Write(*,*) "Total = ",xnn,ynn,znn
      allocate (x(xnn),y(ynn),z(znn)   )
      call dorefine(xs,xn,xsn,x)
      call dorefine(ys,yn,ysn,y)
      call dorefine(zs,zn,zsn,z)

  
      NNODZ = xnn*ynn*znn
      NEL   = sum(xn)*sum(yn)*sum(zn)
      write(*,*) "Mesh = ", NNODZ,NEL 
      allocate( xg(NNODZ,3), IEN(NEL,8))
      call hexahedralize(x,y,z,xnn,ynn,znn,xg,NNODZ,IEN,NEL)

      if (0.eq.1) then
        call writeMesh(xg,NNODZ,IEN,NEL)
        deallocate(x,y,z,xg,IEN)
      else
        allocate( xg2(NNODZ,3), IEN2(NEL,8))
        call removeObstacle(xs(2),xs(3),ys(2),ys(3),zs(1)-1d0,zs(2),    &
                            xg,NNODZ,IEN,NEL,                           &
                            xg2,NNODZ2,IEN2,NEL2)
  
        call writeMesh(xg2(:NNODZ2,:),NNODZ2,IEN2(:NEL2,:),NEL2)
        deallocate(x,y,z,xg,IEN,xg2,IEN2)
      endif

      end program marinMesh
!-------------------------------------------
  
!-------------------------------------------
      subroutine getrefine(xs,dx,xsn,xn)

      implicit none

      integer xsn
      real*8  xs(xsn),dx(xsn-1)
      integer xn(xsn-1)
     
      integer i,j,n

      do i=1,xsn-1 
        xn(i) =int( ceiling((xs(i+1)-xs(i))/dx(i)))
        write(*,*) xs(i),xs(i+1),xn(i)
      enddo

      end subroutine

!-------------------------------------------
  
!-------------------------------------------
      subroutine dorefine(xs,xn,xsn,x)

      implicit none

      integer xsn
      real*8  xs(xsn)
      integer xn(xsn-1)
      real*8 x(sum(xn)+1)
     
      integer i,j,n
    
      n = 0
      do i=1,xsn-1
        write(*,*) i,xn(i)
        do j=1,xn(i)
          n=n+1
          x(n)=xs(i)+(real(j-1)/(real(xn(i))))*(xs(i+1)-xs(i))
          write(*,*) n, x(n)
        enddo
      enddo
      x(n+1) = xs(xsn)
     

      end subroutine

!-------------------------------------------

!-------------------------------------------
      subroutine hexahedralize(x,y,z,xn,yn,zn,xg,NNODZ,IEN,NEL)

      implicit none

      integer xn,yn,zn,NNODZ,NEL
      real*8  x(xn), y(yn), z(zn)
      real*8  xg(NNODZ,3)
      integer IEN(NEL,8)
      integer i,j,k,nn,nxy,jj

      nn=0 
      do k=1,zn
        do j=1,yn
          do i=1,xn
            nn=nn+1
            xg(nn,1) = x(i)
            xg(nn,2) = y(j)
            xg(nn,3) = z(k)    
          enddo
        enddo
      enddo  
      
      nxy = xn*yn
      nn=0 
      do k=0,xn-2
        do j=0,yn-2
          do i=0,zn-2
            nn=nn+1

	    IEN(nn,1) = (k+0) + (j+0)*(xn) + (i+0)*nxy
   	    IEN(nn,2) = (k+1) + (j+0)*(xn) + (i+0)*nxy
            IEN(nn,3) =	(k+1) + (j+1)*(xn) + (i+0)*nxy
            IEN(nn,4) = (k+0) + (j+1)*(xn) + (i+0)*nxy
            IEN(nn,5) =	(k+0) + (j+0)*(xn) + (i+1)*nxy
            IEN(nn,6) = (k+1) + (j+0)*(xn) + (i+1)*nxy
            IEN(nn,7) = (k+1) + (j+1)*(xn) + (i+1)*nxy
            IEN(nn,8) = (k+0) + (j+1)*(xn) + (i+1)*nxy

          enddo
        enddo
      enddo  

      end subroutine

!-------------------------------------------

!-------------------------------------------
      subroutine removeObstacle(x0,x1,y0,y1,z0,z1,   &
                                xg,NNODZ,IEN,NEL,    &
                                xg2,NNODZ2,IEN2,NEL2)

      implicit none

      real*8  x0,x1,y0,y1,z0,z1
      integer NNODZ,NEL
      real*8  xg(NNODZ,3)
      integer IEN(NEL,8)
      integer NNODZ2,NEL2
      real*8  xg2(NNODZ,3)
      integer IEN2(NEL,8)

      integer i,j, o2n(NNODZ)
      real*8  xn(3)

      NNODZ2=0
      do i = 1, NNODZ
  
         if ((xg(i,1).gt.x0).and.(xg(i,1).lt.x1).and. &
             (xg(i,2).gt.y0).and.(xg(i,2).lt.y1).and. &
             (xg(i,3).gt.z0).and.(xg(i,3).lt.z1)) then
           write(*,*) "Removing point = ", xg(i,:)
           o2n(i) = -1
         else
           NNODZ2 = NNODZ2+1           
           xg2(NNODZ2,:) = xg(i,:)
           o2n(i)=NNODZ2-1
         endif    
  
      end do

      NEL2 = 0
      do i = 1, NEL
        if ( (o2n(IEN(i,1)+1).ge.0).and.(o2n(IEN(i,2)+1).ge.0).and. &
             (o2n(IEN(i,3)+1).ge.0).and.(o2n(IEN(i,4)+1).ge.0).and. &
             (o2n(IEN(i,5)+1).ge.0).and.(o2n(IEN(i,6)+1).ge.0).and. &
             (o2n(IEN(i,7)+1).ge.0).and.(o2n(IEN(i,8)+1).ge.0)) then
          xn = (xg(IEN(i,1)+1,:) + xg(IEN(i,2)+1,:)   & 
             +  xg(IEN(i,3)+1,:) + xg(IEN(i,4)+1,:)   &
             +  xg(IEN(i,5)+1,:) + xg(IEN(i,6)+1,:)   &
             +  xg(IEN(i,7)+1,:) + xg(IEN(i,8)+1,:))/8d0   
         
          if ((xn(1).gt.x0).and.(xn(1).lt.x1).and.   &
              (xn(2).gt.y0).and.(xn(2).lt.y1).and.   &
              (xn(3).gt.z0).and.(xn(3).lt.z1)) then
           write(*,*) "  Removing Element ", i, " (based on coord!)"        
          else
            NEL2 = NEL2 + 1
            IEN2(NEL2,1)=o2n(IEN(i,1)+1)
            IEN2(NEL2,2)=o2n(IEN(i,2)+1)
            IEN2(NEL2,3)=o2n(IEN(i,3)+1)
            IEN2(NEL2,4)=o2n(IEN(i,4)+1)
            IEN2(NEL2,5)=o2n(IEN(i,5)+1)
            IEN2(NEL2,6)=o2n(IEN(i,6)+1)
            IEN2(NEL2,7)=o2n(IEN(i,7)+1)
            IEN2(NEL2,8)=o2n(IEN(i,8)+1)
         endif
        !else
        !  write(*,*) "  Removing Element", i
        endif
      end do


      end subroutine

!-------------------------------------------

!-------------------------------------------
      subroutine writeMesh(xg,NNODZ,IEN,NEL)

      implicit none

      integer NNODZ,NEL
      real*8  xg(NNODZ,3)
      integer IEN(NEL,8)
      integer meshf,i

      meshf=11
      open(meshf, file='marinHex.mesh', status='replace')

      write(meshf,*) "HEX"
      write(meshf,*) NNODZ,NEL
      do i = 1, NNODZ
        write(meshf,*)   xg(i,1),   xg(i,2),   xg(i,3)
      end do

      do i = 1, NEL
        write(meshf,'(9I8)') IEN(i,1), IEN(i,2), IEN(i,3), IEN(i,4),IEN(i,5), IEN(i,6), IEN(i,7), IEN(i,8), 0
      end do

      end subroutine

!-------------------------------------------
