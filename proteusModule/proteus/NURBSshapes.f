

      subroutine eval_SHAPE_3D(
     &  NSD,P,Q,R,NSHL,MCP,NCP,OCP,NEL,
     &  U_KNOT,V_KNOT,W_KNOT,B_NET,
     &  ni,nj,nk,u_hat,v_hat,w_hat,nders,
     &  shl,shgradl,shgradg,shhessg,dxidx,DetJ
     &  )

      implicit none

c ----------------------------------------------------------------------
c...  Define stuff that would normally be in common/adjkeep       
      integer  NSD
      
      integer  P,Q,R,NSHL
      integer  MCP,NCP,OCP,NEL
      
      real*8 U_KNOT(MCP+P+1),V_KNOT(NCP+Q+1),W_KNOT(OCP+R+1)
      real*8 B_NET(MCP,NCP,OCP,NSD+1)
      !real*8 B_NET(P+1,Q+1,R+1,NCP,OCP,NSD+1)
          
      real DetJ
  
c     ------------------------------------------------------------------
c...  Element number,hess_flag
c...  u and v coordinates of integration point in parent element
      real*8 u_hat, v_hat, w_hat, du, dv, dw
      integer e,nders
      
c...  Vector of Local basis function vales at (u_hat, v_hat), local and 
c          global gradients.
      
      
      real*8 shl(NSHL),shgradl(NSHL,NSD),shgradg(NSHL,NSD),
     &  shhessg(NSHL,NSD,NSD), tempshl(NSHL),tempshgradl(NSHL,NSD),
     &  shhessl(NSHL,6), tempshhessl(NSHL,6)

      real*8 dxdxi(NSD,NSD), dxidx(NSD,NSD),dxdxixj(NSD,6), locLHS(6,6)
      
c...  Local Variables
c    
c     1D nonrational basis functions and derivs in u and v

      real*8 N(nders+1,P+1), M(nders+1,Q+1),O(nders+1,R+1)    
            
c     u and v coordinates of integration point, denominator and derivative sums
      real*8 u, v, w, denom_sum, derv_sum_U, derv_sum_V,
     &  derv_sum_W, derv_sum_UU,
     &  derv_sum_UV, derv_sum_UW,derv_sum_VV, derv_sum_VW,derv_sum_WW
      
c     NURBS coordinates, counters for loops
      integer ni, nj, nk, i, j, k, icount, aa
      
c     temporary variables
      real*8 tmp

c ----------------------------------------------------------------------
          

c     Get u and v coordinates of integration point
      
      u = ((U_KNOT(ni+1)-U_KNOT(ni))*u_hat + 
     &     U_KNOT(ni+1) + U_KNOT(ni))/2d+0
      v = ((V_KNOT(nj+1)-V_KNOT(nj))*v_hat + 
     &     V_KNOT(nj+1) + V_KNOT(nj))/2d+0
      w = ((W_KNOT(nk+1)-W_KNOT(nk))*w_hat + 
     &     W_KNOT(nk+1) + W_KNOT(nk))/2d+0

c     Get knot span sizes
      
      du = U_KNOT(ni+1)-U_KNOT(ni)
      dv = V_KNOT(nj+1)-V_KNOT(nj)
      dw = W_KNOT(nk+1)-W_KNOT(nk)
      
c     Evalate 1D shape functions and derivatives each direction
      
      
      call dersbasisfuns(ni,P,MCP,u,nders,U_KNOT,N) ! calculate in u dir.
      call dersbasisfuns(nj,Q,NCP,v,nders,V_KNOT,M) ! calculate in v dir.
      call dersbasisfuns(nk,R,OCP,w,nders,W_KNOT,O) ! calculate in v dir.

      
c     Form basis functions and derivatives dR/du and dR/dv
      
      icount = 0
      denom_sum = 0d+0
      derv_sum_U = 0d+0
      derv_sum_V = 0d+0
      derv_sum_W = 0d+0
     
      do k = 0, R
        do j = 0, Q
          do i = 0, P
            icount = icount+1
            
c...        basis functions
            shl(icount) = N(1,P+1-i)*M(1,Q+1-j)*O(1,R+1-k)*
     &        B_NET(ni-i,nj-j,nk-k,NSD+1)
            denom_sum = denom_sum + shl(icount)
            
c...        derivatives
            shgradl(icount,1) = 
     &        N(2,P+1-i)*M(1,Q+1-j)*O(1,R+1-k)*
     &        B_NET(ni-i,nj-j,nk-k,NSD+1) ! u
            derv_sum_U = derv_sum_U + shgradl(icount,1)
            shgradl(icount,2) = 
     &        N(1,P+1-i)*M(2,Q+1-j)*O(1,R+1-k)*
     &        B_NET(ni-i,nj-j,nk-k,NSD+1) ! v
            derv_sum_V = derv_sum_V + shgradl(icount,2)
            shgradl(icount,3) = 
     &        N(1,P+1-i)*M(1,Q+1-j)*O(2,R+1-k)*
     &        B_NET(ni-i,nj-j,nk-k,NSD+1) ! w
            derv_sum_W = derv_sum_W + shgradl(icount,3)
            
          enddo
        enddo
      enddo
      
c     Divide through by denominator
      
      tempshl = shl
      tempshgradl = shgradl
      
      
      do i = 1,NSHl
        shgradl(i,1) = shgradl(i,1)/denom_sum - 
     &    (shl(i)*derv_sum_U)/(denom_sum**2)
        shgradl(i,2) = shgradl(i,2)/denom_sum - 
     &    (shl(i)*derv_sum_V)/(denom_sum**2)
        shgradl(i,3) = shgradl(i,3)/denom_sum - 
     &    (shl(i)*derv_sum_W)/(denom_sum**2)
        shl(i) = shl(i)/denom_sum
      enddo


      
c     
c...  Now calculate gradients.
c

c     calculate dx/dxi
      


      dxdxi = 0d+0
      icount = 0
      
      do k = 0, R
        do j = 0, Q
          do i = 0, P
            icount = icount + 1
            
            dxdxi(1,1) = dxdxi(1,1) +
     &        B_NET(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,1)
            dxdxi(1,2) = dxdxi(1,2) +
     &        B_NET(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,2)
            dxdxi(1,3) = dxdxi(1,3) +
     &        B_NET(ni-i,nj-j,nk-k,1) *
     &        shgradl(icount,3)
            dxdxi(2,1) = dxdxi(2,1) +
     &        B_NET(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,1)
            dxdxi(2,2) = dxdxi(2,2) +
     &        B_NET(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,2)
            dxdxi(2,3) = dxdxi(2,3) +
     &        B_NET(ni-i,nj-j,nk-k,2) *
     &        shgradl(icount,3)
            dxdxi(3,1) = dxdxi(3,1) +
     &        B_NET(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,1)
            dxdxi(3,2) = dxdxi(3,2) +
     &        B_NET(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,2)
            dxdxi(3,3) = dxdxi(3,3) +
     &        B_NET(ni-i,nj-j,nk-k,3) *
     &        shgradl(icount,3)
            
          enddo
        enddo
      enddo

      
c
c.... comPte the inverse of deformation gradient
c
      
      dxidx(1,1) =   dxdxi(2,2) * dxdxi(3,3) 
     &     - dxdxi(3,2) * dxdxi(2,3)
      dxidx(1,2) =   dxdxi(3,2) * dxdxi(1,3) 
     &     - dxdxi(1,2) * dxdxi(3,3)
      dxidx(1,3) =  dxdxi(1,2) * dxdxi(2,3) 
     &     - dxdxi(1,3) * dxdxi(2,2)
      tmp          = 1d+0 / ( dxidx(1,1) * dxdxi(1,1) 
     &     + dxidx(1,2) * dxdxi(2,1)  
     &     + dxidx(1,3) * dxdxi(3,1) )
      dxidx(1,1) = dxidx(1,1) * tmp
      dxidx(1,2) = dxidx(1,2) * tmp
      dxidx(1,3) = dxidx(1,3) * tmp
      dxidx(2,1) = (dxdxi(2,3) * dxdxi(3,1) 
     &     - dxdxi(2,1) * dxdxi(3,3)) * tmp
      dxidx(2,2) = (dxdxi(1,1) * dxdxi(3,3) 
     &     - dxdxi(3,1) * dxdxi(1,3)) * tmp
      dxidx(2,3) = (dxdxi(2,1) * dxdxi(1,3) 
     &     - dxdxi(1,1) * dxdxi(2,3)) * tmp
      dxidx(3,1) = (dxdxi(2,1) * dxdxi(3,2) 
     &  - dxdxi(2,2) * dxdxi(3,1)) * tmp
      dxidx(3,2) = (dxdxi(3,1) * dxdxi(1,2) 
     &     - dxdxi(1,1) * dxdxi(3,2)) * tmp
      dxidx(3,3) = (dxdxi(1,1) * dxdxi(2,2) 
     &     - dxdxi(1,2) * dxdxi(2,1)) * tmp
      
      DetJ = 1d+0/tmp           ! Note that DetJ resides in common


      

      do i = 1, NSHl

        shgradg(i,1) = shgradl(i,1) * dxidx(1,1) + 
     &    shgradl(i,2) * dxidx(2,1) +
     &    shgradl(i,3) * dxidx(3,1)
        shgradg(i,2) = shgradl(i,1) * dxidx(1,2) + 
     &    shgradl(i,2) * dxidx(2,2) +
     &    shgradl(i,3) * dxidx(3,2) 
        shgradg(i,3) = shgradl(i,1) * dxidx(1,3) + 
     &    shgradl(i,2) * dxidx(2,3) +
     &    shgradl(i,3) * dxidx(3,3) 
        
      enddo

      
      if (nders.eq.2) then
        
        icount = 0
        derv_sum_UU = 0d+0
        derv_sum_UV = 0d+0
        derv_sum_UW = 0d+0
        derv_sum_VV = 0d+0
        derv_sum_VW = 0d+0
        derv_sum_WW = 0d+0
        do k = 0, R
          do j = 0,Q
            do i = 0,P
              icount = icount+1
              
c...          2nd derivatives
              tempshhessl(icount,1) = 
     &          N(3,P+1-i)*M(1,Q+1-j)*O(1,R+1-k)*
     &          B_NET(ni-i,nj-j,nk-k,NSD+1) ! ,uu
              derv_sum_UU = derv_sum_UU + tempshhessl(icount,1)
              tempshhessl(icount,2) = 
     &          N(2,P+1-i)*M(2,Q+1-j)*O(1,R+1-k)*
     &          B_NET(ni-i,nj-j,nk-k,NSD+1) ! ,uv
              derv_sum_UV = derv_sum_UV + tempshhessl(icount,2)
              tempshhessl(icount,3) = 
     &          N(2,P+1-i)*M(1,Q+1-j)*O(2,R+1-k)*
     &          B_NET(ni-i,nj-j,nk-k,NSD+1) ! ,uw
              derv_sum_UW = derv_sum_UW + tempshhessl(icount,3)
              tempshhessl(icount,4) = 
     &          N(1,P+1-i)*M(3,Q+1-j)*O(1,R+1-k)*
     &          B_NET(ni-i,nj-j,nk-k,NSD+1) ! ,vv
              derv_sum_VV = derv_sum_VV + tempshhessl(icount,4)
              tempshhessl(icount,5) = 
     &          N(1,P+1-i)*M(2,Q+1-j)*O(2,R+1-k)*
     &          B_NET(ni-i,nj-j,nk-k,NSD+1) ! ,vw
              derv_sum_VW = derv_sum_VW + tempshhessl(icount,5)
              tempshhessl(icount,6) = 
     &          N(1,P+1-i)*M(1,Q+1-j)*O(3,R+1-k)*
     &          B_NET(ni-i,nj-j,nk-k,NSD+1) ! ,ww
              derv_sum_WW = derv_sum_WW + tempshhessl(icount,6)
            enddo
          enddo
        enddo
        

c...    local Hessian
        do i = 1,NSHl
          shhessl(i,1) = tempshhessl(i,1)/denom_sum -
     &      tempshgradl(i,1)*derv_sum_U/(denom_sum**2) -
     &      ((tempshgradl(i,1)*derv_sum_U+tempshl(i)*derv_sum_UU)/
     &      (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_U*derv_sum_U/
     &      (denom_sum**3))
          shhessl(i,2) = tempshhessl(i,2)/denom_sum -
     &      tempshgradl(i,1)*derv_sum_V/(denom_sum**2) -
     &      ((tempshgradl(i,2)*derv_sum_U+tempshl(i)*derv_sum_UV)/
     &      (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_U*derv_sum_V/
     &      (denom_sum**3))
          shhessl(i,3) = tempshhessl(i,3)/denom_sum -
     &      tempshgradl(i,1)*derv_sum_W/(denom_sum**2) -
     &      ((tempshgradl(i,3)*derv_sum_U+tempshl(i)*derv_sum_UW)/
     &      (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_U*derv_sum_W/
     &      (denom_sum**3))         
          shhessl(i,4) = tempshhessl(i,4)/denom_sum -
     &      tempshgradl(i,2)*derv_sum_V/(denom_sum**2) -
     &      ((tempshgradl(i,2)*derv_sum_V+tempshl(i)*derv_sum_VV)/
     &      (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_V*derv_sum_V/
     &      (denom_sum**3))
          shhessl(i,5) = tempshhessl(i,5)/denom_sum -
     &      tempshgradl(i,2)*derv_sum_W/(denom_sum**2) -
     &      ((tempshgradl(i,3)*derv_sum_V+tempshl(i)*derv_sum_VW)/
     &      (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_V*derv_sum_W/
     &      (denom_sum**3))  
          shhessl(i,6) = tempshhessl(i,6)/denom_sum -
     &      tempshgradl(i,3)*derv_sum_W/(denom_sum**2) -
     &      ((tempshgradl(i,3)*derv_sum_W+tempshl(i)*derv_sum_WW)/
     &      (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_W*derv_sum_W/
     &      (denom_sum**3))
        enddo
        
c...    global Hessian

c...    Second derivatives of the geometrical map

        dxdxixj = 0d+0
        icount = 0
        
        do k = 0, R
          do j = 0, Q
            do i = 0, P
              icount = icount + 1
              
              dxdxixj(1,:) = dxdxixj(1,:)+
     &          B_NET(ni-i,nj-j,nk-k,1)*
     &          shhessl(icount,:)
              
              dxdxixj(2,:) = dxdxixj(2,:)+
     &          B_NET(ni-i,nj-j,nk-k,2)*
     &          shhessl(icount,:)
              
              dxdxixj(3,:) = dxdxixj(3,:)+
     &          B_NET(ni-i,nj-j,nk-k,3)*
     &          shhessl(icount,:)
              
              
            enddo
          enddo
        enddo
        
c...    RHS of the matrix eQation for the second derivatives of bases.
c...    Reuse local hess. array        

        shhessl(:,1) = shhessl(:,1) - shgradg(:,1)*dxdxixj(1,1) -
     &    shgradg(:,2)*dxdxixj(2,1) - shgradg(:,3)*dxdxixj(3,1)
        
        shhessl(:,2) = shhessl(:,2) - shgradg(:,1)*dxdxixj(1,2) -
     &    shgradg(:,2)*dxdxixj(2,2) - shgradg(:,3)*dxdxixj(3,2)
        
        shhessl(:,3) = shhessl(:,3) - shgradg(:,1)*dxdxixj(1,3) -
     &    shgradg(:,2)*dxdxixj(2,3) - shgradg(:,3)*dxdxixj(3,3)
        
        shhessl(:,4) = shhessl(:,4) - shgradg(:,1)*dxdxixj(1,4) -
     &    shgradg(:,2)*dxdxixj(2,4) - shgradg(:,3)*dxdxixj(3,4)
        
        shhessl(:,5) = shhessl(:,5) - shgradg(:,1)*dxdxixj(1,5) -
     &    shgradg(:,2)*dxdxixj(2,5) - shgradg(:,3)*dxdxixj(3,5)
        
        shhessl(:,6) = shhessl(:,6) - shgradg(:,1)*dxdxixj(1,6) -
     &    shgradg(:,2)*dxdxixj(2,6) - shgradg(:,3)*dxdxixj(3,6)
        
        
c...    LHS (6x6, same for every basis function)

        locLHS(1,1) = dxdxi(1,1)*dxdxi(1,1)
        locLHS(1,2) = 2d+0*dxdxi(1,1)*dxdxi(2,1)
        locLHS(1,3) = 2d+0*dxdxi(1,1)*dxdxi(3,1)
        locLHS(1,4) = dxdxi(2,1)*dxdxi(2,1)
        locLHS(1,5) = 2d+0*dxdxi(2,1)*dxdxi(3,1)
        locLHS(1,6) = dxdxi(3,1)*dxdxi(3,1)

        locLHS(2,1) = dxdxi(1,1)*dxdxi(1,2)
        locLHS(2,2) = dxdxi(1,1)*dxdxi(2,2) + dxdxi(1,2)*dxdxi(2,1)
        locLHS(2,3) = dxdxi(1,1)*dxdxi(3,2) + dxdxi(1,2)*dxdxi(3,1)
        locLHS(2,4) = dxdxi(2,1)*dxdxi(2,2)
        locLHS(2,5) = dxdxi(2,1)*dxdxi(3,2) + dxdxi(2,2)*dxdxi(3,1)
        locLHS(2,6) = dxdxi(3,1)*dxdxi(3,2)

        locLHS(3,1) = dxdxi(1,1)*dxdxi(1,3)
        locLHS(3,2) = dxdxi(1,1)*dxdxi(2,3) + dxdxi(1,3)*dxdxi(2,1)
        locLHS(3,3) = dxdxi(1,1)*dxdxi(3,3) + dxdxi(1,3)*dxdxi(3,1)
        locLHS(3,4) = dxdxi(2,1)*dxdxi(2,3)
        locLHS(3,5) = dxdxi(2,1)*dxdxi(3,3) + dxdxi(2,3)*dxdxi(3,1)
        locLHS(3,6) = dxdxi(3,1)*dxdxi(3,3)

        locLHS(4,1) = dxdxi(1,2)*dxdxi(1,2)
        locLHS(4,2) = 2d+0*dxdxi(1,2)*dxdxi(2,2)
        locLHS(4,3) = 2d+0*dxdxi(1,2)*dxdxi(3,2)
        locLHS(4,4) = dxdxi(2,2)*dxdxi(2,2)
        locLHS(4,5) = 2d+0*dxdxi(2,2)*dxdxi(3,2)
        locLHS(4,6) = dxdxi(3,2)*dxdxi(3,2)

        locLHS(5,1) = dxdxi(1,2)*dxdxi(1,3)
        locLHS(5,2) = dxdxi(1,2)*dxdxi(2,3) + dxdxi(1,3)*dxdxi(2,2)
        locLHS(5,3) = dxdxi(1,2)*dxdxi(3,3) + dxdxi(1,3)*dxdxi(3,2)
        locLHS(5,4) = dxdxi(2,2)*dxdxi(2,3)
        locLHS(5,5) = dxdxi(2,2)*dxdxi(3,3) + dxdxi(2,3)*dxdxi(3,2)
        locLHS(5,6) = dxdxi(3,2)*dxdxi(3,3)

        locLHS(6,1) = dxdxi(1,3)*dxdxi(1,3)
        locLHS(6,2) = 2d+0*dxdxi(1,3)*dxdxi(2,3)
        locLHS(6,3) = 2d+0*dxdxi(1,3)*dxdxi(3,3)
        locLHS(6,4) = dxdxi(2,3)*dxdxi(2,3)
        locLHS(6,5) = 2d+0*dxdxi(2,3)*dxdxi(3,3)
        locLHS(6,6) = dxdxi(3,3)*dxdxi(3,3)

ccc        write(*,*)"lhs"
ccc        write(*,*) locLHS(:,:)
ccc        write(*,*)"rhs"
ccc        write(*,*) shhessl(1,:)
ccc        write(*,*)"ans"
        
c...    (6x6) - Gaussian elimination
        
        do k = 1,6
          do i = k+1,6
            tmp = locLHS(i,k)/locLHS(k,k)
            do j = k+1,6
              locLHS(i,j) = locLHS(i,j) - tmp*locLHS(k,j)
            enddo
            shhessl(:,i) = shhessl(:,i) - tmp*shhessl(:,k)
          enddo
        enddo
        
        do i = 6,1,-1
          do j = i+1,6
            shhessl(:,i) = shhessl(:,i) - locLHS(i,j)*shhessl(:,j)
          enddo
          shhessl(:,i) = shhessl(:,i)/locLHS(i,i)
        enddo

ccc        write(*,*) shhessl(1,:)
ccc        write(*,*) "end"
        
c...    Assign to global hessian of basis functions
        
        shhessg(:,1,1) = shhessl(:,1)
        shhessg(:,1,2) = shhessl(:,2)
        shhessg(:,1,3) = shhessl(:,3)
        
        shhessg(:,2,1) = shhessl(:,2)
        shhessg(:,2,2) = shhessl(:,4)
        shhessg(:,2,3) = shhessl(:,5)
        
        shhessg(:,3,1) = shhessl(:,3)
        shhessg(:,3,2) = shhessl(:,5)
        shhessg(:,3,3) = shhessl(:,6)
        
      endif     
      
      dxidx(1,:) = dxidx(1,:)*2d+0/du
      dxidx(2,:) = dxidx(2,:)*2d+0/dv
      dxidx(3,:) = dxidx(3,:)*2d+0/dw

      
      DetJ = abs(DetJ)
      
      return
      end
      


      
      
!#######################################################################
      subroutine dersbasisfuns(i,pl,ml,u,nders,u_knotl,ders)
      
      implicit none

c     --------------VARIABLE DECLARATIONS-------------------------------
c...  knot span, degree of curve, number of control points, counters
      integer i, pl, ml, j, r, k, j1, j2, s1, s2, rk, pk,nders
c     parameter vale, vector of knots, derivative matrix
      real*8 u, u_knotl(pl+ml+1), ders(nders+1,pl+1), ndu(pl+1,pl+1),d

      real*8 left(pl+1), right(pl+1), saved, temp, a(2,pl+1)

c     ------------------------------------------------------------------
      
      ndu(1,1) = 1
      do j = 1,pl
         left(j+1) = u - u_knotl(i+1-j)
         right(j+1) = u_knotl(i+j) - u
         saved = 0
         do r = 0,j-1
            ndu(j+1,r+1) = right(r+2) + left(j-r+1)
            temp = ndu(r+1,j)/ndu(j+1,r+1)
            ndu(r+1,j+1) = saved + right(r+2)*temp
            saved = left(j-r+1)*temp
         enddo
         ndu(j+1,j+1) = saved
      enddo
      
                                ! load basis functions
      do j = 0,pl
         ders(1,j+1) = ndu(j+1,pl+1)
      enddo
                                ! comPte derivatives
      do r = 0,pl ! loop over function index
         s1 = 0
         s2 = 1                 ! alternate rows in array a
         a(1,1) = 1
                                ! loop to comPte kth derivative
         do k = 1,nders
            d = 0d+0
            rk = r-k
            pk = pl-k
            if (r >= k) then
               a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1)
               d = a(s2+1,1)*ndu(rk+1,pk+1)
            endif
            if (rk >= -1) then
               j1 = 1
            else 
               j1 = -rk
            endif
            if ((r-1) <= pk) then
               j2 = k-1
            else 
               j2 = pl-r
            endif
            do j = j1,j2
               a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1)
               d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1)
            enddo
            if (r <= pk) then
               a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1)
               d = d + a(s2+1,k+1)*ndu(r+1,pk+1)
            endif
            ders(k+1,r+1) = d
            j = s1
            s1 = s2
            s2 = j              ! switch rows
         enddo
      enddo
      
c     Multiply through by the correct factors
      r = pl
      do k = 1,nders
         do j = 0,pl
            ders(k+1,j+1) = ders(k+1,j+1)*r
         enddo
         r = r*(pl-k)
      enddo

      return
      end
