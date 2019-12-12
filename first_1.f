!  || Shree Geneshay Namah ||
! User Element subroutine for Strain Gradient Elasticity model.
! Dhaval Rasheshkumar Patel
! dhavalrpatel2511@gmail.com
c***********************************************************************
      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,
     1 nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,
     2 jelem,params,ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,  
     4 ddlmag,mdload,pnewdt,jprops,njpro,period)
c     
     include 'aba_param.inc' 
c     
      dimension rhs(mlvarx,*),amatrx(ndofel,ndofel),props(*),svars(*),
     1 energy(*),coords(mcrd,nnode),u(ndofel),du(mlvarx,*),v(ndofel),
     2 a(ndofel),time(2),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*)
c
      parameter (ndim=2, ndof=5, ndi=3, nshr=1, nnodemax=8,
     1     ntens=4, ninpt=9, nsvint=13)
c     
c      ndim  ... number of spatial dimensions
c      ndof  ... number of degrees of freedom per node
c      ndi   ... number of direct stress components
c      nshr  ... number of shear stress component
c      ntens ... total number of stress tensor components (=ndi+nshr)
c      ninpt ... number of integration points
c      nsvint... number of state variables per integration pt 
c     
      dimension stiff(ndof*nnodemax,ndof*nnodemax),force(ndof*nnodemax),
     1 dN(nnodemax),dshape(ndim,nnodemax),xjaci(ndim,ndim),
     2 bmat(nnodemax*ndim),statevLocal(nsvint),stress(ntens), 
     3 ddsdde(ntens,ntens),stran(ntens),dstran(ntens),wght(ninpt)  
c
c
      do i = 1, 18
        ADIS(i) = U(i)
        ADDIS(i) = DU(i)
      end do
c      
      do i = 1, 16
        ASTR(i) = U(i+18)
        ADSTR(i) = DU(i+18)
      end do
c      
      do i = 1, 4
        ALAG(i) = U(i+34)
      end do
c      
      do i = 1, 9
        AX(i) = COORDS(1,i)
        AY(i) = COORDS(2,i)
        AREA = PROPS(1)
        EMOD = PROPS(2)
        ANU  = PROPS(3)
        ADENS= PROPS(4)
        AC_THREE = (3/5)**1/2
        AW_THREE_A = 25/81
        AW_THREE_B = 40/81
        AW_THREE_C = 64/81
        AXI_THREE(1)=-AC_THREE
        AXI_THREE(2)= AC_THREE
        AXI_THREE(3)= AC_THREE
        AXI_THREE(4)=-AC_THREE
        AXI_THREE(5)= 0
        AXI_THREE(6)= AC_THREE
        AXI_THREE(7)= 0
        AXI_THREE(8)=-AC_THREE
        AXI_THREE(9)= 0
        AOMEGA_THREE(1)=-AC_THREE
        AOMEGA_THREE(2)=-AC_THREE
        AOMEGA_THREE(3)= AC_THREE
        AOMEGA_THREE(4)= AC_THREE
        AOMEGA_THREE(5)=-AC_THREE
        AOMEGA_THREE(6)= 0
        AOMEGA_THREE(7)= AC_THREE
        AOMEGA_THREE(8)= 0
        AOMEGA_THREE(9)= 0
      end do
c
      do i = 1, 4
        AW_THREE(i) = AW_THREE_A
        AW_THREE(i+4) = AW_THREE_B
        AW_THREE(9) = AW_THREE_C
      end do
c 
      do kintk = 1, ninpt
c
!       Evaluate shape functions and derivatives
         call shapefcn_U(kintk,ninpt,nnode,ndim,dN_U,dNd_xi)
         
         

c***********************************************************************
c***********************************************************************
      subroutine shapefcn_U(kintk,ninpt,nnode,ndim,dN_U,dNd_xi)
c
      include 'aba_param.inc'
c
      parameter (gausspoint=0.774596669241483d0)
      dimension dN_U(nnode),dNd_xi(ndim,nnode),coord28(2,9)
c     
      data  coord28 /-1.d0, -1.d0,
        2                0.d0, -1.d0,
        3                1.d0, -1.d0,
        4               -1.d0,  0.d0,
        5                0.d0,  0.d0,
        6                1.d0,  0.d0,
        7               -1.d0,  1.d0,
        8                0.d0,  1.d0,      
        9                1.d0,  1.d0/
C     
C  2D 9-nodes
c
c     determine (xi,omega)
        xi = coord28(1,kintk)*gausspoint
        omega = coord28(2,kintk)*gausspoint
c
c     shape functions
        dN_U(1) = 0.25d0*xi*(xi-1.d0)*omega*(omega-1.d0)
        dN_U(2) = 0.25d0*xi*(xi+1.d0)*omega*(omega-1.d0)
        dN_U(3) = 0.25d0*xi*(xi+1.d0)*omega*(omega+1.d0)
        dN_U(4) = 0.25d0*xi*(xi-1.d0)*omega*(omega+1.d0)
        dN_U(5) = 0.5d0*(1.d0-xi*xi)*omega*(omega-1.d0)
        dN_U(6) = 0.5d0*xi*(xi+1.d0)*(1.d0-omega*omega)
        dN_U(7) = 0.5d0*(1.d0-xi*xi)*omega*(omega+1.d0)
        dN_U(8) = 0.5d0*xi*(xi-1.d0)*(1.d0-omega*omega)
        dN_U(9) = (1.d0-xi*xi)*(1.d0-omega*omega)      
c
c     derivative d(Ni)/d(xi)
        dNd_xi(1,1) = 0.5d0*xi*omega-0.25d0*omega*omega+0.25d0*omega
      &  +0.5d0*xi*omega*omega
        dNd_xi(1,2) = -0.5d0*xi*omega+0.5d0*xi*omega*omega
      &  +0.25d0*omega*omega-0.25d0*omega
        dNd_xi(1,3) = 0.5d0*xi*omega*omega+0.25d0*omega*omega
      &  +0.25d0*omega+0.5d0*xi*omega
        dNd_xi(1,4) = -0.25d0*omega*omega-0.25d0*omega
      &  +0.5d0*xi*omega*omega+0.5d0*xi*omega
        dNd_xi(1,5) = -xi*omega*omega+xi*omega
        dNd_xi(1,6) = -xi*omega*omega+0.5d0-0.5d0*omega*omega+xi
        dNd_xi(1,7) = -xi*omega*omega-xi*omega
        dNd_xi(1,8) = -0.5d0+0.5d0*omega*omega+xi-xi*omega*omega
        dNd_xi(1,9) = -2.d0*xi+2.d0*xi*omega*omega
c
c     derivative d(Ni)/d(omega)
        dNd_xi(2,1) = -0.25d0*omega*omega-0.5d0*xi*omega+0.25d0*xi
      &  +0.5d0*xi*xi*omega
        dNd_xi(2,2) = -0.25d0*omega*omega-0.5d0*xi*xi*omega
      &  +0.5d0*xi*omega+0.25d0*xi
        dNd_xi(2,3) = 0.5d0*xi*xi*omega+0.5d0*xi*omega+0.25d0*xi
      &  +0.25d0*xi*xi
        dNd_xi(2,4) = 0.5d0*xi*xi*omega-0.5d0*xi*omega-0.25d0*xi
      &  +0.25d0*xi*xi
        dNd_xi(2,5) = 0.5d0*xi*xi-0.5d0-xi*xi*omega+omega 
        dNd_xi(2,6) = -xi*xi*omega-xi*omega
        dNd_xi(2,7) = -0.5d0*xi*xi+omega-xi*xi*omega+0.5d0
        dNd_xi(2,8) = xi*omega-xi*xi*omega
        dNd_xi(2,9) = -2.d0*omega+2.d0*xi*xi*omega     
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine shapefcn_PSI(kintk,ninpt,nnode,ndim,dN_PSI,dNd_PSI)
c
      include 'aba_param.inc'
c
      parameter (gausspoint=0.774596669241483d0)
      dimension dN_PSI(nnode),dNd_PSI(ndim,nnode),coord28(2,9)
c
      data  coord28 /-1.d0, -1.d0,
            2                0.d0, -1.d0,
            3                1.d0, -1.d0,
            4               -1.d0,  0.d0,
            5                0.d0,  0.d0,
            6                1.d0,  0.d0,
            7               -1.d0,  1.d0,
            8                0.d0,  1.d0,
            9                1.d0,  1.d0/
C
C  2D 9-nodes
c
c     determine (xi,omega)
            xi = coord28(1,kintk)*gausspoint
            omega = coord28(2,kintk)*gausspoint
c     shape functions
            dN_PSI(1) = 0.25d0*(1.d0-xi)*(1.d0-omega)
            dN_PSI(2) = 0.25d0*(1.d0+xi)*(1.d0-omega)
            dN_PSI(3) = 0.25d0*(1.d0+xi)*(1.d0+omega)
            dN_PSI(4) = 0.25d0*(1.d0-xi)*(1.d0+omega)
c     derivative d(Ni)/d(xi)
            dNd_PSI(1,1) = 0.25d0*(-1.d0+omega)
            dNd_PSI(1,2) = 0.25d0*(1.d0-omega)
            dNd_PSI(1,3) = 0.25d0*(1.d0+omega)
            dNd_PSI(1,4) = 0.25d0*(-1.d0-omega)
c     derivative d(Ni)/d(omega)
            dNd_PSI(2,1) = 0.25d0*(-1.d0+xi)
            dNd_PSI(2,2) = 0.25d0*(-1.d0-xi)
            dNd_PSI(2,3) = 0.25d0*(1.d0+xi)
            dNd_PSI(2,4) = 0.25d0*(1.d0-xi)
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine jacobian_xi(jelem,ndim,nnode,coords,dNd_xi,djac,
    1 xjaci,pnewdt)
c
c     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix
c
      include 'aba_param.inc'
      parameter(maxDof=2)
c
      dimension xjac(maxDof,maxDof),xjaci(ndim,*),coords(3,*),
    1 dNd_xi(ndim,*)
c
      do i = 1, ndim
        do j = 1, ndim
          xjac(i,j)  = 0.d0
          xjaci(i,j) = 0.d0
        end do
      end do
c
      do inod= 1, nnode
         do idim = 1, ndim
           do jdim = 1, ndim
             xjac(jdim,idim) = xjac(jdim,idim) +
    1        dNd_xi(jdim,inod)*coords(idim,inod)
           end do
         end do
      end do
c
      djac = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
       if (djac .gt. 0.d0) then
         ! jacobian is positive - o.k.
         xjaci(1,1) =  xjac(2,2)/djac
         xjaci(2,2) =  xjac(1,1)/djac
         xjaci(1,2) = -xjac(1,2)/djac
         xjaci(2,1) = -xjac(2,1)/djac
       else
         ! negative or zero jacobian
         write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
         pnewdt = 0.25d0
       endif
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine jacobian_psi(jelem,ndim,nnode,coords,dNd_PSI,djac,
    1 xjaci_psi,pnewdt)
c
c     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix
c
      include 'aba_param.inc'
      parameter(maxDof=2)
c
      dimension xjac_psi(maxDof,maxDof),xjaci_psi(ndim,*),coords(3,*),
    1 dNd_PSI(ndim,*)
c
      do i = 1, ndim
        do j = 1, ndim
          xjac_psi(i,j)  = 0.d0
          xjaci_psi(i,j) = 0.d0
        end do
      end do
c
      do inod= 1, 4
         do idim = 1, ndim
           do jdim = 1, ndim
             xjac_psi(jdim,idim) = xjac_psi(jdim,idim) +
    1        dNd_PSI(jdim,inod)*coords(idim,inod)
           end do
         end do
      end do
c
      djac_psi = xjac_psi(1,1)*xjac_psi(2,2) - 
    1            xjac_psi(1,2)*xjac_psi(2,1)
       if (djac .gt. 0.d0) then
         ! jacobian is positive - o.k.
         xjaci_psi(1,1) =  xjac_psi(2,2)/djac
         xjaci_psi(2,2) =  xjac_psi(1,1)/djac
         xjaci_psi(1,2) = -xjac_psi(1,2)/djac
         xjaci_psi(2,1) = -xjac_psi(2,1)/djac
       else
         ! negative or zero jacobian
         write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
         pnewdt = 0.25d0
       endif
      return
      end
c***********************************************************************
c***********************************************************************
