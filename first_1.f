!  || Shree Geneshay Namah ||
! User Element subroutine for Strain Gradient Elasticity model.
! Dhaval Rasheshkumar Patel
! dhavalrpatel2511@gmail.com
c**********************************************************************
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
