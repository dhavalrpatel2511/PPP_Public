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
c*********************************************************************** 
      do kintk = 1, ninpt
c
!       Evaluate shape functions and derivatives
         call shapefcn_U(kintk,ninpt,nnode,ndim,dN_U,dNd_xi)
         call shapefcn_PSI(kintk,ninpt,nnode,ndim,dN_PSI,dNd_PSI)
c
!       Evaluate Jacobian of displacement and psi shape function matrix
         call jacobian_xi(jelem,ndim,nnode,coords,dNd_xi,djac,xjaci,
    1         pnewdt)
         call jacobian_omega(jelem,ndim,nnode,coords,dNd_PSI,djac,
    1         xjaci_psi,pnewdt)
c
         if (pnewdt .lt. pnewdtLocal) pnewdtLocal = pnewdt
c
!       Evaluate B_matrix for displcement and psi 
         call bmatrix_U(xjaci,dNd_xi,nnode,ndim,bmat_u)
         call bmatrix_PSI(xjaci_psi,dNd_PSI,nnode,ndim,bmat_psi)
c
            
c
c***********************************************************************
!       assembly of N_matrix
c
         do i = 1, 9
            Nmatrix(1,(2*i-1)) = dN_U(i)
            Nmatrix(2,(2*i)) = dN_U(i)
c             Nmatrix(1,i) = Nmatrix(1,i)+dN_U(i)
c             Nmatrix(2,i) = Nmatrix(2,i)+dN_U(i)
c                
         end do
c
c         
!       assembly of M_matrix    
         do i = 1, 2*nnode
            Mmatrix(1,i) = 0.d0
            Mmatrix(2,i) = 0.d0
            Mmatrix(3,i) = 0.d0
            Mmatrix(4,i) = 0.d0
         end do   
c
         do i = 1, nnode
            Mmatrix(1,2*i-1) = bmat_u(2*i-1)
            Mmatrix(2,2*i)   = bmat_u(2*i)
            Mmatrix(3,2*i-1) = bmat_u(2*i-1)
            Mmatrix(4,2*i)   = bmat_u(2*i)
         end do  
c    
!       assembly of B_matrix    
         do i = 1, 2*nnode
            Bmatrix(1,i) = 0.d0
            Bmatrix(2,i) = 0.d0
            Bmatrix(3,i) = 0.d0
         end do   
c
         do i = 1, nnode
            Bmatrix(1,2*i-1) = bmat_u(2*i-1)
            Bmatrix(2,2*i)   = bmat_u(2*i)
            Bmatrix(3,2*i-1) = bmat_u(2*i)
            Bmatrix(3,2*i)   = bmat_u(2*i-1)
         end do
c
!       assembly of N_PSI_matrix 
         do i = 1, 16
            Npsimatrix(1,i) = 0.d0
            Npsimatrix(2,i) = 0.d0
            Npsimatrix(3,i) = 0.d0
            Npsimatrix(4,i) = 0.d0
         end do
c
         do i = 1, 4
            Npsimatrix(1,4*i-3) = dN_PSI(i)
            Npsimatrix(2,4*i-2) = dN_PSI(i)
            Npsimatrix(3,4*i-1) = dN_PSI(i)
            Npsimatrix(4,4*i)   = dN_PSI(i)
         end do
c
!       assembly of B_PSI_matrix 
         do i = 1, 16
            Bpsimatrix(1,i) = 0.d0
            Bpsimatrix(2,i) = 0.d0
            Bpsimatrix(3,i) = 0.d0
            Bpsimatrix(4,i) = 0.d0
            Bpsimatrix(5,i) = 0.d0
            Bpsimatrix(6,i) = 0.d0
         end do
c
         do i = 1, 4
            Bpsimatrix(1,4*i-3) = bmat_PSI(2*i-1)
            Bpsimatrix(2,4*i-2) = bmat_PSI(2*i)
            Bpsimatrix(3,4*i-1) = bmat_PSI(2*i-1)
            Bpsimatrix(4,4*i)   = bmat_PSI(2*i)
            Bpsimatrix(5,4*i-3) = bmat_PSI(2*i)
            Bpsimatrix(5,4*i-2) = bmat_PSI(2*i-1)
            Bpsimatrix(6,4*i-1) = bmat_PSI(2*i)
            Bpsimatrix(6,4*i)   = bmat_PSI(2*i-1)
         end do
c
!       assembly of N_rho
         do i = 1, 4
            Nrhomatrix(1,i) = 0.d0
            Nrhomatrix(2,i) = 0.d0 
            Nrhomatrix(3,i) = 0.d0
            Nrhomatrix(4,i) = 0.d0
         end do
         do i = 1, 4
            Nrhomatrix(i,i) = 1.d0
         end do
c
********************************************************
c         
      !  update of state variables.
         !GradVariable = 0.d0
         T_Mmatrix = Transpose(Mmatrix)		
         do i = 1, 18
          do j = 1, 4
            GradVariable(j,kintk) = GradVariable(j,kintk) 
     1      +T_Mmatrix(i,j)*displacement(i)
          end do
         end do  
c
c
	     T_Bmatrix = Transpose(Bmatrix)
         !strain = 0.d0
         !Deltastrain = 0.d0
         do i = 1, 18
          do j = 1, 3
            strain(j,kintk) = strain(j,kintk) 
     1      + T_Bmatrix(i,j)*displacement(i)
            Deltastrain(j,kintk) = Deltastrain(j,kintk) 
     1      +T_Bmatrix(i,j)*delta_displacement(i)
          end do
         end do
c
c
	     T_Npsimatrix = Transpose(Npsimatrix)
         !relaxedstrain = 0.d0
         do i = 1, 16
          do j = 1, 4
            relaxedstrain(j,kintk) = relaxedstrain(j,kintk) 
     1      +T_Npsimatrix(i,j)*relaxed_strain(i)
          end do  
         end do
c
c
         !relaxedstraingradient = 0.d0
         !delta_relaxedstraingradient = 0.d0
	     T_Bpsimatrix = Transpose(Bpsimatrix)
         do i = 1, 16
          do j = 1, 6
            relaxedstraingradient(j,kintk)  
     1      =relaxedstraingradient(j,kintk) 
     1      +T_Bpsimatrix(i,j)*relaxed_strain(i)
c 
c       
            delta_relaxedstraingradient(j,kintk)
     1      =delta_relaxedstraingradient(j,kintk) 
     1      +T_Bpsimatrix(i,j)*delta_relaxed_strain(i)
          end do
         end do
c
c
         !Langrangemulti = 0.d0
         do i = 1, 4
          do j = 1, 4
            Langrangemulti(j,kintk) = Langrangemulti(j,kintk) 
     1      +Nrhomatrix(i,j)*lagrange_multi(i)
          end do
         end do
c
c                  
c***********************************************************************
!       call the KUMAT to find the state variables Sigma and Tau
!       and also get U_psilon and Lambda from KUMAT as output.
c
         call KUMAT(Sigma,Tau,U_psilon,Lambda,kintk,lemda,mue,
     1   Deltastrain,delta_relaxedstraingradient)
c     
!       Print the State variable Sigma in .dat file.     
*         write(6,*) "this is Sigma"
*         write(6,*) Sigma(:,1)
*         do i = 1, size(Sigma,kintk)
*            write(6,'(20G12.4)')  Sigma(i,:)
*         end do
!       Print the State variable Tau in .dat file.
*         write(6,*) "this is Tau"
*         write(6,*) Tau(:,kintk)
*         do i = 1, size(Tau,1)
*            write(6,'(20G12.4)')  Tau(i,:)
*         end do
c
c
c
!       store state variables Sigma and Tau in Svars vector. :
         do i = 1, 3
            svars(i+(kintk-1)*9) = Sigma(i,kintk)
         end do
         do i = 4, 9
            svars(i+(kintk-1)*9) = Tau(i-3,kintk)
         end do
!       Print the State variable vector in .dat file.         
*         write(6,*) "this is state vars"
*         do i = 1, 9
*            write(6,'(20G12.4)')  svars(i)
*         end do
c
!-------Find the Element Stiffness Matrices:----------------------------
c
!       Find the UU_stiffness matrix - dimension(18*18)
         do i = 1, 18
          do j = 1, 18
           do k = 1, 3
            do l = 1, 3
               stiffness_UU(i,j) = stiffness_UU(i,j) 
     1         +T_Bmatrix(i,k)*U_psilon(k,l)*T_Bmatrix(j,l)*djac
     1         *Gauss_weight(kintk)
            end do
           end do
          end do
         end do
!       Print the stiffness_UU matrix in .dat file.           
*         write(6,*) "this is sti_UU"
*         do i = 1, size(stiffness_UU,1)
*            write(6,'(20G12.4)')  stiffness_UU(i,:)
*         end do
c
c
!       Find the PsiPsi_stiffness matrix - dimension(16*16)
         do i = 1, 16
          do j = 1, 16
           do k = 1, 6
             do l = 1, 6
               stiffness_PsiPsi(i,j) = stiffness_PsiPsi(i,j) 
     1         +T_Bpsimatrix(i,k)*Lambda(k,l)*T_Bpsimatrix(j,l)
     1          *djac_psi*Gauss_weight(kintk)
             end do
           end do
          end do
         end do
!       Print the stiffness_PsiPsi matrix in .dat file.          
*         write(6,*) "this is sti_psipsi"
*         do i = 1, size(stiffness_PsiPsi,1)
*            write(6,'(20G12.4)')  stiffness_PsiPsi(i,:)
*         end do
c
c
!       Find the Urho_stiffness matrix - dimension(18*4)
         do i = 1, 18
          do j = 1, 4
           do k = 1, 4
               stiffness_Urho(i,j) = stiffness_Urho(i,j) 
     1        +T_Mmatrix(i,k)*Nrhomatrix(j,k)*djac*Gauss_weight(kintk)
           end do
          end do
         end do
!       Print the stiffness_Urho matrix in .dat file.        
*         write(6,*) "this is stiffness_Urho"
*         do i = 1, size(stiffness_Urho,1)
*            write(6,'(20G12.4)')  stiffness_Urho(i,:)
*         end do
c
c
!       Find the psirho_stiffness matrix - dimension(16*4)
         do i = 1, 16
          do j = 1, 4
           do k = 1, 4
               stiffness_psirho(i,j) = stiffness_psirho(i,j) 
     1         +T_Npsimatrix(i,k)*Nrhomatrix(j,k)*djac_psi
     1         *Gauss_weight(kintk)
           end do
          end do
         end do
!       Print the stiffness_psirho matrix in .dat file.         
*         write(6,*) "this is stiffness_psirho"
*         do i = 1, size(stiffness_psirho,1)
*            write(6,'(20G12.4)')  stiffness_psirho(i,:)
*         end do
c
c
!       Find the Mass_U matrix - dimension(18*18)
         T_Nmatrix = Transpose(Nmatrix)
         do i = 1, 18
          do j = 1, 18
           do k = 1, 2
               Mass_U(i,j) = Mass_U(i,j) + Density*T_Nmatrix(i,k)
     1         *T_Nmatrix(j,k)*djac*Gauss_weight(kintk)
           end do
          end do
         end do
!       Print the Mass_U matrix in .dat file.          
*         write(6,*) "this is Mass_U"
*         do i = 1, size(Mass_U,1)
*            write(6,'(20G12.4)')  Mass_U(i,:)
*         end do
c***********************************************************************
c
!       Find the internal force vectors : 
c
!       Find the F_vector which includes the internal forces due to
!       the displcement of all nodes.
         do i = 1, 18 
            do j = 1, 3
               F_vector(i) = F_vector(i) + T_Bmatrix(i,j)*Sigma(j,kintk)
     1         *djac*Gauss_weight(kintk)
*                write(*,*) Bmatrix(i,j)*Sigma(j,kintk)
*     1               *djac*Gauss_weight(kintk)
            end do
            do j = 1, 4
               F_vector(i) = F_vector(i) - T_Mmatrix(i,j) 
     1         *Langrangemulti(j,kintk)*djac*Gauss_weight(kintk)  
            end do
         end do
!       Print the F_vactor in .dat file.         
         write(6,*) "F_vector"
         do i = 1, 18
            write(6,'(20G12.4)')  F_vector(i)
         end do
c
c
!       Find the R_vector which includes the internal forces due to
!       the relaxed strain of nodes(1,2,3,4).
         do i = 1, 16
            do j = 1, 6
              R_vector(i) = R_vector(i) + T_Bpsimatrix(i,j)*Tau(j,kintk)
     1        *djac_psi*Gauss_weight(kintk)
            end do
            do j = 1, 4
               R_vector(i) = R_vector(i) + T_Npsimatrix(i,j)
     1          *Langrangemulti(j,kintk)*djac_psi*Gauss_weight(kintk)
            end do
         end do
!       Print the R_vactor in .dat file.         
         write(6,*) "R_vector"
         do i = 1, 16
            write(6,'(20G12.4)')  R_vector(i)
         end do
c
c
!       Find the S_vector which includes the internal forces due to
!       the langangen multiplier of centre(last) node.
         do i = 1, 4
            do j = 1, 4
               S_vector(i) = S_vector(i) + Nrhomatrix(i,j)
     1        *(relaxedstrain(j,kintk)*djac_psi-GradVariable(j,kintk)
     1        *djac)*Gauss_weight(kintk)
            end do
         end do
!       Print the S_vactor in .dat file.         
         write(6,*) "S_vector"
         do i = 1, 4
            write(6,'(20G12.4)')  S_vector(i)
         end do
c
      end do 
!       loop over integration points is end.
c***********************************************************************
!       Assembaly of AMATRX and RHS acorrding to the arragement of 
!       degrees of freedom of the all 9 nodes of the element.
c      
       tra_stiffness_Urho = TRANSPOSE(stiffness_Urho)   
       tra_stiffness_psirho = TRANSPOSE(stiffness_psirho) 
*      write(*,*) "tra_stiffness_Urho"
*      do i = 1, size(tra_stiffness_Urho,1)
*          write(*,'(20G12.4)')  tra_stiffness_Urho(i,:)
*      end do
!       Print the amat matrix to check the intialization of it.
*      write(*,*) "amat"
*      do i = 1, size(amat,1)
*          write(*,'(20G12.4)')  amat(i,:)
*      end do
c
C           ___                              ____
!      K = |                                     |
!          | K_UU          0          -K_urho    |
!          |  0         K_psipsi       K_psirho  |
!          |-T_K_urho   T_K_psirho        0      |
!          |___                              ____|
c
c
!     amat(1:18,1:18) =  stiffness_UU(1:18,1:18)     
      do i = 1, 18
         do j = 1, 18
            amatrx(i,j) = amatrx(i,j) + stiffness_UU(i,j)
         end do
      end do
c
!     amat(19:34,19:34) =  stiffness_PsiPsi(1:16,1:16)
      do i = 19, 34
         do j = 19, 34
            amatrx(i,j) = amatrx(i,j) + stiffness_PsiPsi(i-18,j-18)   
         end do
      end do
c
!     amat(35:38,19:34) =  tra_stiffness_psirho(1:4,1:16)
      do i = 35, 38
         do j = 19, 34
            amatrx(i,j) = amatrx(i,j) + tra_stiffness_psirho(i-34,j-18)
         end do
      end do
c
!     amat(35:38,1:18) =  tra_stiffness_Urho(1:4,1:18)
      do i = 35, 38
         do j = 1, 18
            amatrx(i,j) = amatrx(i,j) - tra_stiffness_Urho(i-34,j)
         end do
      end do
c
!     amat(1:18,35:38) =  stiffness_Urho(1:18,1:4)
      do i = 1, 18
         do j = 35, 38
            amatrx(i,j) = amatrx(i,j) - stiffness_Urho(i,j-34)
         end do
      end do
c
!     amat(19:34,35:38) =  stiffness_psirho(1:16,1:4)
      do i = 19, 34
         do j = 35, 38
            amatrx(i,j) = amatrx(i,j) + stiffness_psirho(i-18,j-34)
         end do
      end do
c
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
      subroutine bmatrix_U(xjaci,dNd_xi,nnode,ndim,bmat_u)
c
c     Notation: bmat(i) ....dN1/dx, dN1/dy, dN2/dx, dN2/dy...
          include 'aba_param.inc'
          dimension bmat_u(*),dNd_xi(ndim,*),xjaci(ndim,*)
c
          do i = 1, nnode*ndim
              bmat_u(i) = 0.d0
          end do
c
          do inod = 1, nnode
           do ider = 1, ndim
            do idim = 1, ndim
             irow = idim + (inod - 1)*ndim
             bmat_u(irow)=bmat_u(irow)+xjaci(idim,ider)*dNd_xi(ider,inod)
            end do
           end do
          end do
c
          return
          end
c***********************************************************************
c***********************************************************************
      subroutine bmatrix_PSI(xjaci_psi,dNd_PSI,nnode,ndim,bmat_psi)
c
c     Notation: bmat(i) ....dN1/dx, dN1/dy, dN2/dx, dN2/dy...
            include 'aba_param.inc'
            dimension bmat_psi(*),dshape(ndim,*),xjaci(ndim,*)
c
            do i = 1, nnode*ndim
                  bmat_psi(i) = 0.d0
            end do
c
            do inod = 1, 4
              do ider = 1, ndim
                do idim = 1, ndim
                  irow = idim + (inod - 1)*ndim
                  bmat_psi(irow)=bmat_psi(irow)+
    1                            xjaci_psi(idim,ider)*dNd_PSI(ider,inod)
                end do
              end do
            end do
C            
      return
      end
c***********************************************************************
c***********************************************************************
      subroutine KUMAT(Sigma,Tau,U_psilon,Lambda,kintk,lemda,mue,
     1 Deltastrain,delta_relaxedstraingradient)
c
      double precision :: Sigma,Tau,C_matrix,D_matrix
c
      dimension C_matrix(3,3), D_matrix(6,6), sigma(3,9), Tau(6,9),
     1 Lambda(6,6), U_psilon(3,3), Deltastrain(3,9),
     2 delta_relaxedstraingradient(6,9)
*      write(*,*) "this is lemda and mue"
*      write(*,*) lemda
*      write(*,*) mue
c
      C_matrix(1,1) = lemda + 2.d0*mue
      C_matrix(1,2) = lemda
      C_matrix(1,3) = 0.d0
      C_matrix(2,1) = lemda
      C_matrix(2,2) = lemda + 2.d0*mue
      C_matrix(2,3) = 0.d0
      C_matrix(3,1) = 0.d0
      C_matrix(3,2) = 0.d0
      C_matrix(3,3) = 0.5d0*mue
c
*      write(*,*) "this is C_matrix"
*      do i = 1, size(C_matrix,1)
*        write(*,'(20G12.4)')  C_matrix(i,:)
*      end do
      do i = 1, 6
         do j = 1, 6
            D_matrix(j,i) = 0.d0
         end do
      end do
      D_matrix(2,2) = 1.d0
      D_matrix(3,3) = 1.d0
      D_matrix(2,6) = -0.5d0
      D_matrix(3,5) = -0.5d0
      D_matrix(5,3) = -0.5d0
      D_matrix(6,2) = -0.5d0
      D_matrix(5,5) = 0.25d0
      D_matrix(6,6) = 0.25d0
c      
*      write(*,*) "this is D_matrix"
*      do i = 1, size(D_matrix,1)
*        write(*,'(20G12.4)')  D_matrix(i,:)
*      end do
c
      do i = 1, 3
         do j = 1, 3
            Sigma(i,kintk) = Sigma(i,kintk) 
     1      +C_matrix(i,j)*Deltastrain(j,kintk)
c
            U_psilon(i,j) = C_matrix(i,j)
         end do
      end do
c
      do i = 1, 6
         do j = 1, 6
            Tau(i,kintk) = Tau(i,kintk) 
     1       +D_matrix(i,j)*delta_relaxedstraingradient(j,kintk)
            Lambda(i,j) = D_matrix(i,j)
         end do
      end do
c      
      return 
      end subroutine KUMAT
c***********************************************************************
