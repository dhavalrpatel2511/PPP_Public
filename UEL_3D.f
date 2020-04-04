      subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,
     1 nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,
     2 jelem,params,ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,  
     3 ddlmag,mdload,pnewdt,jprops,njpro,period)
c
c
!    Standard format of UEL for abaqus.
!    This statement includes the inbuilt variables/parameters for Abaqus
      include 'aba_param.inc' 
c
      dimension rhs(162),amatrx(162,162),props(4),svars(324),
     1 coords(3,27),u(162),du(1,153),time(*),
     2 a(ndofel),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*),
     4 energy(*),v(ndofel)
c
!    Assign the values for some basics parameters.
      parameter (ndim=3, nnodemax=27, ninpt=27, nsvint=24)
c
c
c    ndim      ... number of spatial dimensions
c    nnodemax  ... number of nodes per element
c    ninpt     ... number of integration points
c    nsvint    ... number of state variables per integration pt     
c
c
!    defining and dimensioning the vectors, matrix.
      double precision, dimension(9)     :: S_vector
      double precision, dimension(27)    :: dN_U  
      double precision, dimension(24)    :: Bmat_psi
      double precision, dimension(9)     :: lagrange_multi         
      double precision, dimension(8)     :: dN_PSI
      double precision, dimension(24)    :: statevLocal
      double precision, dimension(27)    :: Gauss_Weight
      double precision, dimension(81)    :: delta_displacement
      double precision, dimension(72)    :: relaxed_strain
      double precision, dimension(81)    :: F_vector   
      double precision, dimension(72)    :: R_vector      
      double precision, dimension(81)    :: Bmat_u      
      double precision, dimension(81)    :: displacement
      double precision, dimension(72)    :: delta_relaxed_strain   
      double precision, dimension(3,27)   :: dN_U_X
      double precision, dimension(3,8)   :: dN_PSI_X
      double precision, dimension(3,3)   :: xjaci
      double precision, dimension(3,3)   :: xjaci_psi
      double precision, dimension(9,27)   :: GradVariable
      double precision, dimension(6,27)   :: strain
      double precision, dimension(6,27)   :: Deltastrain
      double precision, dimension(9,27)   :: relaxedstrain
      double precision, dimension(18,27)   :: relaxedstraingradient
      double precision, dimension(18,27)  :: delta_relaxedstraingradient
      double precision, dimension(9,27)   :: Langrangemulti
      double precision, dimension(6,27)   :: Sigma
      double precision, dimension(18,27)   :: Tau
      double precision, dimension(9,9)   :: Nrhomatrix
      double precision, dimension(6,6)   :: U_psilon
      double precision, dimension(18,18)   :: Lambda
      double precision, dimension(81,9)  :: stiffness_Urho
      double precision, dimension(72,9)  :: stiffness_psirho
      double precision, dimension(9,81)  :: tra_stiffness_Urho
      double precision, dimension(9,72)  :: tra_stiffness_psirho
      double precision, dimension(81,3)  :: T_Nmatrix
      double precision, dimension(81,9)  :: T_Mmatrix
      double precision, dimension(81,6)  :: T_Bmatrix
      double precision, dimension(72,9)  :: T_Npsimatrix
      double precision, dimension(72,18)  :: T_Bpsimatrix
      double precision, dimension(3,81)  :: Nmatrix
      double precision, dimension(9,81)  :: Mmatrix
      double precision, dimension(6,81)  :: Bmatrix
      double precision, dimension(18,72)  :: Bpsimatrix
      double precision, dimension(9,72)  :: Npsimatrix
      double precision, dimension(81,81) :: stiffness_UU 
      double precision, dimension(72,72) :: stiffness_PsiPsi      
      double precision, dimension(162,162) :: Mass_matrix   
           
c 
      double precision :: E_modulus,Nue,micro_length,Density,
     1  lemda,mue
c    
        S_vector = 0.d0
        dN_U   = 0.d0
        Bmat_psi = 0.d0   
        dN_PSI = 0.d0
        statevLocal = 0.d0
        F_vector = 0.d0
        R_vector  = 0.d0
        Bmat_u = 0.d0
        dN_U_X = 0.d0
        dN_PSI_X = 0.d0
        xjaci = 0.d0
        xjaci_psi = 0.d0
        GradVariable = 0.d0
        strain = 0.d0
        Deltastrain = 0.d0
        relaxedstrain = 0.d0
        relaxedstraingradient = 0.d0
        delta_relaxedstraingradient = 0.d0
        Langrangemulti = 0.d0
        Sigma = 0.d0
        Tau = 0.d0
        Nrhomatrix = 0.d0
        U_psilon = 0.d0
        Lambda = 0.d0
        stiffness_Urho = 0.d0
        stiffness_psirho = 0.d0
        tra_stiffness_Urho = 0.d0
        tra_stiffness_psirho = 0.d0
        T_Nmatrix = 0.d0
        T_Mmatrix = 0.d0
        T_Bmatrix = 0.d0
        T_Npsimatrix = 0.d0
        T_Bpsimatrix = 0.d0
        Nmatrix = 0.d0
        Mmatrix = 0.d0
        Bmatrix = 0.d0
        Bpsimatrix = 0.d0
        Npsimatrix = 0.d0
        stiffness_UU  = 0.d0
        stiffness_PsiPsi = 0.d0
        amat = 0.d0
        Mass_U = 0.d0 
        Mass_matrix = 0.d0
c
c   
!    Guass Weight of Guass integration points in respective order
      data Gauss_Weight /0.1714677641d0,0.1714677641d0,0.1714677641d0,
     1 0.1714677641d0,0.1714677641d0,0.1714677641d0,0.1714677641d0,
     1 0.1714677641d0,0.2743484225d0,0.2743484225d0,0.2743484225d0,
     1 0.2743484225d0,0.2743484225d0,0.2743484225d0,0.2743484225d0,
     1 0.2743484225d0,0.2743484225d0,0.2743484225d0,0.2743484225d0,
     1 0.2743484225d0,0.4389574760d0,0.4389574760d0,0.4389574760d0,
     1 0.4389574760d0,0.4389574760d0,0.4389574760d0,0.7023319616d0/
c
c
!    write the adlmag
      write(*,*) lflags(3)		
c
!    this is to print the U and DU in .dat file 
      write(6,*) "this is U"
      do i = 1, 162
        write(6,*) u(i)
      end do  
      write(6,*) "this is DU"
      do i = 1, 153
        write(6,*) du(1,i)
      end do 
!    this is to print the element number in .dat file      
      write(6,*) jelem				
!    this is to print the coordintes of the nodes 
!    of the respective element in .dat file      
      write(6,*) "this is a coords of nodes of current element"
      do i = 1, size(coords,1)
            write(6,'(20G12.4)') coords(i,:)
      end do
!    this is to print the Element Properties in .dat file 
*      write(*,*) "properties"
*      write(*,*) props(1)
*      write(*,*) props(2)
      E_modulus = props(1)
      Nue = props(2)
      micro_length = props(3)
      Density = props(4)
      lemda = (E_modulus*Nue)/((1.d0+Nue)*(1.d0-2.d0*Nue))
      mue = E_modulus/2.d0/(1.0d0+Nue)
!    this is to print the lame constants in .dat file
*      write(6,*) "this is lemda and mue"	
*      write(6,*) lemda
*      write(6,*) mue
c
!    loop for extracting the displacement and- 
!    delta_displacement value from U ans DU vector.
      do i = 1, 81
         displacement(i) = u(i)
         delta_displacement(i) = du(1,i)
      end do
c
!    loop for extracting the relaxed_strain and-
!    delta_relaxed_strain value from U and DU.
      do i = 1, 72
         relaxed_strain(i) = u(i+81) 
         delta_relaxed_strain(i) = du(1,i+81)
      end do
c 
!    loop for extracting lagrangian multiplier-
!    from the U vector.
      do i = 1, 9
         lagrange_multi(i) = u(i+153) 
      end do
c
c
************************************************************************
!    loop over Gauss Integration Points
!    There are total 9 integration points so it goes from 1 to 9.
c
      do kintk = 1, 27
c
         write(6,*) "The integeration no. is"
         write(6,*) kintk
!     Evaluate shape functions and derivatives
         call shapefcn_U(dN_U,dN_U_X,kintk,ninpt,nnode,ndim)
!     Now we have shape function and its derivates for displacement 
!     degree of freedom. 
*         write(6,*) "shape_function_u"
*         write(6,*) dNd_xi
         call shapefcn_PSI(dN_PSI,dN_PSI_X,kintk,ninpt,nnode,ndim)
*         write(6,*) "shape_function_PSI"
*         write(6,*) dNd_PSI
!     Now we have shape function and its derivates for relaxedstrain 
!     degree of freedom.   
c 
*         do i = 1, size(dN_U,1)
*            write(*,'(20G12.4)') dN_U(:)
*         end do
*
*         do i = 1, size(dNd_xi,1)
*            write(*,'(20G12.4)') dNd_xi(i,:)
*         end do
!       Evaluate Jacobian of displacement and psi shape function matrix
         call Jacobian_UU(djac,xjaci,jelem,ndim,nnode,coords,dN_U_X)
*         write(6,*) "Jacobian of displacement"
*         write(6,*) djac         
*         write(6,*) xjaci
                  
         call Jacobian_psipsi(dj_psi,xjaci_psi,jelem,ndim,nnode,
     1           coords,dN_PSI_X)
*         write(6,*) "Jacobian of PSI"
*         write(6,*) djac_psi
*         write(6,*) xjaci_psi     
c
!       Evaluate B_matrix for displcement and psi 
         call Bmatrix_U(Bmat_u,xjaci,dN_U_X,nnode,ndim)
c~          write(6,*) "bmat_u"         
c~          write(6,*) bmat_u
         call Bmatrix_PSIPSI(Bmat_psi,xjaci_psi,dN_PSI_X,nnode,ndim)
c~          write(6,*) "bmat_psi"         
c~          write(6,*) bmat_psi
c
c
c***********************************************************************
!       Assembly of Nmatrix from shape function for displacement degree
!       freedom.
         do i = 1, 3*nnode
            Nmatrix(1,i) = 0.d0
            Nmatrix(2,i) = 0.d0
            Nmatrix(3,i) = 0.d0
         end do   
c            
         do i = 1, 27
            Nmatrix(1,(3*i-2)) = dN_U(i)
            Nmatrix(2,(3*i-1)) = dN_U(i)
            Nmatrix(3,(3*i))   = dN_U(i)
         end do
!       print the Nmatrix in .dat file            
         write(6,*) "Nmatrix is here"
         do i = 1, size(Nmatrix,1)
            write(6,'(40G12.4)') Nmatrix(i,:)
         end do
c    
c     
!       Assembly of Mmatrix from derivatives of shape functions for 
!       displacement degree freedom.    
         do i = 1, 3*nnode
            Mmatrix(1,i) = 0.d0
            Mmatrix(2,i) = 0.d0
            Mmatrix(3,i) = 0.d0
            Mmatrix(4,i) = 0.d0
            Mmatrix(5,i) = 0.d0
            Mmatrix(6,i) = 0.d0
            Mmatrix(7,i) = 0.d0
            Mmatrix(8,i) = 0.d0
            Mmatrix(9,i) = 0.d0
         end do   
c         
         do i = 1, nnode
            Mmatrix(1,3*i-2) = Bmat_u(3*i-2)
            Mmatrix(2,3*i-1) = Bmat_u(3*i-1)
            Mmatrix(3,3*i)   = Bmat_u(3*i)
            Mmatrix(4,3*i-2) = Bmat_u(3*i-1)
            Mmatrix(5,3*i-1) = Bmat_u(3*i-2)
            Mmatrix(6,3*i-2) = Bmat_u(3*i)
            Mmatrix(7,3*i)   = Bmat_u(3*i-2)
            Mmatrix(8,3*i-1) = Bmat_u(3*i)
            Mmatrix(9,3*i)   = Bmat_u(3*i-1)
         end do  
!       print the Mmatrix in .dat file         
         write(6,*) 'this is mmatrix'
         do i = 1, size(Mmatrix,1)
            write(6,'(85G12.4)') Mmatrix(i,:)
         end do         
c
c
!       Assembly of Bmatrix from derivatives of shape functions for 
!       displacement degree freedom.    
         do i = 1, 3*nnode
            Bmatrix(1,i) = 0.d0
            Bmatrix(2,i) = 0.d0
            Bmatrix(3,i) = 0.d0
            Bmatrix(4,i) = 0.d0
            Bmatrix(5,i) = 0.d0
            Bmatrix(6,i) = 0.d0
         end do   
c         
         do i = 1, nnode
            Bmatrix(1,3*i-2) = Bmat_u(3*i-2)
            Bmatrix(2,3*i-1) = Bmat_u(3*i-1)
            Bmatrix(3,3*i)   = Bmat_u(3*i)
            Bmatrix(4,3*i-2) = Bmat_u(3*i-1)
            Bmatrix(4,3*i-1) = Bmat_u(3*i-2)
            Bmatrix(5,3*i-2) = Bmat_u(3*i)
            Bmatrix(5,3*i)   = Bmat_u(3*i-2)
            Bmatrix(6,3*i-1) = Bmat_u(3*i)
            Bmatrix(6,3*i)   = Bmat_u(3*i-1)
         end do
!       print the Bmatrix in .dat file         
         write(6,*) 'this is Bmatrix'
         do i = 1, size(Bmatrix,1)
            write(6,'(85G12.4)') Bmatrix(i,:)
         end do
c
c
!       Assembly of Npsimatrix from shape functions for 
!       relaxed strain degree freedom. 
         do i = 1, 72
            Npsimatrix(1,i) = 0.d0
            Npsimatrix(2,i) = 0.d0
            Npsimatrix(3,i) = 0.d0
            Npsimatrix(4,i) = 0.d0
            Npsimatrix(5,i) = 0.d0
            Npsimatrix(6,i) = 0.d0
            Npsimatrix(7,i) = 0.d0
            Npsimatrix(8,i) = 0.d0
            Npsimatrix(9,i) = 0.d0
         end do
c         
         do i = 1, 8
            Npsimatrix(1,9*i-8) = dN_PSI(i)
            Npsimatrix(2,9*i-7) = dN_PSI(i)
            Npsimatrix(3,9*i-6) = dN_PSI(i)
            Npsimatrix(4,9*i-5) = dN_PSI(i)
            Npsimatrix(5,9*i-4) = dN_PSI(i)
            Npsimatrix(6,9*i-3) = dN_PSI(i)
            Npsimatrix(7,9*i-2) = dN_PSI(i)
            Npsimatrix(8,9*i-1) = dN_PSI(i)
            Npsimatrix(9,9*i)   = dN_PSI(i)
         end do
!       print the Npsimatrix in .dat file
         write(6,*) 'this is Npsimatrix'
         do i = 1, size(Npsimatrix,1)
            write(6,'(20G12.4)') Npsimatrix(i,:)
         end do
c
c
!       Assembly of Bpsimatrix from derivatives of shape functions for 
!       relaxed strain degree freedom. 
         do i = 1, 72
            Bpsimatrix(1,i) = 0.d0
            Bpsimatrix(2,i) = 0.d0
            Bpsimatrix(3,i) = 0.d0
            Bpsimatrix(4,i) = 0.d0
            Bpsimatrix(5,i) = 0.d0
            Bpsimatrix(6,i) = 0.d0 
            Bpsimatrix(7,i) = 0.d0 
            Bpsimatrix(8,i) = 0.d0 
            Bpsimatrix(9,i) = 0.d0 
            Bpsimatrix(10,i) = 0.d0 
            Bpsimatrix(11,i) = 0.d0 
            Bpsimatrix(12,i) = 0.d0 
            Bpsimatrix(13,i) = 0.d0 
            Bpsimatrix(14,i) = 0.d0 
            Bpsimatrix(15,i) = 0.d0 
            Bpsimatrix(16,i) = 0.d0 
            Bpsimatrix(17,i) = 0.d0 
            Bpsimatrix(18,i) = 0.d0
         end do
c
         do i = 1, 8
            Bpsimatrix(1,9*i-8)  = Bmat_psi(3*i-2)
            Bpsimatrix(2,9*i-7)  = Bmat_psi(3*i-1)
            Bpsimatrix(3,9*i-6)  = Bmat_psi(3*i)
c            
            Bpsimatrix(4,9*i-4)  = Bmat_psi(3*i-2)
            Bpsimatrix(5,9*i-5)  = Bmat_psi(3*i-2)
            Bpsimatrix(5,9*i-8)  = Bmat_psi(3*i-1)
            Bpsimatrix(6,9*i-7)  = Bmat_psi(3*i-2)
            Bpsimatrix(6,9*i-4)  = Bmat_psi(3*i-1)
            Bpsimatrix(7,9*i-5)  = Bmat_psi(3*i-1)
c
            Bpsimatrix(8,9*i-2)  = Bmat_psi(3*i-2)
            Bpsimatrix(9,9*i-3)  = Bmat_psi(3*i-2)
            Bpsimatrix(9,9*i-8)  = Bmat_psi(3*i)
            Bpsimatrix(10,9*i-6) = Bmat_psi(3*i-2)
            Bpsimatrix(10,9*i-2) = Bmat_psi(3*i)
            Bpsimatrix(11,9*i-3) = Bmat_psi(3*i)
c            
            Bpsimatrix(12,9*i)   = Bmat_psi(3*i-1)
            Bpsimatrix(13,9*i-1) = Bmat_psi(3*i-1)
            Bpsimatrix(13,9*i-7) = Bmat_psi(3*i)
            Bpsimatrix(14,9*i-6) = Bmat_psi(3*i-1)
            Bpsimatrix(14,9*i)   = Bmat_psi(3*i)
            Bpsimatrix(15,9*i-1) = Bmat_psi(3*i)
c
            Bpsimatrix(16,9*i)   = Bmat_psi(3*i-2)
            Bpsimatrix(17,9*i-1) = Bmat_psi(3*i-2)
            Bpsimatrix(18,9*i-5) = Bmat_psi(3*i)
c
         end do
!       print the Bpsimatrix in .dat file         
         write(6,*) 'this is Bpsimatrix'
         do i = 1, size(Bpsimatrix,1)
            write(6,'(20G12.4)') Bpsimatrix(i,:)
         end do
c
c
!       Assembly of Nrhomatrix for lagrange_multiplier D.O.F.
         do i = 1, 9
            Nrhomatrix(1,i) = 0.d0
            Nrhomatrix(2,i) = 0.d0 
            Nrhomatrix(3,i) = 0.d0
            Nrhomatrix(4,i) = 0.d0
            Nrhomatrix(5,i) = 0.d0
            Nrhomatrix(6,i) = 0.d0
            Nrhomatrix(7,i) = 0.d0
            Nrhomatrix(8,i) = 0.d0
            Nrhomatrix(9,i) = 0.d0
         end do
c~          do i = 1, 9
c~             Nrhomatrix(i,i) = 1.d0
c~          end do
!       print the Nrhomatrix in .dat file         
*         write(6,*) 'this is Nrhomatrix'
*         do i = 1, size(Nrhomatrix,1)
*            write(6,'(20G12.4)') Nrhomatrix(i,:)
*         end do
c         
c***********************************************************************
!       update of state variables.
!       displacement gradient = GradVariable = 0.d0
!       this is a transpose of Mmatrix
         T_Mmatrix = 0.d0
         T_Mmatrix = Transpose(Mmatrix)		
         do i = 1, 81
          do j = 1, 9
*            write(6,*) T_Mmatrix(i,j)
*            write(6,*) displacement(i)
*            write(6,*) T_Mmatrix(i,j)*displacement(i)
            GradVariable(j,kintk) = GradVariable(j,kintk) 
     1                              +T_Mmatrix(i,j)*displacement(i)
          end do
         end do 
c         
!       Print the displcement gradient in .dat file.           
         write(6,*) "this is gradient variable"
         do i = 1, size(GradVariable,1)
            write(6,'(20G12.4)') GradVariable(i,:)
         end do
c
c
!       strain = 0.d0
!       Deltastrain = 0.d0
!       This is a transpose of Bmatrix
         T_Bmatrix = 0.d0
	     T_Bmatrix = Transpose(Bmatrix)
         do i = 1, 81
          do j = 1, 6
            strain(j,kintk) = strain(j,kintk) 
     1                        +T_Bmatrix(i,j)*displacement(i)
c     
            Deltastrain(j,kintk) = Deltastrain(j,kintk) 
     1                             +T_Bmatrix(i,j)*delta_displacement(i)
          end do
         end do
c
!       Print the Strain in .dat file.
         write(6,*) "this is strain"
         do i = 1, size(strain,1)
            write(6,'(20G12.4)')  strain(i,:)
         end do
c~ !       Print the Delta_Strain in .dat file.         
         write(6,*) "this is Deltastrain"
         do i = 1, size(Deltastrain,1)
            write(6,'(20G12.4)')  Deltastrain(i,:)
         end do
c
c
!       relaxedstrain = 0.d0
!       This is a transpose of Npsimatrix
         T_Npsimatrix = 0.d0
	     T_Npsimatrix = Transpose(Npsimatrix)
         do i = 1, 72
          do j = 1, 9
            relaxedstrain(j,kintk) = relaxedstrain(j,kintk) 
     1                              +T_Npsimatrix(i,j)*relaxed_strain(i)
          end do  
         end do
c
!       Print the relaxedstrain in .dat file.
	     write(6,*) "this is relaxedstrain"
         do i = 1, size(relaxedstrain,1)
            write(6,'(20G12.4)')  relaxedstrain(i,:)
         end do
c
c
!       relaxedstraingradient = 0.d0
!       delta_relaxedstraingradient = 0.d0
!       This is a transpose of Bpsimatrix
         T_Bpsimatrix = 0.d0
	     T_Bpsimatrix = Transpose(Bpsimatrix)
         do i = 1, 72
          do j = 1, 18
            relaxedstraingradient(j,kintk)                        
     1         =relaxedstraingradient(j,kintk)
     1                    +T_Bpsimatrix(i,j)*relaxed_strain(i)
c~             write(6,*) T_Bpsimatrix(i,j)
c~             write(6,*) relaxed_strain(i)  
            delta_relaxedstraingradient(j,kintk)
     1                        =delta_relaxedstraingradient(j,kintk) 
     1                        +T_Bpsimatrix(i,j)*delta_relaxed_strain(i)
          end do
*          write(6,*) relaxedstraingradient(1,kintk)
         end do
c
!       Print the relaxedstraingradient in .dat file.
         write(6,*) "this is relaxedstraingradient"
         do i = 1, size(relaxedstraingradient,1)
            write(6,'(20G12.4)')  relaxedstraingradient(i,:)
         end do
!       Print the delta_relaxedstraingradient in .dat file.         
         write(6,*) "this is delta_relaxedstraingradient"
         do i = 1, size(delta_relaxedstraingradient,1)
            write(6,'(20G12.4)')  delta_relaxedstraingradient(i,:)
         end do
c
c         
         !Langrangemulti = 0.d0
         do i = 1, 9
          do j = 1, 9
            Langrangemulti(j,kintk) = Langrangemulti(j,kintk) 
     1                                +Nrhomatrix(i,j)*lagrange_multi(i)
          end do
         end do
c
!       Print the Langrangemultiplier in .dat file.         
c~          write(6,*) "this is la_multi"
c~          do i = 1, size(Langrangemulti,1)
c~             write(6,'(20G12.4)')  Langrangemulti(i,:)
c~          end do
c***********************************************************************
*          write(6,*) "this is lemda and mue"	
*          write(6,*) lemda
*          write(6,*) mue
*         write(6,*) "this is Deltastrain_1"
*         do i = 1, size(Deltastrain,1)
*            write(6,'(20G12.4)')  Deltastrain(i,:)
*         end do
!       call the KUMAT to find the state variables Sigma and Tau
!       and also get U_psilon and Lambda from KUMAT as output.
         call KUMAT(Sigma,Tau,U_psilon,Lambda,kintk,lemda,mue,
     1        micro_length,strain,relaxedstraingradient)
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
          write(6,*) "this is Upsi_matrix"
          do i = 1, size(U_psilon,1)
            write(6,'(20G12.4)')  U_psilon(i,:)
          end do
          write(6,*) "this is lamda_matrix"
          do i = 1, size(Lambda,1)
            write(6,'(20G12.4)')  Lambda(i,:)
          end do          
c***********************************************************************
c
!       store state variables Sigma and Tau in Svars vector. :
c~          do i = 1, 3
c~             svars(i+(kintk-1)*9) = Sigma(i,kintk)
c~          end do
c~          do i = 4, 9
c~             svars(i+(kintk-1)*9) = Tau(i-3,kintk)
c~          end do
c~ !       Print the State variable vector in .dat file.         
c~          write(6,*) "this is state vars"
c~          do i = 1, 9
c~             write(6,'(20G12.4)')  svars(i+(kintk-1)*9)
c~          end do
c
         do i = 1, 6
           svars(i+(kintk-1)*18) = Sigma(i,kintk)
         end do
         do i = 7, 12
           svars(i+(kintk-1)*18) = strain(i-6,kintk)
         end do
c     
c~          do i = 1, 6
c~             svars(i+(kintk-1)*22) = Sigma(i,kintk)
c~          end do
c~          do i = 4, 9
c~             svars(i+(kintk-1)*22) = Tau(i-3,kintk)
c~          end do
c~          do i = 10, 13
c~             svars(i+(kintk-1)*22) = relaxedstrain(i-9,kintk)
c~          end do
c~          do i = 14, 19
c~             svars(i+(kintk-1)*22) = relaxedstraingradient(i-13,kintk)
c~          end do
c~          do i = 20, 22
c~             svars(i+(kintk-1)*22) = strain(i-19,kintk)
c~          end do
!       Print the State variable vector in .dat file.         
c~          write(6,*) "this is state vars"
c~          do i = 1, 24
c~             write(6,'(20G12.4)')  svars(i+(kintk-1)*24)
c~          end do         
c***********************************************************************
!       Find the Element Stiffness Matrices:
c
*         write(6,*) "Jacobian of displacement_1"
*         write(6,*) djac       
*         write(6,*) "Gausss weight"
*         write(6,*) Gauss_weight(kintk) 
*         write(6,*) "T_Bmatrix"
*         do i = 1, size(T_Bmatrix,1)
*            write(6,'(20G12.4)')  T_Bmatrix(i,:)
*         end do 
!       Find the UU_stiffness matrix - dimension(18*18)
         do i = 1, 81
          do j = 1, 81
           do k = 1, 6
            do l = 1, 6
               stiffness_UU(i,j) = stiffness_UU(i,j) 
     1                 +T_Bmatrix(i,k)*U_psilon(k,l)*T_Bmatrix(j,l)*djac
     1                 *Gauss_weight(kintk)
            end do
           end do
          end do
         end do
!       Print the stiffness_UU matrix in .dat file.           
c~          write(6,*) "this is sti_UU"
c~          do i = 1, size(stiffness_UU,1)
c~             write(6,'(20G12.4)')  stiffness_UU(i,:)
c~          end do
c
c
!       Find the PsiPsi_stiffness matrix - dimension(16*16)
         do i = 1, 72
          do j = 1, 72
           do k = 1, 18
             do l = 1, 18
               stiffness_PsiPsi(i,j) = stiffness_PsiPsi(i,j) 
     1                  +T_Bpsimatrix(i,k)*Lambda(k,l)*T_Bpsimatrix(j,l)
     1                  *dj_psi*Gauss_weight(kintk)
             end do
           end do
          end do
         end do
!       Print the stiffness_PsiPsi matrix in .dat file.          
c~          write(6,*) "this is sti_psipsi"
c~          do i = 1, size(stiffness_PsiPsi,1)
c~             write(6,'(20G12.4)')  stiffness_PsiPsi(i,:)
c~          end do
c
c
!       Find the Urho_stiffness matrix - dimension(18*4)
         do i = 1, 81
          do j = 1, 9
           do k = 1, 9
             stiffness_Urho(i,j) = stiffness_Urho(i,j) 
     1          +T_Mmatrix(i,k)*Nrhomatrix(j,k)*djac*Gauss_weight(kintk)
           end do
          end do
         end do
!       Print the stiffness_Urho matrix in .dat file.        
c~          write(6,*) "this is stiffness_Urho"
c~          do i = 1, size(stiffness_Urho,1)
c~             write(6,'(20G12.4)')  stiffness_Urho(i,:)
c~          end do
c
c
!       Find the psirho_stiffness matrix - dimension(16*4)
         do i = 1, 72
          do j = 1, 9
           do k = 1, 9
               stiffness_psirho(i,j) = stiffness_psirho(i,j) 
     1                       +T_Npsimatrix(i,k)*Nrhomatrix(j,k)*dj_psi
     1                       *Gauss_weight(kintk)
           end do
          end do
         end do
!       Print the stiffness_psirho matrix in .dat file.         
c~          write(6,*) "this is stiffness_psirho"
c~          do i = 1, size(stiffness_psirho,1)
c~             write(6,'(20G12.4)')  stiffness_psirho(i,:)
c~          end do
c
c
!       Find the Mass_U matrix - dimension(18*18)
c~          T_Nmatrix = Transpose(Nmatrix)
c~          do i = 1, 18
c~           do j = 1, 18
c~            do k = 1, 2
c~                Mass_U(i,j) = Mass_U(i,j) + Density*T_Nmatrix(i,k)
c~      1                       *T_Nmatrix(j,k)*djac*Gauss_weight(kintk)
c~            end do
c~           end do
c~          end do
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
         do i = 1, 81 
            do j = 1, 6
               F_vector(i) = F_vector(i) + T_Bmatrix(i,j)*Sigma(j,kintk)
     1                                     *djac*Gauss_weight(kintk)
*               write(6,*) T_Bmatrix(i,j)*Sigma(j,kintk)
*     1                                     *djac*Gauss_weight(kintk)
            end do
            do j = 1, 9
               F_vector(i) = F_vector(i) - T_Mmatrix(i,j) 
     1                 *Langrangemulti(j,kintk)*djac*Gauss_weight(kintk)  
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
         do i = 1, 72
            do j = 1, 18
              R_vector(i) = R_vector(i) + T_Bpsimatrix(i,j)*Tau(j,kintk)
     1                                     *dj_psi*Gauss_weight(kintk)
            end do
            do j = 1, 9
               R_vector(i) = R_vector(i) + T_Npsimatrix(i,j)
     1             *Langrangemulti(j,kintk)*dj_psi*Gauss_weight(kintk)
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
         do i = 1, 9
            do j = 1, 9
               S_vector(i) = S_vector(i) + Nrhomatrix(i,j)
     1           *(relaxedstrain(j,kintk)*dj_psi-GradVariable(j,kintk)
     1           *djac)*Gauss_weight(kintk)
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
      do i = 1, 81
         do j = 1, 81
            amatrx(i,j) = amatrx(i,j) + stiffness_UU(i,j)
         end do
      end do
c
!     amat(19:34,19:34) =  stiffness_PsiPsi(1:16,1:16)
      do i = 82, 153
         do j = 82, 153
            amatrx(i,j) = amatrx(i,j) + stiffness_PsiPsi(i-81,j-81)   
         end do
      end do
c
!     amat(35:38,19:34) =  tra_stiffness_psirho(1:4,1:16)
      do i = 154, 162
         do j = 82, 153 
            amatrx(i,j) = amatrx(i,j) + tra_stiffness_psirho(i-153,j-81)
         end do
      end do
c
!     amat(35:38,1:18) =  tra_stiffness_Urho(1:4,1:18)
      do i = 154, 162
         do j = 1, 81
            amatrx(i,j) = amatrx(i,j) - tra_stiffness_Urho(i-153,j)
         end do
      end do
c
!     amat(1:18,35:38) =  stiffness_Urho(1:18,1:4)
      do i = 1, 81
         do j = 154, 162
            amatrx(i,j) = amatrx(i,j) - stiffness_Urho(i,j-153)
         end do
      end do
c
!     amat(19:34,35:38) =  stiffness_psirho(1:16,1:4)
      do i = 82,153
         do j = 154,162
            amatrx(i,j) = amatrx(i,j) + stiffness_psirho(i-81,j-153)
         end do
      end do
c
c
!     Initialize the mass matrix.
c~       do i = 1, 162
c~             Mass_matrix(i,i) = 1.d0
c~       end do
c
!     assign the amat to the amatrx which is common matrix for
!     stiffness.
!      do i = 1, 38
!         do j = 1, 38
!            amatrx(i,j) = amatrx(i,j) + amat(i,j)
!         end do
!      end do
c
!     assign the 0 value to the diagonal elements(35,36,37,38) of amat.
      do i=153, 162
	     amatrx(i,i) = 1.0e-10
      end do
c
!     if lflags(3) = 1 then Mass_matrix is assigned to the amatrx.
!      if(lflags(3) .eq. 1) then
!            do i = 1, 38
!                do j = 1, 38
!                   amatrx(i,j) = Mass_matrix(i,j)
!                end do
!            end do
!      endif
c      
!     print the amatrx into the .dat file.      
      write(6,*) "amatrx"
      do i = 1, size(amatrx,1)
          write(6,'(40G12.4)')  amatrx(i,:)
      end do
c
c
c
!     Assemble the RHS matrix.
!     rhs(1:18) = F_vector(1:18)
      do i = 1, 81
         rhs(i) = rhs(i) + F_vector(i)
      end do
c
c
!     rhs(19:34) = R_vector(1:16)
      do i = 82, 153
         rhs(i) = rhs(i) + R_vector(i-81)
      end do
c
c
!     rhs(35:38) = S_vector(1:4)
      do i = 154, 162
         rhs(i) = rhs(i) + S_vector(i-153)
      end do
c
c
c~       if (NDLOAD .ge. 1) then
c~           write(*,*) "call KDLOAD"
c~           call KDLOAD(rhs_k,NDLOAD,coords,mdload,ADLMAG,jdltyp)
c~           do i = 1, 16
c~             rhs(i) = rhs_k(i) - rhs(i) 
c~           end do          
c~       end if
c
c
!     print the RHS matrix into the .dat file.      
      write(6,*) "rhs"
      do i = 1, 162
          write(6,'(20G12.4)')  rhs(i)
      end do
c
c
      return
      end subroutine uel       
!     END of the main subroutine UEL.   
c***********************************************************************
      subroutine shapefcn_U(dN_U,dN_U_X,kintk,ninpt,nnode,ndim)
c
      include 'aba_param.inc'
c
      parameter (Gauss_Coord=0.7745966692d0)
      dimension dN_U(nnode),dN_U_X(ndim,nnode),Gauss_coord_X(27),
     1 Gauss_coord_Y(27),Gauss_coord_Z(27)
c     
      data Gauss_coord_X /-1.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0,
     1 0.d0,-1.d0,1.d0,0.d0,-1.d0,1.d0,-1.d0,1.d0,0.d0,-1.d0,1.d0,0.d0,
     1 0.d0,0.d0,-1.d0,1.d0,0.d0,0.d0,0.d0/                  
c      
      data Gauss_coord_Y /-1.d0,-1.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0,1.d0,
     1 -1.d0,0.d0,0.d0,1.d0,-1.d0,-1.d0,1.d0,1.d0,-1.d0,0.d0,0.d0,1.d0,
     1 0.d0,-1.d0,0.d0,0.d0,1.d0,0.d0,0.d0/
c           
      data Gauss_coord_Z /-1.d0,-1.d0,-1.d0,-1.d0,1.d0,1.d0,1.d0,1.d0,
     1 -1.d0,-1.d0,-1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,1.d0,1.d0,
     1 -1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0/     
c     
c  3D 27-nodes
c
c     determine (x,y,z)
        x = Gauss_coord_X(kintk)*Gauss_Coord
        y = Gauss_coord_Y(kintk)*Gauss_Coord
        z = Gauss_coord_Z(kintk)*Gauss_Coord
c
c     shape functions
        dN_U(1) = 0.125d0*x*(x+1.d0)*y*(y-1.d0)*z*(z+1.d0)
        dN_U(2) = 0.125d0*x*(x+1.d0)*y*(y+1.d0)*z*(z+1.d0)
        dN_U(3) = 0.125d0*x*(x+1.d0)*y*(y+1.d0)*z*(z-1.d0)
        dN_U(4) = 0.125d0*x*(x+1.d0)*y*(y-1.d0)*z*(z-1.d0)
        dN_U(5) = 0.125d0*x*(x-1.d0)*y*(y-1.d0)*z*(z+1.d0)
        dN_U(6) = 0.125d0*x*(x-1.d0)*y*(y+1.d0)*z*(z+1.d0) 
        dN_U(7) = 0.125d0*x*(x-1.d0)*y*(y+1.d0)*z*(z-1.d0)               
        dN_U(8) = 0.125d0*x*(x-1.d0)*y*(y-1.d0)*z*(z-1.d0)
c        
        dN_U(9) = 0.25d0*x*(x+1.d0)*(1.d0-y*y)*z*(z+1.d0)
        dN_U(10) = 0.25d0*x*(x+1.d0)*y*(y+1.d0)*(1.d0-z*z)
        dN_U(11) = 0.25d0*x*(x+1.d0)*(1.d0-y*y)*z*(z-1.d0)
        dN_U(12) = 0.25d0*x*(x+1.d0)*y*(y-1.d0)*(1.d0-z*z)
        dN_U(13) = 0.25d0*x*(x-1.d0)*(1.d0-y*y)*z*(z+1.d0)
        dN_U(14) = 0.25d0*x*(x-1.d0)*y*(y+1.d0)*(1.d0-z*z)  
        dN_U(15) = 0.25d0*x*(x-1.d0)*(1.d0-y*y)*z*(z-1.d0)
        dN_U(16) = 0.25d0*x*(x-1.d0)*y*(y-1.d0)*(1.d0-z*z)
        dN_U(17) = 0.25d0*(1.d0-x*x)*y*(y-1.d0)*z*(z+1.d0)
        dN_U(18) = 0.25d0*(1.d0-x*x)*y*(y+1.d0)*z*(z+1.d0)
        dN_U(19) = 0.25d0*(1.d0-x*x)*y*(y+1.d0)*z*(z-1.d0)
        dN_U(20)  = 0.25d0*(1.d0-x*x)*y*(y-1.d0)*z*(z-1.d0)   
c        
        dN_U(21) = 0.5d0*x*(x-1.d0)*(1.d0-y*y)*(1.d0-z*z)
        dN_U(22) = 0.5d0*x*(x+1.d0)*(1.d0-y*y)*(1.d0-z*z)
        dN_U(23) = 0.5d0*(1.d0-x*x)*y*(y+1.d0)*(1.d0-z*z)
        dN_U(24) = 0.5d0*(1.d0-x*x)*y*(y-1.d0)*(1.d0-z*z)
        dN_U(25) = 0.5d0*(1.d0-x*x)*(1.d0-y*y)*z*(z-1.d0)
        dN_U(26) = 0.5d0*(1.d0-x*x)*(1.d0-y*y)*z*(z+1.d0)   
c
        dN_U(27) = (1.d0-x*x)*(1.d0-y*y)*(1.d0-z*z)
c
c     derivative d(Ni)/d(g)
        dN_U_X(1,1) = 0.125d0*(2.d0*x+1.d0)*(y*y-y)*(z*z+z)
        dN_U_X(1,2) = 0.125d0*(2.d0*x+1.d0)*(y*y+y)*(z*z+z)
        dN_U_X(1,3) = 0.125d0*(2.d0*x+1.d0)*(y*y+y)*(z*z-z)
        dN_U_X(1,4) = 0.125d0*(2.d0*x+1.d0)*(y*y-y)*(z*z-z)
        dN_U_X(1,5) = 0.125d0*(2.d0*x-1.d0)*(y*y-y)*(z*z+z)
        dN_U_X(1,6) = 0.125d0*(2.d0*x-1.d0)*(y*y+y)*(z*z+z)
        dN_U_X(1,7) = 0.125d0*(2.d0*x-1.d0)*(y*y+y)*(z*z-z)
        dN_U_X(1,8) = 0.125d0*(2.d0*x-1.d0)*(y*y-y)*(z*z-z)
c        
        dN_U_X(1,9) =0.25d0*(2.d0*x+1.d0)*(1.d0-y*y)*z*(z+1.d0)
        dN_U_X(1,10) =0.25d0*(2.d0*x+1.d0)*y*(y+1.d0)*(1.d0-z*z)
        dN_U_X(1,11) =0.25d0*(2.d0*x+1.d0)*(1.d0-y*y)*z*(z-1.d0)
        dN_U_X(1,12) =0.25d0*(2.d0*x+1.d0)*y*(y-1.d0)*(1.d0-z*z)
        dN_U_X(1,13) =0.25d0*(2.d0*x-1.d0)*(1.d0-y*y)*z*(z+1.d0)
        dN_U_X(1,14) =0.25d0*(2.d0*x-1.d0)*y*(y+1.d0)*(1.d0-z*z)
        dN_U_X(1,15) =0.25d0*(2.d0*x-1.d0)*(1.d0-y*y)*z*(z-1.d0)
        dN_U_X(1,16) =0.25d0*(2.d0*x-1.d0)*y*(y-1.d0)*(1.d0-z*z)
        dN_U_X(1,17) =-0.5d0*x*y*(y-1.d0)*z*(z+1.d0)
        dN_U_X(1,18) =-0.5d0*x*y*(y+1.d0)*z*(z+1.d0)
        dN_U_X(1,19) =-0.5d0*x*y*(y+1.d0)*z*(z-1.d0)                                        
        dN_U_X(1,20) = -0.5d0*x*y*(y-1.d0)*z*(z-1.d0)
c        
        dN_U_X(1,21) = 0.5d0*(2.d0*x-1.d0)*(1.d0-y*y)*(1.d0-z*z)
        dN_U_X(1,22) = 0.5d0*(2.d0*x+1.d0)*(1.d0-y*y)*(1.d0-z*z)
        dN_U_X(1,23) = -0.5d0*(2.d0*x)*y*(y+1.d0)*(1.d0-z*z)
        dN_U_X(1,24) = -0.5d0*(2.d0*x)*y*(y-1.d0)*(1.d0-z*z)
        dN_U_X(1,25) = -0.5d0*(2.d0*x)*(1.d0-y*y)*z*(z-1.d0)
        dN_U_X(1,26) = -0.5d0*(2.d0*x)*(1.d0-y*y)*z*(z+1.d0)
c        
        dN_U_X(1,27) = -(2.d0*x)*(1.d0-y*y)*(1.d0-z*z)
c
c     derivative d(Ni)/d(h)
        dN_U_X(2,1) = 0.125d0*(x*x+x)*(2.d0*y-1.d0)*(z*z+z)
        dN_U_X(2,2) = 0.125d0*(x*x+x)*(2.d0*y+1.d0)*(z*z+z)
        dN_U_X(2,3) = 0.125d0*(x*x+x)*(2.d0*y+1.d0)*(z*z-z)
        dN_U_X(2,4) = 0.125d0*(x*x+x)*(2.d0*y-1.d0)*(z*z-z)
        dN_U_X(2,5) = 0.125d0*(x*x-x)*(2.d0*y-1.d0)*(z*z+z)
        dN_U_X(2,6) = 0.125d0*(x*x-x)*(2.d0*y+1.d0)*(z*z+z)
        dN_U_X(2,7) = 0.125d0*(x*x-x)*(2.d0*y+1.d0)*(z*z-z)                        
        dN_U_X(2,8) = 0.125d0*(x*x-x)*(2.d0*y-1.d0)*(z*z-z)
c        
        dN_U_X(2,9) = -0.25d0*x*(x+1.d0)*(2.d0*y)*z*(z+1.d0)
        dN_U_X(2,10) = 0.25d0*x*(x+1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)
        dN_U_X(2,11) = -0.25d0*x*(x+1.d0)*(2.d0*y)*z*(z-1.d0)                
        dN_U_X(2,12) = 0.25d0*x*(x+1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)
        dN_U_X(2,13) = -0.25d0*x*(x-1.d0)*(2.d0*y)*z*(z+1.d0)
        dN_U_X(2,14) = 0.25d0*x*(x-1.d0)*(2.d0*y+1.d0)*(1.d0-z*z)
        dN_U_X(2,15) = -0.25d0*x*(x-1.d0)*(2.d0*y)*z*(z-1.d0)
        dN_U_X(2,16) = 0.25d0*x*(x-1.d0)*(2.d0*y-1.d0)*(1.d0-z*z)
        dN_U_X(2,17) = 0.25d0*(1.d0-x*x)*(2.d0*y-1.d0)*z*(z+1.d0)
        dN_U_X(2,18) = 0.25d0*(1.d0-x*x)*(2.d0*y+1.d0)*z*(z+1.d0)                                        
        dN_U_X(2,19) = 0.25d0*(1.d0-x*x)*(2.d0*y+1.d0)*z*(z-1.d0)
        dN_U_X(2,20) = 0.25d0*(1.d0-x*x)*(2.d0*y-1.d0)*z*(z-1.d0)
c        
        dN_U_X(2,21) = -0.5d0*x*(x-1.d0)*(2.d0*y)*(1.d0-z*z)
        dN_U_X(2,22) = -0.5d0*x*(x+1.d0)*(2.d0*y)*(1.d0-z*z)
        dN_U_X(2,23) = 0.5d0*(1.d0-x*x)*(2.d0*y+1.d0)*(1.d0-z*z)
        dN_U_X(2,24) = 0.5d0*(1.d0-x*x)*(2.d0*y-1.d0)*(1.d0-z*z)
        dN_U_X(2,25) = -0.5d0*(1.d0-x*x)*(2.d0*y)*z*(z-1.d0)
        dN_U_X(2,26) = -0.5d0*(1.d0-x*x)*(2.d0*y)*z*(z+1.d0)
c        
        dN_U_X(2,27) = -(1.d0-x*x)*(2.d0*y)*(1.d0-z*z)
c   
c     derivative d(Ni)/d(h)
c
        dN_U_X(3,1) = 0.125d0*(x*x+x)*(y*y-y)*(2.d0*z+1.d0)
        dN_U_X(3,2) = 0.125d0*(x*x+x)*(y*y+y)*(2.d0*z+1.d0)
        dN_U_X(3,3) = 0.125d0*(x*x+x)*(y*y+y)*(2.d0*z-1.d0)        
        dN_U_X(3,4) = 0.125d0*(x*x+x)*(y*y-y)*(2.d0*z-1.d0)
        dN_U_X(3,5) = 0.125d0*(x*x-x)*(y*y-y)*(2.d0*z+1.d0)
        dN_U_X(3,6) = 0.125d0*(x*x-x)*(y*y+y)*(2.d0*z+1.d0)
        dN_U_X(3,7) = 0.125d0*(x*x-x)*(y*y+y)*(2.d0*z-1.d0)        
        dN_U_X(3,8) = 0.125d0*(x*x-x)*(y*y-y)*(2.d0*z-1.d0)
c        
        dN_U_X(3,9) = 0.25d0*x*(x+1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)
        dN_U_X(3,10) = -0.25d0*x*(x+1.d0)*y*(y+1.d0)*(2.d0*z)
        dN_U_X(3,11) = 0.25d0*x*(x+1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)        
        dN_U_X(3,12) = -0.25d0*x*(x+1.d0)*y*(y-1.d0)*(2.d0*z)
        dN_U_X(3,13) = 0.25d0*x*(x-1.d0)*(1.d0-y*y)*(2.d0*z+1.d0)
        dN_U_X(3,14) = -0.25d0*x*(x-1.d0)*y*(y+1.d0)*(2.d0*z)
        dN_U_X(3,15) = 0.25d0*x*(x-1.d0)*(1.d0-y*y)*(2.d0*z-1.d0)
        dN_U_X(3,16) = -0.25d0*x*(x-1.d0)*y*(y-1.d0)*(2.d0*z)
        dN_U_X(3,17) = 0.25d0*(1.d0-x*x)*y*(y-1.d0)*(2.d0*z+1.d0)
        dN_U_X(3,18) = 0.25d0*(1.d0-x*x)*y*(y+1.d0)*(2.d0*z+1.d0)
        dN_U_X(3,19) = 0.25d0*(1.d0-x*x)*y*(y+1.d0)*(2.d0*z-1.d0)                                                
        dN_U_X(3,20) = 0.25d0*(1.d0-x*x)*y*(y-1.d0)*(2.d0*z-1.d0)
c        
        dN_U_X(3,21) = -0.5d0*x*(x-1.d0)*(1.d0-y*y)*(2.d0*z)
        dN_U_X(3,22) = -0.5d0*x*(x+1.d0)*(1.d0-y*y)*(2.d0*z)
        dN_U_X(3,23) = -0.5d0*(1.d0-x*x)*y*(y+1.d0)*(2.d0*z)
        dN_U_X(3,24) = -0.5d0*(1.d0-x*x)*y*(y-1.d0)*(2.d0*z)
        dN_U_X(3,25) = 0.5d0*(1.d0-x*x)*(1.d0-y*y)*(2.d0*z-1.d0)
        dN_U_X(3,26) = 0.5d0*(1.d0-x*x)*(1.d0-y*y)*(2.d0*z+1.d0)
c        
        dN_U_X(3,27) = -(1.d0-x*x)*(1.d0-y*y)*(2.d0*z)
c 
      return
      end subroutine shapefcn_U
c***********************************************************************
c***********************************************************************
      subroutine shapefcn_PSI(dN_PSI,dN_PSI_X,kintk,ninpt,nnode,ndim)
c
      include 'aba_param.inc'
c
      parameter (Gauss_Coord=0.7745966692d0)
      dimension dN_PSI(nnode),dN_PSI_X(ndim,nnode),Gauss_coord_X(27),
     1 Gauss_coord_Y(27),Gauss_coord_Z(27)
c     
c     
      data Gauss_coord_X /-1.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0,-1.d0,1.d0,
     1 0.d0,-1.d0,1.d0,0.d0,-1.d0,1.d0,-1.d0,1.d0,0.d0,-1.d0,1.d0,0.d0,
     1 0.d0,0.d0,-1.d0,1.d0,0.d0,0.d0,0.d0/                  
c      
      data Gauss_coord_Y /-1.d0,-1.d0,1.d0,1.d0,-1.d0,-1.d0,1.d0,1.d0,
     1 -1.d0,0.d0,0.d0,1.d0,-1.d0,-1.d0,1.d0,1.d0,-1.d0,0.d0,0.d0,1.d0,
     1 0.d0,-1.d0,0.d0,0.d0,1.d0,0.d0,0.d0/
c           
      data Gauss_coord_Z /-1.d0,-1.d0,-1.d0,-1.d0,1.d0,1.d0,1.d0,1.d0,
     1 -1.d0,-1.d0,-1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,1.d0,1.d0,1.d0,
     1 -1.d0,0.d0,0.d0,0.d0,0.d0,1.d0,0.d0/
C     
C  3D 27-nodes
c
c     determine (x,y,z)
        x = Gauss_coord_X(kintk)*Gauss_Coord
        y = Gauss_coord_Y(kintk)*Gauss_Coord
        z = Gauss_coord_Z(kintk)*Gauss_Coord
c     shape functions
        dN_PSI(1) = 0.125d0*(1.d0-x)*(1.d0-y)*(1.d0-z)
        dN_PSI(2) = 0.125d0*(1.d0+x)*(1.d0-y)*(1.d0-z)
        dN_PSI(3) = 0.125d0*(1.d0+x)*(1.d0+y)*(1.d0-z)
        dN_PSI(4) = 0.125d0*(1.d0-x)*(1.d0+y)*(1.d0-z)
        dN_PSI(5) = 0.125d0*(1.d0-x)*(1.d0-y)*(1.d0+z)
        dN_PSI(6) = 0.125d0*(1.d0+x)*(1.d0-y)*(1.d0+z)
        dN_PSI(7) = 0.125d0*(1.d0+x)*(1.d0+y)*(1.d0+z)
        dN_PSI(8) = 0.125d0*(1.d0-x)*(1.d0+y)*(1.d0+z)
c     derivative d(Ni)/d(xi)
        dN_PSI_X(1,1) = -0.125d0*(1.d0-y)*(1.d0-z)
        dN_PSI_X(1,2) = 0.125d0*(1.d0-y)*(1.d0-z)
        dN_PSI_X(1,3) = 0.125d0*(1.d0+y)*(1.d0-z)
        dN_PSI_X(1,4) = -0.125d0*(1.d0+y)*(1.d0-z)
        dN_PSI_X(1,5) = -0.125d0*(1.d0-y)*(1.d0+z)
        dN_PSI_X(1,6) = 0.125d0*(1.d0-y)*(1.d0+z)
        dN_PSI_X(1,7) = 0.125d0*(1.d0+y)*(1.d0+z)
        dN_PSI_X(1,8) = -0.125d0*(1.d0+y)*(1.d0+z)        
c     derivative d(Ni)/d(omega)
        dN_PSI_X(2,1) = -0.125d0*(1.d0-x)*(1.d0-z)
        dN_PSI_X(2,2) = -0.125d0*(1.d0+x)*(1.d0-z)
        dN_PSI_X(2,3) = 0.125d0*(1.d0+x)*(1.d0-z)
        dN_PSI_X(2,4) = 0.125d0*(1.d0-x)*(1.d0-z)
        dN_PSI_X(2,5) = -0.125d0*(1.d0-x)*(1.d0+z)
        dN_PSI_X(2,6) = -0.125d0*(1.d0+x)*(1.d0+z)
        dN_PSI_X(2,7) = 0.125d0*(1.d0+x)*(1.d0+z)
        dN_PSI_X(2,8) = 0.125d0*(1.d0-x)*(1.d0+z)
c     derivative d(Ni)/d(omega)
        dN_PSI_X(3,1) = -0.125d0*(1.d0-x)*(1.d0-y)
        dN_PSI_X(3,2) = -0.125d0*(1.d0+x)*(1.d0-y)
        dN_PSI_X(3,3) = -0.125d0*(1.d0+x)*(1.d0+y)
        dN_PSI_X(3,4) = -0.125d0*(1.d0-x)*(1.d0+y)
        dN_PSI_X(3,5) = 0.125d0*(1.d0-x)*(1.d0-y)
        dN_PSI_X(3,6) = 0.125d0*(1.d0+x)*(1.d0-y)
        dN_PSI_X(3,7) = 0.125d0*(1.d0+x)*(1.d0+y)
        dN_PSI_X(3,8) = 0.125d0*(1.d0-x)*(1.d0+y)
c     
      return
      end subroutine shapefcn_PSI
c***********************************************************************
      subroutine Jacobian_UU(djac,xjaci,jelem,ndim,nnode,coords,dN_U_X)
c
c     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix 
c                          
      include 'aba_param.inc'
      dimension xjac(3,3),xjaci(3,3),coords(3,*),dN_U_X(3,27)
      djac=0.d0
c
      do i = 1, 3
        do j = 1, 3
          xjac(i,j)  = 0.d0
          xjaci(i,j) = 0.d0
        end do
      end do
c
      do inod= 1, 27
         do i = 1, 3
           do j = 1, 3
              xjac(i,j) = xjac(i,j)+dN_U_X(i,inod)*coords(j,inod) 
           end do
         end do 
      end do
*      write(*,*) "this is xjac"
*      do i = 1, size(xjac,1)
*            write(*,'(20G12.4)') xjac(i,:)
*      end do      
c
      djac =xjac(1,1)*xjac(2,2)*xjac(3,3)+xjac(1,2)*xjac(3,1)*xjac(2,3)
     1     +xjac(1,3)*xjac(2,1)*xjac(3,2)-xjac(1,1)*xjac(3,2)*xjac(2,3)
     1     -xjac(1,2)*xjac(2,1)*xjac(3,3)-xjac(1,3)*xjac(3,1)*xjac(2,2)
*      write(*,*) "det of jacobian"
*      write(*,*) djac

         ! jacobian is positive - o.k.
      xjaci(1,1) =  (xjac(2,2)*xjac(3,3)-xjac(2,3)*xjac(3,2))/djac
      xjaci(1,2) =  -(xjac(1,2)*xjac(3,3)-xjac(3,2)*xjac(1,3))/djac
      xjaci(1,3) =  (xjac(1,2)*xjac(2,3)-xjac(2,2)*xjac(1,3))/djac
      xjaci(2,1) =  -(xjac(2,1)*xjac(3,3)-xjac(3,1)*xjac(2,3))/djac
      xjaci(2,2) =  (xjac(1,1)*xjac(3,3)-xjac(1,3)*xjac(3,1))/djac
      xjaci(2,3) =  -(xjac(1,1)*xjac(2,3)-xjac(2,1)*xjac(1,3))/djac
      xjaci(3,1) =  (xjac(2,1)*xjac(3,2)-xjac(2,2)*xjac(3,1))/djac
      xjaci(3,2) =  -(xjac(1,1)*xjac(3,2)-xjac(1,2)*xjac(3,1))/djac
      xjaci(3,3) =  (xjac(2,2)*xjac(1,1)-xjac(2,1)*xjac(1,2))/djac
*         do i = 1, size(xjaci,1)
*            write(*,'(20G12.4)') xjaci(i,:)
*         end do           
c
      return
      end subroutine Jacobian_UU
c***********************************************************************
c***********************************************************************
      subroutine Jacobian_psipsi(dj_psi,xjaci_psi,jelem,ndim,nnode,
     1 coords,dN_PSI_X)
c
c     Notation: djac - Jac determinant; xjaci - inverse of Jac matrix 
c                          
      include 'aba_param.inc'
c
      dimension xj_p(3,3),xjaci_psi(3,3),coords(3,27),dN_PSI_X(3,8)
      dj_psi = 0.d0
c
      do i = 1, 3
        do j = 1, 3
          xj_p(i,j)  = 0.d0
          xjaci_psi(i,j) = 0.d0
        end do
      end do
c
      do inod= 1, 8
         do i = 1, 3
           do j = 1, 3
             xj_p(i,j)=xj_p(i,j)+dN_PSI_X(i,inod)*coords(j,inod) 
           end do
         end do 
      end do
*      write(*,*) "this is xjac"
*      do i = 1, size(xjac,1)
*            write(*,'(20G12.4)') xjac(i,:)
*      end do      
c     xjac_
      dj_psi=xj_p(1,1)*xj_p(2,2)*xj_p(3,3)+xj_p(1,2)*xj_p(3,1)*xj_p(2,3)
     1     +xj_p(1,3)*xj_p(2,1)*xj_p(3,2)-xj_p(1,1)*xj_p(3,2)*xj_p(2,3)
     1     -xj_p(1,2)*xj_p(2,1)*xj_p(3,3)-xj_p(1,3)*xj_p(3,1)*xj_p(2,2)
*      write(*,*) "det of jacobian"
*      write(*,*) djac

         ! jacobian is positive - o.k.
      xjaci_psi(1,1)=  (xj_p(2,2)*xj_p(3,3)-xj_p(2,3)*xj_p(3,2))/dj_psi
      xjaci_psi(1,2)=  -(xj_p(1,2)*xj_p(3,3)-xj_p(3,2)*xj_p(1,3))/dj_psi
      xjaci_psi(1,3)=  (xj_p(1,2)*xj_p(2,3)-xj_p(2,2)*xj_p(1,3))/dj_psi
      xjaci_psi(2,1)=  -(xj_p(2,1)*xj_p(3,3)-xj_p(3,1)*xj_p(2,3))/dj_psi
      xjaci_psi(2,2)=  (xj_p(1,1)*xj_p(3,3)-xj_p(1,3)*xj_p(3,1))/dj_psi
      xjaci_psi(2,3)=  -(xj_p(1,1)*xj_p(2,3)-xj_p(2,1)*xj_p(1,3))/dj_psi
      xjaci_psi(3,1)=  (xj_p(2,1)*xj_p(3,2)-xj_p(2,2)*xj_p(3,1))/dj_psi
      xjaci_psi(3,2)=  -(xj_p(1,1)*xj_p(3,2)-xj_p(1,2)*xj_p(3,1))/dj_psi
      xjaci_psi(3,3)=  (xj_p(2,2)*xj_p(1,1)-xj_p(2,1)*xj_p(1,2))/dj_psi
      return
      end subroutine Jacobian_psipsi
c***********************************************************************
c***********************************************************************
      subroutine Bmatrix_U(Bmat_u,xjaci,dN_U_X,nnode,ndim)
c
c     Notation: bmat(i) ....dN1/dx, dN1/dy, dN2/dx, dN2/dy...      
         include 'aba_param.inc'
         dimension Bmat_u(*),dN_U_X(ndim,*),xjaci(ndim,*)
c    
         do i = 1, nnode*ndim
            Bmat_u(i) = 0.d0
         end do
c    
         do inod = 1, 27
          do j = 1, ndim
           do i = 1, ndim
            icomp = j + (inod - 1)*ndim
c            
            Bmat_u(icomp)=Bmat_u(icomp)+xjaci(j,i)*dN_U_X(i,inod)
           end do
          end do
         end do 
c    
      return
      end subroutine Bmatrix_U
c    
c***********************************************************************
c***********************************************************************
      subroutine Bmatrix_PSIPSI(Bmat_psi,xjaci_psi,dN_PSI_X,nnode,ndim)
c
c     Notation: Bmat_psi(i) ....dN1/dx, dN1/dy, dN2/dx, dN2/dy...      
         include 'aba_param.inc'
         dimension Bmat_psi(*),dN_PSI_X(ndim,*),xjaci_psi(ndim,*)
c    
         do i = 1, nnode*ndim
            Bmat_psi(i) = 0.d0
         end do
c    
         do inod = 1, 8
          do j = 1, ndim
           do i = 1, ndim
            icomp = j + (inod - 1)*ndim
c            
            Bmat_psi(icomp)=Bmat_psi(icomp)
     1                      +xjaci_psi(j,i)*dN_PSI_X(i,inod)
           end do
          end do
         end do 
c    
      return
      end subroutine Bmatrix_PSIPSI
c    
c***********************************************************************
c***********************************************************************
      subroutine KUMAT(Sigma,Tau,C_matrix,D_matrix,kintk,lemda,mue,
     1 micro_length,strain,relaxedstraingradient)
c
      double precision::Sigma,Tau,C_matrix,D_matrix,mue,micro_length
      double precision::lemda,relaxedstraingradient,strain
c
      dimension C_matrix(6,6), D_matrix(18,18), Sigma(6,27), Tau(18,27),
     1 Lambda(6,6), U_psilon(18,18), strain(6,27),
     2 relaxedstraingradient(18,27)
*      write(*,*) "this is lemda and mue"
      write(*,*) lemda
      write(*,*) mue
      write(*,*) micro_length
c
c
      do i = 1, 6
        do j = 1, 6 
          C_matrix(i,j) = 0.d0
        end do
      end do     
c
      C_matrix(1,1) = lemda + 2.d0*mue
      C_matrix(2,2) = lemda + 2.d0*mue
      C_matrix(3,3) = lemda + 2.d0*mue
      C_matrix(1,2) = lemda
      C_matrix(1,3) = lemda
      C_matrix(2,1) = lemda
      C_matrix(2,3) = lemda
      C_matrix(3,1) = lemda
      C_matrix(3,2) = lemda
      C_matrix(4,4) = 0.5d0*mue
      C_matrix(5,5) = 0.5d0*mue
      C_matrix(6,6) = 0.5d0*mue
      
c
      write(6,*) "this is C_matrix in UMAT"
      do i = 1, size(C_matrix,1)
        write(6,'(20G12.4)')  C_matrix(i,:)
      end do
      do i = 1, 18
         do j = 1, 18
            D_matrix(i,j) = 0.d0
         end do
      end do
c      
c
      D_matrix(4,4) = 1.d0*mue*micro_length*micro_length
      D_matrix(7,7) = 1.d0*mue*micro_length*micro_length
      D_matrix(8,8) = 1.d0*mue*micro_length*micro_length
      D_matrix(11,11) = 1.d0*mue*micro_length*micro_length
      D_matrix(12,12) = 1.d0*mue*micro_length*micro_length
      D_matrix(15,15) = 1.d0*mue*micro_length*micro_length
      D_matrix(4,5) = -0.5d0*mue*micro_length*micro_length
      D_matrix(5,4) = -0.5d0*mue*micro_length*micro_length
      D_matrix(6,7) = -0.5d0*mue*micro_length*micro_length
      D_matrix(7,6) = -0.5d0*mue*micro_length*micro_length
      D_matrix(10,11) = -0.5d0*mue*micro_length*micro_length
      D_matrix(11,10) = -0.5d0*mue*micro_length*micro_length
      D_matrix(12,13) = -0.5d0*mue*micro_length*micro_length
      D_matrix(13,12) = -0.5d0*mue*micro_length*micro_length
      D_matrix(14,15) = -0.5d0*mue*micro_length*micro_length
      D_matrix(15,14) = -0.5d0*mue*micro_length*micro_length
      D_matrix(5,5) = 0.25d0*mue*micro_length*micro_length
      D_matrix(6,6) = 0.25d0*mue*micro_length*micro_length
      D_matrix(9,9) = 0.25d0*mue*micro_length*micro_length
      D_matrix(10,10) = 0.25d0*mue*micro_length*micro_length
      D_matrix(13,13) = 0.25d0*mue*micro_length*micro_length
      D_matrix(14,14) = 0.25d0*mue*micro_length*micro_length
      D_matrix(16,16) = 1.d0*mue*micro_length*micro_length
      D_matrix(17,17) = 1.d0*mue*micro_length*micro_length
      D_matrix(18,18) = 1.d0*mue*micro_length*micro_length
c      
      write(6,*) "this is D_matrix in UMAT"
      do i = 1, size(D_matrix,1)
        write(6,'(20G12.4)')  D_matrix(i,:)
      end do
c
      do i = 1, 6
         do j = 1, 6
            Sigma(i,kintk) = Sigma(i,kintk) 
     1                       +C_matrix(i,j)*strain(j,kintk)
c
c~             U_psilon(i,j) = C_matrix(i,j)
         end do
      end do
      write(6,*) "this is Sigma UMAT" 
      do i = 1, size(Sigma,1)
        write(6,'(40G12.4)')  Sigma(i,:)
      end do
c
      do i = 1, 18
         do j = 1, 18
            Tau(i,kintk) = Tau(i,kintk) 
     1                     +D_matrix(i,j)*relaxedstraingradient(j,kintk)
c~             Lambda(i,j) = D_matrix(i,j)
         end do
      end do
      write(6,*) "this is Tau UMAt" 
      do i = 1, size(Tau,1)
        write(6,'(40G12.4)')  Tau(i,:)
      end do
      return 
      end subroutine KUMAT
c***********************************************************************
c***********************************************************************
      subroutine KDLOAD(rhs_k,NDLOAD,coords,mdload,ADLMAG,jdltyp)
c      
      double precision :: XI,Weight,Knode,KAssembly,BF_gam,KTYP,Bxi,
     1  Bweight,BJ_S,BJ_X,BJ_Y,BN_gam,BN_gam_xi,KNODE_S,KRHS,
     2  ADLMAG,K,rhs_k
c      
      dimension XI(2),Weight(2),Knode(3,4),KAssembly(3,8),BF_gam(3),
     1  BN_gam(3,4),BN_gam_xi(3,4),coords(2,9),
     2  ADLMAG(mdload,1),jdltyp(mdload,1),rhs_k(16)
c      
        XI(1) = -(1.d0/3.d0)**0.5
        XI(2) =  (1.d0/3.d0)**0.5
c      
*        write(*,*) "rhs_inp"
*        write(*,*) rhs
        Weight(1) = 1
        Weight(2) = 1
c
        Knode(1,1) = 1
        Knode(2,1) = 2
        Knode(3,1) = 5
        Knode(1,2) = 2
        Knode(2,2) = 3
        Knode(3,2) = 6
        Knode(1,3) = 3
        Knode(2,3) = 4
        Knode(3,3) = 7
        Knode(1,4) = 1
        Knode(2,4) = 4
        Knode(3,4) = 8
c      
        KAssembly(1,1) = 1
        KAssembly(2,1) = 3
        KAssembly(3,1) = 9
        KAssembly(1,2) = 3
        KAssembly(2,2) = 5
        KAssembly(3,2) = 11
        KAssembly(1,3) = 5
        KAssembly(2,3) = 7
        KAssembly(3,3) = 13
        KAssembly(1,4) = 7
        KAssembly(2,4) = 1
        KAssembly(3,4) = 15
        KAssembly(1,5) = 2
        KAssembly(2,5) = 4
        KAssembly(3,5) = 10
        KAssembly(1,6) = 4
        KAssembly(2,6) = 6
        KAssembly(3,6) = 12
        KAssembly(1,7) = 6
        KAssembly(2,7) = 8
        KAssembly(3,7) = 14
        KAssembly(1,8) = 8
        KAssembly(2,8) = 2
        KAssembly(3,8) = 16
c
c
c
c       do KLOAD = 1, NDLOAD
        KLOAD = NDLOAD 
c
        do i = 1, 3
              BF_gam(i) = 0.d0
        end do  
        JDL = JDLTYP(1,1)
c
        if (JDL .GE. 5) then
          KTYP = JDL - 4
        end if
c
        KTYP = JDL
c
        do kgauss = 1, 2
c          
          Bxi = XI(kgauss)
          Bweight = Weight(kgauss)
c
          BJ_S = 0.d0
          BJ_X = 0.d0
          BJ_Y = 0.d0
c            
          BN_gam(1,1) = 0.5d0*Bxi*(Bxi-1.d0)
          BN_gam(2,1) = 0.5d0*Bxi*(Bxi+1.d0)
          BN_gam(3,1) = 1.d0-Bxi**2
          BN_gam(1,2) = 0.5d0*Bxi*(Bxi-1.d0)
          BN_gam(2,2) = 0.5d0*Bxi*(Bxi+1.d0)
          BN_gam(3,2) = 1.d0-Bxi**2
          BN_gam(1,3) = 0.5d0*Bxi*(Bxi+1.d0)
          BN_gam(2,3) = 0.5d0*Bxi*(Bxi-1.d0)
          BN_gam(3,3) = 1.d0-Bxi**2
          BN_gam(1,4) = 0.5d0*Bxi*(Bxi+1.d0)
          BN_gam(2,4) = 0.5d0*Bxi*(Bxi-1.d0)
          BN_gam(3,4) = 1.d0-Bxi**2
c            
          BN_gam_xi(1,1) = Bxi-0.5d0            
          BN_gam_xi(2,1) = Bxi+0.5d0 
          BN_gam_xi(3,1) = -2.d0*Bxi
          BN_gam_xi(1,2) = Bxi-0.5d0
          BN_gam_xi(2,2) = Bxi+0.5d0
          BN_gam_xi(3,2) = -2.d0*Bxi
          BN_gam_xi(1,3) = Bxi+0.5d0
          BN_gam_xi(2,3) = Bxi-0.5d0
          BN_gam_xi(3,3) = -2.d0*Bxi
          BN_gam_xi(1,4) = Bxi+0.5d0
          BN_gam_xi(2,4) = Bxi-0.5d0
          BN_gam_xi(3,4) = -2.d0*Bxi
c
          do i = 1, 3       
            KNODE_S = Knode(i,KTYP)
            BJ_X = BJ_X  + BN_gam_xi(i,KTYP)*coords(1,KNODE_S)
            BJ_Y = BJ_Y  + BN_gam_xi(i,KTYP)*coords(2,KNODE_S)
            BJ_S = (BJ_X**2+BJ_Y**2)**0.5
          end do           
c
          do i = 1, 3
            BF_gam(i) = BF_gam(i) + BN_gam(i,KTYP)*ADLMAG(KLOAD,1)
     1                    *BJ_S*Bweight
          end do
c
        end do
c      
        do i = 1, 3
  	      K = JDLTYP(1,1)
	      KRHS = KAssembly(i,K)
	      rhs_k(KRHS) = rhs_k(KRHS) + BF_gam(i)
        end do
        write(*,*) "rhs_k"
        write(*,*) rhs_k
c
c         
      return 
      end subroutine KDLOAD
************************************************************************


