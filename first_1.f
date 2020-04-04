c~ This is Ueser Element Subroutine(UEL) for Strain Gradient Elasticity 
c~ Theory(SGET). The main purpose to write this UEL is to implement the 
c~ 2D-FEM Model based on the SGET. This Code is based on the Mixed FEM
c~ formulation derived by the John Y. Shu. The main thing in this UEl is 
c~ to implement the Strain Gradient term into element formulation. 
c~ Therefor, the newly developed QU34L4 element is used.
c~ I code this UEL from my own I dont have any sorurce code reagarding 
c~ theory.In addition to get extra information about this theory and for
c~ basic knowledge of this mixed fem formulation, the diplom thesis of
c~ L.Zybell,which is also refer in Report of this PPP,is provided to me.
c~ I dont have any prior knowledge about the UELcoding so I refer abaqus 
c~ manual for User element routine and also user material routine.
c    
c~ Code written by    : Dhaval Rasheshkumar Patel
c~ Institute          : Technical University Freiberg
c~ Supervisor         : Dr. Sergii Kozinov
c~ Date of Submission : 4th April 2020
c~ Title of PPP       : Implementation of GradientElasticityModel In FEM
c
c***********************************************************************
      subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,
     1 nprops,coords,mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,
     2 jelem,params,ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,  
     3 ddlmag,mdload,pnewdt,jprops,njpro,period)
c
c~      Main Output from UEl : Rhs(right hand side vecor: F_ext - F_int)
c~                           : Amatrx(Stiffness Matrix of Element)
c~                           : Svars(State variables)
                          
c~      Main input to UEl    : U (all degrees of freedom - Displacement,
c~                             Relaxed strain, Langrange_Multiplier) 
c~                           : Du (change in all DOF in every step Disp,
c~                             Relaxed strain, Langrange_Multiplier)  
c
!    Standard format of UEL for abaqus.
!    This statement includes the inbuilt variables/parameters for Abaqus
      include 'aba_param.inc' 
c
      dimension rhs(38),amatrx(38,38),props(4),svars(198),
     1 coords(2,9),u(38),du(1,34),time(2),
     2 a(ndofel),params(*),jdltyp(mdload,*),adlmag(mdload,*),
     3 ddlmag(mdload,*),predef(2,npredf,nnode),lflags(*),jprops(*),
     4 energy(*),v(ndofel)
c
!    Assign the values for some basics parameters.
      parameter (ndim=2, nnodemax=9, ninpt=9, nsvint=22)
c
c
c    ndim      ... number of spatial dimensions
c    nnodemax  ... number of nodes per element
c    ninpt     ... number of integration points
c    nsvint    ... number of state variables per integration pt     
c
c
!    Initialization and dimensioning the vectors, matrix(arrays).
c    internal force vector due to langrange multiplier
      double precision, dimension(4)     :: S_vector
c    shape function matrix related to displcement       
      double precision, dimension(9)     :: dN_U  
c    differential of relaxed strain shape function vector      
      double precision, dimension(8)     :: bmat_psi
c    dof langrange multiplier       
      double precision, dimension(4)     :: lagrange_multi        
c    shape function matrix related to relaxed strain        
      double precision, dimension(4)     :: dN_PSI
c    gauss weight vector       
      double precision, dimension(9)     :: Gauss_Weight
c    output parameter delta displcement for 9 nodes for 2 dimensions       
      double precision, dimension(18)    :: delta_displacement
c    dof relaxed strain      
      double precision, dimension(16)    :: relaxed_strain
c    internal force vector due to displcement      
      double precision, dimension(18)    :: F_vector
c          
      double precision, dimension(18)    :: bmat   
c    internal force vector due to relaxed strain         
      double precision, dimension(16)    :: R_vector 
c    Arranged B_matrix from diffrenti. vecotor of displacement shape fun           
      double precision, dimension(18)    :: bmat_u
c    DOF displacement  
      double precision, dimension(18)    :: displacement
c    output parameter change in relaxed strain in last step. 
      double precision, dimension(16)    :: delta_relaxed_strain   
c    derivation of displacement shape function with respect to xi,omega
      double precision, dimension(2,9)   :: dNd_xi
c    derivation of relaxed_strain shape function wrt xi,omega      
      double precision, dimension(2,4)   :: dNd_PSI
c    Jacobian matrix for DOF displcement      
      double precision, dimension(2,2)   :: xjaci
c    Jacobian matrix for DOf relaxed strain 
      double precision, dimension(2,2)   :: xjaci_psi
c    Four displacement gradient parameter u1,1 u1,2 u2,1 u2,2.
      double precision, dimension(4,9)   :: GradVariable
c    Three strain values for each nodes e11,e22,e12=e21      
      double precision, dimension(3,9)   :: strain
c    change in strain in current step  d-e11, d-e22, d-e12      
      double precision, dimension(3,9)   :: Deltastrain
c    four relaxed strain for gauss point psi11,psi21,psi12,psi22     
      double precision, dimension(4,9)   :: relaxedstrain
c    six relaxed strain gradient values for each guass points       
      double precision, dimension(6,9)   :: relaxedstraingradient
c    six change in relaxed strain gradient values for last step      
      double precision, dimension(6,9)   :: delta_relaxedstraingradient
c    for Langrange values at Centre node of element at each Gauusint pt.      
      double precision, dimension(4,9)   :: Langrangemulti
c    Three Stress values at each Gauss Int. point.      
      double precision, dimension(3,9)   :: Sigma
c    Six Higher Stress values at each Gauss Int. point.      
      double precision, dimension(6,9)   :: Tau
c    Shape Matrix (identtity) related to langrange multiplier DOF.       
      double precision, dimension(4,4)   :: Nrhomatrix
c    Elastic Matrix      
      double precision, dimension(3,3)   :: U_psilon
c    Higher order Elastic Matrix      
      double precision, dimension(6,6)   :: Lambda
c    Stiffness Matrix related to disp. and lagrannge DOF.       
      double precision, dimension(18,4)  :: stiffness_Urho
c    Stiffness Matrix related to relaxed strain and lagrannge DOF.
      double precision, dimension(16,4)  :: stiffness_psirho
c    Transpose of Stiffness_Urho
      double precision, dimension(4,18)  :: tra_stiffness_Urho
c    Transpose of Stiffness_psirho
      double precision, dimension(4,16)  :: tra_stiffness_psirho
c    Transpose of Nmatrix
      double precision, dimension(18,2)  :: T_Nmatrix
c    Transpose of Mmatrix 
      double precision, dimension(18,4)  :: T_Mmatrix
c    Transpose of Bmatrix
      double precision, dimension(18,3)  :: T_Bmatrix
c    Transpose of Npsimatrix 
      double precision, dimension(16,4)  :: T_Npsimatrix
c    Transpose of Bpsimatrix 
      double precision, dimension(16,6)  :: T_Bpsimatrix
c    Shape matrix related to Displacement DOF.
      double precision, dimension(2,18)  :: Nmatrix
c    Displacement Gradient Matrix
      double precision, dimension(4,18)  :: Mmatrix
c    Differential Matrix B-matrix related to Displacement DOF.
      double precision, dimension(3,18)  :: Bmatrix
c    Differential Matrix B-matrix related to relaxed strain DOF.
      double precision, dimension(6,16)  :: Bpsimatrix
c    Shape matrix related to relaxed strain DOF.
      double precision, dimension(4,16)  :: Npsimatrix
c    Stiffness matrix related to only Displacement DOF.      
      double precision, dimension(18,18) :: stiffness_UU 
c    Stiffness matrix related to only Relaxed-Strain DOF.  
      double precision, dimension(16,16) :: stiffness_PsiPsi 
c    Global Striffness Matrix of Element       
      double precision, dimension(38,38) :: amat
c    Mass Matrix      
      double precision, dimension(18,18) :: Mass_U    
c    Mass Matrix (Identity Matrix)        
      double precision, dimension(38,38) :: Mass_matrix   
c    Global Right hand side Vector (F_int- F_ext)      
      double precision, dimension(16)    :: rhs_k            
c 
      double precision :: E_modulus,Nue,micro_length,Density,
     1  lemda,mue
c    
        S_vector = 0.d0
        dN_U   = 0.d0
        bmat_psi = 0.d0   
        dN_PSI = 0.d0
        statevLocal = 0.d0
        F_vector = 0.d0
        bmat  = 0.d0
        R_vector  = 0.d0
        bmat_u = 0.d0
        dNd_xi = 0.d0
        dNd_PSI = 0.d0
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
!    Guass Weight of 9 Guass integration points in respective order
      data Gauss_Weight /0.308641975309d0, 0.308641975309d0, 
     * 0.308641975309d0, 0.308641975309d0, 0.493827160494d0,
     * 0.493827160494d0, 0.493827160494d0, 0.493827160494d0, 
     * 0.79012345679d0/
c
c    This print the lflags which shows the procedure required to solve 
c    the non-linear or linear set of equations.
c
      write(*,*) lflags(3)		
c
!    this is to print the U and DU in .dat file 
      write(6,*) "this is U"
      do i = 1, 38
        write(6,*) u(i)
      end do  
      write(6,*) "this is DU"
      do i = 1, 34
        write(6,*) du(1,i)
      end do 
c      
!    this is to print the element number in .dat file      
      write(6,*) jelem				
c      
!    this is to print the coordintes of the nodes 
!    of the respective element in .dat file      
      write(6,*) "this is a coords of nodes of current element"
      do i = 1, size(coords,1)
            write(6,'(20G12.4)') coords(i,:)
      end do
c      
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
      do i = 1, 18
         displacement(i) = u(i)
         delta_displacement(i) = du(1,i)
      end do
c
!    loop for extracting the relaxed_strain and-
!    delta_relaxed_strain value from U and DU.
      do i = 1, 16
         relaxed_strain(i) = u(i+18) 
         delta_relaxed_strain(i) = du(1,i+18)
      end do
c 
!    loop for extracting lagrangian multiplier-
!    from the U vector.
      do i = 1, 4
         lagrange_multi(i) = u(i+34) 
      end do
c
c
c-----------------------------------------------------------------------
!    loop over Gauss Integration Points
!    There are total 9 integration points so it goes from 1 to 9.
c
      do kintk = 1, 9
c
         write(6,*) "The integeration no. is"
         write(6,*) kintk
!     Evaluate shape functions and derivatives
c
         call shapefcn_U(dN_U,dNd_xi,kintk,ninpt,nnode,ndim)
!     Now we have shape function and its derivates for displacement 
!     degree of freedom. 
*         write(6,*) "shape_function_u"
*         write(6,*) dNd_xi
c
         call shapefcn_PSI(dN_PSI,dNd_PSI,kintk,ninpt,nnode,ndim)
*         write(6,*) "shape_function_PSI"
*         write(6,*) dNd_PSI
!     Now we have shape function and its derivates for relaxedstrain 
!     degree of freedom.   
c 
!       Evaluate Jacobian of displacement and psi shape function matrix
c
         call jacobian_UU(djac,xjaci,jelem,ndim,nnode,coords,dNd_xi,
     1        pnewdt)
*         write(6,*) "Jacobian of displacement"
*         write(6,*) djac         
*         write(6,*) xjaci
c                  
         call jacobian_psipsi(djac_psi,xjaci_psi,jelem,ndim,nnode,
     1        coords,dNd_PSI,pnewdt)
*         write(6,*) "Jacobian of PSI"
*         write(6,*) djac_psi
*         write(6,*) xjaci_psi     
c
!       Evaluate B_matrix for displcement and psi 
c
         call bmatrix_U(bmat_u,xjaci,dNd_xi)
c~          write(6,*) "bmat_u"         
c~          write(6,*) bmat_u
c
         call bmatrix_PSI(bmat_psi,xjaci_psi,dNd_PSI)
c~          write(6,*) "bmat_psi"         
c~          write(6,*) bmat_psi
c
c
c-----------------------------------------------------------------------
c       Assembly of Common Element Matrices :
c
!       Assembly of Nmatrix from shape function for displacement degree
!       freedom.
         do i = 1, 2*nnode
            Nmatrix(1,i) = 0.d0
            Nmatrix(2,i) = 0.d0
         end do   
c            
         do i = 1, 9
            Nmatrix(1,(2*i-1)) = dN_U(i)
            Nmatrix(2,(2*i)) = dN_U(i)
         end do
!       print the Nmatrix in .dat file            
*         write(6,*) "Nmatrix is here"
*         do i = 1, size(Nmatrix,1)
*            write(6,'(20G12.4)') Nmatrix(i,:)
*         end do
c    
c     
!       Assembly of Mmatrix from derivatives of shape functions for 
!       displacement degree freedom.    
         do i = 1, 2*nnode
            Mmatrix(1,i) = 0.d0
            Mmatrix(2,i) = 0.d0
            Mmatrix(3,i) = 0.d0
            Mmatrix(4,i) = 0.d0
         end do   
c         
         do i = 1, nnode
            Mmatrix(1,2*i-1) = bmat_u(2*i-1)
            Mmatrix(2,2*i-1)   = bmat_u(2*i)
            Mmatrix(3,2*i) = bmat_u(2*i-1)
            Mmatrix(4,2*i)   = bmat_u(2*i)
         end do  
!       print the Mmatrix in .dat file         
c~          write(6,*) 'this is mmatrix'
c~          do i = 1, size(Mmatrix,1)
c~             write(6,'(20G12.4)') Mmatrix(i,:)
c~          end do         
c
c
!       Assembly of Bmatrix from derivatives of shape functions for 
!       displacement degree freedom.    
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
!       print the Bmatrix in .dat file         
c~          write(6,*) 'this is Bmatrix'
c~          do i = 1, size(Bmatrix,1)
c~             write(6,'(20G12.4)') Bmatrix(i,:)
c~          end do
c
c
!       Assembly of Npsimatrix from shape functions for 
!       relaxed strain degree freedom. 
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
!       print the Npsimatrix in .dat file
c~          write(6,*) 'this is Npsimatrix'
c~          do i = 1, size(Npsimatrix,1)
c~             write(6,'(20G12.4)') Npsimatrix(i,:)
c~          end do
c
c
!       Assembly of Bpsimatrix from derivatives of shape functions for 
!       relaxed strain degree freedom. 
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
!       print the Bpsimatrix in .dat file         
c~          write(6,*) 'this is Bpsimatrix'
c~          do i = 1, size(Bpsimatrix,1)
c~             write(6,'(20G12.4)') Bpsimatrix(i,:)
c~          end do
c
c
!       Assembly of Nrhomatrix for lagrange_multiplier D.O.F.
         do i = 1, 4
            Nrhomatrix(1,i) = 0.d0
            Nrhomatrix(2,i) = 0.d0 
            Nrhomatrix(3,i) = 0.d0
            Nrhomatrix(4,i) = 0.d0
         end do
         do i = 1, 4
            Nrhomatrix(i,i) = 1.d0
         end do
!       print the Nrhomatrix in .dat file         
*         write(6,*) 'this is Nrhomatrix'
*         do i = 1, size(Nrhomatrix,1)
*            write(6,'(20G12.4)') Nrhomatrix(i,:)
*         end do
c         
c-----------------------------------------------------------------------
c       Calculation of all State Variables.
c
c   ------------------Displacement Gradient--------------------------
!       Displacement gradient = GradVariable = 0.d0
!       this is a transpose of Mmatrix
         T_Mmatrix = 0.d0
         T_Mmatrix = Transpose(Mmatrix)		
c         
         do i = 1, 18
          do j = 1, 4
*            write(6,*) T_Mmatrix(i,j)
*            write(6,*) displacement(i)
*            write(6,*) T_Mmatrix(i,j)*displacement(i)
            GradVariable(j,kintk) = GradVariable(j,kintk) 
     1                              +T_Mmatrix(i,j)*displacement(i)
          end do
         end do 
c         
!       Print the displcement gradient in .dat file.           
c~          write(6,*) "this is gradient variable"
c~          do i = 1, size(GradVariable,1)
c~             write(6,'(20G12.4)') GradVariable(i,:)
c~          end do
c
c   ----------------------Strain, Delta Strain-----------------------
c
!       strain = 0.d0
!       Deltastrain = 0.d0
!       This is a transpose of Bmatrix
         T_Bmatrix = 0.d0
	     T_Bmatrix = Transpose(Bmatrix)
         do i = 1, 18
          do j = 1, 3
            strain(j,kintk) = strain(j,kintk) 
     1                        + T_Bmatrix(i,j)*displacement(i)
c     
            Deltastrain(j,kintk) = Deltastrain(j,kintk) 
     1                            +T_Bmatrix(i,j)*delta_displacement(i)
          end do
         end do
c
!       Print the Strain in .dat file.
c~          write(6,*) "this is strain"
c~          do i = 1, size(strain,1)
c~             write(6,'(20G12.4)')  strain(i,:)
c~          end do
c~ !       Print the Delta_Strain in .dat file.         
c~          write(6,*) "this is Deltastrain"
c~          do i = 1, size(Deltastrain,1)
c~             write(6,'(20G12.4)')  Deltastrain(i,:)
c~          end do
c
c    --------------------Relaxed Strain----------------------------
c
!       relaxedstrain = 0.d0
!       This is a transpose of Npsimatrix
        T_Npsimatrix = 0.d0
	    T_Npsimatrix = Transpose(Npsimatrix)
         do i = 1, 16
          do j = 1, 4
            relaxedstrain(j,kintk) = relaxedstrain(j,kintk) 
     1                              +T_Npsimatrix(i,j)*relaxed_strain(i)
          end do  
         end do
c
!       Print the relaxedstrain in .dat file.
c~ 	     write(6,*) "this is relaxedstrain"
c~          do i = 1, size(relaxedstrain,1)
c~             write(6,'(20G12.4)')  relaxedstrain(i,:)
c~          end do
c
c     ------------------Relaxed Strain Gradient---------------------
c
!       relaxedstraingradient = 0.d0
!       delta_relaxedstraingradient = 0.d0
!       This is a transpose of Bpsimatrix
         T_Bpsimatrix = 0.d0
	     T_Bpsimatrix = Transpose(Bpsimatrix)
         do i = 1, 16
          do j = 1, 6
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
c~          write(6,*) "this is relaxedstraingradient"
c~          do i = 1, size(relaxedstraingradient,1)
c~             write(6,'(20G12.4)')  relaxedstraingradient(i,:)
c~          end do
!       Print the delta_relaxedstraingradient in .dat file.         
c~          write(6,*) "this is delta_relaxedstraingradient"
c~          do i = 1, size(delta_relaxedstraingradient,1)
c~             write(6,'(20G12.4)')  delta_relaxedstraingradient(i,:)
c~          end do
c
c    --------------------Langrange Multiplier----------------------
c         
         !Langrangemulti = 0.d0
         do i = 1, 4
          do j = 1, 4
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
c
c-----------------------------------------------------------------------
!       call the KUMAT to find the state variables Sigma and Tau
c
!       And also get U_psilon(Elastic Matrix) and Lambda(Higehr order
!       -elastic matrix) from KUMAT as output.
c
         call KUMAT(Sigma,Tau,U_psilon,Lambda,kintk,lemda,mue,
     1   micro_length,strain,relaxedstraingradient)
c     
!       Print the State variable Sigma in .dat file.     
*         write(6,*) "this is Sigma"
*         write(6,*) Sigma(:,1)
*         do i = 1, size(Sigma,kintk)
*            write(6,'(20G12.4)')  Sigma(i,:)
*         end do
c
!       Print the State variable Tau in .dat file.
*         write(6,*) "this is Tau"
*         write(6,*) Tau(:,kintk)
*         do i = 1, size(Tau,1)
*            write(6,'(20G12.4)')  Tau(i,:)
*         end do
c
!       Print the Elastic Matrix in .dat file.
*          write(6,*) "this is Upsi_matrix"
*          do i = 1, size(U_psilon,1)
*            write(6,'(20G12.4)')  U_psilon(i,:)
*          end do
c
!       Print the Higher order elastic matrix in .dat file.
*          write(6,*) "this is Lambda"
*          do i = 1, size(Lambda,1)
*            write(6,'(20G12.4)')  Lambda(i,:)
*          end do
c
c-----------------------------------------------------------------------
c
!       store state variables Sigma and Tau in Svars vector. :
c
!       store the 3 stress values
         do i = 1, 3
            svars(i+(kintk-1)*22) = Sigma(i,kintk)
         end do
c  
!       store the 6 higher order stress values       
         do i = 4, 9
            svars(i+(kintk-1)*22) = Tau(i-3,kintk)
         end do
c
!       store the 4 relaxed strain values
         do i = 10, 13
            svars(i+(kintk-1)*22) = relaxedstrain(i-9,kintk)
         end do
c
!       store the 6 Strain Gradient Values
         do i = 14, 19
            svars(i+(kintk-1)*22) = relaxedstraingradient(i-13,kintk)
         end do
c
!       store the 3 normal strain values
         do i = 20, 22
            svars(i+(kintk-1)*22) = strain(i-19,kintk)
         end do
c
!       Print the State variable vector in .dat file.         
c~          write(6,*) "this is state vars"
c~          do i = 1, 22
c~             write(6,'(20G12.4)')  svars(i+(kintk-1)*22)
c~          end do   
c      
c-----------------------------------------------------------------------
!       Calculation of the Element Stiffness Matrices:
c
c
!       Find the UU_stiffness matrix - dimension(18*18)
         do i = 1, 18
          do j = 1, 18
           do k = 1, 3
            do l = 1, 3
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
         do i = 1, 16
          do j = 1, 16
           do k = 1, 6
             do l = 1, 6
               stiffness_PsiPsi(i,j) = stiffness_PsiPsi(i,j) 
     1                  +T_Bpsimatrix(i,k)*Lambda(k,l)*T_Bpsimatrix(j,l)
     1                  *djac_psi*Gauss_weight(kintk)
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
         do i = 1, 18
          do j = 1, 4
           do k = 1, 4
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
         do i = 1, 16
          do j = 1, 4
           do k = 1, 4
               stiffness_psirho(i,j) = stiffness_psirho(i,j) 
     1                       +T_Npsimatrix(i,k)*Nrhomatrix(j,k)*djac_psi
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
c-----------------------------------------------------------------------
c
!       Calculation of the internal force vectors : 
c
!       Find the F_vector which includes the internal forces due to
!       the displcement of all nodes.
         do i = 1, 18 
            do j = 1, 3
               F_vector(i) = F_vector(i) + T_Bmatrix(i,j)*Sigma(j,kintk)
     1                                     *djac*Gauss_weight(kintk)
*               write(6,*) T_Bmatrix(i,j)*Sigma(j,kintk)
*     1                                     *djac*Gauss_weight(kintk)
            end do
            do j = 1, 4
               F_vector(i) = F_vector(i) - T_Mmatrix(i,j) 
     1                 *Langrangemulti(j,kintk)*djac*Gauss_weight(kintk)  
            end do
         end do
!       Print the F_vactor in .dat file.         
c~          write(6,*) "F_vector"
c~          do i = 1, 18
c~             write(6,'(20G12.4)')  F_vector(i)
c~          end do
c
c
!       Find the R_vector which includes the internal forces due to
!       the relaxed strain of nodes(1,2,3,4).
         do i = 1, 16
            do j = 1, 6
              R_vector(i) = R_vector(i) + T_Bpsimatrix(i,j)*Tau(j,kintk)
     1                                     *djac_psi*Gauss_weight(kintk)
            end do
            do j = 1, 4
               R_vector(i) = R_vector(i) + T_Npsimatrix(i,j)
     1             *Langrangemulti(j,kintk)*djac_psi*Gauss_weight(kintk)
            end do
         end do
!       Print the R_vactor in .dat file.         
c~          write(6,*) "R_vector"
c~          do i = 1, 16
c~             write(6,'(20G12.4)')  R_vector(i)
c~          end do
c
c
!       Find the S_vector which includes the internal forces due to
!       the langangen multiplier of centre(last) node.
         do i = 1, 4
            do j = 1, 4
               S_vector(i) = S_vector(i) + Nrhomatrix(i,j)
     1           *(relaxedstrain(j,kintk)*djac_psi-GradVariable(j,kintk)
     1           *djac)*Gauss_weight(kintk)
            end do
         end do
!       Print the S_vactor in .dat file.         
c~          write(6,*) "S_vector"
c~          do i = 1, 4
c~             write(6,'(20G12.4)')  S_vector(i)
c~          end do
c
      end do 
!       loop over integration points is end.
c-----------------------------------------------------------------------
!        Assembaly of AMATRX(Global Stiffness Matrix) and 
!        RHS(Right hand side vector F_int-F_ext) acorrding to the 
!        arragement of degrees of freedom of the 
!        all 9 nodes of the element.
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
c  ------------Assembly of Global Stiffness Matrix-----------------
c           ___                             ____
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
c
c
!     assign the 0 value to the diagonal elements(35,36,37,38) of amat.
      do i=35, 38
	     amatrx(i,i) = 1.0e-10
      end do
c
c
c      
!     print the amatrx into the .dat file.      
      write(6,*) "amatrx"
      do i = 1, size(amatrx,1)
          write(6,'(40G12.4)')  amatrx(i,:)
      end do
c
c
c
c   --------------Assembly of Global RHS-------------------------
c
c                ___     ____
c               |            |
c               |   F_int    |
c       RHS =   |   R_int    |
c               |   S_int    |
c               |___     ____|
c
c
c  
!     Assemble the RHS matrix.
!     rhs(1:18) = F_vector(1:18)
      do i = 1, 18
         rhs(i) = rhs(i) + F_vector(i)
      end do
c
c
!     rhs(19:34) = R_vector(1:16)
      do i = 19, 34
         rhs(i) = rhs(i) + R_vector(i-18)
      end do
c
c
!     rhs(35:38) = S_vector(1:4)
      do i = 35, 38
         rhs(i) = rhs(i) + S_vector(i-34)
      end do
c
c
c     ---This is to add the external force to the internal force.---
c
      if (NDLOAD .ge. 1) then
          write(*,*) "call KDLOAD"
          call KDLOAD(rhs_k,NDLOAD,coords,mdload,ADLMAG,jdltyp)
          do i = 1, 16
            rhs(i) = rhs_k(i) - rhs(i) 
          end do          
      end if
c
c ---------------------------------------------------------------------
c
!     print the RHS matrix into the .dat file.      
      write(6,*) "rhs"
      do i = 1, 38
          write(6,'(20G12.4)')  rhs(i)
      end do
c
c
      return
      end subroutine UEL       
!     END of the main subroutine UEL.   
c***********************************************************************
c
c***********************************************************************
      subroutine shapefcn_U(dN_U,dNd_xi,kintk,ninpt,nnode,ndim)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.774596669241483d0)
      dimension dN_U(nnode),dNd_xi(ndim,nnode),coord28(2,9)
      double precision :: dN_U,dNd_xi
      dN_U = 0.d0
      dNd_xi = 0.d0
c  
c    -------- Coordinates of the Gauss Points in Unit domain --------   
      data  coord28 /-1.d0, -1.d0,
     2                1.d0, -1.d0,
     3                1.d0,  1.d0,
     4               -1.d0,  1.d0,
     5                0.d0, -1.d0,
     6                1.d0,  0.d0,
     7                0.d0,  1.d0,
     8               -1.d0,  0.d0,      
     9                0.d0,  0.d0/
c     
c  2D 9-nodes
c
c     Determine the variable Xi and Omega :
        xi = coord28(1,kintk)*gaussCoord
        omega = coord28(2,kintk)*gaussCoord
c
c     Shape Functions (Quadratic - 3 points in each dirction) :
c
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
c     derivative d(N_U)/d(Xi) :
c
        dNd_xi(1,1) = -0.5d0*xi*omega-0.25d0*omega*omega+0.25d0*omega
     1   +0.5d0*xi*omega*omega
        dNd_xi(1,2) = -0.5d0*xi*omega+0.5d0*xi*omega*omega
     1   +0.25d0*omega*omega-0.25d0*omega
        dNd_xi(1,3) = 0.5d0*xi*omega*omega+0.25d0*omega*omega
     1   +0.25d0*omega+0.5d0*xi*omega
        dNd_xi(1,4) = -0.25d0*omega*omega-0.25d0*omega
     1   +0.5d0*xi*omega*omega+0.5d0*xi*omega
        dNd_xi(1,5) = -xi*omega*omega+xi*omega
        dNd_xi(1,6) = -xi*omega*omega+0.5d0-0.5d0*omega*omega+xi
        dNd_xi(1,7) = -xi*omega*omega-xi*omega
        dNd_xi(1,8) = -0.5d0+0.5d0*omega*omega+xi-xi*omega*omega
        dNd_xi(1,9) = -2.d0*xi+2.d0*xi*omega*omega
c
c     derivative d(N_U)/d(Omega) :
c
        dNd_xi(2,1) = -0.25d0*xi*xi-0.5d0*xi*omega+0.25d0*xi
     1   +0.5d0*xi*xi*omega
        dNd_xi(2,2) = -0.25d0*xi*xi+0.5d0*xi*xi*omega
     1   +0.5d0*xi*omega-0.25d0*xi
        dNd_xi(2,3) = 0.5d0*xi*xi*omega+0.5d0*xi*omega+0.25d0*xi
     1  +0.25d0*xi*xi
        dNd_xi(2,4) = 0.5d0*xi*xi*omega-0.5d0*xi*omega-0.25d0*xi
     1  +0.25d0*xi*xi
        dNd_xi(2,5) = 0.5d0*xi*xi-0.5d0-xi*xi*omega+omega 
        dNd_xi(2,6) = -xi*xi*omega-xi*omega
        dNd_xi(2,7) = -0.5d0*xi*xi+omega-xi*xi*omega+0.5d0
        dNd_xi(2,8) = xi*omega-xi*xi*omega
        dNd_xi(2,9) = -2.d0*omega+2.d0*xi*xi*omega     
c        
      return
      end subroutine shapefcn_U
c***********************************************************************
c
c***********************************************************************
      subroutine shapefcn_PSI(dN_PSI,dNd_PSI,kintk,ninpt,nnode,ndim)
c
      include 'aba_param.inc'
c
      parameter (gaussCoord=0.774596669241483d0)
      dimension dN_PSI(nnode),dNd_PSI(ndim,nnode),coord28(2,9)
      double precision :: dN_PSI,dNd_PSI
      dN_PSI = 0.d0
      dNd_PSI = 0.d0
c
c    -------- Coordinates of the Gauss Points in Unit domain --------      
      data  coord28 /-1.d0, -1.d0,
     2                1.d0, -1.d0,
     3                1.d0,  1.d0,
     4               -1.d0,  1.d0,
     5                0.d0, -1.d0,
     6                1.d0,  0.d0,
     7                0.d0,  1.d0,
     8               -1.d0,  0.d0,      
     9                0.d0,  0.d0/
C     
C  2D 9-nodes
c
c     determine (xi,omega)
         xi = coord28(1,kintk)*gaussCoord
         omega = coord28(2,kintk)*gaussCoord
c         
c     Shape Functions (bi-linear - 2 points in each direction)
c
         dN_PSI(1) = 0.25d0*(1.d0-xi)*(1.d0-omega)
         dN_PSI(2) = 0.25d0*(1.d0+xi)*(1.d0-omega)
         dN_PSI(3) = 0.25d0*(1.d0+xi)*(1.d0+omega)
         dN_PSI(4) = 0.25d0*(1.d0-xi)*(1.d0+omega)
c
c     derivative d(N_U)/d(xi)
c
         dNd_PSI(1,1) = 0.25d0*(-1.d0+omega)
         dNd_PSI(1,2) = 0.25d0*(1.d0-omega)
         dNd_PSI(1,3) = 0.25d0*(1.d0+omega)
         dNd_PSI(1,4) = 0.25d0*(-1.d0-omega)         
c
c     derivative d(N_U)/d(omega)
c
         dNd_PSI(2,1) = 0.25d0*(-1.d0+xi)
         dNd_PSI(2,2) = 0.25d0*(-1.d0-xi) 
         dNd_PSI(2,3) = 0.25d0*(1.d0+xi) 
         dNd_PSI(2,4) = 0.25d0*(1.d0-xi)     
c
      return
      end subroutine shapefcn_PSI
c***********************************************************************
c
c***********************************************************************
      subroutine jacobian_UU(djac,xjaci,jelem,ndim,nnode,coords,dNd_xi,
     1        pnewdt)
c
c     Notation : xjac  - Jacobian Matrix 
c              : djac  - Jacobian determinant 
c              : xjaci - Inverse of Jacobian matrix 
c                          
      include 'aba_param.inc'
      dimension xjac(2,2),xjaci(2,2),coords(2,9),dNd_xi(2,9)
      double precision :: xjac,xjaci,dNd_xi,djac,coords
      djac = 0.d0
c
c     -------Initialization of xjac and xjaci---------
c
      do i = 1, ndim
        do j = 1, ndim
          xjac(i,j)  = 0.d0
          xjaci(i,j) = 0.d0
        end do
      end do
c
c    --------Define Jacobian Matrix xjac----------
c
      do inod= 1, 9
         do idim = 1, 2
           do jdim = 1, 2
              xjac(jdim,idim) = xjac(jdim,idim)
     1                              +dNd_xi(jdim,inod)*coords(idim,inod)
*              write(*,*) xjac(jdim,idim)     
           end do
         end do 
      end do
*      write(*,*) "this is xjac"
*      do i = 1, size(xjac,1)
*            write(*,'(20G12.4)') xjac(i,:)
*      end do      
c
c    --------Determinant of Jacobian Matrix-----------
c
      djac = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
*      write(*,*) "det of jacobian"
*      write(*,*) djac
c
c    -------Define Inverse of Jacobian matrix xjaci---------
c
       if (djac .gt. 0.d0) then
         ! jacobian is positive - o.k.
         xjaci(1,1) =  xjac(2,2)/djac
         xjaci(2,2) =  xjac(1,1)/djac
         xjaci(1,2) = -xjac(1,2)/djac
         xjaci(2,1) = -xjac(2,1)/djac
*         do i = 1, size(xjaci,1)
*            write(*,'(20G12.4)') xjaci(i,:)
*         end do           
       endif
      return
      end subroutine jacobian_UU
c***********************************************************************
c
c***********************************************************************
      subroutine jacobian_psipsi(djac_psi,xjaci_psi,jelem,ndim,nnode,
     1        coords,dNd_PSI,pnewdt)
c
c     Notation : xjac_psi  - Jacobian Matrix 
c              : djac_psi  - Jacobian determinant 
c              : xjaci_psi - Inverse of Jacobian matrix 
c                         
      include 'aba_param.inc'
c
      dimension xjac_psi(2,2),xjaci_psi(ndim,2),coords(2,9),
     1 dNd_PSI(ndim,4)
      double precision :: xjac_psi,xjaci_psi,dNd_PSI,djac_psi,coords
      djac_psi = 0.d0
c
c     -------Initialization of xjac_psi and xjaci_psi---------
c
      do i = 1, ndim
        do j = 1, ndim
          xjac_psi(i,j)  = 0.d0
          xjaci_psi(i,j) = 0.d0
        end do
      end do
c
c    --------Define Jacobian Matrix xjac_psi----------
c
      do inod= 1, 4
         do idim = 1, ndim
           do jdim = 1, ndim
             xjac_psi(jdim,idim) = xjac_psi(jdim,idim) 
     1                             +dNd_PSI(jdim,inod)*coords(idim,inod)                   
           end do
         end do 
      end do
c
c    --------Determinant of Jacobian Matrix-----------
c     
       djac_psi = xjac_psi(1,1)*xjac_psi(2,2) 
     1           -xjac_psi(1,2)*xjac_psi(2,1)
*      write(*,*) djac_psi
c
c    -------Define Inverse of Jacobian matrix xjaci_psi---------
c
       if (djac_psi .gt. 0.d0) then
         ! jacobian is positive - o.k.
         xjaci_psi(1,1) =  xjac_psi(2,2)/djac_psi
         xjaci_psi(2,2) =  xjac_psi(1,1)/djac_psi
         xjaci_psi(1,2) = -xjac_psi(1,2)/djac_psi
         xjaci_psi(2,1) = -xjac_psi(2,1)/djac_psi
       endif
c       
      return
      end subroutine jacobian_psipsi
c***********************************************************************
c
c***********************************************************************
      subroutine bmatrix_U(bmat_u,xjaci,dNd_xi)
c
c     Notation: bmat_u(i) ....dN_U1/dx, dN_U1/dy, dN_U2/dx, dN_U2/dy...
c      
         include 'aba_param.inc'
         dimension bmat_u(18),dNd_xi(2,9),xjaci(2,2)
         double precision :: bmat_u,dNd_xi,xjaci
c  
c  ----------Initialization of the Bmatrix---------------
c  
         do i = 1, nnode*ndim
            bmat_u(i) = 0.d0
         end do
c
c  Find the bmat_u which contains the derivatives of shape functions
c  with respect to cartesian coordinates x and y.
c    
         do i = 1, 9
          do j = 1, 2
           do k = 1, 2
            m = k + (i - 1)*2
c            
            bmat_u(m)=bmat_u(m)+xjaci(k,j)*dNd_xi(j,i)
           end do
          end do
         end do 
c    
      return
      end subroutine bmatrix_U
c    
c***********************************************************************
c
c***********************************************************************
      subroutine bmatrix_PSI(bmat_psi,xjaci_psi,dNd_PSI)
c
c     Notation: bmat_psi(i) ....dN1/dx, dN1/dy, dN2/dx, dN2/dy...      
         include 'aba_param.inc'
         dimension bmat_psi(8),dNd_PSI(2,4),xjaci_psi(2,2)
         double precision :: bmat_psi,dNd_PSI,xjaci_psi
c  
c  ----------Initialization of the Bmatrix---------------
c   
         do i = 1, 4*ndim
               bmat_psi(i) = 0.d0
         end do
c           
c  Find the bmat_psi which contains the derivatives of shape functions
c  with respect to cartesian coordinates x and y.
c   
         do i = 1, 4
          do j = 1, 2
           do k = 1, 2
            n = k + (i - 1)*2
c            
            bmat_psi(n)=bmat_psi(n)+xjaci_psi(k,j)*dNd_PSI(j,i)
           end do
          end do
         end do 
c    
      return
      end subroutine bmatrix_PSI
c    
c***********************************************************************
c
c***********************************************************************
      subroutine KUMAT(Sigma,Tau,C_matrix,D_matrix,kintk,lemda,mue,
     1 micro_length,strain,relaxedstraingradient)
c
      double precision :: Sigma,Tau,C_matrix,D_matrix,mue,micro_length,
     1 strain,lemda,relaxedstraingradient
c
      dimension C_matrix(3,3), D_matrix(6,6), sigma(3,9), Tau(6,9),
     1 Lambda(6,6), U_psilon(3,3), strain(3,9),
     2 relaxedstraingradient(6,9)
c    
c   ------------Intialization of C_matrix and D_matrix--------------- 
c 
      C_matrix = 0.d0
      D_matrix = 0.d0
c
c
*      write(6,*) "this is lemda and mue"	
*      write(6,*) lemda
*      write(6,*) mue
*      write(6,*) "micro"
*      write(6,*) micro_length
c
c   ----------------Elastic Matrix-----------------------
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
      write(6,*) "this is C_matrix"
      do i = 1, size(C_matrix,1)
        write(6,'(20G12.4)')  C_matrix(i,:)
      end do
c
c    ------------Higher order Elastic Matrix----------------
c
c    --------D_matrix for General Strain Gradient solid-----------
c   Introduced by : Amanatidou and Aravas
c
      D_matrix(1,1) = (lemda + 2.d0*mue)*micro_length*micro_length
      D_matrix(1,6) = 0.5d0*lemda*micro_length*micro_length
      D_matrix(2,2) = mue*micro_length*micro_length
      D_matrix(2,6) = 0.5d0*mue*micro_length*micro_length
      D_matrix(3,3) = mue*micro_length*micro_length
      D_matrix(3,5) = 0.5d0*mue*micro_length*micro_length 
      D_matrix(4,4) = (lemda + 2.d0*mue)*micro_length*micro_length
      D_matrix(4,5) = 0.5d0*lemda*micro_length*micro_length
      D_matrix(5,3) = 0.5d0*mue*micro_length*micro_length
      D_matrix(5,4) = 0.5d0*lemda*micro_length*micro_length
      D_matrix(5,5) =0.25d0*(lemda + 3.d0*mue)*micro_length*micro_length
      D_matrix(6,1) = 0.5d0*lemda*micro_length*micro_length
      D_matrix(6,2) = 0.5d0*mue*micro_length*micro_length
      D_matrix(6,6) =0.25d0*(lemda + 3.d0*mue)*micro_length*micro_length
c
c    --------D_matrix for Couple Stress Solid--------------------
c   Introduced by : John Y. shu
c
c~          D_matrix(2,2) = 1.d0*micro_length*micro_length*mue
c~          D_matrix(2,6) = -0.5d0*micro_length*micro_length*mue
c~          D_matrix(3,3) = 1.d0*micro_length*micro_length*mue
c~          D_matrix(3,5) = -0.5d0*micro_length*micro_length*mue
c~          D_matrix(5,3) = -0.5d0*micro_length*micro_length*mue
c~          D_matrix(5,5) = 0.25d0*micro_length*micro_length*mue
c~          D_matrix(6,2) = -0.5d0*micro_length*micro_length*mue
c~          D_matrix(6,6) = 0.25d0*micro_length*micro_length*mue
c
c
      write(6,*) "this is D_matrix"
      do i = 1, size(D_matrix,1)
        write(6,'(20G12.4)')  D_matrix(i,:)
      end do
c
c
c   Determine another state variable Stress(Sigma) :
c
      do i = 1, 3
         do j = 1, 3
            Sigma(i,kintk) = Sigma(i,kintk) 
     1      +C_matrix(i,j)*strain(j,kintk)
c
            write(6,*) C_matrix(i,j)
            write(6,*) strain(j,kintk)
            write(6,*) C_matrix(i,j)*strain(j,kintk)
!           U_psilon(i,j) = C_matrix(i,j)
         end do
      end do
c
!       U_psilon = C_matrix
c
      write(6,*) "this is Sigma" 
      do i = 1, size(Sigma,1)
        write(6,'(20G12.4)')  Sigma(i,:)
      end do
c
c   Determine another state variable Higher order Stress(Tau) :
c
      do i = 1, 6
         do j = 1, 6
            Tau(i,kintk) = Tau(i,kintk) 
     1       +D_matrix(i,j)*relaxedstraingradient(j,kintk)
c     
            write(6,*) D_matrix(i,j)
            write(6,*) relaxedstraingradient(j,kintk)
            write(6,*) D_matrix(i,j)*relaxedstraingradient(j,kintk)
!           Lambda(i,j) = D_matrix(i,j)
         end do
      end do
c  
!       Lambda = D_matrix
c    
      write(6,*) "this is Tau" 
      do i = 1, size(Tau,1)
        write(6,'(20G12.4)')  Tau(i,:)
      end do
c
      return 
      end subroutine KUMAT
c***********************************************************************
c
c***********************************************************************
      subroutine KDLOAD(rhs_k,NDLOAD,coords,mdload,ADLMAG,jdltyp)
c    
c   Subroutine to add the externally applied load to the internal load
c   because in UEL concetrated and distributed load can not directly
c   added to the internal load because of different boundary integrals. 
c  
      double precision :: XI,Weight,Knode,KAssembly,BF_gam,KTYP,Bxi,
     1  Bweight,BJ_S,BJ_X,BJ_Y,BN_gam,BN_gam_xi,KNODE_S,KRHS,
     2  ADLMAG,K,rhs_k
c      
      dimension XI(2),Weight(2),Knode(3,4),KAssembly(3,8),BF_gam(3),
     1  BN_gam(3,4),BN_gam_xi(3,4),coords(2,9),
     2  ADLMAG(mdload,1),jdltyp(mdload,1),rhs_k(16)
c    
c  ----------Two integration points(second order)----------------
c  
        XI(1) = -(1.d0/3.d0)**0.5
        XI(2) =  (1.d0/3.d0)**0.5
c      
c  ----------Gauss weight of two integration points----------
c
        Weight(1) = 1
        Weight(2) = 1
c
c  ----------Nodes on the each edge of the element  
c
c  nodes on bottom edge :
        Knode(1,1) = 1
        Knode(2,1) = 2
        Knode(3,1) = 5
c  nodes on right edge :        
        Knode(1,2) = 2
        Knode(2,2) = 3
        Knode(3,2) = 6
c  nodes on Top edge :        
        Knode(1,3) = 3
        Knode(2,3) = 4
        Knode(3,3) = 7
c  nodes on left edge :        
        Knode(1,4) = 1
        Knode(2,4) = 4
        Knode(3,4) = 8
c      
c
c  Contribution of this all edge nodes into the external force vectors
c  are decided.
c
c  Bottom edge nodes contribution to the x-direction external forces.
        KAssembly(1,1) = 1
        KAssembly(2,1) = 3
        KAssembly(3,1) = 9
c  Right edge nodes contribution to the x-direction external forces.
        KAssembly(1,2) = 3
        KAssembly(2,2) = 5
        KAssembly(3,2) = 11
c  Top edge nodes contribution to the x-direction external forces.        
        KAssembly(1,3) = 5
        KAssembly(2,3) = 7
        KAssembly(3,3) = 13
c  Left edge nodes contribution to the x-direction external forces.        
        KAssembly(1,4) = 7
        KAssembly(2,4) = 1
        KAssembly(3,4) = 15
c  Bottom edge nodes contribution to the y-direction external forces.        
        KAssembly(1,5) = 2
        KAssembly(2,5) = 4
        KAssembly(3,5) = 10
c  Right edge nodes contribution to the y-direction external forces.
        KAssembly(1,6) = 4
        KAssembly(2,6) = 6
        KAssembly(3,6) = 12
c  Top edge nodes contribution to the y-direction external forces.
        KAssembly(1,7) = 6
        KAssembly(2,7) = 8
        KAssembly(3,7) = 14
c  Left edge nodes contribution to the y-direction external forces.
        KAssembly(1,8) = 8
        KAssembly(2,8) = 2
        KAssembly(3,8) = 16
c
c
c
c       do KLOAD = 1, NDLOAD ---for multiple external load 
        KLOAD = NDLOAD ! for single external load
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
c   -----Loop over integration points
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
c   Shape functions for each surface :
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
c   Derivatives of shape function with variable Bxi. :
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
c    end of Loop over integration points
c
        end do
c   
c    Loop to add the external load to the internal load :
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
c***********************************************************************


