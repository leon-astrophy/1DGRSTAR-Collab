!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is the parameter files containing constants necessary for computing MR relations
! Note that we assumes G = Solar Mass = c = 1, all unit conversion are done accordingly
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define double precision !
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit constants !

! Mathematical constants and physical constants !
REAL (DP), PARAMETER :: pi_old = 3.1415926535897932384626433832795E0_DP

! Physical constants to be as one !
REAL (DP), PARAMETER :: gconst = 6.67430D-8
REAL (DP), PARAMETER :: clight = 2.99792458D10
REAL (DP), PARAMETER :: solar = 1.98847D33

! Solar Radius !
REAL (DP), PARAMETER :: rsolar = 6.96342D10

! 1 GeV !
REAL (DP), PARAMETER :: GeV2gram = 1.78266191D-24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unit conversions !

! Conversion between units !
REAL (DP), PARAMETER :: lencgs2code = (clight**2)/(solar*gconst)
REAL (DP), PARAMETER :: masscgs2code = (1.0D0/solar)
REAL (DP), PARAMETER :: tcgs2code = (clight**3)/(solar*gconst)

! Derived conversion !
REAL (DP), PARAMETER :: rhocgs2code = (masscgs2code/lencgs2code**3)
REAL (DP), PARAMETER :: taucgs2code = (masscgs2code*lencgs2code**2/tcgs2code**2)
REAL (DP), PARAMETER :: h_bar = (1.054571817D-27)*(lencgs2code**2*masscgs2code/tcgs2code)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Other unit conversions !

! Conversion to cubic meter !
REAL (DP), PARAMETER :: fm3tocm3 = 1.0D-39

! MeV to erg !
REAL (DP), PARAMETER :: mev2dyncm = 1.6021766339999D-6

! MeV to gram !
REAL (DP), PARAMETER :: mev2gram = 1.78266269594644D-27

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This section is related to the basic parameters governing the simulation !

! Physical length (code unit) of the simulation box !
REAL (DP), PARAMETER :: total_length = 2.0E4_DP

! Value of spatial grid size dx !
REAL (DP), PARAMETER :: dx = 1.0E-1_DP

! The total number of array stored by each variables !
INTEGER, PARAMETER :: length_step = INT (total_length / dx)      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for NM

! Atmospeheric density factor !
REAL (DP), PARAMETER :: rhofac = 1.0D-10

! Number of models per maximum density !
INTEGER, PARAMETER :: n_rho = 101

! Starting log maximum density !
REAL (DP), PARAMETER :: rhostart = 13.0D0

! Ending log maximum density !
REAL (DP), PARAMETER :: rhoend = 15.5D0

! Step size !
REAL (DP), PARAMETER :: drho = (rhoend - rhostart)/(DBLE(n_rho) - 1.0D0)

! Want tabular EOS ? !
INTEGER, PARAMETER :: nmtable_flag = 1

! EOS Table NAME !
character*256 :: nmtable_name = "nuclear_matter_EOS_1.in"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for DM parameters

! Want DM? !
INTEGER, PARAMETER :: dm_flag = 1

! Maximum bisection step !
INTEGER, PARAMETER :: n_max = 1000

! Tolerance !
REAL (DP), PARAMETER :: tor = 1.0D-6

! DM ratio !
REAL (DP), PARAMETER :: dm_ratio = 0.05D0

! Want tabular EOS ? !
INTEGER, PARAMETER :: dmtable_flag = 1

! EOS Table NAME !
character*256 :: dmtable_name = "DM_EOS_3.in"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This section governs the physics of the EOS of NM, ideal degenerate gas

! Baryonic mass for normal matter !			
REAL (DP), PARAMETER :: mb2 = 1.66053906660D-24*masscgs2code

! Fermionic mass (electrons) for normal matter !			
REAL (DP), PARAMETER :: me2 = 9.1093837015D-28*masscgs2code

! Electron fraction for normal matter !
REAL (DP), PARAMETER :: ye2 = 5.0E-1_DP

! Multiplicity factor for normal matter				
REAL (DP), PARAMETER :: gs2 = 2.0E0_DP

! For fermi gas EOS !
REAL (DP), PARAMETER :: a_max2 = (me2**4)/(2.4D1*pi_old**2*h_bar**3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This section governs the physics of the EOS of DM, ideal degenerate gas

! Electron fraction for normal matter !
REAL (DP), PARAMETER :: ye1 = 1.0E0_DP

! Multiplicity factor for normal matter				
REAL (DP), PARAMETER :: gs1 = 2.0E0_DP	

! Baryonic mass for dark matter !			
REAL (DP), PARAMETER :: mb1 = 0.1D0*GeV2gram*masscgs2code

! Fermionic mass (electrons) for dark matter !			
REAL (DP), PARAMETER :: me1 = 0.1D0*GeV2gram*masscgs2code

! For fermi gas EOS !
REAL (DP), PARAMETER :: a_max1 = (me1**4)/(2.4D1*pi_old**2*h_bar**3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!