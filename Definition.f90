!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the definition module contains the key variables used in the hydro code !				  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE DEFINITION
IMPLICIT NONE
SAVE
INCLUDE "Parameter.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for EOS table !

! Eos table !
INTEGER :: eosline1
INTEGER :: eosline2

! Line number of EOS table !
INTEGER :: eoslineno1 = 28800
INTEGER :: eoslineno2 = 28800

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for DM/NM scalar quantity !

! mass !
REAL (DP) :: mass1, mass2

! radius !
REAL (DP) :: rad1, rad2

! For DM Central density !
REAL (DP) :: log10rho1_c

! Central and atmospheric values !
REAL (DP) :: rho1_c, rho1_a, p1_c, p1_a, rhoe1_c, rhoe1_a
REAL (DP) :: rho2_c, rho2_a, p2_c, p2_a, rhoe2_c, rhoe2_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for arrays !

! This is for solving the RK4 equation !
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: y_rk4

! DM primitive variables variables !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: rho1
REAL (DP), ALLOCATABLE, DIMENSION (:) :: epsilon1

! NM primitive variables variables !
REAL (DP), ALLOCATABLE, DIMENSION (:) :: rho2
REAL (DP), ALLOCATABLE, DIMENSION (:) :: epsilon2

! Storing NM and DM EOS table information !
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: eostable1
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: eostable2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE DEFINITION