!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocate arrays related to hydrodynamics variables !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BUILDHYDRO
USE DEFINITION
IMPLICIT NONE

! DM hydrodynamic variables !
ALLOCATE(rho1(-4 : length_step + 5))
ALLOCATE(epsilon1(-4 : length_step + 5))

! NM hydrodynamic variables !
ALLOCATE(rho2(-4 : length_step + 5))
ALLOCATE(epsilon2(-4 : length_step + 5))

END SUBROUTINE