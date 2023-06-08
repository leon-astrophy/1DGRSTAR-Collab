!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the dimensionless fermi momentum used in ideal fermi gas EOS !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FERMIMO (dlfmmo, rho_in, gs, mb, me, ye)
USE DEFINITION
IMPLICIT NONE

! Real Variables !
! the inputed parameters !
REAL (DP), INTENT (IN) :: rho_in, gs, mb, me, ye 

! the output dimensionless fermi momentum !
REAL (DP), INTENT (OUT) :: dlfmmo 

! find the fermi momentum !
dlfmmo = ((6.0E0_DP * pi_old ** (2.0E0_DP) * h_bar ** (3.0E0_DP) * rho_in * ye) &
	/ (gs * mb * me ** (3.0E0_DP))) ** (1.0E0_DP / 3.0E0_DP)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculate the derivative of rho wtih respect to x !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDDXDRHO (dxdrho, rho_in, gs, mb, me, ye)
USE DEFINITION
IMPLICIT NONE

! Real Variables !
! the inputed parameters !
REAL (DP), INTENT (IN) :: rho_in, gs, mb, me, ye 

! the output derivative !
REAL (DP), INTENT (OUT) :: dxdrho

! temporal variable !
REAL (DP) :: x_0

! find the constant !
x_0 = (6.0E0_DP * pi_old ** (2.0E0_DP) * h_bar ** (3.0E0_DP) * ye / (gs * mb * me ** (3.0E0_DP))) ** (1.0E0_DP / 3.0E0_DP)

! Assign derivative !
dxdrho = (x_0*rho_in**(-2.0D0/3.0D0)/3.0D0)

END SUBROUTINE