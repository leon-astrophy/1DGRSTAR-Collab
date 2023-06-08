!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the internal energy !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSINTERNAL (den, pre, internal, type)
USE DEFINITION
IMPLICIT NONE

! Input type !
INTEGER, INTENT (IN) :: type

! Input density !
REAL (DP), INTENT (IN) :: den, pre

! Output value !
REAL (DP), INTENT (OUT) :: internal

! Fermi-momentum !
REAL (DP) :: fermi

! For DM Output !
IF(type == 1) THEN
	CALL FERMIMO (fermi, den, gs1, mb1, me1, ye1)
	IF (fermi<=3.0E-3_DP) THEN
		internal = a_max1*small_energy(fermi)
	ELSE
		internal = a_max1*large_energy(fermi)
	END IF

! For NM !
ELSEIF(type == 2) THEN
	CALL FERMIMO (fermi, den, gs2, mb2, me2, ye2)
	IF (fermi<=3.0E-3_DP) THEN
		internal = a_max2*small_energy(fermi)
	ELSE
		internal = a_max2*large_energy(fermi)
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
	real(DP) function large_energy(x)
	implicit none
	real(DP) :: x
	large_energy = 3.0D0*x*SQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + SQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	real(DP) function small_energy(x)
	implicit none
	real(DP) :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
			+ (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE