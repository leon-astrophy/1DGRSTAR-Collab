!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine aims at finding the density corresponds !
! to a given pressure in the case of ideal degenerate     !
! fermi gas equation of state by interpolation            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO_EOSPTOR (ini_rho, ini_p, eosline, type)
USE DEFINITION
IMPLICIT NONE

! the input pressure !
REAL (DP), INTENT (IN) :: ini_p

! Integer parameter !
INTEGER, INTENT (IN) :: type

! the outputed density !
REAL (DP), INTENT (OUT) :: ini_rho

! distinguish between NM and DM !
INTEGER, INTENT (INOUT) :: eosline

! Integer !
INTEGER :: left, right, m, i
INTEGER :: upper, lower
REAL (DP) :: x1, x2, y1, y2

! We find the density according to the type, 1 is DM, 2 is NM !
IF (type == 1) THEN
  
	! Case by case !
	IF (ini_p < eostable1 (1, 1) .OR. ini_p > eostable1 (eoslineno1, 1)) THEN
		WRITE (*,*) ini_p, eostable1 (1, 1)
		STOP 'Input pressure out of DM EOS table range'
	ELSE

		! Binary search !
		left = 1
		right = eoslineno1
    lower = -1
    upper = -1
    do while (left <= right)
      m = (left + right) / 2
      if (eostable1 (m, 1) < ini_p) then
        left = m + 1
      else
        right = m - 1
      end if
      if (left > right) then
        lower = right
        upper = left
      end if
    end do

		! Interpolate !
		x1 = eostable1 (lower, 1)
		x2 = eostable1 (upper, 1)
		y1 = eostable1 (lower, 2)
		y2 = eostable1 (upper, 2)
		ini_rho = y1*(x2 - ini_p)/(x2 - x1) + y2*(ini_p - x1)/(x2 - x1)

	END IF

ELSEIF (type == 2) THEN

	! Case by case !
	IF (ini_p < eostable2 (1, 1) .OR. ini_p > eostable2 (eoslineno2, 1)) THEN
		WRITE (*,*) ini_p, eostable2 (1, 1)
		STOP 'Input pressure out of NM EOS table range'
	ELSE

		! Binary search !
		left = 1
		right = eoslineno1
    lower = -1
    upper = -1
    do while (left <= right)
      m = (left + right) / 2
      if (eostable2 (m, 1) < ini_p) then
        left = m + 1
      else
        right = m - 1
      end if
      if (left > right) then
        lower = right
        upper = left
      end if
    end do

		! Interpolate !
		x1 = eostable2 (lower, 1)
		x2 = eostable2 (upper, 1)
		y1 = eostable2 (lower, 2)
		y2 = eostable2 (upper, 2)
		ini_rho = y1*(x2 - ini_p)/(x2 - x1) + y2*(ini_p - x1)/(x2 - x1)

	END IF

END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!