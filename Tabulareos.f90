!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads neturon star EOS table extracted form CompOSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READEOSTABLE_NM
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, nlines

! Dummy variables !
REAL (DP) :: dummy

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = nmtable_name) 
DO 
  READ (999,*, END=10) 
  nlines = nlines + 1 
END DO 
10 CLOSE (999) 

! Then, allocate arrays for EOS Table !
ALLOCATE(eostable2(nlines,2))

! Read !
OPEN(UNIT=999, FILE = nmtable_name, ACTION='READ')
DO i = 1, nlines
	READ(999,*) eostable2(i,2), eostable2(i,1)
ENDDO
CLOSE(999)

! assign #
eoslineno2 = nlines

! Convert to code unit !
eostable2(:,1) = eostable2(:,1)*(mev2dyncm/fm3tocm3)*(masscgs2code/lencgs2code/tcgs2code**2)
eostable2(:,2) = eostable2(:,2)*(mev2gram/fm3tocm3)*rhocgs2code

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads DM EOS table 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READEOSTABLE_DM
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, nlines

! Dummy variables !
REAL (DP) :: dummy

! Read the number of lines in the file !
nlines = 0 
OPEN (999, file = dmtable_name) 
DO 
  READ (999,*, END=10) 
  nlines = nlines + 1 
END DO 
10 CLOSE (999) 

! Then, allocate arrays for EOS Table !
ALLOCATE(eostable1(nlines,2))

! Read !
OPEN(UNIT=999, FILE = dmtable_name, ACTION='READ')
DO i = 1, nlines
	READ(999,*) eostable1(i,2), eostable1(i,1)
ENDDO
CLOSE(999)

! assign #
eoslineno1 = nlines

! Convert to code unit !
eostable1(:,1) = eostable1(:,1)*(mev2dyncm/fm3tocm3)*(masscgs2code/lencgs2code/tcgs2code**2)
eostable1(:,2) = eostable1(:,2)*(mev2gram/fm3tocm3)*rhocgs2code

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert density to pressure by using neutron star EOS tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE tabeos_rtop(ini_rho, ini_p, type)
USE DEFINITION
IMPLICIT NONE

! Input and output !
INTEGER, INTENT(IN) :: type
REAL (DP), INTENT(IN) :: ini_rho
REAL (DP), INTENT(OUT) :: ini_p

! Integer !
INTEGER :: left, right, m
INTEGER :: upper, lower
REAL (DP) :: x1, x2, y1, y2

! We find the density according to the type, 1 is DM, 2 is NM !
IF (type == 1) THEN

	! Case by case !
	IF (ini_rho < eostable1 (1, 2) .OR. ini_rho > eostable1 (eoslineno1, 2)) THEN
		WRITE (*,*) ini_rho, eostable1 (1, 2), log10(eostable1 (eoslineno1, 2)/(rhocgs2code))
		STOP 'Input density out of DM EOS table range'
	ELSE

		! Binary search !
		left = 1
		right = eoslineno1
    lower = -1
    upper = -1
    do while (left <= right)
      m = (left + right) / 2
      if (eostable1 (m, 2) < ini_rho) then
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
		x1 = eostable1 (lower, 2)
		x2 = eostable1 (upper, 2)
		y1 = eostable1 (lower, 1)
		y2 = eostable1 (upper, 1)
		ini_p = y1*(x2 - ini_rho)/(x2 - x1) + y2*(ini_rho - x1)/(x2 - x1)

	END IF

ELSEIF (type == 2) THEN

	! Case by case !
	IF (ini_rho < eostable2 (1, 2) .OR. ini_rho > eostable2 (eoslineno2, 2)) THEN
		WRITE (*,*) ini_rho, eostable2 (1, 2)
		STOP 'Input density out of NM EOS table range'
	ELSE

		! Binary search !
		left = 1
		right = eoslineno2
    lower = -1
    upper = -1
    do while (left <= right)
      m = (left + right) / 2
      if (eostable2 (m, 2) < ini_rho) then
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
		x1 = eostable2 (lower, 2)
		x2 = eostable2 (upper, 2)
		y1 = eostable2 (lower, 1)
		y2 = eostable2 (upper, 1)
		ini_p = y1*(x2 - ini_rho)/(x2 - x1) + y2*(ini_rho - x1)/(x2 - x1)
    
	END IF

END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine convert density to pressure by using neutron star EOS tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE tabeos_ptor(ini_rho, ini_p, type)
USE DEFINITION
IMPLICIT NONE

! Input and output !
INTEGER, INTENT(IN) :: type
REAL (DP), INTENT(IN) :: ini_p
REAL (DP), INTENT(OUT) :: ini_rho

! Integer !
INTEGER :: left, right, m
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do 2-point interpolation for a given point of existing 
!datum and output the quantity that the user wished to interpolate 	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LINEAR(x0, x1, y0, y1, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: x0, x1, y0, y1, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: nu0, nu1
REAL (DP) :: de0, de1
REAL (DP) :: l0, l1

! Assign numerator !
nu0 = (x_in - x1)
nu1 = (x_in - x0)

! Assign denominator !
de0 = (x0 - x1)
de1 = (x1 - x0)

! Assign polynominal coefficient !
l0 = nu0/de0
l1 = nu1/de1

! Compute the output !
y_out = l0*y0 + l1*y1

END SUBROUTINE