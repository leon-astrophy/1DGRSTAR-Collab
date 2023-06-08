!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve for the initial hydrostatic equilibrium star !
! assuming a two fluid formalism. We assume the newtonian gravity    !
! and is solving for the initial density profile using RK-5 method   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO1F
USE DEFINITION
IMPLICIT NONE

! The number of hydrostatic equation and the improved accuracy !
INTEGER, PARAMETER :: no_of_eq_ini = 3

! Integer parameters !
INTEGER :: i, j, end

! Distance !
REAL (DP) :: x

! Pressure and density !
REAL (DP) :: ini_p2, ini_rho2, ini_rhoe2

! Variables essential in the RK-5 method !
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: y_zero, y_one, y_two, y_three, y_four, y_five
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: der_one, der_two, der_three, der_four, der_five, der_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Clear array !
y_rk4 = 0.0d0

! We convert density at center and atmosphere to the corresponding pressure !
IF(nmtable_flag == 1) THEN

	! Atmospheric density
	rho2_a = 1.1D0*MINVAL(eostable2(:,2))

	CALL tabeos_rtop(rho2_c, p2_c, 2)
	CALL tabeos_rtop(rho2_a, p2_a, 2)

	! Atmospheric density
	rhoe2_c = 0.0D0
	rhoe2_a = 0.0d0

ELSE

	! Atmospheric density
	rho2_a = rho2_c*rhofac

	CALL GETRHO_EOSRTOP (p2_c, rho2_c, gs2, mb2, me2, ye2, 2)
	CALL GETRHO_EOSRTOP (p2_a, rho2_a, gs2, mb2, me2, ye2, 2)
	CALL EOSINTERNAL (rho2_c, p2_c, rhoe2_c, 2)
	CALL EOSINTERNAL (rho2_a, p2_a, rhoe2_a, 2)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign the value of y at center !
y_rk4 (1, 0) = 0.0E0_DP
y_rk4 (2, 0) = p2_c
y_rk4 (3, 0) = rho2_c + rhoe2_c

! this is the main structure of RK-5 method !
DO j = 0, length_step - 1

	! Update the value of x and y !
	DO i = 1, no_of_eq_ini
		y_zero (i) = y_rk4 (i, j)
	END DO

	! Update distance !
	x = DBLE (j) * dx
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the first step in RK-5 !
	CALL INI_DER_1F (der_one, x, y_zero, no_of_eq_ini)
	
	! We artificially make the value to zero to avoid singularity at center !
	IF (j == 0) THEN
		der_one (2) = 0.0E0_DP
	END IF

	! Update !
	DO i = 1, no_of_eq_ini
		y_one (i) = y_zero (i) + (1.0E0_DP / 4.0E0_DP) * dx * der_one (i)
	END DO

	! Assign pressure !
	ini_p2 = y_one (2)
	
	! We determine whether the pressure reached atmospheric pressure !
	IF (ini_p2 <= p2_a) THEN
		end = j
		EXIT
	END IF
	
	! Convert !
	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF

	! Assign density according to pressure !
	y_one (3) = ini_rho2 + ini_rhoe2

	x = DBLE (j) * dx + (1.0E0_DP / 4.0E0_DP) * dx
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the second step in RK-5 !
	CALL INI_DER_1F (der_two, x, y_one, no_of_eq_ini)

	DO i = 1, no_of_eq_ini
		y_two (i) = y_zero (i) + (1.0E0_DP / 8.0E0_DP) * dx * der_one (i) &
				+ (1.0E0_DP / 8.0E0_DP) * dx * der_two (i)
	END DO

	ini_p2 = y_two (2)
	
	IF (ini_p2 <= p2_a) THEN
		end = j
		EXIT
	END IF

	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF
	
	y_two (3) = ini_rho2 + ini_rhoe2

	x = DBLE (j) * dx + (1.0E0_DP / 4.0E0_DP) * dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the third step in RK-5 !
	CALL INI_DER_1F (der_three, x, y_two, no_of_eq_ini)

	DO i = 1, no_of_eq_ini
		y_three (i) = y_zero (i) - (1.0E0_DP / 2.0E0_DP) * dx * der_two (i) &
			+ dx * der_three (i)
	END DO
	
	ini_p2 = y_three (2)
	
	IF (ini_p2 <= p2_a) THEN
		end = j
		EXIT
	END IF

	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF
	
	y_three (3) = ini_rho2 + ini_rhoe2

	x = DBLE (j) * dx + (1.0E0_DP / 2.0E0_DP) * dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the forth step in RK-5 !
	CALL INI_DER_1F (der_four, x, y_three, no_of_eq_ini)

	DO i = 1, no_of_eq_ini
		y_four (i) = y_zero (i) + (3.0E0_DP / 1.6E1_DP) * dx * der_one (i) &
			+ (9.0E0_DP / 1.6E1_DP) * dx * der_four (i) 
	END DO

	ini_p2 = y_four (2)

	IF (ini_p2 <= p2_a) THEN
		end = j
		EXIT
	END IF

	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF
	
	y_four (3) = ini_rho2 + ini_rhoe2

	x = DBLE (j) * dx + (3.0E0_DP / 4.0E0_DP) * dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the fifth step in RK-5 !
	CALL INI_DER_1F (der_five, x, y_four, no_of_eq_ini)

	DO i = 1, no_of_eq_ini
		y_five (i) = y_zero (i) - (3.0E0_DP / 7.0E0_DP) * dx * der_one (i) &
			+ (2.0E0_DP / 7.0E0_DP) * dx * der_two (i) &
			+ (1.2E1_DP / 7.0E0_DP) * dx * der_three (i) &
			- (1.2E1_DP / 7.0E0_DP) * dx * der_four (i) &
			+ (8.0E0_DP / 7.0E0_DP) * dx * der_five (i)
	END DO

	ini_p2 = y_five (2)

	IF (ini_p2 <= p2_a) THEN
		end = j
		EXIT
	END IF

	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF
	
	y_five (3) = ini_rho2 + ini_rhoe2

	x = DBLE (j) * dx + dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the final step in RK-5 !
	CALL INI_DER_1F (der_new, x, y_five, no_of_eq_ini)

	DO i = 1, no_of_eq_ini
		y_rk4 (i, j + 1) = y_zero (i) + (7.0E0_DP * der_one (i) + 3.2E1_DP * der_three (i) & 
				+ 1.2E1_DP * der_four (i) + 3.2E1_DP * der_five (i) & 
				+ 7.0E0_DP * der_new (i)) * dx / 9.0E1_DP
	END DO

	ini_p2 = y_rk4 (2, j + 1)

	IF (ini_p2 <= p2_a) THEN
		end = j
		EXIT
	END IF
	
	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF

	y_rk4 (3, j + 1) = ini_rho2 + ini_rhoe2

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Print results !
WRITE (101,701) log10(rho2_c/rhocgs2code), y_rk4 (1, end), DBLE (end)*dx/lencgs2code/rsolar

! Format !
701 FORMAT (10ES33.15)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine solve for the initial hydrostatic equilibrium star !
! assuming a two fluid formalism. We assume the newtonian gravity    !
! and is solving for the initial density profile using RK-5 method   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETRHO2F
USE DEFINITION
IMPLICIT NONE

! The number of hydrostatic equation and the improved accuracy !
INTEGER, PARAMETER :: no_of_eq_ini = 6

! Integer parameters !
INTEGER :: i, j, end, atm_1, atm_2

! Distance !
REAL (DP) :: x

! Pressure and density !
REAL (DP) :: ini_p1, ini_rho1, ini_rhoe1
REAL (DP) :: ini_p2, ini_rho2, ini_rhoe2

! Variables essential in the RK-5 method !
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: y_zero, y_one, y_two, y_three, y_four, y_five
REAL (DP), DIMENSION (1 : no_of_eq_ini) :: der_one, der_two, der_three, der_four, der_five, der_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Clear array !
y_rk4 = 0.0d0

! We convert density at center and atmosphere to the corresponding pressure !
IF(dmtable_flag == 1) THEN

	! Atmospheric density
	rho1_a = 1.1D0*MINVAL(eostable1(:,2))

	CALL tabeos_rtop(rho1_c, p1_c, 1)
	CALL tabeos_rtop(rho1_a, p1_a, 1)

	! Atmospheric density
	rhoe1_c = 0.0D0
	rhoe1_a = 0.0d0

ELSE

	! Atmospheric density !
	rho1_a = rho1_c*rhofac

	! We convert density at center and atmosphere to the corresponding pressure !
	CALL GETRHO_EOSRTOP (p1_c, rho1_c, gs1, mb1, me1, ye1, 1)
	CALL GETRHO_EOSRTOP (p1_a, rho1_a, gs1, mb1, me1, ye1, 1)
	CALL EOSINTERNAL (rho1_c, p1_c, rhoe1_c, 1)
	CALL EOSINTERNAL (rho1_a, p1_a, rhoe1_a, 1)

END IF

! We convert density at center and atmosphere to the corresponding pressure !
IF(nmtable_flag == 1) THEN

	! Atmospheric density
	rho2_a = 1.1D0*MINVAL(eostable2(:,2))

	CALL tabeos_rtop(rho2_c, p2_c, 2)
	CALL tabeos_rtop(rho2_a, p2_a, 2)

	! Atmospheric density
	rhoe2_c = 0.0D0
	rhoe2_a = 0.0d0

ELSE

	! Atmospheric density !
	rho2_a = rho2_c*rhofac

	CALL GETRHO_EOSRTOP (p2_c, rho2_c, gs2, mb2, me2, ye2, 2)
	CALL GETRHO_EOSRTOP (p2_a, rho2_a, gs2, mb2, me2, ye2, 2)
	CALL EOSINTERNAL (rho2_c, p2_c, rhoe2_c, 2)
	CALL EOSINTERNAL (rho2_a, p2_a, rhoe2_a, 2)

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! We assign the value of y at center !
y_rk4 (1, 0) = 0.0E0_DP
y_rk4 (2, 0) = p1_c
y_rk4 (3, 0) = rho1_c + rhoe1_c
y_rk4 (4, 0) = 0.0E0_DP
y_rk4 (5, 0) = p2_c
y_rk4 (6, 0) = rho2_c + rhoe2_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
atm_1 = length_step
atm_2 = length_step

! this is the main structure of RK-5 method !
DO j = 0, length_step - 1

	! Exit condition !
	IF (atm_1 < length_step .AND. atm_2 < length_step) THEN
		EXIT
	END IF
	
	! Update the value of x and y !
	DO i = 1, no_of_eq_ini
		y_zero (i) = y_rk4 (i, j)
	END DO

	x = DBLE (j) * dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the first step in RK-5 !
	CALL INI_DER_2F (der_one, x, y_zero, no_of_eq_ini)
	
	! We artificially make the value to zero to avoid singularity at center !
	IF (j == 0) THEN
		der_one (2) = 0.0E0_DP
		der_one (5) = 0.0E0_DP
	END IF

	! If the density reach atmospheric values, no changes in all the quantity !
	IF (y_zero (3) <= rho1_a + rhoe1_a) THEN
		der_one (1) = 0.0E0_DP
		der_one (2) = 0.0E0_DP
	END IF
	IF (y_zero (6) <= rho2_a + rhoe1_a) THEN
		der_one (4) = 0.0E0_DP
		der_one (5) = 0.0E0_DP
	END IF

	! Update !
	DO i = 1, no_of_eq_ini
		y_one (i) = y_zero (i) + (1.0E0_DP / 4.0E0_DP) * dx * der_one (i)
	END DO
	
	! We determine whether the pressure reached atmospheric pressure !
	IF (j < atm_1) THEN
		ini_p1 = y_one (2)
		IF (ini_p1 <= p1_a) THEN
			atm_1 = j
			ini_p1 = p1_a
			y_one (2) = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
		y_one (2) = p1_a
	END IF
	IF (j < atm_2) THEN
		ini_p2 = y_one (5)
		IF (ini_p2 <= p2_a) THEN
			atm_2 = j
			ini_p2 = p2_a
			y_one (5) = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
		y_one (5) = p2_a
	END IF
	
	! Invert !
	IF(dmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho1, ini_p1, 1)
		ini_rhoe1 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eosline1, 1)
		CALL EOSINTERNAL (ini_rho1, ini_p1, ini_rhoe1, 1)
	END IF
	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF

	! Assign density according to pressure !
	If (ini_p1 == p1_a) THEN
		y_one (3) = rho1_a + rhoe1_a
	ELSE
		y_one (3) = ini_rho1 + ini_rhoe1
	END IF
	If (ini_p2 == p2_a) THEN
		y_one (6) = rho2_a + rhoe2_a
	ELSE
		y_one (6) = ini_rho2 + ini_rhoe2
	END IF

	x = DBLE (j) * dx + (1.0E0_DP / 4.0E0_DP) * dx
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the second step in RK-5 !
	CALL INI_DER_2F (der_two, x, y_one, no_of_eq_ini)

	IF (y_one (3) <= rho1_a + rhoe1_a) THEN
		der_two (1) = 0.0E0_DP
		der_two (2) = 0.0E0_DP
	END IF
	IF (y_one (6) <= rho2_a + rhoe2_a) THEN
		der_two (4) = 0.0E0_DP
		der_two (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_two (i) = y_zero (i) + (1.0E0_DP / 8.0E0_DP) * dx * der_one (i) &
				+ (1.0E0_DP / 8.0E0_DP) * dx * der_two (i)
	END DO
	
	IF (j < atm_1) THEN
		ini_p1 = y_two (2)
		IF (ini_p1 <= p1_a) THEN
			atm_1 = j
			ini_p1 = p1_a
			y_two (2) = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
		y_two (2) = p1_a
	END IF
	IF (j < atm_2) THEN
		ini_p2 = y_two (5)
		IF (ini_p2 <= p2_a) THEN
			atm_2 = j
			ini_p2 = p2_a
			y_two (5) = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
		y_two (5) = p2_a
	END IF

	IF(dmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho1, ini_p1, 1)
		ini_rhoe1 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eosline1, 1)
		CALL EOSINTERNAL (ini_rho1, ini_p1, ini_rhoe1, 1)
	END IF
	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF

	If (ini_p1 == p1_a) THEN
		y_two (3) = rho1_a + rhoe1_a
	ELSE
		y_two (3) = ini_rho1 + ini_rhoe1
	END IF
	If (ini_p2 == p2_a) THEN
		y_two (6) = rho2_a + rhoe2_a
	ELSE
		y_two (6) = ini_rho2 + ini_rhoe2
	END IF

	x = DBLE (j) * dx + (1.0E0_DP / 4.0E0_DP) * dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the third step in RK-5 !
	CALL INI_DER_2F (der_three, x, y_two, no_of_eq_ini)
	
	IF (y_two (3) <= rho1_a + rhoe1_a) THEN
		der_three (1) = 0.0E0_DP
		der_three (2) = 0.0E0_DP
	END IF
	IF (y_two (6) <= rho2_a + rhoe2_a) THEN
		der_three (4) = 0.0E0_DP
		der_three (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_three (i) = y_zero (i) - (1.0E0_DP / 2.0E0_DP) * dx * der_two (i) &
			+ dx * der_three (i)
	END DO

	IF (j < atm_1) THEN
		ini_p1 = y_three (2)
		IF (ini_p1 <= p1_a) THEN
			atm_1 = j
			ini_p1 = p1_a
			y_three (2) = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
		y_three (2) = p1_a
	END IF
	IF (j < atm_2) THEN
		ini_p2 = y_three (5)
		IF (ini_p2 <= p2_a) THEN
			atm_2 = j
			ini_p2 = p2_a
			y_three (5) = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
		y_three (5) = p2_a
	END IF

	IF(dmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho1, ini_p1, 1)
		ini_rhoe1 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eosline1, 1)
		CALL EOSINTERNAL (ini_rho1, ini_p1, ini_rhoe1, 1)
	END IF
	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF

	If (ini_p1 == p1_a) THEN
		y_three (3) = rho1_a + rhoe1_a
	ELSE
		y_three (3) = ini_rho1 + ini_rhoe1
	END IF
	If (ini_p2 == p2_a) THEN
		y_three (6) = rho2_a + rhoe2_a
	ELSE
		y_three (6) = ini_rho2 + ini_rhoe2
	END IF

	x = DBLE (j) * dx + (1.0E0_DP / 2.0E0_DP) * dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the forth step in RK-5 !
	CALL INI_DER_2F (der_four, x, y_three, no_of_eq_ini)
	
	IF (y_three (3) <= rho1_a + rhoe1_a) THEN
		der_four (1) = 0.0E0_DP
		der_four (2) = 0.0E0_DP
	END IF
	IF (y_three (6) <= rho2_a + rhoe2_a) THEN
		der_four (4) = 0.0E0_DP
		der_four (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_four (i) = y_zero (i) + (3.0E0_DP / 1.6E1_DP) * dx * der_one (i) &
			+ (9.0E0_DP / 1.6E1_DP) * dx * der_four (i) 
	END DO

	IF (j < atm_1) THEN
		ini_p1 = y_four (2)
		IF (ini_p1 <= p1_a) THEN
			atm_1 = j
			ini_p1 = p1_a
			y_four (2) = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
		y_four (2) = p1_a
	END IF
	IF (j < atm_2) THEN
		ini_p2 = y_four (5)
		IF (ini_p2 <= p2_a) THEN
			atm_2 = j
			ini_p2 = p2_a
			y_four (5) = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
		y_four (5) = p2_a
	END IF

	IF(dmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho1, ini_p1, 1)
		ini_rhoe1 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eosline1, 1)
		CALL EOSINTERNAL (ini_rho1, ini_p1, ini_rhoe1, 1)
	END IF
	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF

	If (ini_p1 == p1_a) THEN
		y_four (3) = rho1_a + rhoe1_a
	ELSE
		y_four (3) = ini_rho1 + ini_rhoe1
	END IF
	If (ini_p2 == p2_a) THEN
		y_four (6) = rho2_a + rhoe2_a
	ELSE
		y_four (6) = ini_rho2 + ini_rhoe2
	END IF

	x = DBLE (j) * dx + (3.0E0_DP / 4.0E0_DP) * dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the fifth step in RK-5 !
	CALL INI_DER_2F (der_five, x, y_four, no_of_eq_ini)
	
	IF (y_four (3) <= rho1_a + rhoe1_a) THEN
		der_five (1) = 0.0E0_DP
		der_five (2) = 0.0E0_DP
	END IF
	IF (y_four (6) <= rho2_a + rhoe2_a) THEN
		der_five (4) = 0.0E0_DP
		der_five (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_five (i) = y_zero (i) - (3.0E0_DP / 7.0E0_DP) * dx * der_one (i) &
			+ (2.0E0_DP / 7.0E0_DP) * dx * der_two (i) &
			+ (1.2E1_DP / 7.0E0_DP) * dx * der_three (i) &
			- (1.2E1_DP / 7.0E0_DP) * dx * der_four (i) &
			+ (8.0E0_DP / 7.0E0_DP) * dx * der_five (i)
	END DO

	IF (j < atm_1) THEN
		ini_p1 = y_five (2)
		IF (ini_p1 <= p1_a) THEN
			atm_1 = j
			ini_p1 = p1_a
			y_five (2) = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
		y_five (2) = p1_a
	END IF

	IF (j < atm_2) THEN
		ini_p2 = y_five (5)
		IF (ini_p2 <= p2_a) THEN
			atm_2 = j
			ini_p2 = p2_a
			y_five (5) = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
		y_five (5) = p2_a
	END IF

	IF(dmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho1, ini_p1, 1)
		ini_rhoe1 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eosline1, 1)
		CALL EOSINTERNAL (ini_rho1, ini_p1, ini_rhoe1, 1)
	END IF
	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF

	If (ini_p1 == p1_a) THEN
		y_five (3) = rho1_a + rhoe1_a
	ELSE
		y_five (3) = ini_rho1 + ini_rhoe1
	END IF
	If (ini_p2 == p2_a) THEN
		y_five (6) = rho2_a + rhoe2_a
	ELSE
		y_five (6) = ini_rho2 + ini_rhoe2
	END IF

	x = DBLE (j) * dx + dx

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! This is the final step in RK-5 !
	CALL INI_DER_2F (der_new, x, y_five, no_of_eq_ini)
	
	IF (y_five (3) <= rho1_a + rhoe1_a) THEN
		der_new (1) = 0.0E0_DP
		der_new (2) = 0.0E0_DP
	END IF
	IF (y_five (6) <= rho2_a + rhoe2_a) THEN
		der_one (4) = 0.0E0_DP
		der_one (5) = 0.0E0_DP
	END IF

	DO i = 1, no_of_eq_ini
		y_rk4 (i, j + 1) = y_zero (i) + (7.0E0_DP * der_one (i) + 3.2E1_DP * der_three (i) & 
				+ 1.2E1_DP * der_four (i) + 3.2E1_DP * der_five (i) & 
				+ 7.0E0_DP * der_new (i)) * dx / 9.0E1_DP
	END DO

	IF (j < atm_1) THEN
		ini_p1 = y_rk4 (2, j + 1)
		IF (ini_p1 <= p1_a) THEN
			atm_1 = j
			ini_p1 = p1_a
			y_rk4 (2, j + 1) = p1_a
		END IF
	ELSE
		ini_p1 = p1_a
		y_rk4 (2, j + 1) = p1_a
	END IF
	IF (j < atm_2) THEN
		ini_p2 = y_rk4 (5, j + 1)
		IF (ini_p2 <= p2_a) THEN
			atm_2 = j
			ini_p2 = p2_a
			y_rk4 (5, j + 1) = p2_a
		END IF
	ELSE
		ini_p2 = p2_a
		y_rk4 (5, j + 1) = p2_a
	END IF

	IF(dmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho1, ini_p1, 1)
		ini_rhoe1 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho1, ini_p1, eosline1, 1)
		CALL EOSINTERNAL (ini_rho1, ini_p1, ini_rhoe1, 1)
	END IF
	IF(nmtable_flag == 1) THEN
		CALL tabeos_ptor(ini_rho2, ini_p2, 2)
		ini_rhoe2 = 0.0D0
	ELSE
		CALL GETRHO_EOSPTOR (ini_rho2, ini_p2, eosline2, 2)
		CALL EOSINTERNAL (ini_rho2, ini_p2, ini_rhoe2, 2)
	END IF
	
	If (ini_p1 == p1_a) THEN
		y_rk4 (3, j + 1) = rho1_a + rhoe1_a
	ELSE
		y_rk4 (3, j + 1) = ini_rho1 + ini_rhoe1
	END IF
	If (ini_p2 == p2_a) THEN
		y_rk4 (6, j + 1) = rho2_a + rhoe2_a
	ELSE
		y_rk4 (6, j + 1) = ini_rho2 + ini_rhoe2
	END IF

END DO	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Look for the maxmimum radius between DM and NM !
end = INT(MAX(DBLE(atm_1), DBLE(atm_2)))

! Assign results !
mass1 = y_rk4 (1, end)
mass2 = y_rk4 (4, end)
rad1 = DBLE (atm_1)*dx/lencgs2code/rsolar
rad2 = DBLE (atm_2)*dx/lencgs2code/rsolar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE