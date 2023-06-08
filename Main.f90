!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve TOV equations for hydrostatic stars !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MASSRADIUS
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j

! Filename !
character(LEN=10) :: mdm_file

! Real !
REAL (DP) :: check_m_last, check_m, drho1c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for initializing and building arrays !

! Build hydro variables !
CALL BUILDHYDRO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! section for constructing EOS table !

! For NM EOS !
IF(nmtable_flag == 1) THEN
	CALL READEOSTABLE_NM
ELSE
	CALL EOSTABLE_NM
END IF

! For DM EOS !
IF(dmtable_flag == 1) THEN
	CALL READEOSTABLE_DM
ELSE
	CALL EOSTABLE_DM
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section for the main functions !

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do for the DM case !
IF(DM_flag == 1) THEN

	! Print !
	WRITE (*,*) 'Solving For 2F TOV ...'
	WRITE (*,*)	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Filename adjustment !
	write(mdm_file,'(f0.4)') dm_ratio

	! the leading zero is missing, add it back !
	if(mdm_file(1:1) == '.') mdm_file = '0' // mdm_file
	if(mdm_file(1:2) == '-.') mdm_file = '-0.' // mdm_file(3:)

	! Openfile !
	OPEN (UNIT = 101, FILE = './outfile/2F-MR-Relation-DMRatio-'// trim(adjustl(mdm_file)) //'.dat', STATUS = 'REPLACE')

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Allocate arrays !
	ALLOCATE(y_rk4(1 : 6, -4 : length_step + 5))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Compute M-Rho relations !
	DO j = 1, n_rho

		! Print step !
		WRITE (*,*) 'step', j

		! Assign NM maximum rhocgs2code !
		rho2_c = rhostart + (DBLE(j) - 1.0D0)*drho
		rho2_c = (1.0D1**(rho2_c)*rhocgs2code)

		! Guess an initial rhocgs2code !
		log10rho1_c = log10(rho2_c)

		! Set step size !
		drho1c = 0.1D0

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Bisection method !
		DO i = 0, n_max

			! Assign central rhocgs2code !
			rho1_c = 10.0D0**(log10rho1_c)			

			! Solve the equilbrium structure !
			CALL GETRHO2F

			! Debug !
			!WRITE (*,*) i, mass1/(mass1+mass2), log10rho1_c

			! Check if the deviation from the target DM mass !
			check_m = mass1/(mass1+mass2) - dm_ratio

			! Make sure you go to the right direction of dtemp
			if(i == 0) then
				if(check_m > 0.0D0) THEN
					drho1c = -drho1c
				end if
			endif

			! Use shooting method    
			if(check_m_last * check_m < 0.0D0 .and. i /= 0) then
				drho1c = -drho1c * 0.5D0
			end if

			! Update !
			log10rho1_c = log10rho1_c + drho1c
			check_m_last = check_m

			! Exit condition !
			IF(ABS(check_m/dm_ratio) < tor) THEN
				EXIT
			END IF

		END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Print results !
		WRITE (101,701) log10(rho2_c/rhocgs2code), mass1, mass2, rad1, rad2

	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Close file !
	CLOSE (101)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now do for the NM case !
ELSE

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Print !
	WRITE (*,*) 'Solving For 1F TOV ...'
	WRITE (*,*)	

	! Allocate arrays !
	ALLOCATE(y_rk4(1 : 6, -4 : length_step + 5))

	! Openfile !
	OPEN (UNIT = 101, FILE = './outfile/1F-MR-Relation-NM.dat', STATUS = 'REPLACE')

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Loop over the densit yrange !
	DO j = 1, n_rho

		! Print step !
		WRITE (*,*) 'step', j

		! Assign NM maximum rhocgs2code !
		rho2_c = rhostart + (DBLE(j) - 1.0D0)*drho
		rho2_c = (1.0D1**(rho2_c)*rhocgs2code)	

		! Solve the equilbrium structure !
		CALL GETRHO1F

	END DO	

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Final section !

! Print !
WRITE (*,*) 'DONE!'

! Format !
701 FORMAT (10ES33.15)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM
