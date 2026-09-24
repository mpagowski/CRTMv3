
!
! test_AerosolScatter_AAOD_TL
!
! AAOD version of test_AerosolScatter_TL: unit test for the aerosol absorption
! optical depth part of CRTM_Compute_AerosolScatter_TL
! (AScat_TL%Absorption_Optical_Depth, AAOD = sum (1-w).ke.rho).
! Same sensor (modis_aqua, channel 1, ~3.8 um, where sea salt absorbs), same
! aerosol (GOCART-GEOS5 sea salt 1) and same three perturbations (relative
! humidity, effective radius, concentration) as the original test.
! Verifies that the AAOD sensitivities to RH and concentration are non-zero,
! that the one to effective radius is exactly zero (no Reff dependence in the
! GOCART-GEOS5 scheme), that AAOD_TL == OD_TL - SSA_TL, and that each non-zero
! tangent-linear matches a central finite difference of the forward model.
! RH is 0.53 (the original uses 0.50, a LUT RH node, where a finite difference
! would straddle an interpolation-stencil switch).
! Not covered here (one layer, one aerosol, scattering on, no clamp active):
! the multi-layer/multi-aerosol accumulation, the scattering-off branch and
! the QC clamps — see the k_matrix regression test test_AAOD.
!

PROGRAM test_AerosolScatter_AAOD_TL

  USE Type_Kinds,             ONLY: fp
  USE Message_Handler,        ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,        ONLY: ZERO, ONE, TWO
  USE CRTM_SpcCoeff,          ONLY: SC, CRTM_SpcCoeff_Load, CRTM_SpcCoeff_Destroy
  USE CRTM_AerosolCoeff,      ONLY: AeroC, CRTM_AerosolCoeff_Load, CRTM_AerosolCoeff_Destroy
  USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type, &
                                    CRTM_Atmosphere_Create, &
                                    CRTM_Atmosphere_Destroy, &
                                    CRTM_Atmosphere_Zero
  USE CRTM_AtmOptics_Define,  ONLY: CRTM_AtmOptics_type, &
                                    CRTM_AtmOptics_Create, &
                                    CRTM_AtmOptics_Destroy, &
                                    CRTM_AtmOptics_Zero
  USE CRTM_AerosolScatter,    ONLY: CRTM_Compute_AerosolScatter, &
                                    CRTM_Compute_AerosolScatter_TL
  USE ASvar_Define,           ONLY: ASvar_type, &
                                    ASvar_Create, ASvar_Destroy
  USE CRTM_ChannelInfo_Define, ONLY: CRTM_ChannelInfo_type

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_AAOD_TL'
  INTEGER :: Error_Status
  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  CHARACTER(256) :: File_Path = './testinput/'

  TYPE(CRTM_Atmosphere_type) :: Atm, Atm_TL, Atm_p, Atm_m
  TYPE(CRTM_AtmOptics_type)  :: AScat, AScat_TL, AScat_p, AScat_m
  TYPE(ASvar_type)           :: ASvar, ASvar_fd
  INTEGER :: SensorIndex = 1
  INTEGER :: ChannelIndex = 1
  INTEGER :: n_Layers = 1
  INTEGER :: n_Aerosols = 1
  INTEGER :: n_Legendre = 4
  INTEGER :: n_Phase = 1
  REAL(fp), PARAMETER :: RH_0      = 0.53_fp    ! off the LUT RH nodes
  REAL(fp), PARAMETER :: H_RH      = 1.0e-4_fp  ! finite-difference step, RH
  REAL(fp), PARAMETER :: EPS_C     = 1.0e-3_fp  ! finite-difference step, relative concentration
  REAL(fp), PARAMETER :: FD_TOL_RH = 1.0e-5_fp  ! cubic RH interpolation: truncation ~1e-7
  REAL(fp), PARAMETER :: FD_TOL_C  = 1.0e-10_fp ! AAOD is linear in the concentration
  REAL(fp), PARAMETER :: ID_TOL    = 1.0e-12_fp
  REAL(fp) :: tl, fd, denom, floor_val

  PRINT *, 'Starting test_AerosolScatter_AAOD_TL...'

  ! Load Coefficients
  Error_Status = CRTM_SpcCoeff_Load( (/Sensor_Id/), File_Path=File_Path )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Using GOCART scheme for RH sensitivity
  Error_Status = CRTM_AerosolCoeff_Load('GOCART-GEOS5', 'AerosolCoeff.GOCART-GEOS5.bin', File_Path=File_Path)
  IF ( Error_Status /= SUCCESS ) STOP 1
  PRINT *, 'Channel wavenumber (cm-1): ', SC(SensorIndex)%Wavenumber(ChannelIndex)

  ! Create structures
  CALL CRTM_Atmosphere_Create( Atm, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_TL, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_AtmOptics_Create( AScat, n_Layers, n_Legendre, n_Phase )
  CALL CRTM_AtmOptics_Create( AScat_TL, n_Layers, n_Legendre, n_Phase )
  CALL CRTM_AtmOptics_Create( AScat_p, n_Layers, n_Legendre, n_Phase )
  CALL CRTM_AtmOptics_Create( AScat_m, n_Layers, n_Legendre, n_Phase )
  CALL ASvar_Create( ASvar, n_Legendre, n_Phase, n_Layers, n_Aerosols )
  CALL ASvar_Create( ASvar_fd, n_Legendre, n_Phase, n_Layers, n_Aerosols )

  ! Populate Input
  Atm%Absorber_ID(1) = 1
  Atm%Aerosol(1)%Type = 6   ! Sea salt 1
  Atm%Aerosol(1)%Concentration(1) = 1.0e-4_fp
  Atm%Aerosol(1)%Effective_Radius(1) = 1.0_fp
  Atm%Relative_Humidity(1) = RH_0

  AScat%Include_Scattering = .TRUE.
  AScat_TL%Include_Scattering = .TRUE.
  AScat_p%Include_Scattering = .TRUE.
  AScat_m%Include_Scattering = .TRUE.

  ! Forward Call to populate ASvar
  CALL CRTM_AtmOptics_Zero( AScat )
  Error_Status = CRTM_Compute_AerosolScatter( Atm, SensorIndex, ChannelIndex, AScat, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1
  PRINT *, 'Forward: OD, AAOD, w (LUT, clamped): ', AScat%Optical_Depth(1), AScat%Absorption_Optical_Depth(1), ASvar%w(1,1)
  IF ( .NOT. ( AScat%Absorption_Optical_Depth(1) > ZERO .AND. &
               AScat%Absorption_Optical_Depth(1) < AScat%Optical_Depth(1) ) ) THEN
    PRINT *, 'FAIL: forward AAOD must be in (0, OD) for this channel (sea salt absorbs in the IR).'
    STOP 1
  END IF
  floor_val = 1.0e-10_fp * AScat%Optical_Depth(1)

  ! --------------------------------------------------------------------------
  ! TL Call for RH
  ! --------------------------------------------------------------------------
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  Atm_TL%Relative_Humidity(1) = 1.0_fp
  CALL CRTM_AtmOptics_Zero( AScat_TL )
  Error_Status = CRTM_Compute_AerosolScatter_TL( Atm, AScat, Atm_TL, SensorIndex, ChannelIndex, AScat_TL, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1
  tl = AScat_TL%Absorption_Optical_Depth(1)
  PRINT *, 'AScat_TL%Absorption_Optical_Depth(1) for RH perturbation: ', tl
  IF ( .NOT. ( ABS(tl) >= 1.0e-12_fp ) ) THEN
     PRINT *, 'FAIL: AAOD sensitivity to RH is ZERO'
     STOP 1
  END IF
  PRINT *, 'SUCCESS: AAOD sensitivity to RH is NON-ZERO'
  CALL Check_Identity( 'RH' )
  ! ...central finite difference in RH
  Atm_p = Atm; Atm_m = Atm
  Atm_p%Relative_Humidity(1) = RH_0 + H_RH
  Atm_m%Relative_Humidity(1) = RH_0 - H_RH
  CALL Forward_Pair()
  fd = ( AScat_p%Absorption_Optical_Depth(1) - AScat_m%Absorption_Optical_Depth(1) ) / ( TWO*H_RH )
  denom = MAX( ABS(tl), ABS(fd), floor_val )
  PRINT *, 'RH: TL = ', tl, ' FD = ', fd, ' rel. diff = ', ABS(tl-fd)/denom
  IF ( .NOT. ( ABS(tl-fd) <= FD_TOL_RH*denom ) ) THEN
     PRINT *, 'FAIL: AAOD tangent-linear vs finite difference (RH)'
     STOP 1
  END IF
  PRINT *, 'SUCCESS: AAOD tangent-linear matches the finite difference (RH)'

  ! --------------------------------------------------------------------------
  ! TL Call for Effective Radius: the GOCART-GEOS5 scheme has no Reff
  ! dependence, so the sensitivity (and the finite difference) must be zero
  ! --------------------------------------------------------------------------
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  Atm_TL%Aerosol(1)%Effective_Radius(1) = 1.0_fp
  CALL CRTM_AtmOptics_Zero( AScat_TL )
  Error_Status = CRTM_Compute_AerosolScatter_TL( Atm, AScat, Atm_TL, SensorIndex, ChannelIndex, AScat_TL, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1
  tl = AScat_TL%Absorption_Optical_Depth(1)
  PRINT *, 'AScat_TL%Absorption_Optical_Depth(1) for Reff perturbation: ', tl
  Atm_p = Atm; Atm_m = Atm
  Atm_p%Aerosol(1)%Effective_Radius(1) = 1.5_fp
  Atm_m%Aerosol(1)%Effective_Radius(1) = 0.5_fp
  CALL Forward_Pair()
  IF ( tl /= ZERO .OR. AScat_p%Absorption_Optical_Depth(1) /= AScat_m%Absorption_Optical_Depth(1) ) THEN
     PRINT *, 'FAIL: AAOD sensitivity to Reff is not exactly zero for the GOCART-GEOS5 scheme'
     STOP 1
  END IF
  PRINT *, 'SUCCESS: AAOD sensitivity to Reff is exactly ZERO (GOCART-GEOS5)'

  ! --------------------------------------------------------------------------
  ! TL Call for Concentration
  ! --------------------------------------------------------------------------
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  Atm_TL%Aerosol(1)%Concentration(1) = 1.0_fp
  CALL CRTM_AtmOptics_Zero( AScat_TL )
  Error_Status = CRTM_Compute_AerosolScatter_TL( Atm, AScat, Atm_TL, SensorIndex, ChannelIndex, AScat_TL, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1
  tl = AScat_TL%Absorption_Optical_Depth(1)
  PRINT *, 'AScat_TL%Absorption_Optical_Depth(1) for Concentration perturbation: ', tl
  IF ( .NOT. ( ABS(tl) >= 1.0e-12_fp ) ) THEN
     PRINT *, 'FAIL: AAOD sensitivity to Concentration is ZERO'
     STOP 1
  END IF
  PRINT *, 'SUCCESS: AAOD sensitivity to Concentration is NON-ZERO'
  CALL Check_Identity( 'Concentration' )
  ! ...dAAOD/dC = (1-w).ke exactly; central finite difference
  Atm_p = Atm; Atm_m = Atm
  Atm_p%Aerosol(1)%Concentration(1) = Atm%Aerosol(1)%Concentration(1) * ( ONE + EPS_C )
  Atm_m%Aerosol(1)%Concentration(1) = Atm%Aerosol(1)%Concentration(1) * ( ONE - EPS_C )
  CALL Forward_Pair()
  fd = ( AScat_p%Absorption_Optical_Depth(1) - AScat_m%Absorption_Optical_Depth(1) ) / &
       ( TWO*EPS_C*Atm%Aerosol(1)%Concentration(1) )
  denom = MAX( ABS(tl), ABS(fd), floor_val )
  PRINT *, 'Concentration: TL = ', tl, ' FD = ', fd, ' (1-w).ke = ', (ONE-ASvar%w(1,1))*ASvar%ke(1,1), &
           ' rel. diff = ', ABS(tl-fd)/denom
  IF ( .NOT. ( ABS(tl-fd) <= FD_TOL_C*denom ) ) THEN
     PRINT *, 'FAIL: AAOD tangent-linear vs finite difference (Concentration)'
     STOP 1
  END IF
  PRINT *, 'SUCCESS: AAOD tangent-linear matches the finite difference (Concentration)'

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_p )
  CALL CRTM_Atmosphere_Destroy( Atm_m )
  CALL CRTM_AtmOptics_Destroy( AScat )
  CALL CRTM_AtmOptics_Destroy( AScat_TL )
  CALL CRTM_AtmOptics_Destroy( AScat_p )
  CALL CRTM_AtmOptics_Destroy( AScat_m )
  CALL ASvar_Destroy( ASvar )
  CALL ASvar_Destroy( ASvar_fd )
  Error_Status = CRTM_SpcCoeff_Destroy()
  Error_Status = CRTM_AerosolCoeff_Destroy()

  PRINT *, 'SUCCESS: test_AerosolScatter_AAOD_TL passed.'

CONTAINS

  ! AAOD_TL == OD_TL - SSA_TL (the SSA slot accumulates the TL of sum rho.ke.w)
  SUBROUTINE Check_Identity( label )
    CHARACTER(*), INTENT(IN) :: label
    REAL(fp) :: err, den
    err = ABS( AScat_TL%Absorption_Optical_Depth(1) - &
               ( AScat_TL%Optical_Depth(1) - AScat_TL%Single_Scatter_Albedo(1) ) )
    den = MAX( ABS(AScat_TL%Optical_Depth(1)), floor_val )
    PRINT *, '  OD_TL = ', AScat_TL%Optical_Depth(1), ' SSA_TL = ', AScat_TL%Single_Scatter_Albedo(1), &
             ' AAOD_TL = ', AScat_TL%Absorption_Optical_Depth(1)
    IF ( .NOT. ( err <= ID_TOL*den ) ) THEN
      PRINT *, 'FAIL: AAOD_TL /= OD_TL - SSA_TL for the '//label//' perturbation: ', err
      STOP 1
    END IF
    PRINT *, 'SUCCESS: AAOD_TL == OD_TL - SSA_TL ('//label//')'
  END SUBROUTINE Check_Identity

  ! Forward model for the perturbed atmospheres Atm_p and Atm_m
  SUBROUTINE Forward_Pair()
    INTEGER :: es
    CALL CRTM_AtmOptics_Zero( AScat_p ); CALL CRTM_AtmOptics_Zero( AScat_m )
    es = CRTM_Compute_AerosolScatter( Atm_p, SensorIndex, ChannelIndex, AScat_p, ASvar_fd ); IF ( es /= SUCCESS ) STOP 1
    es = CRTM_Compute_AerosolScatter( Atm_m, SensorIndex, ChannelIndex, AScat_m, ASvar_fd ); IF ( es /= SUCCESS ) STOP 1
  END SUBROUTINE Forward_Pair

END PROGRAM test_AerosolScatter_AAOD_TL
