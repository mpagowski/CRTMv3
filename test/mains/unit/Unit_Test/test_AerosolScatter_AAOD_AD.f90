
!
! test_AerosolScatter_AAOD_AD
!
! AAOD version of test_AerosolScatter_AD: unit test for the aerosol absorption
! optical depth part of CRTM_Compute_AerosolScatter_AD.
! Verifies adjoint consistency with the tangent-linear model,
!   <H dx, y> == <dx, H^T y>,
! (1) with y = H dx over all AtmOptics outputs INCLUDING Absorption_Optical_Depth
!     (the original test uses Optical_Depth, Single_Scatter_Albedo,
!     Asymmetry_Factor and Phase_Coefficient), and
! (2) with y in the Absorption_Optical_Depth slot only, which isolates the
!     AAOD adjoint terms; the resulting RH adjoint must be non-zero.
! Same sensor (modis_aqua, channel 1), aerosol (GOCART-GEOS5 sea salt 1) and
! perturbation (RH + concentration) as the original test; RH is 0.53 for
! consistency with the TL and K counterparts (which need an off-node RH for
! their finite differences).  Backscat_Coefficient is left out of (1) as in
! the original: kb = 0 in this LUT, and the pre-existing tangent-linear
! accumulates it twice when scattering is on while the adjoint does so once.
! A dot product only proves TL/AD consistency; TL correctness is checked by
! test_AerosolScatter_AAOD_TL, and the clamp branches / scattering-off path
! (not active here: w = 0.16, ke > 0) by the k_matrix regression test test_AAOD.
!

PROGRAM test_AerosolScatter_AAOD_AD

  USE Type_Kinds,             ONLY: fp
  USE Message_Handler,        ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,        ONLY: ZERO, ONE
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
                                    CRTM_Compute_AerosolScatter_TL, &
                                    CRTM_Compute_AerosolScatter_AD
  USE ASvar_Define,           ONLY: ASvar_type, &
                                    ASvar_Create, ASvar_Destroy
  USE CRTM_ChannelInfo_Define, ONLY: CRTM_ChannelInfo_type

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_AAOD_AD'
  INTEGER :: Error_Status
  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  CHARACTER(256) :: File_Path = './testinput/'

  TYPE(CRTM_Atmosphere_type) :: Atm, Atm_TL, Atm_AD
  TYPE(CRTM_AtmOptics_type)  :: AScat, AScat_TL, AScat_AD
  TYPE(ASvar_type)           :: ASvar
  INTEGER :: SensorIndex = 1
  INTEGER :: ChannelIndex = 1
  INTEGER :: n_Layers = 1
  INTEGER :: n_Aerosols = 1
  INTEGER :: n_Legendre = 4
  INTEGER :: n_Phase = 1
  REAL(fp), PARAMETER :: DOT_TOL = 1.0e-12_fp

  REAL(fp) :: sum_tl, sum_ad

  PRINT *, 'Starting test_AerosolScatter_AAOD_AD...'

  ! Load Coefficients
  Error_Status = CRTM_SpcCoeff_Load( (/Sensor_Id/), File_Path=File_Path )
  IF ( Error_Status /= SUCCESS ) STOP 1

  Error_Status = CRTM_AerosolCoeff_Load('GOCART-GEOS5', 'AerosolCoeff.GOCART-GEOS5.bin', File_Path=File_Path)
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Create structures
  CALL CRTM_Atmosphere_Create( Atm, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_TL, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_AD, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_AtmOptics_Create( AScat, n_Layers, n_Legendre, n_Phase )
  CALL CRTM_AtmOptics_Create( AScat_TL, n_Layers, n_Legendre, n_Phase )
  CALL CRTM_AtmOptics_Create( AScat_AD, n_Layers, n_Legendre, n_Phase )
  CALL ASvar_Create( ASvar, n_Legendre, n_Phase, n_Layers, n_Aerosols )

  ! Populate Input
  Atm%Absorber_ID(1) = 1
  Atm%Aerosol(1)%Type = 6 ! Sea Salt
  Atm%Aerosol(1)%Concentration(1) = 1.0e-4_fp
  Atm%Aerosol(1)%Effective_Radius(1) = 1.0_fp
  Atm%Relative_Humidity(1) = 0.53_fp

  ! Populate TL Input: Perturb Relative Humidity and Concentration
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  Atm_TL%Relative_Humidity(1) = 0.1_fp
  Atm_TL%Aerosol(1)%Concentration(1) = 1.0e-5_fp

  AScat%Include_Scattering = .TRUE.
  AScat_TL%Include_Scattering = .TRUE.
  AScat_AD%Include_Scattering = .TRUE.

  ! Forward Call
  CALL CRTM_AtmOptics_Zero( AScat )
  Error_Status = CRTM_Compute_AerosolScatter( Atm, SensorIndex, ChannelIndex, AScat, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1
  PRINT *, 'Forward: OD, AAOD: ', AScat%Optical_Depth(1), AScat%Absorption_Optical_Depth(1)
  IF ( .NOT. ( AScat%Absorption_Optical_Depth(1) > ZERO ) ) THEN
    PRINT *, 'FAIL: forward AAOD is zero, test is invalid.'
    STOP 1
  END IF

  ! TL Call
  CALL CRTM_AtmOptics_Zero( AScat_TL )
  Error_Status = CRTM_Compute_AerosolScatter_TL( Atm, AScat, Atm_TL, SensorIndex, ChannelIndex, AScat_TL, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1
  IF ( .NOT. ( ABS(AScat_TL%Absorption_Optical_Depth(1)) >= 1.0e-12_fp ) ) THEN
    PRINT *, 'FAIL: AAOD tangent-linear is zero, test is invalid.'
    STOP 1
  END IF

  ! --------------------------------------------------------------------------
  ! (1) Adjoint Input: Set AScat_AD = AScat_TL, including Absorption_Optical_Depth
  ! --------------------------------------------------------------------------
  CALL CRTM_AtmOptics_Zero( AScat_AD )
  AScat_AD%Optical_Depth            = AScat_TL%Optical_Depth
  AScat_AD%Absorption_Optical_Depth = AScat_TL%Absorption_Optical_Depth
  AScat_AD%Single_Scatter_Albedo    = AScat_TL%Single_Scatter_Albedo
  AScat_AD%Asymmetry_Factor         = AScat_TL%Asymmetry_Factor
  AScat_AD%Phase_Coefficient        = AScat_TL%Phase_Coefficient

  ! AD Call
  CALL CRTM_Atmosphere_Zero( Atm_AD )
  Error_Status = CRTM_Compute_AerosolScatter_AD( Atm, AScat, AScat_AD, SensorIndex, ChannelIndex, Atm_AD, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Adjoint Test: sum(AScat_TL * AScat_AD_input) == sum(Atm_TL * Atm_AD_output)
  sum_tl = SUM(AScat_TL%Optical_Depth**2) + &
           SUM(AScat_TL%Absorption_Optical_Depth**2) + &
           SUM(AScat_TL%Single_Scatter_Albedo**2) + &
           SUM(AScat_TL%Asymmetry_Factor**2) + &
           SUM(AScat_TL%Phase_Coefficient**2)
  sum_ad = Dot_Atm( Atm_TL, Atm_AD )
  PRINT *, '(1) all outputs incl. AAOD'
  PRINT *, 'LHS (TL sum): ', sum_tl
  PRINT *, 'RHS (AD sum): ', sum_ad
  PRINT *, 'Diff:         ', ABS(sum_tl - sum_ad)
  IF ( .NOT. ( sum_tl > ZERO ) ) THEN
    PRINT *, 'FAIL: LHS is ZERO, test is invalid.'
    STOP 1
  END IF
  PRINT *, 'Rel. Diff:    ', ABS(sum_tl - sum_ad) / sum_tl
  IF ( .NOT. ( ABS(sum_tl - sum_ad) / sum_tl <= DOT_TOL ) ) THEN
     PRINT *, 'FAIL: Adjoint test (1) failed.'
     STOP 1
  END IF
  PRINT *, 'SUCCESS: Adjoint test (1) passed.'

  ! --------------------------------------------------------------------------
  ! (2) Adjoint input in the Absorption_Optical_Depth slot only
  ! --------------------------------------------------------------------------
  CALL CRTM_AtmOptics_Zero( AScat_AD )
  AScat_AD%Absorption_Optical_Depth = AScat_TL%Absorption_Optical_Depth
  CALL CRTM_Atmosphere_Zero( Atm_AD )
  Error_Status = CRTM_Compute_AerosolScatter_AD( Atm, AScat, AScat_AD, SensorIndex, ChannelIndex, Atm_AD, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1
  sum_tl = SUM(AScat_TL%Absorption_Optical_Depth**2)
  sum_ad = Dot_Atm( Atm_TL, Atm_AD )
  PRINT *, '(2) AAOD only'
  PRINT *, 'LHS (TL sum): ', sum_tl
  PRINT *, 'RHS (AD sum): ', sum_ad
  PRINT *, 'Rel. Diff:    ', ABS(sum_tl - sum_ad) / sum_tl
  PRINT *, 'RH adjoint:   ', Atm_AD%Relative_Humidity(1), '  Concentration adjoint: ', Atm_AD%Aerosol(1)%Concentration(1)
  IF ( .NOT. ( ABS(sum_tl - sum_ad) / sum_tl <= DOT_TOL ) ) THEN
     PRINT *, 'FAIL: Adjoint test (2) failed.'
     STOP 1
  END IF
  IF ( .NOT. ( ABS(Atm_AD%Relative_Humidity(1)) > ZERO .AND. ABS(Atm_AD%Aerosol(1)%Concentration(1)) > ZERO ) ) THEN
     PRINT *, 'FAIL: the AAOD adjoint gives no RH or concentration sensitivity.'
     STOP 1
  END IF
  PRINT *, 'SUCCESS: Adjoint test (2) passed.'

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_AD )
  CALL CRTM_AtmOptics_Destroy( AScat )
  CALL CRTM_AtmOptics_Destroy( AScat_TL )
  CALL CRTM_AtmOptics_Destroy( AScat_AD )
  CALL ASvar_Destroy( ASvar )
  Error_Status = CRTM_SpcCoeff_Destroy()
  Error_Status = CRTM_AerosolCoeff_Destroy()

  PRINT *, 'SUCCESS: test_AerosolScatter_AAOD_AD passed.'

CONTAINS

  ! <dx, x_AD> over the atmosphere variables the aerosol operator depends on
  FUNCTION Dot_Atm( a_TL, a_AD ) RESULT( s )
    TYPE(CRTM_Atmosphere_type), INTENT(IN) :: a_TL, a_AD
    REAL(fp) :: s
    s = SUM(a_TL%Relative_Humidity * a_AD%Relative_Humidity) + &
        SUM(a_TL%Aerosol(1)%Concentration * a_AD%Aerosol(1)%Concentration) + &
        SUM(a_TL%Aerosol(1)%Effective_Radius * a_AD%Aerosol(1)%Effective_Radius) + &
        SUM(a_TL%Aerosol(1)%Effective_Variance * a_AD%Aerosol(1)%Effective_Variance)
  END FUNCTION Dot_Atm

END PROGRAM test_AerosolScatter_AAOD_AD
