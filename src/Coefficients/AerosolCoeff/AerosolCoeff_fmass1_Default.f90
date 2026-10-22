!
! AerosolCoeff_fmass1_Default
!
! Module providing default values of fmass1, the fraction of the aerosol
! mass in particles with diameter < 1 um, for the GOCART-GEOS5 aerosol
! types. fmass1 is used to compute the fine-mode AOD and AAOD.
!
! The values are taken from the GMAO aerosol optics tables
!   optics_DU.v2.2.1.legacy.nc4  (dust,           5 bins)
!   optics_SS.v2.2.1.legacy.nc4  (sea salt,       5 bins)
!   optics_NI.v2.2.1.legacy.nc4  (nitrate,        3 bins)
!   optics_OC.v2.2.1.legacy.nc4  (organic carbon, hydrophobic OC1 and hydrophilic OC2)
!   optics_BC.v2.2.1.legacy.nc4  (black carbon,   hydrophobic BC1 and hydrophilic BC2)
!   optics_SU.v2.2.1.legacy.nc4  (sulfate,        a single entry)
! on their 36-point relative humidity grid. The single sulfate entry is used
! for both GOCART-GEOS5 sulfate types. Values that differ from 1 only by
! rounding (within 1e-12) are stored as exactly 1.
!
! The defaults are used until the AerosolCoeff files provide fmass1.
!
!
! CREATION HISTORY:
!       Written by:     Mariusz Pagowski, 2026/10/22
!                       mariusz.pagowski@colorado.edu
!

MODULE AerosolCoeff_fmass1_Default

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use
  USE Type_Kinds         , ONLY: Double
  USE String_Utility     , ONLY: StrUpCase
  USE AerosolCoeff_Define, ONLY: AerosolCoeff_type, &
                                 AerosolCoeff_Associated
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: AerosolCoeff_Set_fmass1_Default


  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(Double), PARAMETER :: ZERO = 0.0_Double
  REAL(Double), PARAMETER :: ONE  = 1.0_Double
  ! Relative humidity within RH_TOL of a table node takes the node value.
  ! The AerosolCoeff RH grid is stored in single precision.
  REAL(Double), PARAMETER :: RH_TOL = 1.0e-6_Double
  ! Table dimensions
  INTEGER, PARAMETER :: N_RH    = 36
  INTEGER, PARAMETER :: N_TYPES = 19
  ! Relative humidity grid of the GMAO tables (fraction)
  REAL(Double), PARAMETER :: FMASS1_RH(N_RH) = (/ &
    0.00_Double, 0.05_Double, 0.10_Double, 0.15_Double, 0.20_Double, 0.25_Double, 0.30_Double, &
    0.35_Double, 0.40_Double, 0.45_Double, 0.50_Double, 0.55_Double, 0.60_Double, 0.65_Double, &
    0.70_Double, 0.75_Double, 0.80_Double, 0.81_Double, 0.82_Double, 0.83_Double, 0.84_Double, &
    0.85_Double, 0.86_Double, 0.87_Double, 0.88_Double, 0.89_Double, 0.90_Double, 0.91_Double, &
    0.92_Double, 0.93_Double, 0.94_Double, 0.95_Double, 0.96_Double, 0.97_Double, 0.98_Double, &
    0.99_Double /)
  ! GOCART-GEOS5 aerosol type names, in AerosolCoeff.GOCART-GEOS5 order
  CHARACTER(*), PARAMETER :: TYPE_NAME(N_TYPES) = (/ CHARACTER(16) :: &
    'Dust 1', 'Dust 2', 'Dust 3', 'Dust 4', 'Dust 5', 'Sea salt 1', 'Sea salt 2', 'Sea salt 3', &
    'Sea salt 4', 'Sea salt 5', 'Organic carbon 1', 'Organic carbon 2', 'Black carbon 1', &
    'Black carbon 2', 'Sulfate 1', 'Sulfate 2', 'Nitrate 1', 'Nitrate 2', 'Nitrate 3' /)
  ! fmass1(RH) for each type
  REAL(Double), PARAMETER :: FMASS1(N_RH, N_TYPES) = RESHAPE( (/ &
    !  1) Dust 1: optics_DU entry 1
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    0.2069999818951323_Double, 0.2069999818951323_Double, 0.2069999818951323_Double, &
    !  2) Dust 2: optics_DU entry 2
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    !  3) Dust 3: optics_DU entry 3
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    !  4) Dust 4: optics_DU entry 4
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    !  5) Dust 5: optics_DU entry 5
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    !  6) Sea salt 1: optics_SS entry 1
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    !  7) Sea salt 2: optics_SS entry 2
    ONE, 0.8852098117395629_Double, 0.8589502323598758_Double, 0.8352797797810071_Double, &
    0.8140995598217247_Double, 0.7935044127734575_Double, 0.7734649480445024_Double, &
    0.7506819498071616_Double, 0.7286332721060674_Double, 0.7057394570915836_Double, &
    0.6805690320159618_Double, 0.6534120418466541_Double, 0.6231167650226753_Double, &
    0.5901694356949353_Double, 0.5537416270552008_Double, 0.5108974827683728_Double, &
    0.46066838679702593_Double, 0.44994590207757723_Double, 0.43831781486930943_Double, &
    0.4259323411976523_Double, 0.4126619528741347_Double, 0.3996446845038238_Double, &
    0.3850251638471669_Double, 0.3706087953027229_Double, 0.35464725642990946_Double, &
    0.3373855874944313_Double, 0.31957511407339734_Double, 0.2998724873513007_Double, &
    0.2784060426916544_Double, 0.2552763617050063_Double, 0.22961840144861645_Double, &
    0.20059365007168792_Double, 0.20059365007168792_Double, 0.20059365007168792_Double, &
    0.20059365007168792_Double, 0.20059365007168792_Double, &
    !  8) Sea salt 3: optics_SS entry 3
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    !  9) Sea salt 4: optics_SS entry 4
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ! 10) Sea salt 5: optics_SS entry 5
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, &
    ! 11) Organic carbon 1: optics_OC entry 1, hydrophobic species
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ! 12) Organic carbon 2: optics_OC entry 2, hydrophilic species
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, 0.9989425521273754_Double, &
    0.9903348720972772_Double, 0.9792306113981765_Double, 0.9677841123308375_Double, &
    0.9522664283272507_Double, 0.9322632644383826_Double, 0.9074088540006063_Double, &
    0.869730565071446_Double, 0.8102859980167078_Double, &
    ! 13) Black carbon 1: optics_BC entry 1, hydrophobic species
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ! 14) Black carbon 2: optics_BC entry 2, hydrophilic species
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    0.9995017838017553_Double, 0.9970127444337061_Double, &
    ! 15) Sulfate 1: optics_SU entry 1
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ONE, 0.9972838757468755_Double, 0.9918411508553346_Double, 0.9837050427082743_Double, &
    0.9701365725662469_Double, 0.96205152527005_Double, 0.9485349260026654_Double, &
    0.9430690112987112_Double, 0.9295879129671007_Double, 0.9161531832515726_Double, &
    0.8973234069648087_Double, 0.884019364340393_Double, 0.8653747504124731_Double, &
    0.841482730378269_Double, 0.8231544620663611_Double, 0.7944365145069403_Double, &
    0.7609785362870125_Double, 0.7031204303629337_Double, &
    ! 16) Sulfate 2: optics_SU entry 1
    ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE, &
    ONE, 0.9972838757468755_Double, 0.9918411508553346_Double, 0.9837050427082743_Double, &
    0.9701365725662469_Double, 0.96205152527005_Double, 0.9485349260026654_Double, &
    0.9430690112987112_Double, 0.9295879129671007_Double, 0.9161531832515726_Double, &
    0.8973234069648087_Double, 0.884019364340393_Double, 0.8653747504124731_Double, &
    0.841482730378269_Double, 0.8231544620663611_Double, 0.7944365145069403_Double, &
    0.7609785362870125_Double, 0.7031204303629337_Double, &
    ! 17) Nitrate 1: optics_NI entry 1
    0.9015300220612703_Double, 0.8754685600732488_Double, 0.8641869769533458_Double, &
    0.8526718231813379_Double, 0.8409572683640271_Double, 0.8290756407630542_Double, &
    0.8200731618073134_Double, 0.8079712574939074_Double, 0.7957814170885158_Double, &
    0.7835288903985432_Double, 0.7712371068514322_Double, 0.7589277333440235_Double, &
    0.7496962367250777_Double, 0.7374031644886134_Double, 0.7251435786076028_Double, &
    0.7098899843362019_Double, 0.6947407853078148_Double, 0.6917256784580202_Double, &
    0.6887159163264547_Double, 0.6857116754717912_Double, 0.6827131279221801_Double, &
    0.6797204409865066_Double, 0.6737532987373842_Double, 0.6707791566942866_Double, &
    0.6648504832872086_Double, 0.661896240702629_Double, 0.656008636000346_Double, &
    0.6501497499582368_Double, 0.6414174434384272_Double, 0.6356346096176347_Double, &
    0.6270205678468646_Double, 0.6156516804851861_Double, 0.6072154947346331_Double, &
    0.593333758578266_Double, 0.5769804820643192_Double, 0.5478982345626923_Double, &
    ! 18) Nitrate 2: optics_NI entry 2
    0.008619187942956411_Double, 0.00587985569176064_Double, 0.005047226871907141_Double, &
    0.0043459415619350175_Double, 0.0037531849875170288_Double, 0.0032504624148518846_Double, &
    0.00292333787989101_Double, 0.002543613815086867_Double, 0.0022187906317615706_Double, &
    0.0019401290567935047_Double, 0.0017004069445423283_Double, 0.001493634423629799_Double, &
    0.0013571470026539404_Double, 0.0011965417196394261_Double, 0.0010570670183908243_Double, &
    0.0009078252612674424_Double, 0.0007819185556116189_Double, 0.0007591712504104095_Double, &
    0.0007371668744851479_Double, 0.0007158784511451548_Double, 0.0006952800834483222_Double, &
    0.000675346904487191_Double, 0.0006373815615488499_Double, 0.0006193044338029723_Double, &
    0.0005848554083575732_Double, 0.0005684436299852488_Double, 0.000537151580957555_Double, &
    0.0005077845942972396_Double, 0.00046705956647962664_Double, 0.0004419529598716734_Double, &
    0.00040708785530401467_Double, 0.00036531136520868025_Double, 0.0003371344239330435_Double, &
    0.0002954324985553737_Double, 0.0002528416545420044_Double, 0.00019148497150574668_Double, &
    ! 19) Nitrate 3: optics_NI entry 3
    8.48582902630505E-06_Double, 4.45512082297915E-06_Double, 3.453946184083825E-06_Double, &
    2.695381880216767E-06_Double, 2.1164888036784053E-06_Double, 1.6716802499167887E-06_Double, &
    1.4057088889281014E-06_Double, 1.1208976183551855E-06_Double, 8.983250807293629E-07_Double, &
    7.234203648335367E-07_Double, 5.852492198010627E-07_Double, 4.7554212723594597E-07_Double, &
    4.0810345932178105E-07_Double, 3.3399511446730565E-07_Double, &
    2.7440137443558985E-07_Double, 2.157469881145747E-07_Double, 1.7056273925616346E-07_Double, &
    1.628349916539681E-07_Double, 1.5548943707059246E-07_Double, 1.4850556899705226E-07_Double, &
    1.418640680344751E-07_Double, 1.3554588689086846E-07_Double, 1.2381609568250826E-07_Double, &
    1.1837216314253048E-07_Double, 1.0825418416687715E-07_Double, &
    1.0355372780347244E-07_Double, 9.480930131677434E-08_Double, 8.686709195537939E-08_Double, &
    7.628690984696194E-08_Double, 7.002075902381564E-08_Double, 6.165260893615018E-08_Double, &
    5.215193668604583E-08_Double, 4.6078772605863654E-08_Double, 3.7605028245175646E-08_Double, &
    2.9614277376815205E-08_Double, 1.936764165554449E-08_Double /), &
    (/ N_RH, N_TYPES /) )


CONTAINS


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       AerosolCoeff_Set_fmass1_Default
!
! PURPOSE:
!       Subroutine to fill the fmass1 component of an AerosolCoeff object
!       with the default values, for coefficient files that do not provide
!       fmass1.
!
!       For the GOCART-GEOS5 scheme each aerosol type is matched by name to
!       the default table, which is linearly interpolated onto the
!       AerosolCoeff RH grid. Types not in the table, and all types for
!       other schemes, are set to one.
!
! CALLING SEQUENCE:
!       CALL AerosolCoeff_Set_fmass1_Default( AerosolCoeff )
!
! OBJECTS:
!       AerosolCoeff:  AerosolCoeff object with the Scheme, Type_Name and
!                      RH components defined.
!                      UNITS:      N/A
!                      TYPE:       TYPE(AerosolCoeff_type)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE AerosolCoeff_Set_fmass1_Default( AerosolCoeff )
    TYPE(AerosolCoeff_type), INTENT(IN OUT) :: AerosolCoeff
    CHARACTER(LEN(AerosolCoeff%Type_Name)) :: Name
    INTEGER :: it, n, ih, i

    IF ( .NOT. AerosolCoeff_Associated( AerosolCoeff ) ) RETURN
    AerosolCoeff%fmass1 = ONE
    IF ( TRIM(AerosolCoeff%Scheme) /= 'GOCART-GEOS5' ) RETURN
    IF ( SIZE(AerosolCoeff%RH) /= AerosolCoeff%n_RH ) RETURN

    Type_Loop: DO it = 1, AerosolCoeff%n_Types
      ! ...Type names in the coefficient files are NUL padded
      Name = AerosolCoeff%Type_Name(it)
      DO i = 1, LEN(Name)
        IF ( IACHAR(Name(i:i)) < 32 ) Name(i:i) = ' '
      END DO
      Name = StrUpCase(ADJUSTL(Name))
      DO n = 1, N_TYPES
        IF ( TRIM(Name) == StrUpCase(TYPE_NAME(n)) ) THEN
          DO ih = 1, AerosolCoeff%n_RH
            AerosolCoeff%fmass1(ih,it) = Linear_Interp( FMASS1_RH, FMASS1(:,n), &
                                                        AerosolCoeff%RH(ih) )
          END DO
          CYCLE Type_Loop
        END IF
      END DO
    END DO Type_Loop

  END SUBROUTINE AerosolCoeff_Set_fmass1_Default


!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

  ! Linear interpolation of y(x) to xi, held constant outside x.
  ! An xi within RH_TOL of a node takes the node value.
  PURE FUNCTION Linear_Interp( x, y, xi ) RESULT( yi )
    REAL(Double), INTENT(IN) :: x(:), y(:), xi
    REAL(Double) :: yi
    INTEGER :: i, n
    n = SIZE(x)
    i = MINLOC(ABS(x - xi), DIM=1)
    IF ( ABS(x(i) - xi) < RH_TOL ) THEN
      yi = y(i)
    ELSE IF ( xi <= x(1) ) THEN
      yi = y(1)
    ELSE IF ( xi >= x(n) ) THEN
      yi = y(n)
    ELSE
      DO i = 1, n-1
        IF ( xi <= x(i+1) ) EXIT
      END DO
      yi = y(i) + (y(i+1)-y(i)) * (xi-x(i)) / (x(i+1)-x(i))
    END IF
  END FUNCTION Linear_Interp

END MODULE AerosolCoeff_fmass1_Default
