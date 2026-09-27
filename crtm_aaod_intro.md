Currently CRTM provides forward, tangent linear (TL), and adjoint (AD)
code for AOD. We want to use the current AOD code to have
corresponding calculations for absorption AOD or AAOD.

The Forward Operator is for AOD
as in:
               aod(k,n,i) = aod(k,n,i) + bext_*qm(iq,k,i)  ! sum over thetracers 

bext is extinction coefficient, qm is (moist) mixing ratio

The Forward Operator is for AAOD
               aaod(k,n,i) = aaod(k,n,i) + (1.-ssa_)*bext_*qm(iq,k,i)  


ssa is single scattering albedo. It is defined in
/work/noaa/wrf-chem/pagowski/jedi/code/jedi-bundle/crtm/src/AtmScatter/AerosolScatter/ASvar_Define.f90

  TYPE :: ASvar_type
  ...
  REAL(fp), ALLOCATABLE :: w(:)

and passed in

/work/noaa/wrf-chem/pagowski/jedi/code/jedi-bundle/crtm/src/AtmScatter/CRTM_AerosolScatter.f90

  USE ASvar_Define, ONLY: ASvar_type, &
                          ASinterp_type


 FUNCTION CRTM_Compute_AerosolScatter
  ...

  TYPE(ASvar_type)          , INTENT(IN OUT) :: ASV

as ASV%w

Also if AerosolScatter%Include_Scattering = .FALSE.

ke  = mass extintion coefficient

is recalculated as

ke = ke * (ONE- w)

so it is becomes Absorption coefficient

ie. by replacing ke (mass extinction coefficient) with absorption
coefficient 

AAOD is calculated in place of AOD.

By adding optional logical argument AAOD in functions


  PUBLIC :: CRTM_AOD
  PUBLIC :: CRTM_AOD_TL
  PUBLIC :: CRTM_AOD_AD
  PUBLIC :: CRTM_AOD_K

in

/work/noaa/wrf-chem/pagowski/jedi/code/jedi-bundle/crtm/src/AtmScatter/CRTM_AOD_Module.f90

Then

  FUNCTION CRTM_AOD( &
    Atmosphere , &  ! Input, M
    ChannelInfo, &  ! Input, M
    RTSolution , &  ! Output, L x M
    Options    ) &  ! Optional input, M

becomes

  FUNCTION CRTM_AOD( &
    Atmosphere , &  ! Input, M
    ChannelInfo, &  ! Input, M
    RTSolution , &  ! Output, L x M
    Options    , &  ! Optional input, M
    AAOD)           ! Optional

    ....
    LOGICAL, OPTIONAL, INTENT(IN)     :: AAOD

and similar for CRTM_AOD_TL, CRTM_AOD_AD, CRTM_AOD_K

I want to set AerosolScatter%Include_Scattering = .FALSE.

and calculate forward, K-matrix, TL, and AD for AAOD.
In this case AAOD would be returned in the same locations/arrays
as AOD in CRTM_AOD, CRTM_AOD_TL, CRTM_AOD_AD, CRTM_AOD_K

Please check is possible and if it will work.

We will be using AeroC%Scheme == 'GOCART-GEOS5'
