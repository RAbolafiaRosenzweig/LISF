!BOP
! !ROUTINE: NoahMP50_convert_to_GEOSexport
! !DESCRIPTION:
!  Converts native Noah-MP state into GEOS/CLSM-ready export variables.
!  Invoked only on the GEOS-coupled path, after NoahMP50_main.
!EOP
subroutine NoahMP50_convert_to_GEOSexport(n)
! !USES:
  use LIS_coreMod
  use LIS_logMod
  use NoahMP50_lsmMod
  use NoahMP50_geosCplMod
  use NoahmpIOVarType        ! for NoahmpIO%dzs(:)

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
! !LOCAL:
  integer :: t, i, nrz
  real    :: fveg, snc, cumdepth, rzdenom, rznumer
  real, parameter :: EPS     = 1.0e-6
  real, parameter :: RZ_TGT  = 1.0     ! root-zone target depth [m]

  do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     associate( nm => NoahMP50_struc(n)%noahmp50(t), &
                gx => NoahMP50_geoscpl(n)%tile(t) )

        fveg = nm%fveg

        ! --- total runoff ---
        gx%totalrunoff = nm%sfcrunoff + nm%udrrunoff

        ! --- sensible heat (SHOUT) ---
        gx%qh_shout = fveg*nm%shg + (1.0-fveg)*nm%shb + nm%shc

        ! --- exchange coefficient (guard fveg=0) ---
        if (fveg > EPS) then
           gx%ch_f = fveg*nm%chv + (1.0-fveg)*nm%chb
        else
           gx%ch_f = nm%chb
        endif

        ! --- root-zone soil moisture: depth-weighted mean over the top
        !     layers whose cumulative thickness first reaches/exceeds 1 m.
        !     Reads thicknesses from NoahmpIO%dzs so it tracks any retuning. ---
        cumdepth = 0.0 ; rznumer = 0.0 ; rzdenom = 0.0 ; nrz = 0
        do i = 1, size(NoahmpIO%dzs)
           rznumer  = rznumer + NoahmpIO%dzs(i)*nm%smc(i)
           rzdenom  = rzdenom + NoahmpIO%dzs(i)
           cumdepth = cumdepth + NoahmpIO%dzs(i)
           nrz = i
           if (cumdepth >= RZ_TGT - EPS) exit    ! include boundary-crossing layer
        enddo
        if (rzdenom > EPS) then
           gx%soilmoist_rz = rznumer / rzdenom
        else
           gx%soilmoist_rz = LIS_rc%udef
        endif

        ! --- 2m specific humidity (fveg-weighted) ---
        gx%qair_t2m = fveg*nm%q2mv + (1.0-fveg)*nm%q2mb

        ! --- 2m air temperature (fveg-weighted) ---
        gx%tair_t2m = fveg*nm%t2mv + (1.0-fveg)*nm%t2mb

        ! --- canopy water content ---
        gx%canopint = nm%canliq + nm%canice

        ! --- snow depth, snow-covered-area form (guard snowc=0) ---
        snc = nm%snowc
        if (snc > EPS) then
           gx%snowdepth_sca = nm%snowh / snc
        else
           gx%snowdepth_sca = nm%snowh
        endif

        ! --- roughness (single source; downstream feeds both Z0 and Z0H) ---
        gx%roughness_f = nm%z0

        ! --- band-resolved albedo: direct SELECT from harvested native
        !     arrays (1=vis, 2=nir). Requires nm%albd(:)/nm%albi(:) to be
        !     populated in NoahMP50_main (see note). ---
        gx%albedo_dir(1) = nm%albsfcdir(1)   ! AlbedoVisDir
        gx%albedo_dir(2) = nm%albsfcdir(2)   ! AlbedoNIRDir
        gx%albedo_dif(1) = nm%albsfcdif(1)   ! AlbedoVisDif
        gx%albedo_dif(2) = nm%albsfcdif(2)   ! AlbedoNIRDif

     end associate
  enddo

end subroutine NoahMP50_convert_to_GEOSexport