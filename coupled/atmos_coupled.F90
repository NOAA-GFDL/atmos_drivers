
module atmos_coupled_mod

!-----------------------------------------------------------------------
!
!          driver for atmospheric gcm using an
!     atmospheric dynamical core and modular physics
!
!-----------------------------------------------------------------------

use atmosphere_mod, only:  atmosphere_up,         &
                           atmosphere_down,       &
                           atmosphere_init,       &
                           atmosphere_end,        &
                           get_bottom_mass,       &
                           get_bottom_wind,       &
                           atmosphere_resolution, &
                           atmosphere_boundary,   &
                           get_atmosphere_axes,   &
                           surf_diff_type

use time_manager_mod, only:  time_type, operator(+), get_time

use    utilities_mod, only:  file_exist, open_file, check_nml_error,  &
                             error_mesg, close_file, get_my_pe,       &
                             read_data, write_data

use  diag_integral_mod, only: diag_integral_init, diag_integral_end, &
                              diag_integral_output

!-----------------------------------------------------------------------

implicit none
private

public  update_atmos_coupled_down, update_atmos_coupled_up,   &
        atmos_coupled_init, atmos_coupled_end, atmos_boundary_data_type

!-----------------------------------------------------------------------

 type atmos_boundary_data_type
     integer                       :: axes(4)
     real, pointer, dimension(:)   :: glon_bnd, glat_bnd,  &
                                       lon_bnd,  lat_bnd
     real, pointer, dimension(:,:) :: t_bot, q_bot, z_bot, p_bot,  &
                                      u_bot, v_bot, p_surf, gust,  &
                                      coszen, flux_sw, flux_lw,    &
                                      lprec, fprec
     type (surf_diff_type)         :: Surf_diff
     type (time_type)              :: Time, Time_step, Time_init
 end type

!-----------------------------------------------------------------------

character(len=256) :: version = '$Id: atmos_coupled.F90,v 1.3 2000/11/22 14:33:18 fms Exp $'
character(len=256) :: tag = '$Name: damascus $'

!-----------------------------------------------------------------------


contains

!#######################################################################

subroutine update_atmos_coupled_down (Atmos, t_surf,  albedo, rough_mom,&
                                            u_star,  b_star, frac_land, &
                                            dtau_dv, tau_x,  tau_y      )

!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs radiation, damping, and vertical diffusion of momentum,
!    tracers, and downward heat/moisture
!
!   in:  t_surf    = surface temperature for radiation calculations
!        albedo    = surface albedo for radiation calculations
!        rough_mom = surface roughness (used for momentum
!        u_star    = friction velocity
!        b_star    = bouyancy scale
!        frac_land = fraction amount of land in a grid box
!        dtau_dv   = derivative of wind stress w.r.t. the 
!                    lowest level wind speed
!
! inout: tau_x     = zonal wind stress
!        tau_y     = meridional wind stress
!
!-----------------------------------------------------------------------

type (atmos_boundary_data_type), intent(inout) :: Atmos
real,  dimension(:,:),  intent(in)    :: t_surf, albedo, rough_mom, &
                                         u_star, b_star, frac_land, &
                                         dtau_dv
real,  dimension(:,:),  intent(inout) :: tau_x,  tau_y
                                      
!-----------------------------------------------------------------------

    call atmosphere_down (Atmos%Time, frac_land,        &
                          t_surf,  albedo, rough_mom,   &
                          u_star,  b_star,              &
                          dtau_dv, tau_x,  tau_y,       &
                          Atmos%gust, Atmos%coszen,     &
                          Atmos%flux_sw, Atmos%flux_lw, &
                          Atmos%Surf_diff               )

!-----------------------------------------------------------------------

 end subroutine update_atmos_coupled_down

!#######################################################################

 subroutine update_atmos_coupled_up (Atmos, frac_land, dt_t_bot, dt_q_bot)

!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs upward vertical diffusion of heat/moisture and
!    moisture processes
!
!   in: frac_land = fraction amount of land in a grid box
!       dt_t_bot  = temperature tendency at the lowest level
!       dt_q_bot  = specific humidity tendency at the lowest level
!
!-----------------------------------------------------------------------

type (atmos_boundary_data_type), intent(inout) :: Atmos
real,  dimension(:,:),  intent(in)    :: frac_land, dt_t_bot, dt_q_bot
                                      
!-----------------------------------------------------------------------

    Atmos%Surf_diff%delta_t_n = dt_t_bot
    Atmos%Surf_diff%delta_q_n = dt_q_bot

    call atmosphere_up (Atmos%Time,  frac_land, Atmos%Surf_diff, &
                        Atmos%lprec, Atmos%fprec)

!   --- advance time ---

    Atmos % Time = Atmos % Time + Atmos % Time_step


    call get_bottom_mass (Atmos % t_bot, Atmos % q_bot,  &
                          Atmos % p_bot, Atmos % z_bot,  &
                          Atmos % p_surf                 )

    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)


!------ global integrals ------

    call diag_integral_output (Atmos % Time)

!-----------------------------------------------------------------------

end subroutine update_atmos_coupled_up

!#######################################################################

subroutine atmos_coupled_init (Atmos, Time_init, Time, Time_step)

type (atmos_boundary_data_type), intent(inout) :: Atmos
type (time_type), intent(in) :: Time_init, Time, Time_step

  integer :: unit, mlon, mlat, nlon, nlat, sec, day
  real    :: dt
!-----------------------------------------------------------------------
!---- print version number to logfile ----

   unit = open_file ('logfile.out', action='append')
   if (get_my_pe() == 0) &
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
   call close_file (unit)

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step

!-----------------------------------------------------------------------
!  ----- initialize atmospheric model -----

    call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                          Atmos%Surf_diff )
                           
!-----------------------------------------------------------------------
!---- allocate space ----

    call atmosphere_resolution (mlon, mlat, global=.true.)
    call atmosphere_resolution (nlon, nlat, global=.false.)

    allocate ( Atmos % glon_bnd (mlon+1),    &
               Atmos % glat_bnd (mlat+1),    &
               Atmos %  lon_bnd (nlon+1),    &
               Atmos %  lat_bnd (nlat+1),    &
               Atmos % t_bot    (nlon,nlat), &
               Atmos % q_bot    (nlon,nlat), &
               Atmos % z_bot    (nlon,nlat), &
               Atmos % p_bot    (nlon,nlat), &
               Atmos % u_bot    (nlon,nlat), &
               Atmos % v_bot    (nlon,nlat), &
               Atmos % p_surf   (nlon,nlat), &
               Atmos % gust     (nlon,nlat), &
               Atmos % flux_sw  (nlon,nlat), &
               Atmos % flux_lw  (nlon,nlat), &
               Atmos % coszen   (nlon,nlat), &
               Atmos % lprec    (nlon,nlat), &
               Atmos % fprec    (nlon,nlat)  )

!-----------------------------------------------------------------------
!------ get initial state for dynamics -------

    call get_atmosphere_axes ( Atmos % axes )

    call atmosphere_boundary ( Atmos % glon_bnd, Atmos % glat_bnd, &
                               global=.true. )
    call atmosphere_boundary ( Atmos %  lon_bnd, Atmos %  lat_bnd, &
                               global=.false. )

    call get_bottom_mass (Atmos % t_bot, Atmos % q_bot,  &
                          Atmos % p_bot, Atmos % z_bot,  &
                          Atmos % p_surf                 )

    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)

!------ read initial state for several atmospheric fields ------

    if (file_exist('INPUT/atmos_coupled.res')) then
        unit = open_file ('INPUT/atmos_coupled.res',  &
                          form='native', action='read')
        call read_data ( unit, Atmos % lprec )
        call read_data ( unit, Atmos % fprec )
        call read_data ( unit, Atmos % gust  )
        call close_file (unit)

!       ---- convert precip mass to precip rate ----
        call get_time (Atmos % Time_step, sec, day)
        dt = float(sec + 86400*day)
        Atmos % lprec = Atmos % lprec / dt
        Atmos % fprec = Atmos % fprec / dt

    else
        Atmos % lprec = 0.0
        Atmos % fprec = 0.0
        Atmos % gust  = 1.0
    endif

!------ initialize global integral package ------

    call diag_integral_init (Atmos % Time_init, Atmos % Time,  &
                             Atmos % lon_bnd,   Atmos % lat_bnd)

!-----------------------------------------------------------------------

end subroutine atmos_coupled_init

!#######################################################################

subroutine atmos_coupled_end (Atmos)

type (atmos_boundary_data_type), intent(in) :: Atmos
integer :: unit, sec, day
real    :: dt

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----
                                              
  call atmosphere_end (Atmos % Time)

!------ global integrals ------

  call diag_integral_end (Atmos % Time)

!---- compute real time step ----
  call get_time (Atmos % Time_step, sec, day)
  dt = float(sec + 86400*day)

!------ write several atmospheric fields ------
!     (convert precip rate to precip mass)

  unit = open_file ('RESTART/atmos_coupled.res',   &
                    form='native', action='write')
  call write_data ( unit, Atmos % lprec * dt )
  call write_data ( unit, Atmos % fprec * dt )
  call write_data ( unit, Atmos % gust       )
  call close_file (unit)

!-------- deallocate space --------

  deallocate ( Atmos % glon_bnd , &
               Atmos % glat_bnd , &
               Atmos %  lon_bnd , &
               Atmos %  lat_bnd , &
               Atmos % t_bot    , &
               Atmos % q_bot    , &
               Atmos % z_bot    , &
               Atmos % p_bot    , &
               Atmos % u_bot    , &
               Atmos % v_bot    , &
               Atmos % p_surf   , &
               Atmos % gust     , &
               Atmos % flux_sw  , &
               Atmos % flux_lw  , &
               Atmos % coszen   , &
               Atmos % lprec    , &
               Atmos % fprec      )

!-----------------------------------------------------------------------

end subroutine atmos_coupled_end

!#######################################################################

end module atmos_coupled_mod

