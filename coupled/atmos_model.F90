
module atmos_model_mod

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
                           atmosphere_domain,     &
                           atmosphere_boundary,   &
                           get_atmosphere_axes,   &
                           surf_diff_type

use time_manager_mod, only:  time_type, operator(+), get_time

use mpp_domains_mod,  only: domain2d

use    utilities_mod, only:  file_exist, open_file, check_nml_error,  &
                             error_mesg, close_file, get_my_pe,       &
                             read_data, write_data, FATAL

use  diag_integral_mod, only: diag_integral_init, diag_integral_end, &
                              diag_integral_output


!-----------------------------------------------------------------------

implicit none
private

public  update_atmos_model_down, update_atmos_model_up,   &
        atmos_model_init, atmos_model_end, atmos_boundary_data_type

!-----------------------------------------------------------------------

 type atmos_boundary_data_type
     type (domain2d)               :: domain
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

!quantities going from land+ice to atmos
type, public :: land_ice_atmos_boundary_type
!defined here but declared by coupler_main, allocated by flux_exchange_init
   real, dimension(:,:), pointer :: t, albedo, land_frac
   real, dimension(:,:), pointer :: dt_t, dt_q
   real, dimension(:,:), pointer :: u_flux, v_flux, dtaudv, u_star, b_star, rough_mom
   real, dimension(:,:,:), pointer :: data !collective field for "named" fields above
   integer :: xtype             !REGRID, REDIST or DIRECT
end type land_ice_atmos_boundary_type
!quantities going from land alone to atmos (none at present)
type, public :: land_atmos_boundary_type
!   private
   real, dimension(:,:), pointer :: data
end type land_atmos_boundary_type
!quantities going from ice alone to atmos (none at present)
type, public :: ice_atmos_boundary_type
!   private
   real, dimension(:,:), pointer :: data
end type ice_atmos_boundary_type

!-----------------------------------------------------------------------

character(len=256) :: version = '$Id: atmos_model.F90,v 1.1 2002/01/30 18:53:47 fms Exp $'
character(len=256) :: tag = '$Name: galway $'

character(len=80) :: restart_format = 'atmos_coupled_mod restart format 01'

!-----------------------------------------------------------------------


contains

!#######################################################################

subroutine update_atmos_model_down( Surface_boundary, Atmos )
!Balaji: this routine is now a driver for any atmosphere
!this comment only applies to the current atmosphere.f90 and should be moved there
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
!        land_frac = fraction amount of land in a grid box
!        dtau_dv   = derivative of wind stress w.r.t. the 
!                    lowest level wind speed
!
! inout: tau_x     = zonal wind stress
!        tau_y     = meridional wind stress
!
!-----------------------------------------------------------------------

  type(land_ice_atmos_boundary_type), intent(inout) :: Surface_boundary
type (atmos_boundary_data_type), intent(inout) :: Atmos
!real,  dimension(:,:),  intent(in)    :: t_surf, albedo, rough_mom, &
!                                         u_star, b_star, land_frac, &
!                                         dtau_dv
!real,  dimension(:,:),  intent(inout) :: tau_x,  tau_y
                                      
!-----------------------------------------------------------------------

    call atmosphere_down (Atmos%Time, Surface_boundary%land_frac,        &
                          Surface_boundary%t,  Surface_boundary%albedo, Surface_boundary%rough_mom,   &
                          Surface_boundary%u_star,  Surface_boundary%b_star,              &
                          Surface_boundary%dtaudv, Surface_boundary%u_flux,  Surface_boundary%v_flux,       &
                          Atmos%gust, Atmos%coszen,     &
                          Atmos%flux_sw, Atmos%flux_lw, &
                          Atmos%Surf_diff               )

!-----------------------------------------------------------------------

 end subroutine update_atmos_model_down

!#######################################################################

 subroutine update_atmos_model_up( Surface_boundary, Atmos )

!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs upward vertical diffusion of heat/moisture and
!    moisture processes
!
!   in: land_frac = fraction amount of land in a grid box
!       dt_t_bot  = temperature tendency at the lowest level
!       dt_q_bot  = specific humidity tendency at the lowest level
!
!-----------------------------------------------------------------------

  type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
type (atmos_boundary_data_type), intent(inout) :: Atmos
!real,  dimension(:,:),  intent(in)    :: land_frac, dt_t_bot, dt_q_bot
                                      
!-----------------------------------------------------------------------


    Atmos%Surf_diff%delta_t = Surface_boundary%dt_t
    Atmos%Surf_diff%delta_q = Surface_boundary%dt_q

    call atmosphere_up (Atmos%Time,  Surface_boundary%land_frac, Atmos%Surf_diff, &
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

end subroutine update_atmos_model_up

!#######################################################################

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

type (atmos_boundary_data_type), intent(inout) :: Atmos
type (time_type), intent(in) :: Time_init, Time, Time_step

  integer :: unit, log_unit
  integer :: mlon, mlat, nlon, nlat, sec, day, ipts, jpts, dt, dto
  character(len=80) :: control
!-----------------------------------------------------------------------

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
    call atmosphere_domain     (Atmos%domain)

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

!-----------------------------------------------------------------------
!---- print version number to logfile ----

   log_unit = open_file ('logfile.out', action='append')
   if (get_my_pe() == 0) &
       write (log_unit,'(/,80("="),/(a))') trim(version), trim(tag)

!------ read initial state for several atmospheric fields ------

    if (file_exist('INPUT/atmos_coupled.res')) then
        unit = open_file ('INPUT/atmos_coupled.res',  &
                          form='native', action='read')
       !--- check version number (format) of restart file ---
        read  (unit) control
        if (trim(control) /= trim(restart_format)) call error_mesg &
             ('coupled_atmos_init', 'invalid restart format', FATAL)
       !--- check resolution and time step ---
        read  (unit) ipts,jpts,dto
        if (ipts /= mlon .or. jpts /= mlat) call error_mesg &
        ('coupled_atmos_init', 'incorrect resolution on restart file', FATAL)
       !--- read data ---
        call read_data ( unit, Atmos % lprec )
        call read_data ( unit, Atmos % fprec )
        call read_data ( unit, Atmos % gust  )
        call close_file (unit)

       !---- if the time step has changed then convert ----
       !        tendency to conserve mass of water
        call get_time (Atmos % Time_step, sec, day)
        dt = sec + 86400*day  ! integer seconds
        if (dto /= dt) then
            Atmos % lprec = Atmos % lprec * real(dto)/real(dt)
            Atmos % fprec = Atmos % fprec * real(dto)/real(dt)
            if (get_my_pe() == 0) write (log_unit,50)
 50         format (/,'The model time step changed .... &
                      &modifying precipitation tendencies')
        endif

    else
        Atmos % lprec = 0.0
        Atmos % fprec = 0.0
        Atmos % gust  = 1.0
    endif

!------ close log file unit -----

   call close_file (log_unit)

!------ initialize global integral package ------

    call diag_integral_init (Atmos % Time_init, Atmos % Time,  &
                             Atmos % lon_bnd,   Atmos % lat_bnd)

!-----------------------------------------------------------------------

end subroutine atmos_model_init

!#######################################################################

subroutine atmos_model_end (Atmos)

type (atmos_boundary_data_type), intent(in) :: Atmos
integer :: unit, sec, day, dt

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----
                                              
  call atmosphere_end (Atmos % Time)

!------ global integrals ------

  call diag_integral_end (Atmos % Time)

!---- compute integer time step (in seconds) ----
  call get_time (Atmos % Time_step, sec, day)
  dt = sec + 86400*day

!------ write several atmospheric fields ------
!        also resolution and time step

  unit = open_file ('RESTART/atmos_coupled.res',   &
                    form='native', action='write')
  if (get_my_pe() == 0) then
      write (unit) restart_format
      write (unit) size(Atmos%glon_bnd)-1, size(Atmos%glat_bnd)-1, dt
  endif
  call write_data ( unit, Atmos % lprec )
  call write_data ( unit, Atmos % fprec )
  call write_data ( unit, Atmos % gust  )
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

end subroutine atmos_model_end

!#######################################################################

end module atmos_model_mod

