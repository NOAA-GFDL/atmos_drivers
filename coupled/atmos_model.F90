
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

use          fms_mod, only:  file_exist, open_restart_file, error_mesg,  &
                             FATAL, close_file, mpp_pe, mpp_root_pe,     &
                             write_version_number, stdlog,               &
                             read_data, write_data, &
                             mpp_clock_id, mpp_clock_begin, mpp_clock_end, &
                             clock_flag_default, CLOCK_COMPONENT

use  diag_integral_mod, only: diag_integral_init, diag_integral_end, &
                              diag_integral_output

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: register_tracers

!-----------------------------------------------------------------------

implicit none
private

public  update_atmos_model_down, update_atmos_model_up,   &
        atmos_model_init, atmos_model_end, atmos_data_type

!-----------------------------------------------------------------------

 type atmos_data_type
     type (domain2d)               :: domain
     integer                       :: axes(4)
     real, pointer, dimension(:)   :: glon_bnd =>NULL(), &
                                      glat_bnd =>NULL(), &
                                       lon_bnd =>NULL(), &
                                       lat_bnd =>NULL()
     real, pointer, dimension(:,:) :: t_bot =>NULL(), &
                                      q_bot =>NULL(), &
                                      z_bot =>NULL(), &
                                      p_bot =>NULL(), &
                                      u_bot =>NULL(), &
                                      v_bot =>NULL(), &
                                      p_surf =>NULL(), &
                                      gust  =>NULL(),  &
                                      coszen =>NULL(), &
                                      flux_sw =>NULL(), &
                                      flux_lw =>NULL(), &
                                      lprec =>NULL(),   &
                                      fprec =>NULL()
     type (surf_diff_type)         :: Surf_diff
     type (time_type)              :: Time, Time_step, Time_init
     integer, pointer              :: pelist(:) =>NULL()
     logical                       :: pe
 end type

type, public :: land_ice_atmos_boundary_type
!variables of this type are declared by coupler_main, allocated by flux_exchange_init
!quantities going from land+ice to atmos
!       t         = surface temperature for radiation calculations
!       albedo    = surface albedo for radiation calculations
!       land_frac = fraction amount of land in a grid box
!       dt_t      = temperature tendency at the lowest level
!       dt_q      = specific humidity tendency at the lowest level
!       u_flux    = zonal wind stress
!       v_flux    = meridional wind stress
!       dtaudv    = derivative of wind stress w.r.t. the lowest level wind speed
!       u_star    = friction velocity
!       b_star    = bouyancy scale
!       q_star    = moisture scale
!       rough_mom = surface roughness (used for momentum
   real, dimension(:,:), pointer :: t =>NULL(), &
                                    albedo =>NULL(), &
                                    land_frac =>NULL()
   real, dimension(:,:), pointer :: dt_t =>NULL(), &
                                    dt_q =>NULL()
   real, dimension(:,:), pointer :: u_flux =>NULL(), &
                                    v_flux =>NULL(), &
                                    dtaudv =>NULL(), &
                                    u_star =>NULL(), &
                                    b_star =>NULL(), &
                                    q_star =>NULL(), &
                                    rough_mom =>NULL()
   real, dimension(:,:,:), pointer :: data =>NULL() !collective field for "named" fields above
   integer :: xtype             !REGRID, REDIST or DIRECT
end type land_ice_atmos_boundary_type
!quantities going from land alone to atmos (none at present)
type, public :: land_atmos_boundary_type
!   private
   real, dimension(:,:), pointer :: data =>NULL()
end type land_atmos_boundary_type
!quantities going from ice alone to atmos (none at present)
type, public :: ice_atmos_boundary_type
!   private
   real, dimension(:,:), pointer :: data =>NULL()
end type ice_atmos_boundary_type
!Balaji
integer, private :: atmClock
!-----------------------------------------------------------------------

character(len=256) :: version = '$Id: atmos_model.F90,v 10.0 2003/10/24 22:00:21 fms Exp $'
character(len=256) :: tag = '$Name: jakarta $'

character(len=80) :: restart_format = 'atmos_coupled_mod restart format 01'

!-----------------------------------------------------------------------


contains

!#######################################################################

subroutine update_atmos_model_down( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs radiation, damping, and vertical diffusion of momentum,
!    tracers, and downward heat/moisture
!
!-----------------------------------------------------------------------

  type(land_ice_atmos_boundary_type), intent(inout) :: Surface_boundary
  type (atmos_data_type), intent(inout) :: Atmos
                                      
!-----------------------------------------------------------------------
  call mpp_clock_begin(atmClock)

    call atmosphere_down (Atmos%Time, Surface_boundary%land_frac,        &
                          Surface_boundary%t,  Surface_boundary%albedo, Surface_boundary%rough_mom,   &
                          Surface_boundary%u_star,  Surface_boundary%b_star,              &
                          Surface_boundary%q_star, &
                          Surface_boundary%dtaudv, Surface_boundary%u_flux,  Surface_boundary%v_flux,       &
                          Atmos%gust, Atmos%coszen,     &
                          Atmos%flux_sw, Atmos%flux_lw, &
                          Atmos%Surf_diff               )

!-----------------------------------------------------------------------

  call mpp_clock_end(atmClock)
 end subroutine update_atmos_model_down

!#######################################################################

 subroutine update_atmos_model_up( Surface_boundary, Atmos )

!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs upward vertical diffusion of heat/moisture and
!    moisture processes
!
!-----------------------------------------------------------------------

type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
type (atmos_data_type), intent(inout) :: Atmos
                                      
!-----------------------------------------------------------------------
  call mpp_clock_begin(atmClock)


    Atmos%Surf_diff%delta_t = Surface_boundary%dt_t
    Atmos%Surf_diff%delta_q = Surface_boundary%dt_q

    call atmosphere_up (Atmos%Time,  Surface_boundary%land_frac, Atmos%Surf_diff, &
                        Atmos%lprec, Atmos%fprec, Atmos%gust)

!   --- advance time ---

    Atmos % Time = Atmos % Time + Atmos % Time_step


    call get_bottom_mass (Atmos % t_bot, Atmos % q_bot,  &
                          Atmos % p_bot, Atmos % z_bot,  &
                          Atmos % p_surf                 )

    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)


!------ global integrals ------

    call diag_integral_output (Atmos % Time)

!-----------------------------------------------------------------------
  call mpp_clock_end(atmClock)

end subroutine update_atmos_model_up

!#######################################################################

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

type (atmos_data_type), intent(inout) :: Atmos
type (time_type), intent(in) :: Time_init, Time, Time_step

  integer :: unit, ntrace, ntprog, ntdiag, ntfamily
  integer :: mlon, mlat, nlon, nlat, sec, day, ipts, jpts, dt, dto
  character(len=80) :: control
!-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step

!-----------------------------------------------------------------------
! how many tracers have been registered?
!  (will print number below)
   call register_tracers ( MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfamily )
   if ( ntfamily > 0 ) call error_mesg ('atmos_model', 'ntfamily > 0', FATAL)

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

   call write_version_number ( version, tag )

!  number of tracers
   if (mpp_pe() == mpp_root_pe()) then
        write (stdlog(), '(a,i3)') 'Number of tracers =', ntrace
        write (stdlog(), '(a,i3)') 'Number of prognostic tracers =', ntprog
        write (stdlog(), '(a,i3)') 'Number of diagnostic tracers =', ntdiag
   endif

!------ read initial state for several atmospheric fields ------

    if (file_exist('INPUT/atmos_coupled.res')) then
        unit = open_restart_file ('INPUT/atmos_coupled.res', 'read')
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
            if (mpp_pe() == mpp_root_pe()) write (stdlog(),50)
 50         format (/,'The model time step changed .... &
                      &modifying precipitation tendencies')
        endif

    else
        Atmos % lprec = 0.0
        Atmos % fprec = 0.0
        Atmos % gust  = 1.0
    endif

!------ initialize global integral package ------

    call diag_integral_init (Atmos % Time_init, Atmos % Time,  &
                             Atmos % lon_bnd,   Atmos % lat_bnd)

!-----------------------------------------------------------------------
atmClock = mpp_clock_id( 'Atmosphere', flags=clock_flag_default, grain=CLOCK_COMPONENT )
end subroutine atmos_model_init

!#######################################################################

subroutine atmos_model_end (Atmos)

type (atmos_data_type), intent(inout) :: Atmos
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

  unit = open_restart_file ('RESTART/atmos_coupled.res', 'write')
  if (mpp_pe() == mpp_root_pe()) then
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

