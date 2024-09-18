!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Atmos Drivers project.
!*
!* This is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* It is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
module atmos_model_mod
!-----------------------------------------------------------------------
!<OVERVIEW>
!  Driver for the atmospheric model, contains routines to advance the
!  atmospheric model state by one time step.
!</OVERVIEW>

!<DESCRIPTION>
!     This version of atmos_model_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the atmospheric model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Most atmospheric processes (dynamics,radiation,etc.) are performed
!     in the down routine. The up routine finishes the vertical diffusion
!     and computes moisture related terms (convection,large-scale condensation,
!     and precipitation).

!     The boundary variables needed by other component models for coupling
!     are contained in a derived data type. A variable of this derived type
!     is returned when initializing the atmospheric model. It is used by other
!     routines in this module and by coupling routines. The contents of
!     this derived type should only be modified by the atmospheric model.

!</DESCRIPTION>

use mpp_mod,            only: mpp_pe, mpp_root_pe, mpp_clock_id, mpp_clock_begin
use mpp_mod,            only: mpp_clock_end, CLOCK_COMPONENT, MPP_CLOCK_SYNC
use mpp_mod,            only: mpp_min, mpp_max, mpp_error, mpp_chksum
use mpp_domains_mod,   only : mpp_get_compute_domain
use mpp_domains_mod,    only: domain2d
use mpp_mod,            only: mpp_get_current_pelist_name, mpp_set_current_pelist
use mpp_mod,            only: input_nml_file, stdlog, stdout
use fms2_io_mod,        only: file_exists
use fms_mod,            only: write_version_number
use fms_mod,            only: clock_flag_default, error_mesg
use fms_mod,            only: check_nml_error
use diag_manager_mod,   only: diag_send_complete_instant
use time_manager_mod,   only: time_type, get_time, get_date, &
                              operator(+), operator(-)
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_number_tracers, get_tracer_names
use xgrid_mod,          only: grid_box_type
use atmosphere_mod,     only: atmosphere_init
use atmosphere_mod,     only: atmosphere_restart
use atmosphere_mod,     only: atmosphere_end, get_bottom_mass, get_bottom_wind
use atmosphere_mod,     only: atmosphere_state_update
use atmosphere_mod,     only: atmos_phys_driver_statein
use atmosphere_mod,     only: atmosphere_control_data
use atmosphere_mod,     only: atmosphere_resolution, atmosphere_domain
use atmosphere_mod,     only: atmosphere_grid_bdry, atmosphere_grid_ctr
use atmosphere_mod,     only: atmosphere_dynamics, atmosphere_diag_axes
use atmosphere_mod,     only: atmosphere_etalvls, atmosphere_hgt
use atmosphere_mod,     only: atmosphere_diss_est, atmosphere_nggps_diag
use atmosphere_mod,     only: atmosphere_scalar_field_halo
use atmosphere_mod,     only: set_atmosphere_pelist
use atmosphere_mod,     only: atmosphere_coarse_graining_parameters
use atmosphere_mod,     only: atmosphere_coarse_diag_axes
use atmosphere_mod,     only: atmosphere_coarsening_strategy
use atmosphere_mod,     only: Atm, mygrid
use block_control_mod,  only: block_control_type, define_blocks_packed
use IPD_typedefs,       only: IPD_init_type, IPD_control_type, &
                              IPD_data_type, IPD_diag_type,    &
                              IPD_restart_type, kind_phys
use IPD_driver,         only: IPD_initialize, IPD_setup_step, &
                              IPD_radiation_step,             &
                              IPD_physics_step1,              &
                              IPD_physics_step2, IPD_physics_end

use coupler_types_mod, only : coupler_2d_bc_type
use diag_integral_mod, only : diag_integral_init

#ifdef STOCHY
use stochastic_physics, only: init_stochastic_physics,         &
                              run_stochastic_physics
use stochastic_physics_sfc, only: run_stochastic_physics_sfc
#endif
use FV3GFS_io_mod,      only: FV3GFS_restart_read, FV3GFS_restart_write, &
                              FV3GFS_IPD_checksum,                       &
                              gfdl_diag_register, gfdl_diag_output, &
                              FV3GFS_restart_write_coarse, FV3GFS_diag_register_coarse, &
                              sfc_data_override
use FV3GFS_io_mod,      only: register_diag_manager_controlled_diagnostics, register_coarse_diag_manager_controlled_diagnostics
use FV3GFS_io_mod,      only: send_diag_manager_controlled_diagnostic_data
use fv_iau_mod,         only: iau_external_data_type,getiauforcing,iau_initialize
use module_ocean,       only: ocean_init
#if defined (USE_COSP)
use cosp2_test,         only: cosp2_driver
#endif
#if defined (COSP_OFFLINE)
use cosp2_test,         only: cosp2_offline
#endif
!-----------------------------------------------------------------------

implicit none
private

public update_atmos_model_radiation
public update_atmos_model_state
public update_atmos_model_dynamics
public atmos_model_init, atmos_model_end, atmos_data_type
public atmos_model_restart
public land_ice_atmos_boundary_type
public update_atmos_model_up, update_atmos_model_down, lnd_ice_atm_bnd_type_chksum
public atmos_data_type_chksum, lnd_atm_bnd_type_chksum, ice_atm_bnd_type_chksum
public atm_stock_pe
!-----------------------------------------------------------------------

!<PUBLICTYPE >
! FROM AM4/src/atmos_phys/atmos_param/vert_diff/vert_diff.F90
! Surface diffusion derived type data
type surf_diff_type

  real, pointer, dimension(:,:) :: dtmass  => NULL(),   & !dt / mass (mass of lower atmospheric layer)
                                   dflux_t => NULL(),   & ! temperature flux derivative at the top of the lowest atmospheric
                                                          ! layer with respect to the temperature of that layer  (J/(m2 K))
                                   delta_t => NULL(),   & ! the increment in temperature in the lowest atmospheric layer
                                   delta_u => NULL(),   & ! same for u
                                   delta_v => NULL()      ! same for v
  real, pointer, dimension(:,:,:) :: tdt_dyn => NULL(), & ! no explanation found in the vert_diff.F90
                                     qdt_dyn => NULL(), &! no explanation found in the vert_diff.F90
                                     dgz_dyn => NULL(), &! no explanation found in the vert_diff.F90
                                     ddp_dyn => NULL(), &! no explanation found in the vert_diff.F90

                                     tdt_rad => NULL()   !miz

  real, pointer, dimension(:,:,:) :: dflux_tr => NULL(),& ! tracer flux tendency
                                     delta_tr => NULL()   ! tracer tendency
end type surf_diff_type
!<PUBLICTYPE >

!<PUBLICTYPE >
 type atmos_data_type
! Atmosphere related variables
     type (domain2d)               :: domain             ! domain decomposition
     type (domain2d)               :: domain_for_read    ! domain decomposition for reads
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     real,                 pointer, dimension(:,:) :: lon_bnd  => null() ! local longitude axis grid box corners in radians.
     real,                 pointer, dimension(:,:) :: lat_bnd  => null() ! local latitude axis grid box corners in radians.
     real(kind=kind_phys), pointer, dimension(:,:) :: lon      => null() ! local longitude axis grid box centers in radians.
     real(kind=kind_phys), pointer, dimension(:,:) :: lat      => null() ! local latitude axis grid box centers in radians.
     real, pointer, dimension(:,:) :: t_bot    => null() ! temperature at lowest model level
     real, pointer, dimension(:,:,:) :: tr_bot => null() ! tracers at lowest model level
     real, pointer, dimension(:,:) :: z_bot    => null() ! height above the surface for the lowest model level
     real, pointer, dimension(:,:) :: p_bot    => null() ! pressure at lowest model level
     real, pointer, dimension(:,:) :: u_bot    => null() ! zonal wind component at lowest model level
     real, pointer, dimension(:,:) :: v_bot    => null() ! meridional wind component at lowest model level
     real, pointer, dimension(:,:) :: p_surf   => null() ! surface pressure
     real, pointer, dimension(:,:) :: slp      => null() ! sea level pressure
     real, pointer, dimension(:,:) :: gust     => null() ! gustiness factor
     real, pointer, dimension(:,:) :: coszen   => null() ! cosine of the zenith angle
     real, pointer, dimension(:,:) :: flux_sw  => null() ! net shortwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: flux_sw_dir            =>null()
     real, pointer, dimension(:,:) :: flux_sw_dif            =>null()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dir   =>null()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dif   =>null()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dir =>null()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dif =>null()
     real, pointer, dimension(:,:) :: flux_sw_vis            =>null()
     real, pointer, dimension(:,:) :: flux_sw_vis_dir        =>null()
     real, pointer, dimension(:,:) :: flux_sw_vis_dif        =>null()
     real, pointer, dimension(:,:) :: flux_lw  => null() ! net longwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: lprec    => null() ! mass of liquid precipitation since last time step (Kg/m2)
     real, pointer, dimension(:,:) :: fprec    => null() ! ass of frozen precipitation since last time step (Kg/m2)
     logical, pointer, dimension(:,:) :: maskmap =>null()! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. The other logical
                                                         ! processors would be all land points and
                                                         ! are not assigned to actual processors.
                                                         ! This need not be assigned if all logical
                                                         ! processors are used. This variable is dummy and need
                                                         ! not to be set, but it is needed to pass compilation.
     type (surf_diff_type)         :: Surf_diff          ! store data needed by the multi-step version of the diffusion algorithm
     type (time_type)              :: Time               ! current time
     type (time_type)              :: Time_step          ! atmospheric time step.
     type (time_type)              :: Time_init          ! reference time.
     integer                       :: iau_offset         ! iau running window length
     integer, pointer              :: pelist(:) =>null() ! pelist where atmosphere is running.
     logical                       :: pe                 ! current pe.
     type(grid_box_type)           :: grid               ! hold grid information needed for 2nd order conservative flux exchange
                                                         ! to calculate gradient on cubic sphere grid.
     integer                       :: layout(2)          ! computer task laytout
     logical                       :: regional           ! true if domain is regional
     logical                       :: bounded_domain     ! true if domain is bounded
     real(kind=8), pointer, dimension(:) :: ak
     real(kind=8), pointer, dimension(:) :: bk
     real(kind=8), pointer, dimension(:,:,:) :: layer_hgt
     real(kind=8), pointer, dimension(:,:,:) :: level_hgt
     real(kind=kind_phys), pointer, dimension(:,:) :: dx
     real(kind=kind_phys), pointer, dimension(:,:) :: dy
     real(kind=8), pointer, dimension(:,:) :: area
     type(domain2d)                :: coarse_domain      ! domain decomposition of the coarse grid
     logical                       :: write_coarse_restart_files  ! whether to write coarse restart files
     logical                       :: write_only_coarse_intermediate_restarts  ! whether to write only coarse intermediate restart files
     character(len=64)             :: coarsening_strategy  ! Strategy for coarse-graining diagnostics and restart files
     type(coupler_2d_bc_type)      :: fields   ! array of fields used for additional tracers
end type atmos_data_type
!</PUBLICTYPE >

!<PUBLICTYPE >
type land_ice_atmos_boundary_type
   ! variables of this type are declared by coupler_main, allocated by flux_exchange_init.
!quantities going from land+ice to atmos
   real, dimension(:,:),   pointer :: t              =>null() ! surface temperature for radiation calculations
   real, dimension(:,:),   pointer :: t_ocean        =>null() ! ocean surface temperature for shield physics coming
                                                              ! from Ice%t_surf through xgrid !joseph
   real, dimension(:,:),   pointer :: u_ref          =>null() ! surface zonal wind (cjg: PBL depth mods) !bqx
   real, dimension(:,:),   pointer :: v_ref          =>null() ! surface meridional wind (cjg: PBL depth mods) !bqx
   real, dimension(:,:),   pointer :: t_ref          =>null() ! surface air temperature (cjg: PBL depth mods)
   real, dimension(:,:),   pointer :: q_ref          =>null() ! surface air specific humidity (cjg: PBL depth mods)
   real, dimension(:,:),   pointer :: albedo         =>null() ! surface albedo for radiation calculations
   real, dimension(:,:),   pointer :: albedo_vis_dir =>null()
   real, dimension(:,:),   pointer :: albedo_nir_dir =>null()
   real, dimension(:,:),   pointer :: albedo_vis_dif =>null()
   real, dimension(:,:),   pointer :: albedo_nir_dif =>null()
   real, dimension(:,:),   pointer :: land_frac      =>null() ! fraction amount of land in a grid box
   real, dimension(:,:),   pointer :: dt_t           =>null() ! temperature tendency at the lowest level
   real, dimension(:,:,:), pointer :: dt_tr          =>null() ! tracer tendency at the lowest level
   real, dimension(:,:),   pointer :: u_flux         =>null() ! zonal wind stress
   real, dimension(:,:),   pointer :: v_flux         =>null() ! meridional wind stress
   real, dimension(:,:),   pointer :: wind           =>null()
   real, dimension(:,:),   pointer :: thv_atm        =>null()
   real, dimension(:,:),   pointer :: thv_surf       =>null()
   real, dimension(:,:),   pointer :: dtaudu         =>null() ! derivative of zonal wind stress w.r.t. the lowest zonal level wind speed
   real, dimension(:,:),   pointer :: dtaudv         =>null() ! derivative of meridional wind stress w.r.t. the lowest meridional level wind speed
   real, dimension(:,:),   pointer :: u_star         =>null() ! friction velocity
   real, dimension(:,:),   pointer :: b_star         =>null() ! bouyancy scale
   real, dimension(:,:),   pointer :: q_star         =>null() ! moisture scale
   real, dimension(:,:),   pointer :: shflx          =>null() ! sensible heat flux !miz
   real, dimension(:,:),   pointer :: lhflx          =>null() ! latent heat flux   !miz
   real, dimension(:,:),   pointer :: rough_mom      =>null() ! surface roughness (used for momentum)
   real, dimension(:,:),   pointer :: rough_heat     =>null() ! surface roughness (used for heat) ! kgao
   real, dimension(:,:),   pointer :: frac_open_sea  =>null() ! non-seaice fraction (%)
   real, dimension(:,:,:), pointer :: data           =>null() !collective field for "named" fields above
   integer                         :: xtype                   !REGRID, REDIST or DIRECT
end type land_ice_atmos_boundary_type

!<PUBLICTYPE >
type :: land_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>NULL() ! quantities going from land alone to atmos (none at present)
end type land_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
!quantities going from ice alone to atmos (none at present)
type :: ice_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>NULL() ! quantities going from ice alone to atmos (none at present)
end type ice_atmos_boundary_type
!</PUBLICTYPE >

integer :: fv3Clock, getClock, overrideClock, setupClock, radClock, physClock, diagClock, shieldClock

!-----------------------------------------------------------------------
integer :: blocksize       = 1
logical :: chksum_debug    = .false.
logical :: dycore_only     = .false.
logical :: debug           = .false.
logical :: sync            = .false.
logical :: first_time_step = .false.
logical :: fprint          = .true.
logical :: ignore_rst_cksum = .false. ! enforce (.false.) or override (.true.) data integrity restart checksums
logical :: fullcoupler_fluxes = .false. ! get surface fluxes from the full coupler ! for mom6 coupling
real, dimension(4096) :: fdiag = 0. ! xic: TODO: this is hard coded, space can run out in some cases. Should make it allocatable.
logical :: fdiag_override = .false. ! lmh: if true overrides fdiag and fhzer: all quantities are zeroed out
                                    ! after every calcluation, output interval and accumulation/avg/max/min
                                    ! are controlled by diag_manager, fdiag controls output interval only
logical :: do_netcdf_restart = .true.
logical :: restart_tbot_qbot = .false.
integer :: nxblocks = 1
integer :: nyblocks = 1
namelist /atmos_model_nml/ do_netcdf_restart, restart_tbot_qbot, nxblocks, nyblocks, &
                           blocksize, chksum_debug, dycore_only, debug, sync, first_time_step, fdiag, fprint, &
                           fdiag_override, ignore_rst_cksum, fullcoupler_fluxes

type (time_type) :: diag_time, diag_time_fhzero
logical :: fdiag_fix = .false.

!--- concurrent and decoupled radiation and physics variables
!----------------
!  IPD containers ! For SHiELD physics
!----------------
type(IPD_init_type)                 :: Init_parm
type(IPD_control_type)              :: IPD_Control
type(IPD_data_type),    allocatable :: IPD_Data(:)  ! number of blocks
type(IPD_diag_type)                 :: IPD_Diag(250)
type(IPD_restart_type)              :: IPD_Restart

!--------------
! IAU container
!--------------
type(iau_external_data_type)        :: IAU_Data

!-----------------
!  Block container
!-----------------
type (block_control_type), target   :: Atm_block

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

real(kind=kind_phys), parameter :: zero = 0.0_kind_phys

contains

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_down">
!
! <OVERVIEW>
!   compute the atmospheric tendencies for dynamics, radiation,
!   vertical diffusion of momentum, tracers, and heat/moisture.
! !!!IMPORTANT!!! taken from atmos_null, no op for shield runs
! </OVERVIEW>
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down".
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_model_down( Surface_boundary, Atmos )
!   </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>

subroutine update_atmos_model_down( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs radiation, damping, and vertical diffusion of momentum,
!    tracers, and downward heat/moisture
!
!-----------------------------------------------------------------------

  type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
  type (atmos_data_type), intent(inout) :: Atmos

  return
! !!!IMPORTANT!!!no op for shield runs
end subroutine update_atmos_model_down
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_up">
!
!-----------------------------------------------------------------------
! <OVERVIEW>
!   upward vertical diffusion of heat/moisture and moisture processes
! !!!IMPORTANT!!! taken from atmos_null, no op for shield runs
! </OVERVIEW>

!<DESCRIPTION>
!   Called every time step as the atmospheric driver to finish the upward
!   sweep of the tridiagonal elimination for heat/moisture and compute the
!   convective and large-scale tendencies.  The atmospheric variables are
!   advanced one time step and tendencies set back to zero. 
!</DESCRIPTION>

! <TEMPLATE>
!     call  update_atmos_model_up( Surface_boundary, Atmos )
! </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <IN NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </IN>
subroutine update_atmos_model_up( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs upward vertical diffusion of heat/moisture and
!    moisture processes
!
!-----------------------------------------------------------------------

   type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
   type (atmos_data_type), intent(in) :: Atmos
 
   return
! !!!IMPORTANT!!!no op for shield runs

end subroutine update_atmos_model_up


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_radiation">
!
! <OVERVIEW>
!   compute the atmospheric tendencies for dynamics, radiation,
!   vertical diffusion of momentum, tracers, and heat/moisture.
!!!IMPORTANT!!! SHiELD physics are run here
! </OVERVIEW>
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down".
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_model_radiation( Surface_boundary, Atmos)
!   </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </INOUT>
subroutine update_atmos_model_radiation (Surface_boundary, Atmos) ! name change to match the full coupler call
! subroutine update_atmos_radiation_physics (Atmos) !original
! SHiELD physics are run here
!-----------------------------------------------------------------------
  type (atmos_data_type), intent(inout) :: Atmos
  type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
  ! dont know if Surface_boundary is needed here for the fluxes in shield, seems to be used as follow in am4
  ! AM4/src/atmos_phys/atmos_param/radiation/driver/radiation_driver.F90
!--- local variables---
    integer :: nb, jdat(8)
    integer :: nthrds

    call set_atmosphere_pelist() ! should be called before local clocks since they are defined on local atm(n)%pelist
    call mpp_clock_begin(shieldClock)

#ifdef OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif

    if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "statein driver"
!--- get atmospheric state from the dynamic core
    call mpp_clock_begin(getClock)
    if (IPD_control%do_skeb) call atmosphere_diss_est (IPD_control%skeb_npass) !  do smoothing for SKEB
    call atmos_phys_driver_statein (IPD_data, Atm_block)
    call mpp_clock_end(getClock)

!--- get varied surface data
    call mpp_clock_begin(overrideClock)
    call sfc_data_override (Atmos%Time, IPD_data, Atm_block, IPD_Control)
    call mpp_clock_end(overrideClock)

!--- if dycore only run, set up the dummy physics output state as the input state
    if (dycore_only) then
      do nb = 1,Atm_block%nblks
        IPD_Data(nb)%Stateout%gu0 = IPD_Data(nb)%Statein%ugrs
        IPD_Data(nb)%Stateout%gv0 = IPD_Data(nb)%Statein%vgrs
        IPD_Data(nb)%Stateout%gt0 = IPD_Data(nb)%Statein%tgrs
        IPD_Data(nb)%Stateout%gq0 = IPD_Data(nb)%Statein%qgrs
      enddo
    else
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "setup step"

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--- for atmos-ocean coupling: pass surface fluxes from coupler to SHiELD (by joseph and kun)
      if (fullcoupler_fluxes) then
        if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "call apply_sfc_data_to_IPD"
        call apply_sfc_data_to_IPD (Surface_boundary)
      endif
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

!--- update IPD_Control%jdat(8)
      jdat(:) = 0
      call get_date (Atmos%Time, jdat(1), jdat(2), jdat(3),  &
                                 jdat(5), jdat(6), jdat(7))
      IPD_Control%jdat(:) = jdat(:)
!--- execute the IPD atmospheric setup step
      call mpp_clock_begin(setupClock)
      call IPD_setup_step (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart)

#ifdef STOCHY
!--- call stochastic physics pattern generation / cellular automata
      if (IPD_Control%do_sppt .OR. IPD_Control%do_shum .OR. IPD_Control%do_skeb .OR. IPD_Control%do_sfcperts) then
         call run_stochastic_physics(IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Coupling, nthrds)
      end if
#endif

      call mpp_clock_end(setupClock)

!below are the old routines update_atmos_radiation and update_atmos_physics
!no need to continue for dycore_only runs, check logic in SHiELD/atmos_model.F90
      if (dycore_only) return

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "radiation driver"
!--- execute the IPD atmospheric radiation subcomponent (RRTM)
      call mpp_clock_begin(radClock)
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_radiation_step (IPD_Control, IPD_Data(nb), IPD_Diag, IPD_Restart)
      enddo
      call mpp_clock_end(radClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RADIATION STEP  ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "physics driver"
!--- execute the IPD atmospheric physics step1 subcomponent (main physics driver)
      call mpp_clock_begin(physClock)
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_physics_step1 (IPD_Control, IPD_Data(nb), IPD_Diag, IPD_Restart)
      enddo
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP1   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--- for atmos-ocean coupling: pass rad and prec fluxes from IPD to Atmos structure (by kun)
      if (fullcoupler_fluxes) then
        if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "call apply_fluxes_from_IPD_to_Atmos"
        call apply_fluxes_from_IPD_to_Atmos (Atmos)
      endif
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "stochastic physics driver"
!--- execute the IPD atmospheric physics step2 subcomponent (stochastic physics driver)
      call mpp_clock_begin(physClock)
!$OMP parallel do default (none) &
!$OMP            schedule (dynamic,1), &
!$OMP            shared   (Atm_block, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart) &
!$OMP            private  (nb)
      do nb = 1,Atm_block%nblks
        call IPD_physics_step2 (IPD_Control, IPD_Data(nb), IPD_Diag, IPD_Restart)
      enddo
      call mpp_clock_end(physClock)

      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'PHYSICS STEP2   ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif
      call getiauforcing(IPD_Control,IAU_data)
      if (mpp_pe() == mpp_root_pe() .and. debug) write(6,*) "end of radiation and physics step"
    endif

    call mpp_clock_end(shieldClock)
    call mpp_set_current_pelist() !should exit with global pelist to accomodate the full coupler atmos clock

!-----------------------------------------------------------------------
 !end subroutine update_atmos_radiation_physics
end subroutine update_atmos_model_radiation ! name change to match the full coupler call
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atm_stock_pe">
!
! <OVERVIEW>
!  returns the total stock in atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Called to compute and return the total stock (e.g., water, heat, etc.)
! in the atmospheric on the current PE.
! </DESCRIPTION>

! <TEMPLATE>
!   call atm_stock_pe (Atmos, index, value)
! </TEMPLATE>

! <INOUT NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>
!
! <IN NAME="index" TYPE="integer">
!   Index of stock to be computed.
! </IN>
!
! <OUT NAME="value" TYPE="real">
!   Value of stock on the current processor.
! </OUT>

! from coupled atmos_model
subroutine atm_stock_pe (Atm, index, value)

type (atmos_data_type), intent(inout) :: Atm
integer,                intent(in)    :: index
real,                   intent(out)   :: value

   value = 0.0

end subroutine atm_stock_pe

!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!     This routine allocates storage and returns a variable of type
!     atmos_boundary_data_type, and also reads a namelist input and restart file.
! </DESCRIPTION>

! <TEMPLATE>
!     call atmos_model_init (Atmos, Time_init, Time, Time_step, &
!                             do_concurrent_radiation_in)
! </TEMPLATE>

! <IN NAME="Time_init" TYPE="type(time_type)" >
!   The base (or initial) time of the experiment.
! </IN>

! <IN NAME="Time" TYPE="type(time_type)" >
!   The current time.
! </IN>

! <IN NAME="Time_step" TYPE="type(time_type)" >
!   The atmospheric model/physics time step.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step, do_concurrent_radiation) !argument change to match the full coupler
!subroutine atmos_model_init (Atmos, Time_init, Time, Time_step, iau_offset) !check how iau_offset should work

#ifdef OPENMP
  use omp_lib
#endif
  use mpp_mod, only: mpp_npes

  type (atmos_data_type), intent(inout) :: Atmos
  type (time_type), intent(in) :: Time_init, Time, Time_step
  logical, intent(in) :: do_concurrent_radiation
  !integer, intent(inout) :: iau_offset = 0
!--- local variables ---
  integer :: iau_offset = 0
  integer :: unit, ntdiag, ntfamily, i, j, k
  integer :: mlon, mlat, nlon, nlat, nlev, sec, dt, sec_prev
  integer :: ierr, io, logunit
  integer :: idx, tile_num
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: blk, ibs, ibe, jbs, jbe
  real(kind=kind_phys) :: dt_phys
  real, allocatable :: q(:,:,:,:), p_half(:,:,:)
  character(len=80) :: control
  character(len=64) :: filename, filename2, pelist_name
  character(len=132) :: text
  logical :: p_hydro, hydro, fexist
  logical :: do_inline_mp, do_cosp
  logical, save :: block_message = .true.
  integer :: bdat(8), cdat(8)
  integer :: ntracers
  integer :: kdt_prev
  character(len=32), allocatable, target :: tracer_names(:)
  integer :: coarse_diagnostic_axes(4)
  integer :: nthrds
  !-----------------------------------------------------------------------

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
   Atmos % iau_offset = iau_offset
   call get_time (Atmos % Time_step, sec)
   call get_time (Atmos%Time - Atmos%Time_init, sec_prev)
   dt_phys = real(sec)      ! integer seconds
   kdt_prev = int(sec_prev / dt_phys)

   logunit = stdlog()

!-----------------------------------------------------------------------
! initialize atmospheric model -----

!---------- initialize atmospheric dynamics -------
   call atmosphere_init (Atmos%Time_init, Atmos%Time, Atmos%Time_step,&
                         Atmos%grid, Atmos%area, IAU_Data)

   IF ( file_exists('input.nml')) THEN
      read(input_nml_file, nml=atmos_model_nml, iostat=io)
      ierr = check_nml_error(io, 'atmos_model_nml')
   endif
!-----------------------------------------------------------------------
   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   call atmosphere_resolution (nlon, nlat, global=.false.)
   call atmosphere_resolution (mlon, mlat, global=.true.)

   call alloc_atmos_data_type (nlon, nlat, ntracers, Atmos)
   call atmosphere_domain (Atmos%domain, Atmos%domain_for_read, Atmos%layout, Atmos%regional, &
                           Atmos%bounded_domain)
   call alloc_atmos_data_surfdiff_type (ntracers, Atmos)
   call atmosphere_diag_axes (Atmos%axes)
   call atmosphere_etalvls (Atmos%ak, Atmos%bk, flip=.true.)
   call atmosphere_grid_bdry (Atmos%lon_bnd, Atmos%lat_bnd, global=.false.)
   call atmosphere_grid_ctr (Atmos%lon, Atmos%lat)
   call atmosphere_hgt (Atmos%layer_hgt, 'layer', relative=.false., flip=.true.)
   call atmosphere_hgt (Atmos%level_hgt, 'level', relative=.false., flip=.true.)
   call atmosphere_coarse_graining_parameters(Atmos%coarse_domain, Atmos%write_coarse_restart_files, &
                                              Atmos%write_only_coarse_intermediate_restarts)
   call atmosphere_coarsening_strategy(Atmos%coarsening_strategy)

!-----------------------------------------------------------------------
!--- before going any further check definitions for 'blocks'
!-----------------------------------------------------------------------
   call atmosphere_control_data (isc, iec, jsc, jec, nlev, p_hydro, hydro, tile_num, &
                                 do_inline_mp, do_cosp)
   call define_blocks_packed ('atmos_model', Atm_block, isc, iec, jsc, jec, nlev, &
                              blocksize, block_message)

   allocate(IPD_Data(Atm_block%nblks))

#ifdef OPENMP
   nthrds = omp_get_max_threads()
#else
   nthrds = 1
#endif

!--- update IPD_Control%jdat(8)
   bdat(:) = 0
   call get_date (Time_init, bdat(1), bdat(2), bdat(3),  &
                             bdat(5), bdat(6), bdat(7))
   cdat(:) = 0
   call get_date (Time,      cdat(1), cdat(2), cdat(3),  &
                             cdat(5), cdat(6), cdat(7))
   call get_number_tracers(MODEL_ATMOS, num_tracers=ntracers)
   allocate (tracer_names(ntracers))
   do i = 1, ntracers
     call get_tracer_names(MODEL_ATMOS, i, tracer_names(i))
   enddo
!--- setup IPD Init_parm
   Init_parm%me              =  mpp_pe()
   Init_parm%master          =  mpp_root_pe()
   Init_parm%tile_num        =  tile_num
   Init_parm%isc             =  isc
   Init_parm%jsc             =  jsc
   Init_parm%nx              =  nlon
   Init_parm%ny              =  nlat
   Init_parm%levs            =  nlev
   Init_parm%cnx             =  mlon
   Init_parm%cny             =  mlat
   Init_parm%gnx             =  Init_parm%cnx*4
   Init_parm%gny             =  Init_parm%cny*2
   Init_parm%nlunit          =  9999
   Init_parm%logunit         =  logunit
   Init_parm%bdat(:)         =  bdat(:)
   Init_parm%cdat(:)         =  cdat(:)
   Init_parm%dt_dycore       =  dt_phys
   Init_parm%dt_phys         =  dt_phys
   Init_parm%iau_offset      =  Atmos%iau_offset
   Init_parm%blksz           => Atm_block%blksz
   Init_parm%ak              => Atmos%ak
   Init_parm%bk              => Atmos%bk
   Init_parm%xlon            => Atmos%lon
   Init_parm%xlat            => Atmos%lat
   Init_parm%area            => Atmos%area
   Init_parm%tracer_names    => tracer_names
   Init_parm%hydro           = hydro
   Init_parm%do_inline_mp    = do_inline_mp
   Init_parm%do_cosp         = do_cosp

   allocate(Init_parm%input_nml_file, mold=input_nml_file)
   Init_parm%input_nml_file  => input_nml_file
   Init_parm%fn_nml='using internal file'

   call IPD_initialize (IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, Init_parm)

! ensure sfc_coupled is properly set (needs to be true when using ocean coupling)
   if (fullcoupler_fluxes) then
        if (mpp_pe() == mpp_root_pe()) print *, "using ocean coupling - enforce sfc_coupled in SHiELD phys to be true"
        IPD_Control%sfc_coupled = .true.
   else
        if (mpp_pe() == mpp_root_pe()) print *, "no ocean coupling - enforce sfc_coupled in SHiELD phys to be false"
        IPD_Control%sfc_coupled = .false.
   endif

#ifdef STOCHY
   if (IPD_Control%do_sppt .OR. IPD_Control%do_shum .OR. IPD_Control%do_skeb .OR. IPD_Control%do_sfcperts) then
      ! Initialize stochastic physics
      call init_stochastic_physics(IPD_Control, Init_parm, mpp_npes(), nthrds)
      if (mpp_pe() == mpp_root_pe()) print *,'do_skeb=',IPD_Control%do_skeb
   end if

   if (IPD_Control%do_sfcperts) then
      ! Get land surface perturbations here (move to GFS_time_vary
      ! step if wanting to update each time-step)
      call run_stochastic_physics_sfc(IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Coupling)
   end if
#endif

   Atm(mygrid)%flagstruct%do_diss_est = IPD_Control%do_skeb

!  initialize the IAU module
   call iau_initialize (IPD_Control,IAU_data,Init_parm)

   IPD_Control%kdt_prev = kdt_prev

!--- initialize slab ocean model or mixed layer ocean model
   if (IPD_Control%do_ocean) call ocean_init (IPD_Control, Init_parm%logunit, input_nml_file)

   Init_parm%blksz           => null()
   Init_parm%ak              => null()
   Init_parm%bk              => null()
   Init_parm%xlon            => null()
   Init_parm%xlat            => null()
   Init_parm%area            => null()
   Init_parm%tracer_names    => null()
   deallocate (tracer_names)

   !--- update tracers in FV3 with any initialized during the physics/radiation init phase
!rab   call atmosphere_tracer_postinit (IPD_Data, Atm_block)

   call atmosphere_nggps_diag (Time, init=.true.)
!   call gfdl_diag_register (Time, IPD_Data(:)%Sfcprop, IPD_Data(:)%IntDiag, IPD_Data%Cldprop, &
!        Atm_block, Atmos%axes, IPD_Control%nfxr, IPD_Control%ldiag3d, &
!        IPD_Control%nkld, IPD_Control%levs)
   call gfdl_diag_register (Time, IPD_Data(:)%Sfcprop, IPD_Data(:)%IntDiag, IPD_Control, IPD_Data%Cldprop, &
        Atm_block, Atmos%axes)
!   call register_diag_manager_controlled_diagnostics(Time, IPD_Data(:)%IntDiag, Atm_block%nblks, Atmos%axes)
   call register_diag_manager_controlled_diagnostics(Time, IPD_Data(:)%Sfcprop, IPD_Data(:)%IntDiag, IPD_Control, Atm_block%nblks, Atmos%axes)
   if (Atm(mygrid)%coarse_graining%write_coarse_diagnostics) then
       call atmosphere_coarse_diag_axes(coarse_diagnostic_axes)
       call FV3GFS_diag_register_coarse(Time, coarse_diagnostic_axes)
       call register_coarse_diag_manager_controlled_diagnostics(Time, coarse_diagnostic_axes)
    endif
   if (.not. dycore_only) &
      call FV3GFS_restart_read (IPD_Data, IPD_Restart, Atm_block, IPD_Control, Atmos%domain_for_read, ignore_rst_cksum)
      if (chksum_debug) then
        if (mpp_pe() == mpp_root_pe()) print *,'RESTART READ  ', IPD_Control%kdt, IPD_Control%fhour
        call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
      endif

   !--- set the initial diagnostic timestamp
   diag_time = Time
   if (Atmos%iau_offset > zero) then
     call get_time (Atmos%Time - Atmos%Time_init, sec)
     if (sec < Atmos%iau_offset*3600) then
       diag_time = Atmos%Time_init
       diag_time_fhzero = Atmos%Time
     endif
   endif



!-----------------------------------------------------------------------
!------ get initial state for dynamics -------
    call get_bottom_mass (Atmos % t_bot,  Atmos % tr_bot, &
                          Atmos % p_bot,  Atmos % z_bot,  &
                          Atmos % p_surf, Atmos % slp     )

    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)


   !--- print version number to logfile
   call write_version_number ( version, tagname )

   !--- write the namelist to a log file
   write (stdlog(), nml=atmos_model_nml)

   !--- get fdiag
#ifdef GFS_PHYS
!--- check fdiag to see if it is an interval or a list
   if (fdiag_override) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---OVERRIDING fdiag: USING SETTINGS IN diag_table for GFS PHYSICS DIAGS"
      IPD_Control%fhzero = dt_phys / 3600.
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---fhzero IS SET TO dt_atmos: ALL DIAGNOSTICS ARE SINGLE-STEP"
   else
      if (nint(fdiag(2)) == 0) then
         fdiag_fix = .true.
         do i = 2, size(fdiag,1)
            fdiag(i) = fdiag(1) * real(i)
         enddo
      endif
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---fdiag",fdiag(1:40)
   endif
#endif

   setupClock = mpp_clock_id( '---GFS Step Setup     ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   overrideClock = mpp_clock_id( '---GFS Override       ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   radClock   = mpp_clock_id( '---GFS Radiation      ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   physClock  = mpp_clock_id( '---GFS Physics        ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   diagClock  = mpp_clock_id( '---GFS Diag           ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   getClock   = mpp_clock_id( '---GFS Get State      ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   if (sync) then
     fv3Clock = mpp_clock_id( '--FV3 Dycore          ', flags=clock_flag_default+MPP_CLOCK_SYNC, grain=CLOCK_COMPONENT )
   else
     fv3Clock = mpp_clock_id( '--FV3 Dycore          ', flags=clock_flag_default, grain=CLOCK_COMPONENT )
   endif
   shieldClock= mpp_clock_id( '--SHiELD Physics      ', flags=clock_flag_default, grain=CLOCK_COMPONENT )

   call diag_integral_init (Atmos % Time_init, Atmos % Time,  &
                            Atmos % lon_bnd(:,:),  &
                            Atmos % lat_bnd(:,:), Atmos % area)

   call mpp_set_current_pelist() !should exit with global pelist to accomodate the full coupler atmos clock

end subroutine atmos_model_init
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_dynamics"
!
! <OVERVIEW>
! update atmosphere dynamics 
! </OVERVIEW>
! <DESCRIPTION>
!     This routine advances prognostic atmosphere variables in time
! </DESCRIPTION>
! <TEMPLATE>
!     call update_atmos_model_dynamics (Atmos)
! </TEMPLATE>
subroutine update_atmos_model_dynamics (Atmos)
! run the atmospheric dynamics to advect the properties
  type (atmos_data_type), intent(in) :: Atmos ! should this be 'inout'

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call atmosphere_dynamics (Atmos%Time)
    call mpp_clock_end(fv3Clock)
    call mpp_set_current_pelist() !should exit with global pelist to accomodate the full coupler atmos clock

end subroutine update_atmos_model_dynamics
! </SUBROUTINE>


!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_state"
!
! <OVERVIEW>
! update atmosphere physics 
! </OVERVIEW>
! <DESCRIPTION>
!     This routine update atmosphere physics variables
! </DESCRIPTION>
! <TEMPLATE>
!     call update_atmos_model_state (Atmos)
! </TEMPLATE>
subroutine update_atmos_model_state (Atmos)
! to update the model state after all concurrency is completed
  type (atmos_data_type), intent(inout) :: Atmos
!--- local variables
  integer :: isec,seconds,isec_fhzero
  real(kind=kind_phys) :: time_int, time_intfull
  integer :: is, ie, js, je, kt

    call set_atmosphere_pelist()
    call mpp_clock_begin(fv3Clock)
    call atmosphere_state_update (Atmos%Time, IPD_Data, IAU_Data, Atm_block)
    call mpp_clock_end(fv3Clock)

    call mpp_clock_begin(shieldClock)
    call mpp_clock_begin(diagClock)

    if (chksum_debug) then
      if (mpp_pe() == mpp_root_pe()) print *,'UPDATE STATE    ', IPD_Control%kdt, IPD_Control%fhour
      call FV3GFS_IPD_checksum(IPD_Control, IPD_Data, Atm_block)
    endif

!------ advance time ------
    Atmos % Time = Atmos % Time + Atmos % Time_step

    call get_bottom_mass (Atmos % t_bot,  Atmos % tr_bot, &
                          Atmos % p_bot,  Atmos % z_bot,  &
                          Atmos % p_surf, Atmos % slp     )
    call get_bottom_wind (Atmos % u_bot, Atmos % v_bot)
    call atmosphere_control_data(is, ie, js, je, kt)
    call send_diag_manager_controlled_diagnostic_data(Atmos%Time, &
       Atm_block, IPD_Data, IPD_Control%nx, IPD_Control%ny, IPD_Control%levs, &
       Atm(mygrid)%coarse_graining%write_coarse_diagnostics, &
       real(Atm(mygrid)%delp(is:ie,js:je,:), kind=kind_phys), &
       Atmos%coarsening_strategy, real(Atm(mygrid)%ptop, kind=kind_phys))

    call get_time (Atmos%Time - diag_time, isec)
    call get_time (Atmos%Time - Atmos%Time_init, seconds)

    time_int = real(isec)
    if (ANY(nint(fdiag(:)*3600.0) == seconds) .or. (fdiag_fix .and. mod(seconds, nint(fdiag(1)*3600.0)) .eq. 0) .or. (IPD_Control%kdt == 1 .and. first_time_step) ) then
      if (mpp_pe() == mpp_root_pe()) write(6,*) "---isec,seconds",isec,seconds
      if (mpp_pe() == mpp_root_pe()) write(6,*) ' gfs diags time since last bucket empty: ',time_int/3600.,'hrs'
      call atmosphere_nggps_diag(Atmos%Time)
    endif
    if (ANY(nint(fdiag(:)*3600.0) == seconds) .or. (fdiag_fix .and. mod(seconds, nint(fdiag(1)*3600.0)) .eq. 0) .or. (IPD_Control%kdt == 1 .and. first_time_step)) then
      if(Atmos%iau_offset > zero) then
        if( time_int - Atmos%iau_offset*3600. > zero ) then
          time_int = time_int - Atmos%iau_offset*3600.
        else if(seconds == Atmos%iau_offset*3600) then
          call get_time (Atmos%Time - diag_time_fhzero, isec_fhzero)
          time_int = real(isec_fhzero)
          if (mpp_pe() == mpp_root_pe()) write(6,*) "---iseczero",isec_fhzero
        endif
      endif
      time_intfull = real(seconds)
      if(Atmos%iau_offset > zero) then
        if( time_intfull - Atmos%iau_offset*3600. > zero) then
          time_intfull = time_intfull - Atmos%iau_offset*3600.
        endif
      endif
#if defined (USE_COSP)
      !-----------------------------------------------------------------------
      ! The CFMIP Observation Simulator Package (COSP)
      ! Added by Linjiong Zhou
      ! May 2021
      !-----------------------------------------------------------------------

      if (IPD_Control%do_cosp) then

        call cosp2_driver (IPD_Control, IPD_Data(:)%Grid, IPD_Data(:)%Statein, &
                           IPD_Data(:)%Stateout, IPD_Data(:)%Sfcprop, &
                           IPD_Data(:)%Radtend, IPD_Data(:)%Intdiag, Init_parm)

      endif
#endif
#if defined (COSP_OFFLINE)
      !-----------------------------------------------------------------------
      ! Output Variables for the Offline CFMIP Observation Simulator Package (COSP)
      ! Added by Linjiong Zhou
      ! Nov 2022
      !-----------------------------------------------------------------------

      if (IPD_Control%do_cosp) then

        call cosp2_offline (IPD_Control, IPD_Data(:)%Statein, &
                           IPD_Data(:)%Stateout, IPD_Data(:)%Sfcprop, &
                           IPD_Data(:)%Radtend, IPD_Data(:)%Intdiag, Init_parm)

      endif
#endif
      call gfdl_diag_output(Atmos%Time, Atm_block, IPD_Data, IPD_Control%nx, IPD_Control%ny, fprint, &
                            IPD_Control%levs, 1, 1, 1.d0, time_int, time_intfull, &
                            IPD_Control%fhswr, IPD_Control%fhlwr, &
                            mod(seconds, nint(fdiag(1)*3600.0)) .eq. 0, &
                            Atm(mygrid)%coarse_graining%write_coarse_diagnostics,&
                            real(Atm(mygrid)%delp(is:ie,js:je,:), kind=kind_phys), &
                            Atmos%coarsening_strategy, real(Atm(mygrid)%ptop, kind=kind_phys))
      call diag_send_complete_instant (Atmos%Time)
      if (mod(isec,nint(3600*IPD_Control%fhzero)) == 0) diag_time = Atmos%Time
    endif

    call mpp_clock_end(diagClock)
    call mpp_clock_end(shieldClock)
    call mpp_set_current_pelist() !should exit with global pelist to accomodate the full coupler atmos clock

 end subroutine update_atmos_model_state
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="apply_sfc_data_to_IPD">
!
! <OVERVIEW>
! update IPD layer variables
! </OVERVIEW>

! <DESCRIPTION>
! apply coupler variables from the exchange grid as input for SHiELD physics
! through the IPD layer
! </DESCRIPTION>

! <TEMPLATE>
!     call apply_sfc_data_to_IPD (Surface_boundary)
! </TEMPLATE>

subroutine apply_sfc_data_to_IPD (Surface_boundary)
!
!By Joseph and Kun: 
!Here we use sfc-layer variables/fluxes over ocean points from 
!the coupler code (through xgrid and saved in atmos%surface_boundary)
!to update variables in SHiELD physics 
!
  type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
  integer :: nb, blen, ix, i, j

  do nb = 1,Atm_block%nblks
     blen = Atm_block%blksz(nb)
     do ix = 1, blen
        i = Atm_block%index(nb)%ii(ix)
        j = Atm_block%index(nb)%jj(ix)
        ! sensible heat flux (rho*cp_air*t_flux)
        IPD_Data(nb)%Sfcprop%shflx(ix)  = Surface_boundary%shflx(i,j)
        ! moisture flux (rho*q_flux)
        IPD_Data(nb)%Sfcprop%lhflx(ix)  = Surface_boundary%lhflx(i,j)
        ! only do ocean points for the fields below 
        if (nint(IPD_Data(nb)%Sfcprop%slmsk(ix)) == 0) then
          ! sea surface temp 
          IPD_Data(nb)%Sfcprop%tsfc(ix)   = Surface_boundary%t_ocean(i,j)
          ! roughness length for momentum in cm
          IPD_Data(nb)%Sfcprop%zorl(ix)   = 100.* Surface_boundary%rough_mom(i,j)
          ! roughness length for heat in cm
          IPD_Data(nb)%Sfcprop%ztrl(ix)   = 100.* Surface_boundary%rough_heat(i,j)
          ! ustar
          IPD_Data(nb)%Sfcprop%uustar(ix) = Surface_boundary%u_star(i,j)
        endif
     enddo
  enddo

end subroutine apply_sfc_data_to_IPD
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="apply_fluxes_from_IPD_to_Atmos">
!
! <OVERVIEW>
! update Atmos variables
! </OVERVIEW>

! <DESCRIPTION>
! apply surface shortwave, longwave, precip fluxes from SHiELD physics to Atmos data structure, 
! which will be put on exchange grid in the coupler for other component models to use 
! 
! created by Kun Gao (kun.gao@noaa.gov)
! 
! </DESCRIPTION>

! <TEMPLATE>
!     call apply_fluxes_from_IPD_to_Atmos ( Atmos )
! </TEMPLATE>

subroutine apply_fluxes_from_IPD_to_Atmos ( Atmos )

  type (atmos_data_type), intent(inout) :: Atmos
  integer :: nb, blen, ix, i, j

  real :: nirbmdi, visbmdi, nirbmui, visbmui
  real :: nirdfdi, visdfdi, nirdfui, visdfui

  do nb = 1,Atm_block%nblks
     blen = Atm_block%blksz(nb)
     do ix = 1, blen
        i = Atm_block%index(nb)%ii(ix) - Atm_block%isc + 1
        j = Atm_block%index(nb)%jj(ix) - Atm_block%jsc + 1

        ! --- shortwave

        ! eight sw components: dir/dif; nir/vis; down/up
        nirbmdi = IPD_Data(nb)%Coupling%nirbmdi(ix)
        nirdfdi = IPD_Data(nb)%Coupling%nirdfdi(ix)
        visbmdi = IPD_Data(nb)%Coupling%visbmdi(ix)
        visdfdi = IPD_Data(nb)%Coupling%visdfdi(ix)
        nirbmui = IPD_Data(nb)%Coupling%nirbmui(ix)
        nirdfui = IPD_Data(nb)%Coupling%nirdfui(ix)
        visbmui = IPD_Data(nb)%Coupling%visbmui(ix)
        visdfui = IPD_Data(nb)%Coupling%visdfui(ix)

        ! cosine of zenith angle
        Atmos%coszen(i,j) = IPD_Data(nb)%Radtend%coszen(ix)

        ! visible down
        Atmos%flux_sw_down_vis_dir(i,j) = visbmdi                         ! downward visible sw flux at surface - direct
        Atmos%flux_sw_down_vis_dif(i,j) = visdfdi                         ! downward visible sw flux at surface - diffused

        ! visible net
        Atmos%flux_sw_vis_dir(i,j) = visbmdi - visbmui                    ! net (down-up) visible sw flux at surface - direct
        Atmos%flux_sw_vis_dif(i,j) = visdfdi - visdfui                    ! net (down-up) visible sw flux at surface - diffused

        ! total down
        Atmos%flux_sw_down_total_dir(i,j) = nirbmdi + visbmdi             ! downward total sw flux at surface - direct
        Atmos%flux_sw_down_total_dif(i,j) = nirdfdi + visdfdi             ! downward total sw flux at surface - diffused

        ! total net
        Atmos%flux_sw_dir(i,j) = (nirbmdi + visbmdi) - (nirbmui + visbmui)! net (down-up) sw flux at surface - direct
        Atmos%flux_sw_dif(i,j) = (nirdfdi + visdfdi) - (nirdfui + visbmui)! net (down-up) sw flux at surface - diffused

        ! total net and visible net; not used on exchange grid (not important)
        Atmos%flux_sw(i,j)     = Atmos%flux_sw_dir(i,j) + Atmos%flux_sw_dif(i,j)
        Atmos%flux_sw_vis(i,j) = Atmos%flux_sw_vis_dir(i,j) + Atmos%flux_sw_vis_dir(i,j)

        ! --- longwave
        ! total downward lw flux at sfc
        Atmos%flux_lw(i,j) = IPD_Data(nb)%Radtend%sfcflw(ix)%dnfxc

        ! --- precip rate (kg/m**2/s)
        if ( IPD_Data(nb)%Sfcprop%srflag(ix) .lt. 0.5) then  ! rain (srflag = 0)
          Atmos%lprec(i,j) = 1./IPD_Control%dtp * IPD_Data(nb)%Sfcprop%tprcp(ix)
          Atmos%fprec(i,j) = 0.
        else                                                 ! snow (srflag = 1)
          Atmos%lprec(i,j) = 0.
          Atmos%fprec(i,j) = 1./IPD_Control%dtp * IPD_Data(nb)%Sfcprop%tprcp(ix)
        endif
     enddo
  enddo

end subroutine apply_fluxes_from_IPD_to_Atmos
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_end (Atmos)
  type (atmos_data_type), intent(inout) :: Atmos
!---local variables
  integer :: idx

    call IPD_physics_end (IPD_Control)

!-----------------------------------------------------------------------
!---- termination routine for atmospheric model ----

    call atmosphere_end (Atmos % Time, Atmos%grid)
    if (.not. dycore_only) then
       call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
            IPD_Control, Atmos%domain)
       if (Atmos%write_coarse_restart_files) then
          call FV3GFS_restart_write_coarse(IPD_Data, IPD_Restart, Atm_block, &
            IPD_Control, Atmos%coarse_domain)
       endif
    endif

end subroutine atmos_model_end
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_restart">
! <DESCRIPTION>
!  Write out restart files registered through register_restart_file
! </DESCRIPTION>
subroutine atmos_model_restart(Atmos, timestamp)
  type (atmos_data_type),   intent(inout) :: Atmos
  character(len=*),  intent(in)           :: timestamp

    call atmosphere_restart(timestamp)
    if (.not. dycore_only) then
       if (.not. Atmos%write_only_coarse_intermediate_restarts) then
          call FV3GFS_restart_write (IPD_Data, IPD_Restart, Atm_block, &
               IPD_Control, Atmos%domain, timestamp)
       endif
       if (Atmos%write_coarse_restart_files) then
          call FV3GFS_restart_write_coarse(IPD_Data, IPD_Restart, Atm_block, &
               IPD_Control, Atmos%coarse_domain, timestamp)
       endif
    endif
end subroutine atmos_model_restart
! </SUBROUTINE>

!#######################################################################
!#######################################################################
! <SUBROUTINE NAME="atmos_data_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the atmos_data_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the atmos_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, atm)
! </TEMPLATE>

! <IN NAME="Atm" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields in the atmos_data_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine atmos_data_type_chksum(id, timestep, atm)
type(atmos_data_type), intent(in) :: atm
    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    integer :: n, outunit

100 format("CHECKSUM::",A32," = ",Z20)
101 format("CHECKSUM::",A16,a,'%',a," = ",Z20)

  outunit = stdout()
  write(outunit,*) 'BEGIN CHECKSUM(Atmos_data_type):: ', id, timestep
  write(outunit,100) ' atm%lon_bnd                ', mpp_chksum(atm%lon_bnd               )
  write(outunit,100) ' atm%lat_bnd                ', mpp_chksum(atm%lat_bnd               )
  write(outunit,100) ' atm%lon                    ', mpp_chksum(atm%lon                   )
  write(outunit,100) ' atm%lat                    ', mpp_chksum(atm%lat                   )
  write(outunit,100) ' atm%t_bot                  ', mpp_chksum(atm%t_bot                 )
  do n = 1, size(atm%tr_bot,3)
  write(outunit,100) ' atm%tr_bot(:,:,n)          ', mpp_chksum(atm%tr_bot(:,:,n)         )
  enddo
  write(outunit,100) ' atm%z_bot                  ', mpp_chksum(atm%z_bot                 )
  write(outunit,100) ' atm%p_bot                  ', mpp_chksum(atm%p_bot                 )
  write(outunit,100) ' atm%u_bot                  ', mpp_chksum(atm%u_bot                 )
  write(outunit,100) ' atm%v_bot                  ', mpp_chksum(atm%v_bot                 )
  write(outunit,100) ' atm%p_surf                 ', mpp_chksum(atm%p_surf                )
  write(outunit,100) ' atm%slp                    ', mpp_chksum(atm%slp                   )
  write(outunit,100) ' atm%gust                   ', mpp_chksum(atm%gust                  )
  write(outunit,100) ' atm%coszen                 ', mpp_chksum(atm%coszen                )
  write(outunit,100) ' atm%flux_sw                ', mpp_chksum(atm%flux_sw               )
  write(outunit,100) ' atm%flux_sw_dir            ', mpp_chksum(atm%flux_sw_dir           )
  write(outunit,100) ' atm%flux_sw_dif            ', mpp_chksum(atm%flux_sw_dif           )
  write(outunit,100) ' atm%flux_sw_down_vis_dir   ', mpp_chksum(atm%flux_sw_down_vis_dir  )
  write(outunit,100) ' atm%flux_sw_down_vis_dif   ', mpp_chksum(atm%flux_sw_down_vis_dif  )
  write(outunit,100) ' atm%flux_sw_down_total_dir ', mpp_chksum(atm%flux_sw_down_total_dir)
  write(outunit,100) ' atm%flux_sw_down_total_dif ', mpp_chksum(atm%flux_sw_down_total_dif)
  write(outunit,100) ' atm%flux_sw_vis            ', mpp_chksum(atm%flux_sw_vis           )
  write(outunit,100) ' atm%flux_sw_vis_dir        ', mpp_chksum(atm%flux_sw_vis_dir       )
  write(outunit,100) ' atm%flux_sw_vis_dif        ', mpp_chksum(atm%flux_sw_vis_dif       )
  write(outunit,100) ' atm%flux_lw                ', mpp_chksum(atm%flux_lw               )
  write(outunit,100) ' atm%lprec                  ', mpp_chksum(atm%lprec                 )
  write(outunit,100) ' atm%fprec                  ', mpp_chksum(atm%fprec                 )
!  call surf_diff_type_chksum(id, timestep, atm%surf_diff)

end subroutine atmos_data_type_chksum

! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="lnd_ice_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the land_ice_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the land_ice_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_data_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains fields in the land_ice_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine lnd_ice_atm_bnd_type_chksum(id, timestep, bnd_type)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_ice_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(lnd_ice_Atm_bnd_type):: ', id, timestep
100 format("CHECKSUM::",A32," = ",Z20)
    write(outunit,100) 'lnd_ice_atm_bnd_type%t             ',mpp_chksum(bnd_type%t              )
    write(outunit,100) 'lnd_ice_atm_bnd_type%t_ocean       ',mpp_chksum(bnd_type%t_ocean        )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo        ',mpp_chksum(bnd_type%albedo         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_vis_dir',mpp_chksum(bnd_type%albedo_vis_dir )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_nir_dir',mpp_chksum(bnd_type%albedo_nir_dir )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_vis_dif',mpp_chksum(bnd_type%albedo_vis_dif )
    write(outunit,100) 'lnd_ice_atm_bnd_type%albedo_nir_dif',mpp_chksum(bnd_type%albedo_nir_dif )
    write(outunit,100) 'lnd_ice_atm_bnd_type%land_frac     ',mpp_chksum(bnd_type%land_frac      )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dt_t          ',mpp_chksum(bnd_type%dt_t           )
    do n = 1, size(bnd_type%dt_tr,3)
    write(outunit,100) 'lnd_ice_atm_bnd_type%dt_tr(:,:,n)  ',mpp_chksum(bnd_type%dt_tr(:,:,n)   )
    enddo
    write(outunit,100) 'lnd_ice_atm_bnd_type%u_flux        ',mpp_chksum(bnd_type%u_flux         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%v_flux        ',mpp_chksum(bnd_type%v_flux         )

    write(outunit,100) 'lnd_ice_atm_bnd_type%wind          ',mpp_chksum(bnd_type%wind           )
    write(outunit,100) 'lnd_ice_atm_bnd_type%thv_atm       ',mpp_chksum(bnd_type%thv_atm        )
    write(outunit,100) 'lnd_ice_atm_bnd_type%thv_surf      ',mpp_chksum(bnd_type%thv_surf       )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dtaudu        ',mpp_chksum(bnd_type%dtaudu         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%dtaudv        ',mpp_chksum(bnd_type%dtaudv         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%u_star        ',mpp_chksum(bnd_type%u_star         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%b_star        ',mpp_chksum(bnd_type%b_star         )
    write(outunit,100) 'lnd_ice_atm_bnd_type%q_star        ',mpp_chksum(bnd_type%q_star         )
#ifndef use_AM3_physics
    write(outunit,100) 'lnd_ice_atm_bnd_type%shflx         ',mpp_chksum(bnd_type%shflx          )!miz
    write(outunit,100) 'lnd_ice_atm_bnd_type%lhflx         ',mpp_chksum(bnd_type%lhflx          )!miz
#endif
    write(outunit,100) 'lnd_ice_atm_bnd_type%rough_mom     ',mpp_chksum(bnd_type%rough_mom      )
    write(outunit,100) 'lnd_ice_atm_bnd_type%rough_heat    ',mpp_chksum(bnd_type%rough_heat     )!kgao
!    write(outunit,100) 'lnd_ice_atm_bnd_type%data          ',mpp_chksum(bnd_type%data           )

end subroutine lnd_ice_atm_bnd_type_chksum
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="lnd_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the land_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the land_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call lnd_atm_bnd_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(land_atmos_boundary_type)">
!   Derived-type variable that contains fields in the land_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine lnd_atm_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(lnd_atmos_boundary_type):: ', id, timestep
!    write(outunit,100) 'lnd_atm_bnd_type%data',mpp_chksum(bnd_type%data)

100 format("CHECKSUM::",A32," = ",Z20)

end subroutine lnd_atm_bnd_type_chksum
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="ice_atm_bnd_type_chksum">
!
! <OVERVIEW>
!  Print checksums of the various fields in the ice_atmos_boundary_type.
! </OVERVIEW>

! <DESCRIPTION>
!  Routine to print checksums of the various fields in the ice_atmos_boundary_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call ice_atm_bnd_type_chksum(id, timestep, bnd_type)
! </TEMPLATE>

! <IN NAME="bnd_type" TYPE="type(ice_atmos_boundary_type)">
!   Derived-type variable that contains fields in the ice_atmos_boundary_type.
! </INOUT>
!
! <IN NAME="id" TYPE="character">
!   Label to differentiate where this routine in being called from.
! </IN>
!
! <IN NAME="timestep" TYPE="integer">
!   An integer to indicate which timestep this routine is being called for.
! </IN>
!
subroutine ice_atm_bnd_type_chksum(id, timestep, bnd_type)
  use fms_mod,                 only: stdout
  use mpp_mod,                 only: mpp_chksum

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(ice_atmos_boundary_type), intent(in) :: bnd_type
 integer ::   n, outunit

    outunit = stdout()
    write(outunit,*) 'BEGIN CHECKSUM(ice_atmos_boundary_type):: ', id, timestep
!    write(outunit,100) 'ice_atm_bnd_type%data',mpp_chksum(data_type%data)

100 format("CHECKSUM::",A32," = ",Z20)


end subroutine ice_atm_bnd_type_chksum
! </SUBROUTINE>


  subroutine alloc_atmos_data_type (nlon, nlat, ntprog, Atmos)
   integer, intent(in) :: nlon, nlat, ntprog
   type(atmos_data_type), intent(inout) :: Atmos

    allocate ( Atmos % lon_bnd  (nlon+1,nlat+1), &
               Atmos % lat_bnd  (nlon+1,nlat+1), &
               Atmos % lon      (nlon,nlat), &
               Atmos % lat      (nlon,nlat), &
               Atmos % t_bot    (nlon,nlat), &
               Atmos % tr_bot   (nlon,nlat, ntprog), &
               Atmos % z_bot    (nlon,nlat), &
               Atmos % p_bot    (nlon,nlat), &
               Atmos % u_bot    (nlon,nlat), &
               Atmos % v_bot    (nlon,nlat), &
               Atmos % p_surf   (nlon,nlat), &
               Atmos % slp      (nlon,nlat), &
               Atmos % gust     (nlon,nlat), &
               Atmos % flux_sw  (nlon,nlat), &
               Atmos % flux_sw_dir (nlon,nlat), &
               Atmos % flux_sw_dif (nlon,nlat), &
               Atmos % flux_sw_down_vis_dir (nlon,nlat), &
               Atmos % flux_sw_down_vis_dif (nlon,nlat), &
               Atmos % flux_sw_down_total_dir (nlon,nlat), &
               Atmos % flux_sw_down_total_dif (nlon,nlat), &
               Atmos % flux_sw_vis (nlon,nlat), &
               Atmos % flux_sw_vis_dir (nlon,nlat), &
               Atmos % flux_sw_vis_dif(nlon,nlat), &
               Atmos % flux_lw  (nlon,nlat), &
               Atmos % coszen   (nlon,nlat), &
               Atmos % lprec    (nlon,nlat), &
               Atmos % fprec    (nlon,nlat)  )

    Atmos % flux_sw                 = 0.0
    Atmos % flux_lw                 = 0.0
    Atmos % flux_sw_dir             = 0.0
    Atmos % flux_sw_dif             = 0.0
    Atmos % flux_sw_down_vis_dir    = 0.0
    Atmos % flux_sw_down_vis_dif    = 0.0
    Atmos % flux_sw_down_total_dir  = 0.0
    Atmos % flux_sw_down_total_dif  = 0.0
    Atmos % flux_sw_vis             = 0.0
    Atmos % flux_sw_vis_dir         = 0.0
    Atmos % flux_sw_vis_dif         = 0.0
    Atmos % coszen                  = 0.0
    Atmos % tr_bot                  = 0.0 

  end subroutine alloc_atmos_data_type


  subroutine alloc_atmos_data_surfdiff_type (ntprog, Atmos)
   integer, intent(in) :: ntprog
   type(atmos_data_type), intent(inout) :: Atmos
   integer :: is, ie, js, je

    call mpp_get_compute_domain(Atmos%domain,is,ie,js,je)

    allocate ( Atmos % Surf_diff % dtmass(is:ie, js:je) )
    allocate ( Atmos % Surf_diff % dflux_t(is:ie, js:je) )
    allocate ( Atmos % Surf_diff % delta_t(is:ie, js:je) )
    allocate ( Atmos % Surf_diff % delta_u(is:ie, js:je) )
    allocate ( Atmos % Surf_diff % delta_v(is:ie, js:je) )
    allocate ( Atmos % Surf_diff % dflux_tr(is:ie, js:je, ntprog) )
    allocate ( Atmos % Surf_diff % delta_tr(is:ie, js:je, ntprog) )

    Atmos % Surf_diff % dtmass   = 0.0
    Atmos % Surf_diff % dflux_t  = 0.0
    Atmos % Surf_diff % delta_t  = 0.0
    Atmos % Surf_diff % delta_u  = 0.0
    Atmos % Surf_diff % delta_v  = 0.0
    Atmos % Surf_diff % dflux_tr = 0.0
    Atmos % Surf_diff % delta_tr = 0.0
  end subroutine alloc_atmos_data_surfdiff_type

  subroutine dealloc_atmos_data_type (Atmos)
   type(atmos_data_type), intent(inout) :: Atmos
    deallocate (Atmos%lon_bnd,                &
                Atmos%lat_bnd,                &
                Atmos%lon,                    &
                Atmos%lat,                    &
                Atmos%t_bot,                  &
                Atmos%tr_bot,                 &
                Atmos%z_bot,                  &
                Atmos%p_bot,                  &
                Atmos%u_bot,                  &
                Atmos%v_bot,                  &
                Atmos%p_surf,                 &
                Atmos%slp,                    &
                Atmos%gust,                   &
                Atmos%flux_sw,                &
                Atmos%flux_sw_dir,            &
                Atmos%flux_sw_dif,            &
                Atmos%flux_sw_down_vis_dir,   &
                Atmos%flux_sw_down_vis_dif,   &
                Atmos%flux_sw_down_total_dir, &
                Atmos%flux_sw_down_total_dif, &
                Atmos%flux_sw_vis,            &
                Atmos%flux_sw_vis_dir,        &
                Atmos%flux_sw_vis_dif,        &
                Atmos%flux_lw,                &
                Atmos%coszen,                 &
                Atmos%lprec,                  &
                Atmos%fprec  )
  end subroutine dealloc_atmos_data_type

end module atmos_model_mod
