
program atmos_model

!-----------------------------------------------------------------------
!
!  Main program for running a stand-alone atmospheric dynamical core.
!
!-----------------------------------------------------------------------

use   atmosphere_mod, only: atmosphere_init, atmosphere_end, atmosphere, atmosphere_domain

use time_manager_mod, only: time_type, set_time, get_time,  &
                            operator(+), operator (<), operator (>), &
                            operator (/=), operator (/), operator (*)

use          fms_mod, only: file_exist, check_nml_error,                &
                            error_mesg, FATAL, WARNING,                 &
                            mpp_pe, mpp_root_pe, fms_init, fms_end,     &
                            stdlog, write_version_number,               &
                            open_namelist_file, open_restart_file,      &
                            mpp_clock_id, mpp_clock_begin,              &
                            mpp_clock_end, CLOCK_COMPONENT, set_domain, nullify_domain
use       fms_io_mod, only: fms_io_exit

use  mpp_domains_mod, only: domain2d
use       mpp_io_mod, only: mpp_open, mpp_close, MPP_ASCII, MPP_OVERWR, &
                            MPP_SEQUENTIAL, MPP_SINGLE, MPP_RDONLY, MPP_DELETE

use diag_manager_mod, only: diag_manager_init, diag_manager_end, get_base_date

use  field_manager_mod, only: MODEL_ATMOS
use tracer_manager_mod, only: register_tracers
use       memutils_mod, only: print_memuse_stats

implicit none

!-----------------------------------------------------------------------

character(len=128), parameter :: version = &
'$Id: atmos_model.F90,v 14.0 2007/03/15 21:59:43 fms Exp $'

character(len=128), parameter :: tag = &
'$Name: nalanda_2007_04 $'

!-----------------------------------------------------------------------
!       ----- model time -----
! there is no calendar associated with model of this type
! therefore, year=0, month=0 are assumed

   type (time_type) :: Time, Time_init, Time_end, Time_step_atmos
   integer :: num_atmos_calls, na

! ----- model initial date -----

   integer :: date_init(6) ! note: year=month=0

! ----- timing flags -----

   integer :: id_init, id_loop, id_end
   integer, parameter :: timing_level = 1

!-----------------------------------------------------------------------
   character(len=80) :: text
!-----------------------------------------------------------------------
   type(domain2d), save :: atmos_domain  ! This variable must be treated as read-only
!-----------------------------------------------------------------------

      integer, dimension(4) :: current_time = (/ 0, 0, 0, 0 /)
      integer :: days=0, hours=0, minutes=0, seconds=0
      integer :: dt_atmos = 0
      integer :: memuse_interval = 72

      namelist /main_nml/ current_time, dt_atmos,  &
                          days, hours, minutes, seconds, memuse_interval

!#######################################################################

 call fms_init ( )
 call atmos_model_init 

!   ------ atmosphere integration loop -------

    call mpp_clock_begin (id_loop)

    do na = 1, num_atmos_calls

       call atmosphere (Time)

       Time = Time + Time_step_atmos

       if(modulo(na,memuse_interval) == 0) then
         write( text,'(a,i4)' )'Main loop at timestep=',na
         call print_memuse_stats(text)
       endif

    enddo

    call mpp_clock_end (id_loop)

!   ------ end of atmospheric time step loop -----

 call atmos_model_end
 call fms_io_exit
 call fms_end

contains

!#######################################################################

   subroutine atmos_model_init

!-----------------------------------------------------------------------
    integer :: unit, ierr, io
    integer :: ntrace, ntprog, ntdiag, ntfamily
    integer :: date(6)
    type (time_type) :: Run_length
!-----------------------------------------------------------------------
!----- initialization timing identifiers ----

 id_init = mpp_clock_id ('MAIN: initialization', grain=CLOCK_COMPONENT)
 id_loop = mpp_clock_id ('MAIN: time loop'     , grain=CLOCK_COMPONENT)
 id_end  = mpp_clock_id ('MAIN: termination'   , grain=CLOCK_COMPONENT)

 call mpp_clock_begin (id_init)

!-------------------------------------------
! how many tracers have been registered?
!  (will print number below)
   call register_tracers ( MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfamily )


!----- read namelist -------

   unit = open_namelist_file ( )
   ierr=1; do while (ierr /= 0)
          read  (unit, nml=main_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'main_nml')
   enddo
10 call mpp_close (unit)

!----- write namelist to logfile -----

   call write_version_number (version,tag)
   if ( mpp_pe() == mpp_root_pe() ) write (stdlog(), nml=main_nml)

   if (dt_atmos == 0) then
     call error_mesg ('program atmos_model', 'dt_atmos has not been specified', FATAL)
   endif

!----- read restart file -----

   if (file_exist('INPUT/atmos_model.res')) then
       call mpp_open (unit, 'INPUT/atmos_model.res', action=MPP_RDONLY, nohdrs=.true.)
       read  (unit,*) date
       call mpp_close (unit)
   else
    ! use namelist time if restart file does not exist
      date(1:2) = 0
      date(3:6) = current_time
   endif

!----- write current/initial date actually used to logfile file -----

    if ( mpp_pe() == mpp_root_pe() ) then
      write (stdlog(),16) date(3:6)
    endif

 16 format ('  current time used = day',i5,' hour',i3,2(':',i2.2)) 

!  print number of tracers to logfile
   if (mpp_pe() == mpp_root_pe()) then
        write (stdlog(), '(a,i3)') 'Number of tracers =', ntrace
        write (stdlog(), '(a,i3)') 'Number of prognostic tracers =', ntprog
        write (stdlog(), '(a,i3)') 'Number of diagnostic tracers =', ntdiag
   endif

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

    call diag_manager_init

!----- always override initial/base date with diag_manager value -----

!----- get the base date in the diag_table from the diag_manager ----
!      this base date is typically the starting date for the
!      experiment and is subtracted from the current date

    call get_base_date ( date_init(1), date_init(2), date_init(3), &
                         date_init(4), date_init(5), date_init(6)  )

  ! make sure base date does not have a year or month specified
    if ( date_init(1)+date_init(2) /= 0 ) then
         call error_mesg ('program atmos_model', 'invalid base base - &
                          &must have year = month = 0', FATAL)
    endif

!----- set initial and current time types ------
!----- set run length and compute ending time -----

    Time_init  = set_time (date_init(4)*3600 + date_init(5)*60 + date_init(6), date_init(3))
    Time       = set_time (date     (4)*3600 + date     (5)*60 + date     (6), date     (3))
    Run_length = set_time (       hours*3600 +      minutes*60 +      seconds, days        )
    Time_end   = Time + Run_length

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

      call mpp_open (unit, 'time_stamp.out', form=MPP_ASCII, action=MPP_OVERWR, &
                     access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date

!     compute ending time in days,hours,minutes,seconds
      call get_time ( Time_end, date(6), date(3) )  ! gets sec,days
      date(4) = date(6)/3600; date(6) = date(6) - date(4)*3600
      date(5) = date(6)/60  ; date(6) = date(6) - date(5)*60

      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date

      call mpp_close (unit)

  20  format (6i7,2x,'day')   ! can handle day <= 999999

!-----------------------------------------------------------------------
!--- compute the time steps ---
!    determine number of iterations through the time integration loop 
!    must be evenly divisible

      Time_step_atmos = set_time (dt_atmos,0)
      num_atmos_calls = Run_length / Time_step_atmos

!-----------------------------------------------------------------------
!----- initial (base) time must not be greater than current time -----

   if ( Time_init > Time ) call error_mesg ('program atmos_model',  &
                   'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of atmos time step ------

   if ( num_atmos_calls * Time_step_atmos /= Run_length )  &
        call error_mesg ('program atmos_model',  &
           'run length must be multiple of atmosphere time step', FATAL)
   
!-----------------------------------------------------------------------
!------ initialize atmospheric model ------

      call atmosphere_init (Time_init, Time, Time_step_atmos)
      call atmosphere_domain(atmos_domain)

!-----------------------------------------------------------------------
!   open and close dummy file in restart dir to check if dir exists

      call mpp_open  (unit, 'RESTART/file' )
      call mpp_close (unit, action=MPP_DELETE)

!  ---- terminate timing ----
   call mpp_clock_end (id_init)

!-----------------------------------------------------------------------

   call print_memuse_stats('atmos_model_init')
   end subroutine atmos_model_init

!#######################################################################

   subroutine atmos_model_end

   integer :: unit, date(6)
!-----------------------------------------------------------------------
   call mpp_clock_begin (id_end)

      call atmosphere_end

!----- compute current time in days,hours,minutes,seconds -----

      date(1:2) = 0
      call get_time ( Time, date(6), date(3) )
      date(4) = date(6)/3600; date(6) = date(6) - date(4)*3600
      date(5) = date(6)/60  ; date(6) = date(6) - date(5)*60

!----- check time versus expected ending time ----

      if (Time /= Time_end) call error_mesg ('program atmos_model',  &
              'final time does not match expected ending time', WARNING)

!----- write restart file ------

      if ( mpp_pe() == mpp_root_pe() ) then
           call mpp_open (unit, 'RESTART/atmos_model.res', form=MPP_ASCII, action=MPP_OVERWR, &
                          access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )
           write (unit,'(6i6,8x,a)') date, &
                 'Current model time: year, month, day, hour, minute, second'
           call mpp_close (unit)
      endif

!----- final output of diagnostic fields ----
      call set_domain(atmos_domain)  ! This assumes all output fields are on the atmos domain

      call diag_manager_end (Time)

      call nullify_domain()

      call mpp_clock_end (id_end)
!-----------------------------------------------------------------------

   end subroutine atmos_model_end

!#######################################################################
! routines to set/get date when no calendar is set (i.e., yr=0 and mo=0)
!#######################################################################

end program atmos_model

