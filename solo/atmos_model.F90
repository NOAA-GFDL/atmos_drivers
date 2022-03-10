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

program atmos_model

!-----------------------------------------------------------------------
!
!  Main program for running a stand-alone atmospheric dynamical core.
!
!-----------------------------------------------------------------------

use FMS
use FMSconstants
use atmosphere_mod, only: atmosphere_init, atmosphere_end, atmosphere, atmosphere_domain

!--- FMS old io
use fms_io_mod, only: fms_io_exit!< This can't be removed until fms_io is not used at all

implicit none

!-----------------------------------------------------------------------

!---- version number
! Include variable "version" to be written to log file.
#include <file_version.h>

!-----------------------------------------------------------------------
!       ----- model time -----

   integer :: calendartype = INVALID_CALENDAR
   type (time_type) :: Time, Time_init, Time_end, Time_step_atmos
   integer :: num_atmos_calls, na
   type (time_type) :: Time_tmp ! used to facilitate some operations on time_type variables.
   integer :: days_tmp          ! used together with Time_tmp

! ----- model initial date -----

   integer :: date_init(6)

! ----- timing flags -----

   integer :: id_init, id_loop, id_end
   integer, parameter :: timing_level = 1

!-----------------------------------------------------------------------
   character(len=80) :: text
!-----------------------------------------------------------------------
   type(domain2d), save :: atmos_domain  ! This variable must be treated as read-only
!-----------------------------------------------------------------------
   character(len=17) :: calendar = 'no_calendar      '  !< The calendar type used by the current integration.  Valid values are
                                                        !! consistent with the time_manager module: 'gregorian', 'julian',
                                                        !! 'noleap', or 'thirty_day'. All values must be lower case.
   integer, dimension(4) :: current_time = (/ 0, 0, 0, 0/) !< The current time integration starts with (DD,HH,MM,SS)
   integer :: years=0    !< Number of years the current integration will be run
   integer :: months=0   !< Number of months the current integration will be run
   integer :: days=0     !< Number of days the current integration will be run
   integer :: hours=0    !< Number of hours the current integration will be run
   integer :: minutes=0  !< Number of minutes the current integration will be run
   integer :: seconds=0  !< Number of seconds the current integration will be run
   integer :: dt_atmos = 0  !< Atmospheric model time step in seconds
   integer :: memuse_interval = 72  !< Output memory statistics every <N> time steps
   integer :: atmos_nthreads = 1  !< Number of OpenMP threads to use in the atmosphere
   logical :: use_hyper_thread = .false.  !< If .TRUE>, affinity placement (if activated) will consider virtual cores
                                          !! in the placement algorithm

   namelist /main_nml/ calendar, current_time, dt_atmos,  &
                       years, months, days, hours, minutes, seconds, &
                       memuse_interval, atmos_nthreads, use_hyper_thread

!#######################################################################

 call fms_init ()
 call atmos_model_init ()

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
    integer :: time_stamp_unit !< Unif of the time_stamp file
    integer :: stdout_unit !< Unit of stdout file
    integer :: ascii_unit  !< Unit of a dummy ascii file
    character(len=:), dimension(:), allocatable :: restart_file !< Restart file saved as a string
    integer :: ierr, io, logunit
    integer :: ntrace, ntprog, ntdiag, ntfamily
    integer :: date(6)
    type (time_type) :: Run_length
!$  integer :: omp_get_thread_num
    integer :: yr, mo, total_days, total_seconds
!-----------------------------------------------------------------------
!----- initialization timing identifiers ----

 id_init = mpp_clock_id ('MAIN: initialization', grain=CLOCK_COMPONENT)
 id_loop = mpp_clock_id ('MAIN: time loop'     , grain=CLOCK_COMPONENT)
 id_end  = mpp_clock_id ('MAIN: termination'   , grain=CLOCK_COMPONENT)

 logunit = stdlog()

 call mpp_clock_begin (id_init)

!-------------------------------------------
! how many tracers have been registered?
!  (will print number below)
   call register_tracers ( MODEL_ATMOS, ntrace, ntprog, ntdiag, ntfamily )


!----- read namelist -------

   read (input_nml_file, nml=main_nml, iostat=io)
   ierr = check_nml_error(io, 'main_nml')

!----- write namelist to logfile -----

   call write_version_number ("program atmos_model", version)
   if ( mpp_pe() == mpp_root_pe() ) write (logunit, nml=main_nml)

   if (dt_atmos == 0) then
     call error_mesg ('program atmos_model', 'dt_atmos has not been specified', FATAL)
   endif

   if(lowercase(calendar) == 'no_calendar') then
     calendartype = NO_CALENDAR
   else if(lowercase(calendar) == 'thirty_day_months') then
     calendartype = THIRTY_DAY_MONTHS
   else if(lowercase(calendar) == 'noleap') then
     calendartype = NOLEAP
   else if(lowercase(calendar) == 'julian') then
     calendartype = JULIAN
   else if(lowercase(calendar) == 'gregorian') then
     calendartype = GREGORIAN
   else
     call error_mesg ('program atmos_model', trim(calendar)//' is an invalid value for calendar', FATAL)
   endif
   call set_calendar_type(calendartype)

!----- read restart file -----

   if (file_exists('INPUT/atmos_model.res')) then
       call ascii_read('INPUT/atmos_model.res', restart_file)
       read  (restart_file(1),*) date
       deallocate(restart_file)
   else
    ! use namelist time if restart file does not exist
      if(calendartype == NO_CALENDAR) then
        date(1:2) = 0
        date(3:6) = current_time
      else
        date(1:3) = 1
        date(4:6) = current_time(2:4)
      endif
   endif

!----- write current/initial date actually used to logfile file -----

    if ( mpp_pe() == mpp_root_pe() ) then
      write (logunit,16) date
    endif

 16 format ('  current time used = ',i4,'-',i2,'-',i2,1x,i3,2(':',i2.2))

!  print number of tracers to logfile
   if (mpp_pe() == mpp_root_pe()) then
        write (logunit, '(a,i3)') 'Number of tracers =', ntrace
        write (logunit, '(a,i3)') 'Number of prognostic tracers =', ntprog
        write (logunit, '(a,i3)') 'Number of diagnostic tracers =', ntdiag
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

!----- set initial and current time types ------
!----- set run length and compute ending time -----
    if(calendartype == NO_CALENDAR) then
       Time_init  = set_time(date_init(4)*int(SECONDS_PER_HOUR)+date_init(5)*int(SECONDS_PER_MINUTE)+date_init(6),date_init(3))
       Time       = set_time(date     (4)*int(SECONDS_PER_HOUR)+date     (5)*int(SECONDS_PER_MINUTE)+date     (6),date     (3))
       Run_length = set_time(       hours*int(SECONDS_PER_HOUR)+     minutes*int(SECONDS_PER_MINUTE)+     seconds,days        )
    else
       Time_init  = set_date(date_init(1),date_init(2), date_init(3),date_init(4),date_init(5),date_init(6))
       Time       = set_date(date(1),date(2),date(3),date(4),date(5),date(6))
       Time_tmp = Time
       total_days = 0
       do yr=1,years
         days_tmp = days_in_year(Time_tmp)
         total_days = total_days + days_tmp
         Time_tmp = Time_tmp + set_time(0,days_tmp)
       enddo
       do mo=1,months
         days_tmp = days_in_month(Time_tmp)
         total_days = total_days + days_tmp
         Time_tmp = Time_tmp + set_time(0,days_tmp)
       enddo
       total_days = total_days + days
       total_seconds = hours*3600 + minutes*60 + seconds
       Run_length    = set_time (total_seconds,total_days)
    endif
    Time_end   = Time + Run_length
!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

      if ( mpp_pe() == mpp_root_pe() ) open (newunit=time_stamp_unit, file='time_stamp.out', status='replace', form='formatted')
      if ( mpp_pe() == mpp_root_pe() ) write (time_stamp_unit,20) date

!     compute ending time in days,hours,minutes,seconds
      if(calendartype == NO_CALENDAR) then
        call get_time ( Time_end, date(6), date(3) )  ! gets sec,days
        date(4) = date(6)/int(SECONDS_PER_HOUR); date(6) = date(6) - date(4)*int(SECONDS_PER_HOUR)
        date(5) = date(6)/int(SECONDS_PER_MINUTE)  ; date(6) = date(6) - date(5)*int(SECONDS_PER_MINUTE)
      else
        call get_date(Time_end, date(1), date(2), date(3), date(4), date(5), date(6))
      endif
      if ( mpp_pe() == mpp_root_pe() ) write (time_stamp_unit,20) date

      if ( mpp_pe() == mpp_root_pe() ) close (time_stamp_unit)

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

!$ call omp_set_num_threads(atmos_nthreads)
   call fms_affinity_set('Atmos Program', use_hyper_thread, atmos_nthreads)
   if (mpp_pe() .eq. mpp_root_pe()) then
     stdout_unit=stdout()
     write(stdout_unit,*) ' starting ',atmos_nthreads,' OpenMP threads per MPI-task'
     call flush(stdout_unit)
   endif

   call atmosphere_init (Time_init, Time, Time_step_atmos)
   call atmosphere_domain(atmos_domain)

   !----- add back in if needed -----
   !call data_override_init(Atm_domain_in = atmos_domain)

!-----------------------------------------------------------------------
!   open and close dummy file in restart dir to check if dir exists
   call mpp_set_current_pelist()
   if ( mpp_pe().EQ.mpp_root_pe() ) then
     open(newunit = ascii_unit, file='RESTART/file', status='replace', form='formatted')
     close(ascii_unit,status="delete")
   endif

!  ---- terminate timing ----
   call mpp_clock_end (id_init)

!-----------------------------------------------------------------------

   call print_memuse_stats('atmos_model_init')
   end subroutine atmos_model_init

!#######################################################################

   subroutine atmos_model_end

   integer :: restart_file_unit, date(6)
!-----------------------------------------------------------------------
      call mpp_clock_begin (id_end)

      call atmosphere_end

!----- compute current time in days,hours,minutes,seconds -----

      if(calendartype == NO_CALENDAR) then
        date(1:2) = 0
        call get_time ( Time, date(6), date(3) )
        date(4) = date(6)/int(SECONDS_PER_HOUR); date(6) = date(6) - date(4)*int(SECONDS_PER_HOUR)
        date(5) = date(6)/int(SECONDS_PER_MINUTE); date(6) = date(6) - date(5)*int(SECONDS_PER_MINUTE)
      else
        call get_date(Time, date(1), date(2), date(3), date(4), date(5), date(6))
      endif

!----- check time versus expected ending time ----

      if (Time /= Time_end) call error_mesg ('program atmos_model',  &
              'final time does not match expected ending time', WARNING)

!----- write restart file ------

      if ( mpp_pe() == mpp_root_pe() ) then
           open(newunit = restart_file_unit, file='RESTART/atmos_model.res', status='replace', form='formatted')
           write (restart_file_unit,'(6i6,8x,a)') date, &
                 'Current model time: year, month, day, hour, minute, second'
           close (restart_file_unit)
      endif

!----- final output of diagnostic fields ----
      call diag_manager_end (Time)

      call mpp_clock_end (id_end)
!-----------------------------------------------------------------------

   end subroutine atmos_model_end

!#######################################################################
! routines to set/get date when no calendar is set (i.e., yr=0 and mo=0)
!#######################################################################

end program atmos_model
