
program atmos_model

!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------

use   atmosphere_mod, only: atmosphere_init, atmosphere_end, atmosphere

use time_manager_mod, only: time_type, set_time, get_time,  &
                            operator(+), operator (<), operator (>), &
                            operator (/=), operator (/), operator (*)

use    utilities_mod, only: open_file, file_exist, check_nml_error,   &
                            error_mesg, FATAL, WARNING,               &
                            get_my_pe, utilities_init, utilities_end, &
                            close_file, check_system_clock,           &
                            mpp_clock_init, mpp_clock_begin,          &
                            mpp_clock_end

use diag_manager_mod, only: diag_manager_init, diag_manager_end, get_base_date


implicit none

!-----------------------------------------------------------------------

character(len=128), parameter :: version = &
'$Id: atmos_model.F90,v 1.4 2001/04/13 15:31:11 fms Exp $'

character(len=128), parameter :: tag = &
'$Name: eugene $'

!-----------------------------------------------------------------------
! ----- model time -----

   type (time_type) :: Time, Time_init, Time_end, Time_step_atmos
   integer :: num_atmos_calls, na

! ----- coupled model initial date -----

   integer :: date_init(6)

! ----- timing flags -----

   integer :: id_init, id_loop, id_end
   integer, parameter :: timing_level = 1

!-----------------------------------------------------------------------

      integer, dimension(4) :: current_time = (/ 0, 0, 0, 0 /)
      logical :: override = .false.  ! override restart values for date
      integer :: days=0, hours=0, minutes=0, seconds=0
      integer :: dt_atmos = 0

      namelist /main_nml/ current_time, override, dt_atmos, &
                          days, hours, minutes, seconds

!#######################################################################

 call utilities_init ( )

 call atmos_model_init 

 call check_system_clock ('END OF INITIALIZATION')
 call mpp_clock_begin (id_loop)

!   ------ atmosphere integration loop -------

    do na = 1, num_atmos_calls

       call atmosphere (Time)

       Time = Time + Time_step_atmos

    enddo

!   ------ end of atmospheric time step loop -----

 call mpp_clock_end (id_loop)
 call check_system_clock ('END OF TIME LOOP')

 call atmos_model_end
 call utilities_end

 stop

contains

!#######################################################################

   subroutine atmos_model_init

!-----------------------------------------------------------------------
    integer :: total_days, total_seconds, unit, log_unit, ierr, io, id, jd, kd
    integer :: date(6)
    type (time_type) :: Run_length
    logical :: use_namelist
!-----------------------------------------------------------------------
!----- initialization timing identifiers ----

 id_init = mpp_clock_init ('atmos_model', 'atmos_model_init', timing_level)
 id_loop = mpp_clock_init ('atmos_model', 'time_loop'       , timing_level)
 id_end  = mpp_clock_init ('atmos_model', 'atmos_model_end' , timing_level)

 call mpp_clock_begin (id_init)

!----- read namelist -------

   unit = open_file ('input.nml', action='read')
   ierr=1; do while (ierr /= 0)
          read  (unit, nml=main_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'main_nml')
   enddo
10 call close_file (unit)

!----- write namelist to logfile (close log_unit later) -----

   log_unit = open_file ('logfile.out', action='append')
   if ( get_my_pe() == 0 ) then
     write (log_unit,'(/,80("="),/(a))') trim(version), trim(tag)
     write (log_unit, nml=main_nml)
   endif

   if(dt_atmos == 0) then
     call error_mesg ('program atmos_model', 'dt_atmos has not been specified', FATAL)
   endif

!----- read restart file -----

   if (file_exist('INPUT/atmos_model.res')) then
       unit = open_file ('INPUT/atmos_model.res', form='native', action='read')
       read  (unit) date
       call close_file (unit)
       use_namelist = .false.
   else
       use_namelist = .true.
   endif

!----- override date with namelist values ------
!----- (either no restart or override flag on) ---

 if ( use_namelist .or. override ) then
      date(1:2) = 0
      date(3:6) = current_time
 endif

!----- write current/initial date actually used to logfile file -----

    if ( get_my_pe() == 0 ) then
      write (log_unit,16) date(3:6)
    endif
      call close_file (log_unit)

 16 format ('  current time used = day',i5,' hour',i3,2(':',i2.2)) 

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

    call diag_manager_init

!----- always override initial/base date with diag_manager value -----

    call get_base_date ( date_init(1), date_init(2), date_init(3), &
                         date_init(4), date_init(5), date_init(6)  )

    if ( date_init(1)+date_init(2) /= 0 ) then
         call error_mesg ('program atmos_model', 'invalid base base - &
                          &must have year = month = 0', FATAL)
    endif

!----- set initial and current time types ------

    total_seconds = date_init(4)*3600 + date_init(5)*60 + date_init(6)
    Time_init = set_time ( total_seconds, date_init(3) )

    total_seconds = date(4)*3600 + date(5)*60 + date(6)
    Time      = set_time ( total_seconds, date(3) )

!----- compute the ending time -----

      total_seconds = hours*3600 + minutes*60 + seconds
      Run_length    = set_time ( total_seconds, days )
      Time_end      = Time + Run_length

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

      unit = open_file ('time_stamp.out', action='write')

      if ( get_my_pe() == 0 ) write (unit,20) date

!     ---- compute ending time in days,hours,minutes,seconds ----
      call get_time ( Time_end, date(6), date(3) )
      date(4) = date(6)/3600; date(6) = date(6) - date(4)*3600
      date(5) = date(6)/60  ; date(6) = date(6) - date(5)*60

      if ( get_my_pe() == 0 ) write (unit,20) date

      call close_file (unit)

  20  format (6i4,2x,'day')

!-----------------------------------------------------------------------
!----- compute the time steps ------
!----- determine maximum number of iterations per loop ------

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

!-----------------------------------------------------------------------
!---- open and close output restart to make sure directory is there ----

      unit = open_file ('RESTART/atmos_model.res',  &
                        form='native', action='write')
      call close_file (unit, status='delete')


!  ---- terminate timing ----
   call mpp_clock_end (id_init)

!-----------------------------------------------------------------------

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

      if ( get_my_pe() == 0 ) then
           unit = open_file ('RESTART/atmos_model.res',  &
                             form='native', action='write')
           write (unit) date
           call close_file (unit)
      endif

!----- final output of diagnostic fields ----

      call diag_manager_end (Time)


      call mpp_clock_end (id_end)
!-----------------------------------------------------------------------

   end subroutine atmos_model_end

!#######################################################################

end program atmos_model

