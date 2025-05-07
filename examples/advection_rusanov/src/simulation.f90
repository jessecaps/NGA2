!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,            only: WP
  use geometry,             only: cfg
  use rusanov_class,        only: rusanov
  use hyperbolic_advection, only: make_advec_rusanov
  use timetracker_class,    only: timetracker
  use ensight_class,        only: ensight
  use partmesh_class,       only: partmesh
  use event_class,          only: event
  use monitor_class,        only: monitor
  implicit none
  private

  !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
  type(rusanov), public :: fs
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> simulation monitor file
  type(monitor) :: mfile, cflfile

  !> simulation control functions
  public :: simulation_init, simulation_run, simulation_final

contains

  !> circles and donuts shape
  ! TODO fill between zero and one on boundaries
  subroutine init_donut(iR, oR, U)
    implicit none
    real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(out) :: U
    real(WP), intent(in) :: iR, oR
    real(WP), dimension(3) :: c, x
    real(WP) :: r2
    integer :: i, j, k

    U(:,:,:) = 0.0_WP

    c = 0.5 * (/ cfg%xL, cfg%yL, cfg%zL /)

    do k = cfg%kmin_, cfg%kmax_
      x(3) = cfg%zm(k)
      do j = cfg%jmin_, cfg%jmax_
        x(2) = cfg%ym(j)
        do i = cfg%imin_, cfg%imax_
          x(1) = cfg%xm(i)
          r2 = sum((c - x)**2)
          if (iR**2 .lt. r2 .and. r2 .lt. oR**2) U(i,j,k) = 1.0_WP
        end do
      end do
    end do

  end subroutine init_donut

  !> bar so we can test 1d behavior
  subroutine init_bar(xmin, xmax, U)
    implicit none
    real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(out) :: U
    real(WP), intent(in) :: xmin, xmax
    integer :: i

    U(:,:,:) = 0.0_WP

    do i = cfg%imin_, cfg%imax_
      if (xmin .lt. cfg%xm(i) .and. cfg%xm(i) .lt. xmax) then
        U(i,:,:) = 1.0_WP
      end if
    end do

  end subroutine

  !> initialization of problem solver
  subroutine simulation_init()
    use param, only: param_read
    use string, only: str_medium
    use messager, only: die
    implicit none

    integer :: numfields
    character, dimension(:), allocatable :: fields

    ! read shapes
    read_shapes: block
      character(len=str_medium) :: fieldbuffer
      integer :: i, j

      call param_read('Shapes', fieldbuffer)
      fieldbuffer = adjustl(fieldbuffer)

      ! check for duplicates
      do i = 1, str_medium
        if (fieldbuffer(i:i) .ne. ' ') then
          do j = i+1, str_medium
            if (fieldbuffer(i:i) .eq. fieldbuffer(j:j)) then
              call die("field name '" // fieldbuffer(i:i) // "' appears twice")
            end if
          end do
        end if
      end do

      ! count valid fields and warn user
      numfields = 0
      do i = 1, str_medium
        select case (fieldbuffer(i:i))
        case ('c')  ! circle
          numfields = numfields + 1
        case ('d')  ! donut
          numfields = numfields + 1
        case ('b')  ! bar
          numfields = numfields + 1
        case (' ')  ! skip spaces
        case default
          write(*, *) "unrecognized shape: '", fieldbuffer(i:i), "'"
        end select
      end do
      if (numfields .eq. 0) call die("must advect at least one shape")

      ! read fields into array
      allocate(fields(numfields))
      numfields = 0
      do i = 1, str_medium
        select case (fieldbuffer(i:i))
        case ('c')  ! circle
          numfields = numfields + 1
          fields(numfields) = 'c'
        case ('d')  ! donut
          numfields = numfields + 1
          fields(numfields) = 'd'
        case ('b')  ! bar
          numfields = numfields + 1
          fields(numfields) = 'b'
        case default
        end select
      end do

    end block read_shapes

    !TODO need to figure out what the right interaction with this is
    initialize_timetracker: block

      time = timetracker(amRoot=cfg%amRoot)
      call param_read('Max cfl number', time%cflmax)
      time%dt = 1e-3_WP
      time%itmax = 1

    end block initialize_timetracker

    ! create a single-phase flow solver without bconds
    create_and_initialize_flow_solver: block
      real(WP), dimension(3) :: vel
      real(WP) :: traversal_length

      ! read velocity direction, find velocity
      call param_read('Advection direction', vel, 'Advec dir')
      traversal_length = sqrt(sum(vel**2)) * cfg%xL
      vel = vel * traversal_length

      ! call constructor
      fs = make_advec_rusanov(cfg, numfields, vel)

    end block create_and_initialize_flow_solver

    ! prepare initial fields
    initialize_fields: block
      integer :: i
      real(WP), dimension(:,:,:), pointer :: fieldptr

      do i = 1, numfields
        fieldptr => fs%Uc(i,:,:,:)
        select case (fields(i))
        case ('c')  ! circle
          call init_donut(0.0_WP, 0.375_WP * cfg%xL, fieldptr)
        case ('d')  ! donut
          call init_donut(0.250_WP * cfg%xL, 0.375_WP * cfg%xL, fieldptr)
        case ('b')  ! bar
          call init_bar(0.333_WP * cfg%xL, 0.666_WP * cfg%xL, fieldptr)
        end select
      end do

    end block initialize_fields

    ! Add Ensight output
    create_ensight: block
      use string, only: str_short
      integer :: i
      character(len=4) :: istr
      real(WP), dimension(:,:,:), pointer :: ptr1, ptr2, ptr3

      ! Create Ensight output from cfg
      ens_out=ensight(cfg=cfg, name='AdvectionTest')

      ! Create event for Ensight output
      ens_evt=event(time=time, name='Ensight output')
      call param_read('Ensight output period', ens_evt%tper)

      ! Add variables to output
      do i = 1, numfields
        ptr1 => fs%Uc(i,:,:,:)
        write(istr,'(I2.2)') i
        call ens_out%add_scalar('U'//istr, ptr1)
      end do

      ptr1 => fs%params(1,:,:,:)
      ptr2 => fs%params(2,:,:,:)
      ptr3 => fs%params(3,:,:,:)
      call ens_out%add_vector('velocity', ptr1, ptr2, ptr3)

      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)

    end block create_ensight

    ! Create a monitor file
    create_monitor: block
      use string, only: str_short
      character(len=str_short) :: fieldname
      integer :: i
      real(WP), pointer :: real_ptr

      ! Prepare some info about fields
      call fs%get_cfl(time%dt,time%cfl)
      call fs%get_range()

      ! Create simulation monitor
      mfile = monitor(fs%cfg%amRoot, 'simulation')
      call mfile%add_column(time%n, 'Timestep number')
      call mfile%add_column(time%t, 'Time')
      call mfile%add_column(time%dt, 'Timestep size')
      call mfile%add_column(time%cfl, 'Maximum CFL')
      do i = 1, numfields
        write(fieldname,'("[",a,"]min")') fields(i)
        real_ptr => fs%Umin(i)
        call mfile%add_column(real_ptr, fieldname)
        write(fieldname,'("[",a,"]max")') fields(i)
        real_ptr => fs%Umax(i)
        call mfile%add_column(real_ptr, fieldname)
      end do
      call mfile%write()

      ! Create CFL monitor
      cflfile = monitor(fs%cfg%amRoot, 'cfl')
      call cflfile%add_column(time%n, 'Timestep number')
      call cflfile%add_column(time%t, 'Time')
      call cflfile%add_column(fs%CFL_x, 'CFLx')
      call cflfile%add_column(fs%CFL_y, 'CFLy')
      call cflfile%add_column(fs%CFL_z, 'CFLz')
      call cflfile%write()

    end block create_monitor

  end subroutine simulation_init

  !> do the thing
  subroutine simulation_run
    implicit none

    ! Perform time integration
    do while (.not.time%done())

      ! Increment time
      call fs%get_cfl(time%dt, time%cfl)
      call time%adjust_dt()
      call time%increment()

      ! take step (Strang)
      fs%dU(:, :, :, :) = 0.0_WP
      !call fs%apply_bcond(time%t, time%dt)
      call fs%calc_dU_x(0.5 * time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:, :, :, :) = 0.0_WP
      !call fs%apply_bcond(time%t, time%dt)
      call fs%calc_dU_y(time%dt)
      fs%Uc = fs%Uc + fs%dU
      fs%dU(:, :, :, :) = 0.0_WP
      !call fs%apply_bcond(time%t, time%dt)
      call fs%calc_dU_x(0.5 * time%dt)
      fs%Uc = fs%Uc + fs%dU

      ! Output to ensight
      if (ens_evt%occurs()) call ens_out%write_data(time%t)

      ! Perform and output monitoring
      call fs%get_range()
      call mfile%write()
      call cflfile%write()

    end do

  end subroutine simulation_run

  !> Finalize the NGA2 simulation
  subroutine simulation_final
    implicit none

    ! Get rid of all objects - need destructors
    !TODO
    ! monitor
    ! ensight
    ! bcond
    ! timetracker


  end subroutine simulation_final

end module simulation

