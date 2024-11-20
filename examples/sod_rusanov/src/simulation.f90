!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,            only: WP
  use geometry,             only: cfg
  use rusanov_class,        only: rusanov, EXTRAPOLATE
  use hyperbolic_euler,     only: make_euler_rusanov, euler_tocons,           &
    & euler_tophys, SOD_PHYS_L, SOD_PHYS_R, DIATOMIC_GAMMA
  use timetracker_class,    only: timetracker
  use ensight_class,        only: ensight
  use partmesh_class,       only: partmesh
  use event_class,          only: event
  use monitor_class,        only: monitor
  use string, only: str_medium
  implicit none
  private

  !> Flow solver and a time tracker
  type(rusanov), public :: fs
  type(timetracker), public :: time

  !> Ensight postprocessing
  character(len=str_medium) :: ens_out_name
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> physical value arrays for ensight output
  real(WP), dimension(:,:,:,:), pointer :: phys_out

  !> simulation monitor file
  type(monitor) :: mfile, cflfile, consfile

  !> simulation control functions
  public :: simulation_init, simulation_run, simulation_final

contains

  !> update phys
  subroutine update_phys_out()
    implicit none
    integer :: i, j, k

    do k = cfg%kmin_, cfg%kmax_
      do j = cfg%jmin_, cfg%jmax_
        do i = cfg%imin_, cfg%imax_
          call euler_tophys(DIATOMIC_GAMMA, fs%Uc(:,i,j,k), phys_out(1:5,i,j,k))
          phys_out(6,i,j,k) = sqrt(sum(phys_out(2:4,i,j,k)**2) / (DIATOMIC_GAMMA * phys_out(5,i,j,k) / phys_out(1,i,j,k)))
        end do
      end do
    end do

    call cfg%sync(phys_out)

  end subroutine

  !> functions localize the sides of the domain
  function side_locator_x_l(pg, i, j, k) result(loc)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i, j, k
    logical :: loc
    loc = .false.
    if (i.lt.pg%imin) loc = .true.
  end function side_locator_x_l
  function side_locator_x_r(pg, i, j, k) result(loc)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i, j, k
    logical :: loc
    loc = .false.
    if (i.gt.pg%imax) loc = .true.
  end function side_locator_x_r
  function side_locator_y_l(pg, i, j, k) result(loc)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i, j, k
    logical :: loc
    loc = .false.
    if (j.lt.pg%jmin) loc = .true.
  end function side_locator_y_l
  function side_locator_y_r(pg, i, j, k) result(loc)
    use pgrid_class, only: pgrid
    implicit none
    class(pgrid), intent(in) :: pg
    integer, intent(in) :: i, j, k
    logical :: loc
    loc = .false.
    if (j.gt.pg%jmax) loc = .true.
  end function side_locator_y_r

  !> initialization of problem solver
  subroutine simulation_init()
    use param, only: param_read
    use messager, only: die
    implicit none

    !TODO need to figure out what the right interaction with this is
    initialize_timetracker: block

      time = timetracker(amRoot=cfg%amRoot)
      call param_read('Max cfl number', time%cflmax)
      time%dt = 1e-3_WP
      time%itmax = 1

    end block initialize_timetracker

    ! create a single-phase flow solver
    create_and_initialize_flow_solver: block

      ! call constructor
      fs = make_euler_rusanov(cfg)

      ! set bcs
      call fs%add_bcond(name='openxl' , type=EXTRAPOLATE,                     &
        & locator=side_locator_x_l, dir='xl')
      call fs%add_bcond(name='openxr' , type=EXTRAPOLATE,                     &
        & locator=side_locator_x_r, dir='xr')
      call fs%add_bcond(name='openyl' , type=EXTRAPOLATE,                     &
        & locator=side_locator_y_l, dir='yl')
      call fs%add_bcond(name='openyr' , type=EXTRAPOLATE,                     &
        & locator=side_locator_y_r, dir='yr')

    end block create_and_initialize_flow_solver

    ! prepare initial fields
    initialize_fields: block
      real(WP), dimension(5) :: lval, rval, val
      real(WP) :: angle, x0, y0, x1, y1, w
      integer :: i, j, k

      call param_read('Mesh angle', angle)

      write(ens_out_name,'(A,F0.1)') 'SodAngle', angle

      angle = angle * 2.0_WP * 4.0_WP * atan(1.0_WP) / 360

      call euler_tocons(DIATOMIC_GAMMA, SOD_PHYS_L, lval)
      call euler_tocons(DIATOMIC_GAMMA, SOD_PHYS_R, rval)

        do j = cfg%jmino_, cfg%jmaxo_
          y0 = cfg%y(j); y1 = cfg%y(j+1);
          do i = cfg%imino_, cfg%imaxo_
            x0 = cfg%x(i); x1 = cfg%x(i+1);
            if (cos(angle) * x1 + sin(angle) * y1 .le. 0.0_WP) then
              val = lval
            else if (cos(angle) * x0 + sin(angle) * y0 .ge. 0.0_WP) then
              val = rval
            else
              w = 0.5_WP * (cotan(angle) * x0  + y0) * (tan(angle) * y0 + x0)
              if (-cotan(angle) * x0 .gt. y1) then
                w = w - 0.5_WP * (tan(angle) * y1 + x0) * (cotan(angle) * x0 + y1)
              end if
              if (-tan(angle) * y0 .gt. x1) then
                w = w - 0.5_WP * (tan(angle) * y0 + x1) * (cotan(angle) * x1 + y0)
              end if
              w = w / (cfg%dx(i) * cfg%dy(j))
              val = lval * w + rval * (1.0_WP - w)
            end if
            do k = cfg%kmino_, cfg%kmaxo_
              fs%Uc(:,i,j,k) = val
            end do
          end do
        end do

    end block initialize_fields

    ! Add Ensight output
    create_ensight: block
      use string, only: str_short
      real(WP), dimension(:,:,:), pointer :: ptr1, ptr2, ptr3

      ! create array to hold physical coordinates
      allocate(phys_out(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,        &
        & cfg%kmino_:cfg%kmaxo_))

      ! Create Ensight output from cfg
      ens_out = ensight(cfg=cfg, name=ens_out_name)

      ! Create event for Ensight output
      ens_evt = event(time=time, name='Ensight output')
      call param_read('Ensight output period', ens_evt%tper)

      ! Add variables to output
      ptr1 => phys_out(1,:,:,:)
      call ens_out%add_scalar('density', ptr1)
      ptr1 => phys_out(2,:,:,:)
      ptr2 => phys_out(3,:,:,:)
      ptr3 => phys_out(4,:,:,:)
      call ens_out%add_vector('velocity', ptr1, ptr2, ptr3)
      ptr1 => phys_out(5,:,:,:)
      call ens_out%add_scalar('pressure', ptr1)
      ptr1 => fs%params(1,:,:,:)
      call ens_out%add_scalar('gamma', ptr1)
      ptr1 => phys_out(6,:,:,:)
      call ens_out%add_scalar('Ma', ptr1)

      ! Output to ensight
      if (ens_evt%occurs()) then
        call update_phys_out()
        call ens_out%write_data(time%t)
      end if

    end block create_ensight

    ! Create a monitor file
    create_monitor: block
      use string, only: str_short
      real(WP), pointer :: real_ptr

      ! Prepare some info about fields
      call fs%get_cfl(time%dt, time%cfl)
      call fs%get_range()

      ! Create simulation monitor
      mfile = monitor(fs%cfg%amRoot, 'simulation')
      call mfile%add_column(time%n, 'Timestep number')
      call mfile%add_column(time%t, 'Time')
      call mfile%add_column(time%dt, 'Timestep size')
      call mfile%add_column(time%cfl, 'Maximum CFL')
      ! not sure why this doesn't work?
      !call mfile%add_column(real_ptr, fields(i:i)//'min')
      !call mfile%add_column(real_ptr, fields(i:i)//'max')
      real_ptr => fs%Umin(1)
      call mfile%add_column(real_ptr, 'dens_min')
      real_ptr => fs%Umax(1)
      call mfile%add_column(real_ptr, 'dens_max')
      real_ptr => fs%Umin(2)
      call mfile%add_column(real_ptr, 'momx_min')
      real_ptr => fs%Umax(2)
      call mfile%add_column(real_ptr, 'momx_max')
      real_ptr => fs%Umin(3)
      call mfile%add_column(real_ptr, 'momy_min')
      real_ptr => fs%Umax(3)
      call mfile%add_column(real_ptr, 'momy_max')
      real_ptr => fs%Umin(4)
      call mfile%add_column(real_ptr, 'momz_min')
      real_ptr => fs%Umax(4)
      call mfile%add_column(real_ptr, 'momz_max')
      real_ptr => fs%Umin(5)
      call mfile%add_column(real_ptr, 'totE_min')
      real_ptr => fs%Umax(5)
      call mfile%add_column(real_ptr, 'totE_max')
      call mfile%write()

      ! Create CFL monitor
      cflfile = monitor(fs%cfg%amRoot, 'cfl')
      call cflfile%add_column(time%n, 'Timestep number')
      call cflfile%add_column(time%t, 'Time')
      call cflfile%add_column(fs%CFL_x, 'CFLx')
      call cflfile%add_column(fs%CFL_y, 'CFLy')
      call cflfile%add_column(fs%CFL_z, 'CFLz')
      call cflfile%write()

      ! Create conservation monitor
      consfile = monitor(fs%cfg%amRoot, 'conservation')
      call consfile%add_column(time%n, 'Timestep number')
      call consfile%add_column(time%t, 'Time')
      real_ptr => fs%Uint(1)
      call consfile%add_column(real_ptr, 'dens_int')
      real_ptr => fs%Uint(2)
      call consfile%add_column(real_ptr, 'momx_int')
      real_ptr => fs%Uint(3)
      call consfile%add_column(real_ptr, 'momy_int')
      real_ptr => fs%Uint(4)
      call consfile%add_column(real_ptr, 'momz_int')
      real_ptr => fs%Uint(5)
      call consfile%add_column(real_ptr, 'totE_int')
      call consfile%write()

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
      call fs%apply_bcond(time%t, time%dt)
      fs%dU(:, :, :, :) = 0.0_WP
      call fs%calc_dU_x(0.5 * time%dt)
      fs%Uc = fs%Uc + fs%dU
      call fs%apply_bcond(time%t, time%dt)
      fs%dU(:, :, :, :) = 0.0_WP
      call fs%calc_dU_y(time%dt)
      fs%Uc = fs%Uc + fs%dU
      call fs%apply_bcond(time%t, time%dt)
      fs%dU(:, :, :, :) = 0.0_WP
      call fs%calc_dU_x(0.5 * time%dt)
      fs%Uc = fs%Uc + fs%dU

      ! Output to ensight
      if (ens_evt%occurs()) then
        call update_phys_out()
        call ens_out%write_data(time%t)
      end if

      ! Perform and output monitoring
      call fs%get_range()
      call mfile%write()
      call cflfile%write()
      call consfile%write()

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

