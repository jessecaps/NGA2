!>  Rusanov solver class
!>
!>  This first-order Rusanov scheme is extremely diffusive, but foolproof (as
!>  long as a flux function is known.
!>
!>  Originally written by John P Wakefield, Feb 2023.
!>
!>
module rusanov_class
  use precision,      only: WP
  use string,         only: str_medium
  use config_class,   only: config
  use iterator_class, only: iterator
  use hyperbolic,     only: eigenvals_ftype, flux_ftype
  implicit none
  private

  ! expose type/constructor/methods
  public :: rusanov, rus_bc

  ! expose interfaces
  public :: eigenvals_ftype, flux_ftype

  ! List of known available bcond types for this solver
  integer(1), parameter, public :: EXTRAPOLATE = 1_1          !< 0th order extrapolation
  integer(1), parameter, public :: REFLECT = 2_1              !< reflection by velmask

  !> boundary conditions for the hyperbolic solver
  type :: rus_bc
    type(rus_bc), pointer :: next                             !< linked list of bconds
    character(len=str_medium) :: name = 'UNNAMED_BCOND'       !< bcond name (default UNNAMED_BCOND)
    integer(1) :: type                                        !< bcond type
    type(iterator) :: itr                                     !< this is the iterator for the bcond
    integer(1) :: dir                                         !< bcond direction 0-6
  end type rus_bc

  !> solver object definition
  type :: rusanov

    ! config / pgrid object
    class(config), pointer :: cfg                             !< config / pgrid information

    ! name
    character(len=str_medium) :: name = 'UNNAMED_RUSANOV'     !< solver name

    ! system information
    integer :: N                                              !< system dimension
    integer :: P                                              !< number of parameters
    procedure(eigenvals_ftype), pointer, nopass :: evals_x, evals_y, evals_z
    procedure(flux_ftype), pointer, nopass :: flux_x, flux_y, flux_z
    logical, dimension(:), allocatable :: vel_mask_x, vel_mask_y, vel_mask_z

    ! boundary condition list
    integer :: nbc                                            !< number of bcond for our solver
    type(rus_bc), pointer :: first_bc                         !< list of bcond for our solver

    ! flow variables, parameter arrays
    real(WP), dimension(:,:,:,:), pointer :: Uc, dU           !< state variables
    real(WP), dimension(:,:,:,:), pointer :: params           !< params for evals/flux

    ! CFL numbers
    real(WP) :: CFL_x, CFL_y, CFL_z                           !< global CFL numbers (over dt)
    logical :: have_CFL_x, have_CFL_y, have_CFL_z

    ! monitoring quantities
    real(WP), dimension(:), pointer :: Umin, Umax             !< state variable range
    real(WP), dimension(:), pointer :: Uint                   !< integral of state vars over domain
    logical :: have_Urange

  contains

    ! standard interface
    procedure :: print => rusanov_print                       !< output solver to the screen
    procedure :: add_bcond                                    !< add a boundary condition
    procedure :: get_bcond                                    !< get a boundary condition
    procedure :: apply_bcond                                  !< apply all boundary conditions
    procedure :: get_cfl                                      !< get maximum CFL
    procedure :: get_range                                    !< calculate min/max field values
    procedure :: get_min => get_range                         !< compatibility
    procedure :: get_max => get_range                         !< compatibility
    procedure :: calc_dU_x, calc_dU_y, calc_dU_z              !< take step

  end type rusanov

  !> declare constructors
  interface rusanov; procedure rusanov_from_args; end interface

contains

  !! standard interface class functions

  !> default constructor for rusanov-type flow solver
  function rusanov_from_args(cfg, name, N, P, evals_x, evals_y, evals_z,      &
      flux_x, flux_y, flux_z) result(this)
    use messager, only: die
    implicit none
    type(rusanov) :: this
    class(config), target, intent(in) :: cfg
    character(len=*), intent(in), optional :: name
    integer, intent(in) :: N, P
    procedure(eigenvals_ftype), pointer, intent(in) :: evals_x, evals_y, evals_z
    procedure(flux_ftype), pointer, intent(in) :: flux_x, flux_y, flux_z
    integer :: imino, imaxo, jmino, jmaxo, kmino, kmaxo, minmino, maxmaxo

    ! check number of ghost cells
    if (cfg%no .lt. 1) call die("must have at least one ghost cell")

    ! point to pgrid object
    this%cfg => cfg

    ! set the name for the solver
    if (present(name)) this%name = trim(adjustl(name))

    ! system information
    this%N = N; this%P = P;
    this%evals_x => evals_x; this%evals_y => evals_y; this%evals_z => evals_z;
    this%flux_x => flux_x; this%flux_y => flux_y; this%flux_z => flux_z;

    ! boundary condition list
    this%nbc = 0
    this%first_bc => NULL()

    ! get array sizes
    imino = this%cfg%imino_; imaxo = this%cfg%imaxo_;
    jmino = this%cfg%jmino_; jmaxo = this%cfg%jmaxo_;
    kmino = this%cfg%kmino_; kmaxo = this%cfg%kmaxo_;
    minmino = min(imino, jmino, kmino); maxmaxo = max(imaxo, jmaxo, kmaxo)

    ! allocate flow variables and parameter arrays
    allocate(this%params(P,imino:imaxo,jmino:jmaxo,kmino:kmaxo))
    allocate(this%Uc(N,imino:imaxo,jmino:jmaxo,kmino:kmaxo))
    allocate(this%dU(N,imino:imaxo,jmino:jmaxo,kmino:kmaxo))

    ! allocate velocity masks
    allocate(this%vel_mask_x(N), this%vel_mask_y(N), this%vel_mask_z(N))

    ! initialize CFL computations
    this%have_CFL_x = .false.
    this%have_CFL_y = .false.
    this%have_CFL_z = .false.

    ! initialize monitoring
    allocate(this%Umin(N), this%Umax(N), this%Uint(N))
    this%have_Urange = .false.

  end function

  !> solver print function
  subroutine rusanov_print(this)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    class(rusanov), intent(in) :: this

    if (this%cfg%amRoot) then
      write(output_unit,'("rusanov-type solver [",a,"] for config [",a,"]")')   &
        & trim(this%name), trim(this%cfg%name)
    end if

  end subroutine

  !> Add a boundary condition
  subroutine add_bcond(this, name, type, locator, dir)
    use iterator_class, only: locator_gen_ftype
    use string,         only: lowercase
    use messager,       only: die
    implicit none
    class(rusanov), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(1), intent(in) :: type
    procedure(locator_gen_ftype) :: locator
    character(len=2), intent(in) :: dir
    type(rus_bc), pointer :: new_bc

    ! prepare new bcond
    allocate(new_bc)
    new_bc%name = trim(adjustl(name))
    new_bc%type = type
    new_bc%itr = iterator(pg=this%cfg, name=new_bc%name, locator=locator)
    select case (lowercase(dir))
    case ('c');              new_bc%dir = 0_1
    case ('-x', 'x-', 'xl'); new_bc%dir = 2_1
    case ('+x', 'x+', 'xr'); new_bc%dir = 1_1
    case ('-y', 'y-', 'yl'); new_bc%dir = 4_1
    case ('+y', 'y+', 'yr'); new_bc%dir = 3_1
    case ('-z', 'z-', 'zl'); new_bc%dir = 6_1
    case ('+z', 'z+', 'zr'); new_bc%dir = 5_1
    case default; call die('[rusanov add_bcond] unknown bcond direction')
    end select

    ! insert it in front
    new_bc%next => this%first_bc
    this%first_bc => new_bc
    this%nbc = this%nbc + 1

  end subroutine add_bcond

  !> get a boundary condition
  subroutine get_bcond(this, name, my_bc)
    use messager, only: die
    implicit none
    class(rusanov), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(rus_bc), pointer, intent(out) :: my_bc

    my_bc => this%first_bc
    do while (associated(my_bc))
      if (trim(my_bc%name).eq.trim(name)) return
      my_bc => my_bc%next
    end do

    call die('[rusanov get_bcond] Boundary condition was not found')

  end subroutine get_bcond

  !> enforce boundary condition
  !> note that this interpolates (ensuring ghost cells have the right value),
  !> but all wave-related constraints are imposed during the step
  subroutine apply_bcond(this, t, dt)
    use messager, only: die
    implicit none
    class(rusanov), intent(inout) :: this
    real(WP), intent(in) :: t, dt
    logical, dimension(this%N) :: vel_mask
    integer :: i, j, k, m, n, iref, jref, kref
    type(rus_bc), pointer :: my_bc

    ! Traverse bcond list
    my_bc => this%first_bc
    do while (associated(my_bc))

      ! Only processes inside the bcond work here
      if (my_bc%itr%amIn) then

        ! Select appropriate action based on the bcond type
        select case (my_bc%type)
          case (0_1)                        !< do nothing
          case (1_1)                        !< extrapolate
            do m = 1, my_bc%itr%no_
              i = my_bc%itr%map(1,m)
              j = my_bc%itr%map(2,m)
              k = my_bc%itr%map(3,m)
              iref = min(this%cfg%imax, max(this%cfg%imin, i))
              jref = min(this%cfg%jmax, max(this%cfg%jmin, j))
              kref = min(this%cfg%kmax, max(this%cfg%kmin, k))
              this%Uc(:,i,j,k) = this%Uc(:,iref,jref,kref)
            end do
          case (2_1)                        !< reflect
            do m = 1, my_bc%itr%no_
              i = my_bc%itr%map(1,m)
              j = my_bc%itr%map(2,m)
              k = my_bc%itr%map(3,m)
              iref = min(this%cfg%imax, max(this%cfg%imin, i))
              jref = min(this%cfg%jmax, max(this%cfg%jmin, j))
              kref = min(this%cfg%kmax, max(this%cfg%kmin, k))
              if (i .ne. iref) then
                if (i .lt. iref) then
                  iref = i + 2 * (iref - i) - 1
                else
                  iref = i - 2 * (i - iref) + 1
                end if
                vel_mask(:) = this%vel_mask_x(:)
              end if
              if (j .ne. jref) then
                if (j .lt. jref) then
                  jref = j + 2 * (jref - j) - 1
                else
                  jref = j - 2 * (j - jref) + 1
                end if
                vel_mask(:) = this%vel_mask_y(:)
              end if
              if (k .ne. kref) then
                if (k .lt. kref) then
                  kref = k + 2 * (kref - k) - 1
                else
                  kref = k - 2 * (k - kref) + 1
                end if
                vel_mask(:) = this%vel_mask_z(:)
              end if
              do n = 1, this%N
                if (vel_mask(n)) then
                  this%Uc(n,i,j,k) = - this%Uc(n,iref,jref,kref)
                else
                  this%Uc(n,i,j,k) = + this%Uc(n,iref,jref,kref)
                end if
              end do
            end do
          case default
            call die('[rusanov apply_bcond] Unknown bcond type')
        end select

        ! zero normal velocity?
        !masked_type = iand(my_bc%type, 4_1)
        !if (masked_type .eq. 4_1) then
        !  select case (my_bc%dir)
        !  case (0_1)
        !    vel_mask(:) = this%vel_mask_x(:) .and. this%vel_mask_y(:) .and.   &
        !      & this%vel_mask_z(:)
        !  case (1_1, 2_1)
        !    vel_mask(:) = this%vel_mask_x(:)
        !  case (3_1, 4_1)
        !    vel_mask(:) = this%vel_mask_y(:)
        !  case (5_1, 6_1)
        !    vel_mask(:) = this%vel_mask_z(:)
        !  end select
        !  do m = 1, my_bc%itr%n_
        !    i = my_bc%itr%map(1,m)
        !    j = my_bc%itr%map(2,m)
        !    k = my_bc%itr%map(3,m)
        !    do n = 1, this%N
        !      if (vel_mask(n)) this%Uc(n,i,j,k) = 0.0_WP
        !    end do
        !  end do
        !end if

      end if

      ! Sync full fields after each bcond - this should be optimized
      call this%cfg%sync(this%Uc)
      call this%cfg%sync(this%dU)

      ! Move on to the next bcond
      my_bc => my_bc%next

    end do

  end subroutine apply_bcond

  !> get CFL number
  subroutine get_cfl(this, dt, cfl)
    use messager, only: die
    use mpi_f08,  only: MPI_ALLREDUCE, MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(rusanov), intent(inout) :: this
    real(WP), intent(in) :: dt
    real(WP), intent(out) :: cfl
    real(WP), dimension(this%N) :: evals
    real(WP) :: task_CFL_x, task_CFL_y, task_CFL_z
    logical :: have_CFL
    integer :: i, j, k, ierrx, ierry, ierrz

    have_CFL = this%have_CFL_x .and. this%have_CFL_y .and. this%have_CFL_z

    if (.not. have_CFL) then

      call this%cfg%sync(this%Uc)

      task_CFL_x = 0.0_WP; task_CFL_y = 0.0_WP; task_CFL_z = 0.0_WP;
      do k = this%cfg%kmin_, this%cfg%kmax_
        do j = this%cfg%jmin_, this%cfg%jmax_
          do i = this%cfg%imin_, this%cfg%imax_
            call this%evals_x(this%P, this%N, this%params(:,i,j,k), this%Uc(:,i,j,k), evals)
            task_CFL_x = max(task_CFL_x, this%cfg%dxi(i) * maxval(abs(evals)))
            call this%evals_y(this%P, this%N, this%params(:,i,j,k), this%Uc(:,i,j,k), evals)
            task_CFL_y = max(task_CFL_y, this%cfg%dyi(i) * maxval(abs(evals)))
            call this%evals_z(this%P, this%N, this%params(:,i,j,k), this%Uc(:,i,j,k), evals)
            task_CFL_z = max(task_CFL_z, this%cfg%dzi(i) * maxval(abs(evals)))
          end do
        end do
      end do

      call MPI_ALLREDUCE(task_CFL_x, this%CFL_x, 1, MPI_REAL_WP, MPI_MAX,       &
        & this%cfg%comm, ierrx)
      call MPI_ALLREDUCE(task_CFL_y, this%CFL_y, 1, MPI_REAL_WP, MPI_MAX,       &
        & this%cfg%comm, ierry)
      call MPI_ALLREDUCE(task_CFL_z, this%CFL_z, 1, MPI_REAL_WP, MPI_MAX,       &
        & this%cfg%comm, ierrz)

      if (ierrx .ne. 0) call die("error reducing CFL")
      if (ierry .ne. 0) call die("error reducing CFL")
      if (ierrz .ne. 0) call die("error reducing CFL")

      this%have_CFL_x = .true.
      this%have_CFL_y = .true.
      this%have_CFL_z = .true.

    end if

    cfl = dt * max(this%CFL_x, this%CFL_y, this%CFL_z)

  end subroutine get_cfl

  !> get min/max field values
  subroutine get_range(this)
    use messager, only: die
    use mpi_f08,  only: MPI_ALLREDUCE, MPI_MIN, MPI_MAX, MPI_SUM
    use parallel, only: MPI_REAL_WP
    implicit none
    class(rusanov), intent(inout) :: this
    real(WP), dimension(this%N) :: task_Umin, task_Umax, task_Uint
    integer :: i, j, k, ierr

    if (.not. this%have_Urange) then

      task_Umin(:) = + huge(1.0_WP)
      task_Umax(:) = - huge(1.0_WP)
      task_Uint(:) = 0.0_WP

      do k=this%cfg%kmin_,this%cfg%kmax_
        do j=this%cfg%jmin_,this%cfg%jmax_
          do i=this%cfg%imin_,this%cfg%imax_
            task_Umin = min(task_Umin, this%Uc(:,i,j,k))
            task_Umax = max(task_Umax, this%Uc(:,i,j,k))
            task_Uint = task_Uint + this%Uc(:,i,j,k) * this%cfg%dx(i) *       &
              & this%cfg%dy(j) * this%cfg%dz(k)
          end do
        end do
      end do

      call MPI_ALLREDUCE(task_Umin, this%Umin, this%N, MPI_REAL_WP, MPI_MIN,  &
        & this%cfg%comm, ierr)
      if (ierr .ne. 0) call die("error reducing Umin")
      call MPI_ALLREDUCE(task_Umax, this%Umax, this%N, MPI_REAL_WP, MPI_MAX,  &
        & this%cfg%comm, ierr)
      if (ierr .ne. 0) call die("error reducing Umax")
      call MPI_ALLREDUCE(task_Uint, this%Uint, this%N, MPI_REAL_WP, MPI_SUM,  &
        & this%cfg%comm, ierr)
      if (ierr .ne. 0) call die("error reducing Uint")

      this%have_Urange = .true.

    end if

  end subroutine get_range

  !> compute dU in x direction
  subroutine calc_dU_x(this, dt)
    implicit none
    class(rusanov), intent(inout) :: this
    real(WP), intent(in) :: dt
    integer :: j, k

    call this%cfg%sync(this%Uc)

    do k = this%cfg%kmin_, this%cfg%kmax_
      do j = this%cfg%jmin_, this%cfg%jmax_
        call rusanov_1d(this%N, this%P, size(this%Uc, 2), this%flux_x,        &
          this%evals_x, this%cfg%dx, dt, this%params(:,:,j,k),                &
          this%Uc(:,:,j,k), this%dU(:,:,j,k))
      end do
    end do

    call this%cfg%sync(this%dU)

    this%have_CFL_x = .false.; this%have_Urange = .false.;

  end subroutine calc_dU_x

  !> compute dU in y direction
  subroutine calc_dU_y(this, dt)
    implicit none
    class(rusanov), intent(inout) :: this
    real(WP), intent(in) :: dt
    integer :: i, k, jmino_, jmaxo_
    real(WP), dimension(:,:), allocatable :: U1d, dU1d, p1d

    if (this%cfg%ny .lt. 2) return

    call this%cfg%sync(this%Uc)

    jmino_ = this%cfg%jmino_; jmaxo_ = this%cfg%jmaxo_;

    allocate(U1d(this%N,jmino_:jmaxo_), p1d(this%P,jmino_:jmaxo_), dU1d(this%N,jmino_:jmaxo_))

    do k = this%cfg%kmin_, this%cfg%kmax_
      do i = this%cfg%imin_, this%cfg%imax_
        U1d(:,:) = this%Uc(:,i,jmino_:jmaxo_,k)
        dU1d(:,:) = this%dU(:,i,jmino_:jmaxo_,k)
        p1d(:,:) = this%params(:,i,jmino_:jmaxo_,k)
        call rusanov_1d(this%N, this%P, size(U1d, 2), this%flux_y,            &
          this%evals_y, this%cfg%dy, dt, p1d, U1d, dU1d)
        this%dU(:,i,jmino_:jmaxo_,k) = dU1d
      end do
    end do

    deallocate(U1d, dU1d, p1d)

    call this%cfg%sync(this%dU)

    this%have_CFL_y = .false.; this%have_Urange = .false.;

  end subroutine calc_dU_y

  !> compute dU in z direction
  subroutine calc_dU_z(this, dt)
    implicit none
    class(rusanov), intent(inout) :: this
    real(WP), intent(in) :: dt
    integer :: i, j, kmino_, kmaxo_
    real(WP), dimension(:,:), allocatable :: U1d, dU1d, p1d

    if (this%cfg%nz .lt. 2) return

    call this%cfg%sync(this%Uc)

    kmino_ = this%cfg%kmino_; kmaxo_ = this%cfg%kmaxo_;

    allocate(U1d(this%N,kmino_:kmaxo_), p1d(this%P,kmino_:kmaxo_), dU1d(this%N,kmino_:kmaxo_))

    do j = this%cfg%jmin_, this%cfg%jmax_
      do i = this%cfg%imin_, this%cfg%imax_
        U1d(:,:) = this%Uc(:,i,j,kmino_:kmaxo_)
        dU1d(:,:) = this%dU(:,i,j,kmino_:kmaxo_)
        p1d(:,:) = this%params(:,i,j,kmino_:kmaxo_)
        call rusanov_1d(this%N, this%P, size(U1d, 2), this%flux_z,            &
          this%evals_z, this%cfg%dz, dt, p1d, U1d, dU1d)
        this%dU(:,i,j,kmino_:kmaxo_) = dU1d
      end do
    end do

    deallocate(U1d, dU1d, p1d)

    call this%cfg%sync(this%dU)

    this%have_CFL_z = .false.; this%have_Urange = .false.;

  end subroutine calc_dU_z

  ! requires one ghost cell
  ! left as an independent function, but called is by class
  pure subroutine rusanov_1d(N, P, M, flux, evals, dxs, dt, params, U, dU)
    implicit none
    integer, intent(in) :: N, P, M
    procedure(eigenvals_ftype), pointer, intent(in) :: evals
    procedure(flux_ftype), pointer, intent(in) :: flux
    real(WP), dimension(M), intent(in) :: dxs                   ! size M
    real(WP), intent(in) :: dt
    real(WP), dimension(P,M), intent(in) :: params              ! size (P, M)
    real(WP), dimension(N,M), intent(in) :: U                   ! size (N, M)
    real(WP), dimension(N,M), intent(inout) :: dU               ! size (N, M)
    real(WP), dimension(N) :: F, Fl, Fr, nudiff, laml, lamr
    integer :: j

    dU(:,:) = 0.0_WP

    call flux(P, N, params(:,1), U(:,1), Fl)
    call evals(P, N, params(:,1), U(:,1), laml)
    laml(:) = abs(laml(:))

    do j = 2, M

      call flux(P, N, params(:,j), U(:,j), Fr)
      call evals(P, N, params(:,j), U(:,j), lamr)
      lamr(:) = abs(lamr(:))

      nudiff = (U(:,j) - U(:,j-1)) * max(maxval(laml), maxval(lamr))

      ! not the right division by dx, but whatever
      F = dt * 0.5_WP / (0.5_WP * (dxs(j-1) + dxs(j))) * (Fr + Fl - nudiff)

      dU(:,j-1) = dU(:,j-1) - F
      dU(:,j  ) = dU(:,j  ) + F

      Fl(:) = Fr(:); laml(:) = lamr(:);

    end do

  end subroutine rusanov_1d

end module rusanov_class

