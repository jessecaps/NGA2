
!> contains necessary functions for using the muscl class (or other general hyperbolic
!> solvers) to solve basic advection problems

!> Originally written by John P Wakefield in December 2022.

module hyperbolic_advection
  use precision,     only: WP
  use string,        only: str_medium
  use config_class,  only: config
  use hyperbolic,    only: eigenvals_ftype, rsolver_ftype, limiter_ftype, flux_ftype
  use muscl_class,   only: muscl
  use rusanov_class, only: rusanov
  implicit none

  real(WP), parameter :: advec_muscl_cflsafety = 0.99_WP
  real(WP), parameter :: advec_muscl_divzero_eps = 1e-12_WP
  character(len=str_medium), parameter :: ADVEC_MUSCL_NAME = 'MUSCL_CONST_ADVEC'
  character(len=str_medium), parameter :: ADVEC_RUS_NAME = 'RUSANOV_CONST_ADVEC'

contains

  !> muscl factory
  function make_advec_muscl(cfg, N, limiter, velocity) result(solver)
    implicit none
    type(muscl) :: solver
    class(config), target, intent(in) :: cfg
    integer, intent(in) :: N
    integer(1), intent(in) :: limiter
    real(WP), dimension(3), intent(in) :: velocity
    character(len=str_medium) :: name_actual
    procedure(eigenvals_ftype), pointer :: evals_x_ptr, evals_y_ptr, evals_z_ptr
    procedure(rsolver_ftype), pointer :: rsolv_x_ptr, rsolv_y_ptr, rsolv_z_ptr
    integer :: i

    name_actual = ADVEC_MUSCL_NAME

    evals_x_ptr => advec_evals_x; rsolv_x_ptr => advec_rsolv_x;
    evals_y_ptr => advec_evals_y; rsolv_y_ptr => advec_rsolv_y;
    evals_z_ptr => advec_evals_z; rsolv_z_ptr => advec_rsolv_z;

    ! build solver
    solver = muscl(cfg, name_actual, N, 3, evals_x_ptr, evals_y_ptr,          &
      & evals_z_ptr, rsolv_x_ptr, rsolv_y_ptr, rsolv_z_ptr, limiter,          &
      & advec_muscl_divzero_eps, advec_muscl_cflsafety)

    ! set velocity
    do i = 1, 3
      solver%params(i,:,:,:) = velocity(i)
    end do

    ! set velocity mask
    solver%vel_mask_x(:) = .false.
    solver%vel_mask_y(:) = .false.
    solver%vel_mask_z(:) = .false.

  end function make_advec_muscl

  !> rusanov factory
  function make_advec_rusanov(cfg, N, velocity) result(solver)
    implicit none
    type(rusanov) :: solver
    class(config), target, intent(in) :: cfg
    integer, intent(in) :: N
    real(WP), dimension(3), intent(in) :: velocity
    procedure(eigenvals_ftype), pointer :: evals_x_ptr, evals_y_ptr, evals_z_ptr
    procedure(flux_ftype), pointer :: flux_x_ptr, flux_y_ptr, flux_z_ptr
    integer :: i

    evals_x_ptr => advec_evals_x; flux_x_ptr => advec_flux_x;
    evals_y_ptr => advec_evals_y; flux_y_ptr => advec_flux_y;
    evals_z_ptr => advec_evals_z; flux_z_ptr => advec_flux_z;

    ! build solver
    solver = rusanov(cfg, ADVEC_RUS_NAME, N, 3, evals_x_ptr, evals_y_ptr,     &
      & evals_z_ptr, flux_x_ptr, flux_y_ptr, flux_z_ptr)

    ! put velocity in param array
    do i = 1, 3
      solver%params(i,:,:,:) = velocity(i)
    end do

    ! set velocity mask
    solver%vel_mask_x(:) = .false.
    solver%vel_mask_y(:) = .false.
    solver%vel_mask_z(:) = .false.

  end function make_advec_rusanov

  pure subroutine advec_evals_x(P, N, params, U, evals)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: U
    real(WP), dimension(N), intent(out) :: evals

    ! silence unused variable warning
    evals(1) = U(1)

    evals(:) = params(1)

  end subroutine advec_evals_x

  pure subroutine advec_evals_y(P, N, params, U, evals)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: U
    real(WP), dimension(N), intent(out) :: evals

    ! silence unused variable warning
    evals(1) = U(1)

    evals(:) = params(2)

  end subroutine advec_evals_y

  pure subroutine advec_evals_z(P, N, params, U, evals)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: U
    real(WP), dimension(N), intent(out) :: evals

    ! silence unused variable warning
    evals(1) = U(1)

    evals(:) = params(3)

  end subroutine advec_evals_z

  pure subroutine advec_flux_x(P, N, params, U, flux)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: flux

    flux(:) = params(1) * U(:)

  end subroutine advec_flux_x

  pure subroutine advec_flux_y(P, N, params, U, flux)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: flux

    flux(:) = params(2) * U(:)

  end subroutine advec_flux_y

  pure subroutine advec_flux_z(P, N, params, U, flux)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: flux

    flux(:) = params(3) * U(:)

  end subroutine advec_flux_z

  pure subroutine advec_rsolv_simple(v, Ul, Ur, rs)
    real(WP), intent(in) :: v
    real(WP), dimension(:), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs
    integer :: i

    rs(:,:) = 0.0_WP
    rs(:,1) = min(v, 0.0_WP)
    rs(:,2) = max(v, 0.0_WP)
    rs(:,3) = Ur(:) - Ul(:)

    rs(1:size(Ul,1),5:(size(Ul,1)+4)) = 1.0_WP

  end subroutine

  pure subroutine advec_rsolv_x(P, N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs

    call advec_rsolv_simple(0.5_WP * (pl(1) + pr(1)), Ul, Ur, rs)

  end subroutine advec_rsolv_x

  pure subroutine advec_rsolv_y(P, N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs

    call advec_rsolv_simple(0.5_WP * (pl(2) + pr(2)), Ul, Ur, rs)

  end subroutine advec_rsolv_y

  pure subroutine advec_rsolv_z(P, N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs

    call advec_rsolv_simple(0.5_WP * (pl(3) + pr(3)), Ul, Ur, rs)

  end subroutine advec_rsolv_z

end module hyperbolic_advection

