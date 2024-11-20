
!> contains necessary functions for using the muscl class (or other general
!> hyperbolic solvers) to solve the Euler equations for an ideal (gamma) gas

!> Originally written by John P Wakefield in December 2022.

module hyperbolic_euler
  use precision,     only: WP
  use string,        only: str_medium
  use config_class,  only: config
  use hyperbolic,    only: VANLEER, eigenvals_ftype, rsolver_ftype, limiter_ftype, flux_ftype
  use muscl_class,   only: muscl
  use rusanov_class, only: rusanov
  implicit none

  real(WP), parameter :: euler_muscl_cflsafety = 0.92_WP
  real(WP), parameter :: euler_muscl_divzero_eps = 1e-9_WP
  character(len=str_medium), parameter :: euler_muscl_name = 'MUSCL_EULER'

  real(WP), parameter :: DIATOMIC_GAMMA = 1.4_WP
  real(WP), dimension(5), parameter :: SOD_PHYS_L = (/ 1.0_WP, 0.0_WP,        &
    & 0.0_WP, 0.0_WP, 1.0_WP /)
  real(WP), dimension(5), parameter :: SOD_PHYS_R = (/ 0.125_WP, 0.0_WP,      &
    & 0.0_WP, 0.0_WP, 0.1_WP /)

contains

  !> muscl factory
  function make_euler_muscl(cfg, limiter, gma) result(solver)
    implicit none
    type(muscl) :: solver
    class(config), target, intent(in) :: cfg
    integer(1), optional, intent(in) :: limiter
    real(WP), optional, intent(in) :: gma
    integer(1) :: limiter_actual
    real(WP) :: gma_actual
    character(len=str_medium) :: name_actual
    procedure(eigenvals_ftype), pointer :: evals_x_ptr, evals_y_ptr, evals_z_ptr
    procedure(rsolver_ftype), pointer :: rsolv_x_ptr, rsolv_y_ptr, rsolv_z_ptr

    if (present(limiter)) then
      limiter_actual = limiter
    else
      limiter_actual = VANLEER
    end if

    if (present(gma)) then
      gma_actual = gma
    else
      gma_actual = DIATOMIC_GAMMA
    end if

    name_actual = euler_muscl_name

    evals_x_ptr => euler_evals_x; rsolv_x_ptr => euler_rsolv_x;
    evals_y_ptr => euler_evals_y; rsolv_y_ptr => euler_rsolv_y;
    evals_z_ptr => euler_evals_z; rsolv_z_ptr => euler_rsolv_z;

    ! build solver
    solver = muscl(cfg, name_actual, 5, 1, evals_x_ptr, evals_y_ptr,          &
      & evals_z_ptr, rsolv_x_ptr, rsolv_y_ptr, rsolv_z_ptr, limiter_actual,   &
      & euler_muscl_divzero_eps, euler_muscl_cflsafety)

    ! set param array to hold gamma
    solver%params(1,:,:,:) = gma_actual

    ! set velocity mask
    solver%vel_mask_x(:) = (/ .false., .true., .false., .false., .false. /)
    solver%vel_mask_y(:) = (/ .false., .false., .true., .false., .false. /)
    solver%vel_mask_z(:) = (/ .false., .false., .false., .true., .false. /)

  end function make_euler_muscl

  !> rusanov factory
  function make_euler_rusanov(cfg, gma) result(solver)
    implicit none
    type(rusanov) :: solver
    class(config), target, intent(in) :: cfg
    real(WP), optional, intent(in) :: gma
    real(WP) :: gma_actual
    procedure(eigenvals_ftype), pointer :: evals_x_ptr, evals_y_ptr, evals_z_ptr
    procedure(flux_ftype), pointer :: flux_x_ptr, flux_y_ptr, flux_z_ptr

    if (present(gma)) then
      gma_actual = gma
    else
      gma_actual = DIATOMIC_GAMMA
    end if

    evals_x_ptr => euler_evals_x; flux_x_ptr => euler_flux_x;
    evals_y_ptr => euler_evals_y; flux_y_ptr => euler_flux_y;
    evals_z_ptr => euler_evals_z; flux_z_ptr => euler_flux_z;

    ! build solver
    solver = rusanov(cfg, 'EULER_RUSANOV', 5, 1, evals_x_ptr, evals_y_ptr,    &
      & evals_z_ptr, flux_x_ptr, flux_y_ptr, flux_z_ptr)

    ! set param array to hold gamma
    solver%params(1,:,:,:) = gma_actual

    ! set velocity mask
    solver%vel_mask_x(:) = (/ .false., .true., .false., .false., .false. /)
    solver%vel_mask_y(:) = (/ .false., .false., .true., .false., .false. /)
    solver%vel_mask_z(:) = (/ .false., .false., .false., .true., .false. /)

  end function make_euler_rusanov

  !> Convert to Physical Coordinates (velocity and pressure to momentum and energy)
  pure subroutine euler_tophys(gma, cons, phys)
    implicit none
    real(WP), intent(in) :: gma
    real(WP), dimension(5), intent(in) :: cons
    real(WP), dimension(5), intent(out) :: phys
    real(WP) :: KE

    phys(1) = cons(1)
    phys(2) = cons(2) / cons(1)
    phys(3) = cons(3) / cons(1)
    phys(4) = cons(4) / cons(1)
    KE = 0.5_WP * sum(cons(2:4)**2) / cons(1)
    phys(5) = (gma - 1.0_WP) * (cons(5) - KE)

  end subroutine

  !> Convert to Conserved Coordinates (velocity and pressure to momentum and energy)
  pure subroutine euler_tocons(gma, phys, cons)
    implicit none
    real(WP), intent(in) :: gma
    real(WP), dimension(5), intent(in) :: phys
    real(WP), dimension(5), intent(out) :: cons
    real(WP) :: KE

    cons(1) = phys(1)
    cons(2) = phys(1) * phys(2)
    cons(3) = phys(1) * phys(3)
    cons(4) = phys(1) * phys(4)
    KE = 0.5_WP * phys(1) * sum(phys(2:4)**2)
    cons(5) = phys(5) / (gma - 1.0_WP) + KE

  end subroutine

  pure subroutine euler_flux_x(P, N, params, U, flux)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: flux
    real(WP) :: pressure

    pressure = (params(1) - 1.0_WP) * (U(5) - 0.5_WP * (U(2)**2 + U(3)**2 +   &
      U(4)**2)/ U(1))

    flux(1) = U(2)
    flux(2) = U(2)**2 / U(1) + pressure
    flux(3) = U(2) * U(3) / U(1)
    flux(4) = U(2) * U(4) / U(1)
    flux(5) = U(2) / U(1) * (U(5) + pressure)

  end subroutine euler_flux_x

  pure subroutine euler_flux_y(P, N, params, U, flux)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: flux
    real(WP), dimension(5) :: u_permute

    u_permute(:) = u(:); u_permute(2) = u(3); u_permute(3) = u(2);

    call euler_flux_x(P, N, params, u_permute, flux)

  end subroutine euler_flux_y

  pure subroutine euler_flux_z(P, N, params, U, flux)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: flux
    real(WP), dimension(N) :: u_permute

    u_permute(:) = u(:); u_permute(2) = u(4); u_permute(4) = u(2);

    call euler_flux_x(P, N, params, u_permute, flux)

  end subroutine euler_flux_z

  pure subroutine euler_evals_x(P, N, params, U, evals)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: evals

    call euler_evals_1d(params(1), U, evals)

  end subroutine euler_evals_x

  pure subroutine euler_evals_y(P, N, params, U, evals)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: evals
    real(WP), dimension(5) :: u_permute

    u_permute(:) = u(:); u_permute(2) = u(3); u_permute(3) = u(2);

    call euler_evals_1d(params(1), u_permute, evals)

  end subroutine euler_evals_y

  pure subroutine euler_evals_z(P, N, params, U, evals)
    implicit none
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: params
    real(WP), dimension(N), intent(in) :: u
    real(WP), dimension(N), intent(out) :: evals
    real(WP), dimension(N) :: u_permute

    u_permute(:) = u(:); u_permute(2) = u(4); u_permute(4) = u(2);

    call euler_evals_1d(params(1), u_permute, evals)

  end subroutine euler_evals_z

  !> eigenvalues
  pure subroutine euler_evals_1d(gma, U, lambda)
    implicit none
    real(WP), intent(in) :: gma
    real(WP), dimension(5), intent(in) :: U
    real(WP), dimension(5), intent(out) :: lambda
    real(WP) :: v, p, c, KE

    KE = 0.5_WP * sum(U(2:4)**2) / U(1)
    v = U(2) / U(1)
    p = (gma - 1) * (U(5) - KE)
    c = sqrt(gma * p / U(1))

    lambda = (/ v - c, v, v, v, v + c /)

  end subroutine euler_evals_1d

  pure subroutine euler_rsolv_x(P, N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs

    call euler_rsolv_roe_1d(0.5_WP * (pl(1) + pr(1)), Ul, Ur, rs)

  end subroutine euler_rsolv_x

  pure subroutine euler_rsolv_y(P, N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs
    real(WP), dimension(N) :: Uln, Urn
    real(WP), dimension(5) :: rs_row

    Uln(:) = Ul(:); Uln(2) = Ul(3); Uln(3) = Ul(2);
    Urn(:) = Ur(:); Urn(2) = Ur(3); Urn(3) = Ur(2);

    call euler_rsolv_roe_1d(0.5_WP * (pl(1) + pr(1)), Uln, Urn, rs)

    rs_row(:) = rs(2,5:9); rs(2,5:9) = rs(3,5:9); rs(3,5:9) = rs_row(:);

  end subroutine euler_rsolv_y

  pure subroutine euler_rsolv_z(P, N, pl, Ul, pr, Ur, rs)
    integer, intent(in) :: P, N
    real(WP), dimension(P), intent(in) :: pl, pr
    real(WP), dimension(N), intent(in) :: Ul, Ur
    real(WP), dimension(:,:), intent(out) :: rs
    real(WP), dimension(N) :: Uln, Urn
    real(WP), dimension(5) :: rs_row

    Uln(:) = Ul(:); Uln(2) = Ul(4); Uln(4) = Ul(2);
    Urn(:) = Ur(:); Urn(2) = Ur(4); Urn(4) = Ur(2);

    call euler_rsolv_roe_1d(0.5_WP * (pl(1) + pr(1)), Uln, Urn, rs)

    rs_row(:) = rs(2,5:9); rs(2,5:9) = rs(4,5:9); rs(4,5:9) = rs_row(:);

  end subroutine euler_rsolv_z

  !> Roe Solver
  ! see page 354-ish of Toro's book
  ! here `Ul' and `Ur' are the whole state, (u, v, w) is the velocity vector
  pure subroutine euler_rsolv_roe_1d(gma, Ul, Ur, rs)
    implicit none
    real(WP), intent(in) :: gma
    real(WP), dimension(5), intent(in) :: Ul, Ur
    real(WP), dimension(5,9), intent(out) :: rs
    real(WP) :: rrhol, rrhor, drho, rho, pl, pr, Hl, Hr, u, v, w, H, a, du5b

    ! Roe averages
    rrhol = sqrt(Ul(1)); rrhor = sqrt(Ur(1));
    rho = rrhol * rrhor
    pl = (gma - 1.0_WP) * (Ul(5) - 0.5_WP * sum(Ul(2:4)**2) / Ul(1))
    pr = (gma - 1.0_WP) * (Ur(5) - 0.5_WP * sum(Ur(2:4)**2) / Ur(1))
    Hl = (Ul(5) + pl) / Ul(1); Hr = (Ur(5) + pr) / Ur(1);
    u = (rrhol * Ul(2) / Ul(1) + rrhor * Ur(2) / Ur(1)) / (rrhol + rrhor)
    v = (rrhol * Ul(3) / Ul(1) + rrhor * Ur(3) / Ur(1)) / (rrhol + rrhor)
    w = (rrhol * Ul(4) / Ul(1) + rrhor * Ur(4) / Ur(1)) / (rrhol + rrhor)
    H = (rrhol * Hl + rrhor * Hr) / (rrhol + rrhor)
    a = sqrt((gma - 1.0_WP) * (H - 0.5_WP * (u**2 + v**2 + w**2)))

    ! lambdal and lambdar
    rs(1,2) = u - a; rs(2:4,2) = u; rs(5,2) = u + a;
    rs(:,1) = min(rs(:,2), 0.0_WP)
    rs(:,2) = rs(:,2) - rs(:,1)

    ! α2
    drho = Ur(1) - Ul(1)
    du5b = Ur(5) - Ul(5) - (Ur(3) - Ul(3) - v * drho) * v - (Ur(4) - Ul(4) -  &
      & w * drho) * w
    rs(2,3) = (gma - 1.0_WP) / a**2 * (                                       &
      & (Ur(1) - Ul(1)) * (H - u**2) + u * (Ur(2) - Ul(2)) - du5b             &
      & )

    ! α1
    rs(1,3) = (drho * (u + a) - (Ur(2) - Ul(2)) - a * rs(2,3)) / (2 * a)

    ! α5
    rs(5,3) = drho - rs(1,3) - rs(2,3)

    ! α3 and α4
    rs(3,3) = Ur(3) - Ul(3) - v * drho
    rs(4,3) = Ur(4) - Ul(4) - w * drho

    ! βs
    rs(:,4) = 0.0_WP

    ! eigenvectors
    rs(1,5) = 1.0_WP; rs(2,5) = u - a; rs(3,5) = v; rs(4,5) = w; rs(5,5) = H - u * a;
    rs(1,9) = 1.0_WP; rs(2,9) = u + a; rs(3,9) = v; rs(4,9) = w; rs(5,9) = H + u * a;
    rs(1,6) = 1.0_WP; rs(2,6) = u; rs(3,6) = v; rs(4,6) = w; rs(5,6) = 0.5_WP * (u**2 + v**2 + w**2);
    rs(1,7) = 0.0_WP; rs(2,7) = 0.0_WP; rs(3,7) = 1.0_WP; rs(4,7) = 0.0_WP; rs(5,7) = v;
    rs(1,8) = 0.0_WP; rs(2,8) = 0.0_WP; rs(3,8) = 0.0_WP; rs(4,8) = 1.0_WP; rs(5,8) = w;

  end subroutine euler_rsolv_roe_1d

  !! subroutines for Solving Riemann Problems Exactly

  !> From Toro, returns fK and fKprime
  pure subroutine euler_fK(gma, rhoK, pK, cK, p, fK, fKprime)
    implicit none
    real(WP), intent(in) :: gma
    real(WP), intent(in) :: rhoK, pK, cK, p
    real(WP), intent(out) :: fK, fKprime
    real(WP) :: AK, BK

    AK = 2 * ((gma + 1) * rhoK)**(-1); BK = (gma - 1) / (gma + 1) * pK;
    if (p > pK) then                ! shock
      fK  = (p - pK) * sqrt(AK / (p + BK))
      fKprime = sqrt(AK / (BK + p)) * (1 - (p - pK) / (2 * (BK + p)))
    else                            ! rarefaction
      fK = 2 * cK / (gma - 1) * ((p / pK)**((gma - 1) / (2 * gma)) - 1)
      fKprime = (rhoK * cK)**(-1) * (p / pK)**(- (gma + 1) / (2 * gma))
    end if

  end subroutine

  !> returns f and fprime
  pure subroutine euler_f(gma, rhoL, vL, pL, cL, rhoR, vR, pR, cR, p, f, fprime)
    implicit none
    real(WP), intent(in) :: gma
    real(WP), intent(in) :: rhoL, vL, pL, cL, rhoR, vR, pR, cR, p
    real(WP), intent(out) :: f, fprime
    real(WP) :: fL, fLp, fR, fRp

    call euler_fK(gma, rhoL, pL, cL, p, fL, fLp)
    call euler_fK(gma, rhoR, pR, cR, p, fR, fRp)

    f = fL + fR + vR - vL
    fprime = fLp + fRp

  end subroutine

  !> returns pressure and velocity (in that order) of the mid state
  pure subroutine euler_findmid(gma, rhoL, vL, pL, rhoR, vR, pR, delta, maxiters, pM, vM)
    implicit none
    real(WP), intent(in) :: gma
    real(WP), intent(in) :: rhoL, vL, pL, rhoR, vR, pR, delta
    integer, intent(in) :: maxiters
    real(WP), intent(out) :: pM, vM
    real(WP) :: z, f, fp, cL, cR, fL, fLp, fR, fRp
    integer(WP) :: n

    cL = sqrt(gma * pL / rhoL); cR = sqrt(gma * pR / rhoR);

    ! start with two rarefaction solution
    z = (gma - 1) / (2 * gma)
    pM = ((cL + cR - 0.5_WP * (gma - 1) * (vR - vL)) / (cL / pL**z + cR /       &
      & pR**z))**(z**(-1))

    ! do rootfinding
    call euler_f(gma, rhoL, vL, pL, cL, rhoR, vR, pR, cR, pM, f, fp)
    n = 0
    do while (abs(f) < delta .and. n < maxiters)
      pM = pM - f / fp
      call euler_f(gma, rhoL, vL, pL, cL, rhoR, vR, pR, cR, pM, f, fp)
      n = n + 1
    end do

    call euler_fK(gma, rhoL, pL, cL, pM, fL, fLp)
    call euler_fK(gma, rhoR, pR, cR, pM, fR, fRp)
    vM = 0.5_WP * (vL + vR + fR - fL)

  end subroutine

end module hyperbolic_euler

