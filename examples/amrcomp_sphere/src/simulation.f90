!> AMR compressible sphere test case with shock initialization
module simulation
   use precision,         only: WP
   use amrgrid_class,     only: amrgrid
   use amrcomp_class,     only: amrcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible solver
   type(timetracker) :: time
   type(amrcomp), target :: fs
   type(amrdata) :: dQdt,k1,k2,k3,k4
   type(amrdata) :: Umag,Mach

   !> IBs
   type(amrdata) :: VF
   
   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: Rec_tag=huge(1.0_WP)
   real(WP) :: Res_tag=huge(1.0_WP)
   
   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile
   
   !> Stiffened gas EOS parameters
   real(WP) :: Gamma,Pinf,Cv

   !> Flow parameters
   real(WP) :: M2,Xs                  !< Post-shock Mach and shock location
   real(WP) :: Ms                     !< Shock Mach number
   real(WP) :: rho1,p1,u1             !< Pre-shock state
   real(WP) :: rho2,p2,u2             !< Post-shock state
   real(WP) :: Reynolds,Prandtl       !< Viscous parameters

   !> Sponge parameters
   real(WP) :: R_spg=3.0_WP
   real(WP) :: L_spg=1.0_WP
   real(WP) :: nu_spg=0.01_WP
   
contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for sphere
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=sqrt(xyz(1)**2+xyz(2)**2+xyz(3)**2)-0.5_WP
   end function sphere_levelset

   !> P=EOS(RHO,I) - Stiffened gas
   pure real(WP) function get_P(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_P=RHO*I*(Gamma-1.0_WP)-Gamma*Pinf
   end function get_P
   
   !> T=f(RHO,P)
   pure real(WP) function get_T(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_T=(P+Pinf)/(Cv*RHO*(Gamma-1.0_WP))
   end function get_T
   
   !> C=f(RHO,P) - Speed of sound
   pure real(WP) function get_C(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_C=sqrt(Gamma*(P+Pinf)/RHO)
   end function get_C

   !> I=EOS(RHO,P)
   pure real(WP) function get_I(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_I=(P+Gamma*Pinf)/(RHO*(Gamma-1.0_WP))
   end function get_I

   !> Compute viscosity using Sutherland's law, zero bulk viscosity, and set diffusivity based on Prandtl number
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pT,pQ,pVisc,pBeta,pDiff
      real(WP) :: r_cyl,blend
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pT=>fs%T%mf(lvl)%dataptr(mfi)
            pQ=>fs%Q%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Sutherland's law
               pVisc(i,j,k,1)=Reynolds**(-1.0_WP)*(1.4042_WP*pT(i,j,k,1)**1.5_WP)/(pT(i,j,k,1)+0.4042_WP)
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Heat diffusivity: k = Cp*mu/Pr = Cv*Gamma*mu/Pr
               pDiff(i,j,k,1)=Gamma*Cv*pVisc(i,j,k,1)/Prandtl
               ! Apply sponge layer viscosity
               r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2+(amr%zlo+(real(k,WP)+0.5_WP)*amr%dz(lvl))**2)
               if (r_cyl.gt.R_spg) then
                  blend=min((r_cyl-R_spg)/L_spg,1.0_WP)**2
                  pVisc(i,j,k,1)=max(pVisc(i,j,k,1),blend*nu_spg*pQ(i,j,k,1))
                  pDiff(i,j,k,1)=max(pDiff(i,j,k,1),blend*nu_spg*pQ(i,j,k,1))
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities
   
   !> User init callback - set normal shock profile
   subroutine shock_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: rho,P,U,IE,H
      integer :: i
      call amrex_mfiter_build(mfi,ba,dm,tiling=.true.)
      do while (mfi%next())
         ! Get pointer to data
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         ! Get tilebox with overlap
         bx=mfi%growntilebox(solver%nover)
         do i=bx%lo(1),bx%hi(1)
            ! Evaluate Heaviside function
            H=Hshock(x=Xs-(solver%amr%xlo+(real(i,WP)+0.5_WP)*solver%amr%dx(lvl)),delta=0.5_WP*solver%amr%dx(lvl))
            ! Interpolate between post-shock (H=1, right of shock) and pre-shock (H=0, left of shock)
            rho=rho1+(rho2-rho1)*H
            U=u1+(u2-u1)*H
            P=p1+(p2-p1)*H
            IE=get_I(rho,P)
            ! Set conserved variables
            pQ(i,:,:,1)=rho
            pQ(i,:,:,2)=rho*U
            pQ(i,:,:,3)=0.0_WP
            pQ(i,:,:,4)=0.0_WP
            pQ(i,:,:,5)=rho*IE
         end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shock_init

   !> Apply inflow BC at low-x (face=1)
   subroutine shock_dirichlet(solver,lvl,time,face,bx,pQ)
      use amrex_amr_module, only: amrex_box
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      integer :: i,j,k
      select case (face)
       case (1)  ! X-LOW: Dirichlet inflow with pre-shock (stationary) values
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            pQ(i,j,k,1)=rho2
            pQ(i,j,k,2)=rho2*u2
            pQ(i,j,k,3)=0.0_WP
            pQ(i,j,k,4)=0.0_WP
            pQ(i,j,k,5)=rho2*get_I(rho2,p2)
         end do; end do; end do
      end select
   end subroutine shock_dirichlet

   !> Initialize fluid volume fraction
   subroutine init_VF(data,lvl,time,ba,dm)
      use mms_geom, only: initialize_volume_moments
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrdata), intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz
      integer :: i,j,k
      real(WP), parameter :: VFlo=1.0e-12_WP
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pVF=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            call initialize_volume_moments(lo=[data%amr%xlo+real(i  ,WP)*dx,data%amr%ylo+real(j  ,WP)*dy,data%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[data%amr%xlo+real(i+1,WP)*dx,data%amr%ylo+real(j+1,WP)*dy,data%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=3,VFlo=VFlo,VF=pVF(i,j,k,1),BL=BL,BG=BG)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_VF

   !> Tagger based on vorticity and divergence
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVisc
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
      real(WP) :: vort_mag,div_neg,rho,mu,Rec,Res,dist
      integer :: i,j,k
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         pVisc=>solver%visc%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get rho and mu
            rho=pQ(i,j,k,1)
            mu=pVisc(i,j,k,1)
            if (mu.le.0.0_WP) mu=1.0_WP/Reynolds
            ! Get velocity gradient
            dudx=0.5_WP*dxi*(pQ(i+1,j,k,2)/max(pQ(i+1,j,k,1),solver%rho_floor)-pQ(i-1,j,k,2)/max(pQ(i-1,j,k,1),solver%rho_floor))
            dudy=0.5_WP*dyi*(pQ(i,j+1,k,2)/max(pQ(i,j+1,k,1),solver%rho_floor)-pQ(i,j-1,k,2)/max(pQ(i,j-1,k,1),solver%rho_floor))
            dudz=0.5_WP*dzi*(pQ(i,j,k+1,2)/max(pQ(i,j,k+1,1),solver%rho_floor)-pQ(i,j,k-1,2)/max(pQ(i,j,k-1,1),solver%rho_floor))
            dvdx=0.5_WP*dxi*(pQ(i+1,j,k,3)/max(pQ(i+1,j,k,1),solver%rho_floor)-pQ(i-1,j,k,3)/max(pQ(i-1,j,k,1),solver%rho_floor))
            dvdy=0.5_WP*dyi*(pQ(i,j+1,k,3)/max(pQ(i,j+1,k,1),solver%rho_floor)-pQ(i,j-1,k,3)/max(pQ(i,j-1,k,1),solver%rho_floor))
            dvdz=0.5_WP*dzi*(pQ(i,j,k+1,3)/max(pQ(i,j,k+1,1),solver%rho_floor)-pQ(i,j,k-1,3)/max(pQ(i,j,k-1,1),solver%rho_floor))
            dwdx=0.5_WP*dxi*(pQ(i+1,j,k,4)/max(pQ(i+1,j,k,1),solver%rho_floor)-pQ(i-1,j,k,4)/max(pQ(i-1,j,k,1),solver%rho_floor))
            dwdy=0.5_WP*dyi*(pQ(i,j+1,k,4)/max(pQ(i,j+1,k,1),solver%rho_floor)-pQ(i,j-1,k,4)/max(pQ(i,j-1,k,1),solver%rho_floor))
            dwdz=0.5_WP*dzi*(pQ(i,j,k+1,4)/max(pQ(i,j,k+1,1),solver%rho_floor)-pQ(i,j,k-1,4)/max(pQ(i,j,k-1,1),solver%rho_floor))
            ! Get vorticity magnitude
            vort_mag=sqrt((dwdy-dvdz)**2+(dudz-dwdx)**2+(dvdx-dudy)**2)
            ! Get dilatation
            div_neg=min(dudx+dvdy+dwdz,0.0_WP)
            ! Tag based on cell Reynolds numbers
            Rec=rho*vort_mag*solver%amr%min_meshsize(lvl)**2/mu
            if (Rec.gt.Rec_tag) tagarr(i,j,k,1)=SETtag
            ! Also tag based on cell shock Reynolds number
            Res=rho*abs(div_neg)*solver%amr%min_meshsize(lvl)**2/mu
            if (Res.gt.Res_tag) tagarr(i,j,k,1)=SETtag
            ! Tag near sphere surface
            dist=sphere_levelset([solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz],time)
            if (dist.lt.5.0_WP*dx.and.dist.gt.-dx) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Apply IB forcing - zero Q inside solid
   subroutine apply_ibm()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF
      real(WP) :: sum_VF,sum_VFQ1,sum_VFQ5,myVF
      integer :: i,j,k,lvl,ii,jj,kk
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pQ=>fs%Q%mf(lvl)%dataptr(mfi)
            pVF=>VF%mf(lvl)%dataptr(mfi)
            ! Get interior tilebox
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Skip pure fluid cells
               if (pVF(i,j,k,1).eq.1.0_WP) cycle
               ! Scale Q(2-4) by VF
               pQ(i,j,k,2)=pVF(i,j,k,1)*pQ(i,j,k,2)
               pQ(i,j,k,3)=pVF(i,j,k,1)*pQ(i,j,k,3)
               pQ(i,j,k,4)=pVF(i,j,k,1)*pQ(i,j,k,4)
               ! VF-weighted neighbor average for Q(1) and Q(5)
               sum_VF=0.0_WP; sum_VFQ1=0.0_WP; sum_VFQ5=0.0_WP
               do kk=-1,1; do jj=-1,1; do ii=-1,1
                  if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                  sum_VF  =sum_VF  +pVF(i+ii,j+jj,k+kk,1)
                  sum_VFQ1=sum_VFQ1+pVF(i+ii,j+jj,k+kk,1)*pQ(i+ii,j+jj,k+kk,1)
                  sum_VFQ5=sum_VFQ5+pVF(i+ii,j+jj,k+kk,1)*pQ(i+ii,j+jj,k+kk,5)
               end do; end do; end do
               if (sum_VF.gt.0.0_WP) then
                  pQ(i,j,k,1)=pVF(i,j,k,1)*pQ(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ1/sum_VF
                  pQ(i,j,k,5)=pVF(i,j,k,1)*pQ(i,j,k,5)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ5/sum_VF
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine apply_ibm
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: A,B,C
         ! EoS parameters
         call param_read('Gamma',Gamma)
         Pinf=0.0_WP
         ! Shock parameters (input is M2, post-shock lab Mach)
         call param_read('Mach number',M2)
         call param_read('Shock location',Xs)
         ! Post-shock normalization: rho2=1, u2=1, u1=0, T2=1
         rho2=1.0_WP
         p2=1.0_WP/(Gamma*M2**2)
         ! Quadratic for rho1: A*rho1^2 - B*rho1 + C = 0
         A=2.0_WP*Gamma*p2+(Gamma-1.0_WP)
         B=4.0_WP*Gamma*p2+(Gamma+1.0_WP)
         C=2.0_WP*Gamma*p2
         rho1=(B-sqrt(B**2-4.0_WP*A*C))/(2.0_WP*A)  ! smaller root for compression
         ! Shock-fixed frame velocities and pressure
         u1=1.0_WP/(1.0_WP-rho1)
         u2=u1-1.0_WP
         p1=p2-rho1/(1.0_WP-rho1)
         if (p1.le.0.0_WP) call die('[simulation_init] Cannot achieve requested Mach number - negative pre-shock pressure')
         ! Shock Mach number
         Ms=u1/sqrt(Gamma*p1/rho1)
         ! Shift to lab frame: pre-shock stationary
         u2=1.0_WP
         u1=0.0_WP
         ! Cv from T2=1
         Cv=p2/(rho2*(Gamma-1.0_WP))
         ! Viscous parameters
         call param_read('Reynolds number',Reynolds)
         call param_read('Prandtl number',Prandtl)
         ! Log shock conditions
         write(message,'("[Post-shock Mach] M2=",es12.5)') M2; call log(message)
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  rho1=",es12.5," p1=",es12.5)') rho1,p1; call log(message)
         write(message,'("[Post-shock] rho2=",es12.5," p2=",es12.5)') rho2,p2; call log(message)
         write(message,'("[Cv=",es12.5,"]")') Cv; call log(message)
      end block init_eos_and_flow
      
      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='amrcomp_sphere'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-05.0_WP; amr%xhi=+15.0_WP
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         call amr%initialize()
      end block create_amrgrid
      
      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
      end block initialize_timetracker

      ! Initialize compressible solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         ! Create flow solver
         call fs%initialize(amr=amr)
         ! Provide thermodynamic model
         fs%getP=>get_P
         fs%getC=>get_C
         fs%getT=>get_T
         ! Set initial conditions
         fs%Q%lo_bc(1,:)=amrex_bc_ext_dir
         fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%user_init=>shock_init
         ! Set boundary conditions
         fs%user_bc=>shock_dirichlet
      end block create_solver

      ! Create VF for IB
      create_VF: block
         use amrdata_class, only: amrex_interp_reinit
         call VF%initialize(amr,name='VF',ncomp=1,ng=fs%nover,interp=amrex_interp_reinit)
         VF%user_init=>init_VF
         call VF%register()
      end block create_VF
      
      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=5,ng=0,interp=amrex_interp_none); call dQdt%register()
         call k1%initialize(amr,name='k1',ncomp=5,ng=0,interp=amrex_interp_none); call k1%register()
         call k2%initialize(amr,name='k2',ncomp=5,ng=0,interp=amrex_interp_none); call k2%register()
         call k3%initialize(amr,name='k3',ncomp=5,ng=0,interp=amrex_interp_none); call k3%register()
         call k4%initialize(amr,name='k4',ncomp=5,ng=0,interp=amrex_interp_none); call k4%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=amrex_interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=amrex_interp_none); call Mach%register()
      end block create_workspace
      
      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Rec',Rec_tag)
         call param_read('Tagging Res',Res_tag)
         ! Create initial grid
         call amr%init_from_scratch(time=time%t)
         ! Compute viscosities
         call get_viscosities()
         ! Add artificial bulk viscosity
         call fs%add_viscartif(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding
      
      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'sphere')
         call viz%add_scalar(fs%Q,1,'RHO')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%U,1,'U')
         call viz%add_scalar(fs%V,1,'V')
         call viz%add_scalar(fs%W,1,'W')
         call viz%add_scalar(fs%I,1,'I')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(VF,1,'VF')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_viz
      
      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmin,'Pmin')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Qmin(1),'RHOmin')
         call mfile%add_column(fs%Qmax(1),'RHOmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Mass')
         call consfile%add_column(fs%Qint(2),'U Momentum')
         call consfile%add_column(fs%Qint(3),'V Momentum')
         call consfile%add_column(fs%Qint(4),'W Momentum')
         call consfile%add_column(fs%Qint(5),'Internal energy')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block create_monitors

   end subroutine simulation_init
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old conserved variables
         call fs%Qold%copy(src=fs%Q)
         
         ! ===== RK2/RK4 Stage 1: k1 = f(t, Q) =====
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,time=time%t); call k1%copy(src=dQdt)
         
         ! ===== RK2/RK4 Stage 2: k2 = f(t+dt/2, Q+k1/2) =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=0.5_WP*time%dt,src=k1)  ! Q=Qold+dt/2*k1
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         call apply_ibm(); call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt) ! IB forcing
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,time=time%t+0.5_WP*time%dt); call k2%copy(src=dQdt)

         ! ===== FOR RK2, RUN THIS =====
         call fs%Q%copy(src=fs%Qold)
         call fs%Q%saxpy(a=time%dt,src=k2)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call apply_ibm(); call fs%Q%average_down(); call fs%Q%fill(time=time%t) ! IB forcing
         ! =============================
         
         ! ===== RK4 Stage 3: k3 = f(t+dt/2, Q+k2/2) =====
         !call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=0.5_WP*time%dt,src=k2)  ! Q=Qold+dt/2*k2
         !call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         !call apply_ibm(); call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt) ! IB forcing
         !call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,time=time%t+0.5_WP*time%dt); call k3%copy(src=dQdt)
         
         ! ===== RK4 Stage 4: k4 = f(t+dt, Q+k3) =====
         !call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=time%dt,src=k3)  ! Q=Qold+dt*k3
         !call fs%Q%average_down(); call fs%Q%fill(time=time%t+time%dt)
         !call apply_ibm(); call fs%Q%average_down(); call fs%Q%fill(time=time%t+time%dt) ! IB forcing
         !call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,time=time%t+time%dt); call k4%copy(src=dQdt)
         
         ! ===== RK4 Combination: Q = Qold + (k1 + 2*k2 + 2*k3 + k4)/6 =====
         !call fs%Q%copy(src=fs%Qold)
         !call fs%Q%saxpy(a=time%dt/6.0_WP,src=k1)
         !call fs%Q%saxpy(a=time%dt/3.0_WP,src=k2)
         !call fs%Q%saxpy(a=time%dt/3.0_WP,src=k3)
         !call fs%Q%saxpy(a=time%dt/6.0_WP,src=k4)
         !call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         !call apply_ibm(); call fs%Q%average_down(); call fs%Q%fill(time=time%t) ! IB forcing

         ! Recompute primitive variables
         call fs%get_primitive(fs%Q)
         
         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add artificial bulk viscosity
         call fs%add_viscartif(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time%t)

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         
      end do
      
   end subroutine simulation_run
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call k1%finalize()
      call k2%finalize()
      call k3%finalize()
      call k4%finalize()
      call VF%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final
   
end module simulation
