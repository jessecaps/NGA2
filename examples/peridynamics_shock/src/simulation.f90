!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP,SP
   use geometry,          only: cfg
   use spcomp_class,      only: spcomp
   use lss_class,         only: lss
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, an incompressible flow solver and corresponding time tracker
   type(spcomp),      public :: fs
   type(lss),         public :: ls
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile,sfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: dQdt
   real(WP), dimension(:,:,:)    , allocatable :: Ui,Vi,Wi,Ma,beta,visc,visc_t,div

   !> Post-shock viscosity and temperature
   real(WP) :: visc0,T0

   !> Equations of state
   real(WP) :: Pinf,Gamma,Cv,Prandtl

   !> Flow parameters
   real(WP) :: Ms,Xs,Rcyl
   real(WP) :: rho1,p1,u1,M1
   real(WP) :: rho2,p2,u2,M2
   real(WP) :: Re

   !> Max timestep size for solid solver
   real(WP) :: ls_dt,ls_dt_max
   
 contains


   !> Function that returns a smooth Heaviside of thickness delta
   real(WP) function Hshock(x,delta)
     real(WP), intent(in) :: x,delta
     ! Goes from 0 to 1 as x goes from begative to positive
     Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock


   !> P=EOS(RHO,I)
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
   !> RHO=f(T,P)
   pure real(WP) function get_RHO(T,P)
     implicit none
     real(WP), intent(in) :: T,P
     get_RHO=(P+Pinf)/(Cv*T*(Gamma-1.0_WP))
   end function get_RHO
   !> I=EOS(RHO,P)
   pure real(WP) function get_I(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_I=(P+Gamma*Pinf)/(RHO*(Gamma-1.0_WP))
   end function get_I
   !> C=f(RHO,P)
   pure real(WP) function get_C(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_C=sqrt(Gamma*(P+Pinf)/RHO)
   end function get_C
   !> S=f(RHO,P)
   pure real(WP) function get_S(RHO,P)
     implicit none
     real(WP), intent(in) :: RHO,P
     get_S=Cv*log((P+Pinf)/RHO**Gamma)
   end function get_S


   !> Calculate viscosities
   subroutine prepare_viscosities()
     implicit none
     integer :: i,j,k
     real(WP) :: S
     ! Get viscosity from Sutherland's law
     S=110.4_WP/273.15_WP*T0
     do k=fs%cfg%kmino_,fs%cfg%kmaxo_
        do j=fs%cfg%jmino_,fs%cfg%jmaxo_
           do i=fs%cfg%imino_,fs%cfg%imaxo_
              visc(i,j,k)=visc0*(T0+S)/(fs%T(i,j,k)+S)*(fs%T(i,j,k)/T0)**1.5_WP
           end do
        end do
     end do
     ! Get LAD
     call fs%get_viscartif(dt=time%dt,beta=beta); fs%BETA=fs%Q(:,:,:,1)*beta
     ! Get eddy viscosity
     call fs%get_vreman   (dt=time%dt,visc=visc_t); fs%VISC=fs%Q(:,:,:,1)*visc_t+visc
     ! Recompute thermal conductivity
     fs%diff=Gamma*Cv*fs%visc/Prandtl
     ! Add LAD
     fs%VISC=fs%VISC+0.002_WP*fs%BETA
   end subroutine prepare_viscosities


   !> Calculate velocity divergence
   subroutine get_div()
     implicit none
     integer :: i,j,k
     do k=fs%cfg%kmino_,fs%cfg%kmaxo_-1; do j=fs%cfg%jmino_,fs%cfg%jmaxo_-1; do i=fs%cfg%imino_,fs%cfg%imaxo_-1
        div(i,j,k)=fs%dxi*(fs%U(i+1,j,k)-fs%U(i,j,k))+fs%dyi*(fs%V(i,j+1,k)-fs%V(i,j,k))+fs%dzi*(fs%W(i,j,k+1)-fs%W(i,j,k))
     end do; end do; end do
     call fs%cfg%sync(div)
     if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) div(fs%cfg%imaxo,:,:)=div(fs%cfg%imaxo-1,:,:)
     if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.fs%cfg%npy) div(:,fs%cfg%jmaxo,:)=div(:,fs%cfg%jmaxo-1,:)
     if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.fs%cfg%npz) div(:,:,fs%cfg%kmaxo)=div(:,:,fs%cfg%kmaxo-1)
   end subroutine get_div


   !> Overwrite conserved variables using volume-of-solid IBM
   subroutine apply_ibm()
     implicit none
     integer :: i,j,k,ii,jj,kk
     real(WP) :: sum_VF,sum_VFQ1,sum_VFQ2
     real(WP), dimension(:,:,:), allocatable :: Q1old,Q2old
     allocate(Q1old(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_)); Q1old=fs%Q(:,:,:,1)
     allocate(Q2old(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_)); Q2old=fs%Q(:,:,:,2)
     do k=cfg%kmin_,cfg%kmax_
        do j=cfg%jmin_,cfg%jmax_
           do i=cfg%imin_,cfg%imax_
              if (ls%VF(i,j,k).eq.0.0_WP) cycle
              ! Neumann: VF-weighted neighbor average for Q(1) and Q(2)
              sum_VF=0.0_WP; sum_VFQ1=0.0_WP; sum_VFQ2=0.0_WP
              do kk=-1,1; do jj=-1,1; do ii=-1,1
                 if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                 sum_VF  =sum_VF  +(1.0_WP-ls%VF(i+ii,j+jj,k+kk))
                 sum_VFQ1=sum_VFQ1+(1.0_WP-ls%VF(i+ii,j+jj,k+kk))*Q1old(i+ii,j+jj,k+kk)
                 sum_VFQ2=sum_VFQ2+(1.0_WP-ls%VF(i+ii,j+jj,k+kk))*Q2old(i+ii,j+jj,k+kk)
              end do; end do; end do
              if (sum_VF.gt.0.0_WP) then
                 fs%Q(i,j,k,1)=(1.0_WP-ls%VF(i,j,k))*fs%Q(i,j,k,1)+ls%VF(i,j,k)*sum_VFQ1/sum_VF
                 fs%Q(i,j,k,2)=(1.0_WP-ls%VF(i,j,k))*fs%Q(i,j,k,2)+ls%VF(i,j,k)*sum_VFQ2/sum_VF
              end if
              ! No-slip now that density is determined
              fs%Q(i,j,k,3)=(1.0_WP-0.5_WP*(ls%VF(i-1,j,k)+ls%VF(i,j,k)))*fs%Q(i,j,k,3)+0.5_WP*(fs%Q(i-1,j,k,1)+fs%Q(i,j,k,1))*0.5_WP*(ls%VFU(i-1,j,k)+ls%VFU(i,j,k))
              fs%Q(i,j,k,4)=(1.0_WP-0.5_WP*(ls%VF(i,j-1,k)+ls%VF(i,j,k)))*fs%Q(i,j,k,4)+0.5_WP*(fs%Q(i,j-1,k,1)+fs%Q(i,j,k,1))*0.5_WP*(ls%VFV(i,j-1,k)+ls%VFV(i,j,k))
              fs%Q(i,j,k,5)=(1.0_WP-0.5_WP*(ls%VF(i,j,k-1)+ls%VF(i,j,k)))*fs%Q(i,j,k,5)+0.5_WP*(fs%Q(i,j,k-1,1)+fs%Q(i,j,k,1))*0.5_WP*(ls%VFW(i,j,k-1)+ls%VFW(i,j,k))
           end do
        end do
     end do
     ! Communicate
     call fs%cfg%sync(fs%Q(:,:,:,1))
     call fs%cfg%sync(fs%Q(:,:,:,2))
     call fs%cfg%sync(fs%Q(:,:,:,3))
     call fs%cfg%sync(fs%Q(:,:,:,4))
     call fs%cfg%sync(fs%Q(:,:,:,5))
     ! Rebuild primitive variables
     call fs%get_primitive()
   end subroutine apply_ibm


   !> Apply boundary conditions
   subroutine apply_bconds()
     implicit none
     integer :: i,j,k

     ! Apply clipped Neumann on primitive variables in x+
     if (.not.fs%cfg%xper.and.fs%cfg%iproc.eq.fs%cfg%npx) then
        do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do j=fs%cfg%jmino_,fs%cfg%jmaxo_
           ! Copy over from imax to imax+1 and above
           do i=fs%cfg%imax+1,fs%cfg%imaxo
              ! Copy primitive variables
              ls%VF(i,j,k)=ls%VF(fs%cfg%imax,j,k)
              fs%Q(i,j,k,1)=fs%Q(fs%cfg%imax,j,k,1)
              fs%P(i,j,k)=fs%P(fs%cfg%imax,j,k)
              fs%I(i,j,k)=fs%I(fs%cfg%imax,j,k)
              fs%U(i,j,k)=max(fs%U(fs%cfg%imax,j,k),0.0_WP)
              fs%V(i,j,k)=fs%V(fs%cfg%imax,j,k)
              fs%W(i,j,k)=fs%W(fs%cfg%imax,j,k)
           end do
        end do; end do
     end if

     ! Apply clipped Neumann on primitive variables in y+
     if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.fs%cfg%npy) then
        do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
           ! Copy over from jmax to jmax+1 and above
           do j=fs%cfg%jmax+1,fs%cfg%jmaxo
              ! Copy primitive variables
              ls%VF(i,j,k)=ls%VF(i,fs%cfg%jmax,k)
              fs%Q(i,j,k,1)=fs%Q(i,fs%cfg%jmax,k,1)
              fs%P(i,j,k)=fs%P(i,fs%cfg%jmax,k)
              fs%I(i,j,k)=fs%I(i,fs%cfg%jmax,k)
              fs%U(i,j,k)=fs%U(i,fs%cfg%jmax,k)
              fs%V(i,j,k)=max(fs%V(i,fs%cfg%jmax,k),0.0_WP)
              fs%W(i,j,k)=fs%W(i,fs%cfg%jmax,k)
           end do
        end do; end do
     end if

     ! Apply clipped Neumann on primitive variables in y-
     if (.not.fs%cfg%yper.and.fs%cfg%jproc.eq.1) then
        do k=fs%cfg%kmino_,fs%cfg%kmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
           ! First copy over V from jmin+1 to jmin
           fs%V(i,fs%cfg%jmin,k)=min(fs%V(i,fs%cfg%jmin+1,k),0.0_WP)
           ! Then copy over from jmin to jmin-1 and below
           do j=fs%cfg%jmino,fs%cfg%jmin-1
              ! Copy primitive variables
              ls%VF(i,j,k)=ls%VF(i,fs%cfg%jmin,k)
              fs%Q(i,j,k,1)=fs%Q(i,fs%cfg%jmin,k,1)
              fs%P(i,j,k)=fs%P(i,fs%cfg%jmin,k)
              fs%I(i,j,k)=fs%I(i,fs%cfg%jmin,k)
              fs%U(i,j,k)=fs%U(i,fs%cfg%jmin,k)
              fs%V(i,j,k)=min(fs%V(i,fs%cfg%jmin,k),0.0_WP)
              fs%W(i,j,k)=fs%W(i,fs%cfg%jmin,k)
           end do
        end do; end do
     end if

      ! Apply clipped Neumann on primitive variables in z+
     if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.fs%cfg%npz) then
        do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
           ! Copy over from kmax to kmax+1 and above
           do k=fs%cfg%kmax+1,fs%cfg%kmaxo
              ! Copy primitive variables
              ls%VF(i,j,k)=ls%VF(i,j,fs%cfg%kmax)
              fs%Q(i,j,k,1)=fs%Q(i,j,fs%cfg%kmax,1)
              fs%P(i,j,k)=fs%P(i,j,fs%cfg%kmax)
              fs%I(i,j,k)=fs%I(i,j,fs%cfg%kmax)
              fs%U(i,j,k)=fs%U(i,j,fs%cfg%kmax)
              fs%V(i,j,k)=fs%V(i,j,fs%cfg%kmax)
              fs%W(i,j,k)=max(fs%W(i,j,fs%cfg%kmax),0.0_WP)
           end do
        end do; end do
     end if

     ! Apply clipped Neumann on primitive variables in z-
     if (.not.fs%cfg%zper.and.fs%cfg%kproc.eq.1) then
        do j=fs%cfg%jmino_,fs%cfg%jmaxo_; do i=fs%cfg%imino_,fs%cfg%imaxo_
           ! First copy over W from kmin+1 to kmin
           fs%W(i,j,fs%cfg%kmin)=min(fs%W(i,j,fs%cfg%kmin+1),0.0_WP)
           ! Then copy over from kmin to kmin-1 and below
           do k=fs%cfg%kmino,fs%cfg%kmin-1
              ! Copy primitive variables
              ls%VF(i,j,k)=ls%VF(i,j,fs%cfg%kmin)
              fs%Q(i,j,k,1)=fs%Q(i,j,fs%cfg%kmin,1)
              fs%P(i,j,k)=fs%P(i,j,fs%cfg%kmin)
              fs%I(i,j,k)=fs%I(i,j,fs%cfg%kmin)
              fs%U(i,j,k)=fs%U(i,j,fs%cfg%kmin)
              fs%V(i,j,k)=fs%V(i,j,fs%cfg%kmin)
              fs%W(i,j,k)=min(fs%W(i,j,fs%cfg%kmin),0.0_WP)
           end do
        end do; end do
     end if

     ! Rebuild conserved quantities
     fs%Q(:,:,:,2)=fs%Q(:,:,:,1)*fs%I
     call fs%get_momentum()

   end subroutine apply_bconds


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists
      implicit none

      
      ! Create compressible flow solver
      create_flow_solver: block
        call fs%initialize(cfg=cfg,name='Compressible NS')
      end block create_flow_solver


      ! Allocate work arrays
      allocate_work_arrays: block
        allocate(dQdt  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:fs%nQ,1:4))
        allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(Ma    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(beta  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(visc  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(visc_t(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
        allocate(div   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays


      ! Initialize eos and flow parameters
      initialize_parameters: block
        use string,   only: str_long
        use messager, only: log
        use param,    only: param_read
        character(str_long) :: message
        ! Set Pinf to zero
        Pinf=0.0_WP
        ! Read in Gamma
        call param_read('Gamma',Gamma)
        ! Read in Prandtl number
        call param_read('Prandtl number',Prandtl)
        ! Read in shock Mach number and location
        call param_read('Shock Mach number',Ms)
        call param_read('Shock location',Xs)
        ! First generate static shock with normalized pre-shock conditions
        M1=Ms
        rho1=1.0_WP
        rho2=rho1*(Gamma+1.0_WP)*M1**2/((Gamma-1.0_WP)*M1**2+2.0_WP)
        p1=0.25_WP*rho1/Gamma*((Gamma+1.0_WP)*M1/(M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
        p2=p1*(2.0_WP*Gamma/(Gamma+1.0_WP)*(M1**2-1.0_WP)+1.0_WP)
        u1=M1*sqrt(Gamma*p1/rho1)
        u2=u1*rho1/rho2
        ! Now shift frame of reference to obtain moving shock
        u2=abs(u2-u1); M2=u2/sqrt(Gamma*p2/rho2); u1=0.0_WP; M1=u1/sqrt(Gamma*p1/rho1)
        ! Set heat capacities corresponding to a normalized pre-shock
        Cv=(p1+Pinf)/(rho1*(Gamma-1.0_WP))
        ! Get reference temperature based on post-shock conditions
        T0=get_T(rho2,p2)
        ! Define viscosity based on post-shock Reynolds number
        call param_read('Cylinder radius',Rcyl)
        call param_read('Reynolds number',Re); visc0=rho2*2.0_WP*Rcyl*u2/Re
        ! Output case info
        if (cfg%amRoot) then
           write(message,'("[Gas EOS]               =>  Gamma=",es12.5)')    Gamma; call log(message)
           write(message,'("[Gas EOS]               =>     Cv=",es12.5)')       Cv; call log(message)
           write(message,'("[Shock Mach number]     =>     Ms=",es12.5)')       Ms; call log(message)
           write(message,'("[Pre -shock conditions] =>   rho1=",es12.5)')     rho1; call log(message)
           write(message,'("[Pre -shock conditions] =>     p1=",es12.5)')       p1; call log(message)
           write(message,'("[Pre -shock conditions] =>     u1=",es12.5)')       u1; call log(message)
           write(message,'("[Pre -shock conditions] =>     M1=",es12.5)')       M1; call log(message)
           write(message,'("[Post-shock conditions] =>   rho2=",es12.5)')     rho2; call log(message)
           write(message,'("[Post-shock conditions] =>     p2=",es12.5)')       p2; call log(message)
           write(message,'("[Post-shock conditions] =>     u2=",es12.5)')       u2; call log(message)
           write(message,'("[Post-shock conditions] =>     M2=",es12.5)')       M2; call log(message)
           write(message,'("[Gas Reynolds]          =>     Re=",es12.5)')       Re; call log(message)
           write(message,'("[Gas viscosity]         =>     mu=",es12.5)')    visc0; call log(message)
        end if
      end block initialize_parameters


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
        time=timetracker(amRoot=cfg%amRoot)
        call param_read('Max timestep size',time%dtmax)
        call param_read('Max cfl number',time%cflmax)
        call param_read('Max time',time%tmax)
        time%dt=time%dtmax
        time%itmax=2
      end block initialize_timetracker


      ! Initialize Lagrangian solid solver
      initialize_lss: block
         real(WP) :: dx,mu,kk,max_stretch,Lx,Ly,Lz
         real(WP) :: xmin,xmax,ymin,ymax,zmin,zmax
         integer :: np,nt
         type triangle_type
            real(WP), dimension(3) :: norm
            real(WP), dimension(3) :: v1
            real(WP), dimension(3) :: v2
            real(WP), dimension(3) :: v3
         end type triangle_type
         type(triangle_type), dimension(:), allocatable :: t
         
         ! Create solver
         ls=lss(cfg=cfg,name='solid')
         
         ! Set material properties
         call param_read('Elastic Modulus',ls%elastic_modulus)
         call param_read('Poisson Ratio',ls%poisson_ratio)
         call param_read('Solid density',ls%rho)
         call param_read('Critical Energy Release Rate',ls%crit_energy)

         ! Maximum timestep size used for particles
         call param_read('Particle timestep size',ls_dt_max,default=huge(1.0_WP))
         ls_dt=min(ls_dt_max,time%dtmax)
         
         ! Discretization
         ls%delta=fs%cfg%min_meshsize

         ! Output some info on stretch
         mu=ls%elastic_modulus/(2.0_WP+2.0_WP*ls%poisson_ratio)
         kk=ls%elastic_modulus/(3.0_WP-6.0_WP*ls%poisson_ratio)
         max_stretch=sqrt(ls%crit_energy/((3.0_WP*mu+(kk-5.0_WP*mu/3.0_WP)*0.75_WP**4)*ls%delta))
         
         ! Only root process initializes solid particles
         if (ls%cfg%amRoot) then
            ! Read the STL file and get domain extents and levelset
            read_bin: block
              use messager, only: die
              integer :: p,iunit,ierr
              character(len=80) :: partfile
              call param_read('Particle file',partfile)
              open(newunit=iunit,file=trim(partfile),access="stream",form="unformatted",action="read",status="old",iostat=ierr)
              if(ierr.ne.0) call die('[read_stl] Could not open file: '//trim(partfile))
              read(iunit) np
              call ls%resize(np)
              do p=1,np
                 read(iunit) ls%p(p)%pos(1), ls%p(p)%pos(2), ls%p(p)%pos(3), ls%p(p)%vol
                 ! Set object id and velocity
                 ls%p(p)%id=-2
                 ls%p(p)%vel=0.0_WP
                 ! Zero out force
                 ls%p(p)%Abond=0.0_WP
                 ls%p(p)%Afluid=0.0_WP
                 ! Locate the particle on the mesh
                 ls%p(p)%ind=ls%cfg%get_ijk_global(ls%p(p)%pos,[ls%cfg%imin,ls%cfg%jmin,ls%cfg%kmin])
                 ! Assign a unique integer to particle
                 ls%p(p)%i=p
                 ! Activate the particle
                 ls%p(p)%flag=0
              end do
              close(iunit)
            end block read_bin
         end if
      
         ! Communicate particles
         call ls%sync()

         ! Get initial volume fraction
         call ls%update_VF()
         
         ! Initalize bonds
         call ls%bond_init()

         if (ls%cfg%amRoot) then
            print*,"===== Solid Setup Description ====="
            print*,'Number of particles', np
            print*,'Maximum stretching =',max_stretch
         end if
         
      end block initialize_lss


     ! Create partmesh object for visualizing Lagrangian particles
      create_pmesh: block
         use lss_class, only: max_bond
         integer :: i,n,nbond
         pmesh=partmesh(nvar=3,nvec=2,name='solid')
         pmesh%varname(1)='failfrac'
         pmesh%varname(2)='dilatation'
         pmesh%varname(3)='flag'
         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='bond_force'
         call ls%update_partmesh(pmesh)
         do i=1,ls%np_
            pmesh%var(1,i)=0.0_WP
            nbond=0
            do n=1,max_bond
               if (ls%p(i)%ibond(n).gt.0) nbond=nbond+1
            end do
            if (ls%p(i)%nbond.gt.0) then
               pmesh%var(1,i)=1.0_WP-real(nbond,WP)/real(ls%p(i)%nbond,WP)
            else
               pmesh%var(1,i)=0.0_WP
            end if
            pmesh%var(2,i)  =ls%p(i)%dil
            pmesh%var(3,i)  =ls%p(i)%flag
            pmesh%vec(:,1,i)=ls%p(i)%vel
            pmesh%vec(:,2,i)=ls%p(i)%Abond
         end do
      end block create_pmesh


      ! Initialize variables
      initialize_variables: block
        integer :: i,j,k
        ! Provide thermodynamic model
        fs%getP=>get_P; fs%getC=>get_C; fs%getS=>get_S; fs%getT=>get_T
        ! Initialize primary variables to normal shock
        do k=cfg%kmino_,cfg%kmaxo_
           do j=cfg%jmino_,cfg%jmaxo_
              do i=cfg%imino_,cfg%imaxo_
                 fs%U(i,j,k)  =u2*Hshock(Xs-fs%cfg%x(i),delta=0.5_WP*fs%dx)
                 fs%V(i,j,k)  =0.0_WP
                 fs%W(i,j,k)  =0.0_WP
                 fs%Q(i,j,k,1)=rho1+(rho2-rho1)*Hshock(Xs-fs%cfg%xm(i),delta=0.5_WP*fs%dx)
                 fs%P(i,j,k)  =p1  +(p2  -p1  )*Hshock(Xs-fs%cfg%xm(i),delta=0.5_WP*fs%dx)
                 fs%I(i,j,k)  =get_I(fs%Q(i,j,k,1),fs%P(i,j,k))
              end do
           end do
        end do
        ! Initialize conserved variables
        fs%Q(:,:,:,2)=fs%Q(:,:,:,1)*fs%I
        call fs%get_momentum()
        ! Rebuild primitive variables
        call fs%get_primitive()
        ! Interpolate velocity
        call fs%interp_vel(Ui,Vi,Wi)
        ! Compute local Mach number
        Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C
        ! Compute dilatation
        call get_div()
      end block initialize_variables


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='shock')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('P',fs%P)
         call ens_out%add_scalar('T',fs%T)
         call ens_out%add_scalar('Mach',Ma)
         call ens_out%add_scalar('beta',beta)
         call ens_out%add_scalar('visc',visc)
         call ens_out%add_scalar('visc_t',visc_t)
         call ens_out%add_scalar('div',div) 
         call ens_out%add_scalar('VFs',ls%VF)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create monitor files
      create_monitor: block
        real(WP) :: cfl
        ! Prepare some info about fields
        call ls%get_cfl(time%dt,time%cfl)
        call fs%get_cfl(time%dt,cfl); time%cfl=max(cfl,time%cfl)
        call fs%get_info()
        call ls%get_max()
        ! Create simulation monitor
        mfile=monitor(fs%cfg%amRoot,'simulation')
        call mfile%add_column(time%n,'Timestep number')
        call mfile%add_column(time%t,'Time')
        call mfile%add_column(time%dt,'Timestep size')
        call mfile%add_column(time%cfl,'Maximum CFL')
        call mfile%add_column(fs%Umax,'Umax')
        call mfile%add_column(fs%Vmax,'Vmax')
        call mfile%add_column(fs%Wmax,'Wmax')
        call mfile%add_column(fs%RHOmax,'max(RHO)')
        call mfile%add_column(fs%RHOmin,'min(RHO)')
        call mfile%add_column(fs%Pmax  ,'max(P)'  )
        call mfile%add_column(fs%Pmin  ,'min(P)'  )
        call mfile%add_column(fs%Tmax  ,'max(T)'  )
        call mfile%add_column(fs%Tmin  ,'min(T)'  )
        call mfile%write()
        ! Create CFL monitor
        cflfile=monitor(fs%cfg%amRoot,'cfl')
        call cflfile%add_column(time%n,'Timestep number')
        call cflfile%add_column(time%t,'Time')
        call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
        call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
        call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
        call cflfile%add_column(fs%CFLa_x,'Acoustic xCFL')
        call cflfile%add_column(fs%CFLa_y,'Acoustic yCFL')
        call cflfile%add_column(fs%CFLa_z,'Acoustic zCFL')
        call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
        call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
        call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
        call cflfile%write()
        ! Create conservation monitor
        consfile=monitor(fs%cfg%amRoot,'conservation')
        call consfile%add_column(time%n,'Timestep number')
        call consfile%add_column(time%t,'Time')
        call consfile%add_column(fs%Qint(1),'Mass')
        call consfile%add_column(fs%Qint(2),'Energy')
        call consfile%add_column(fs%Qint(3),'U Momentum')
        call consfile%add_column(fs%Qint(4),'V Momentum')
        call consfile%add_column(fs%Qint(5),'W Momentum')
        call consfile%add_column(fs%RHOKint,'Kinetic Energy')
        call consfile%add_column(fs%RHOSint,'Entropy')
        call consfile%write()
        ! Create solid monitor
        sfile=monitor(ls%cfg%amRoot,'solid')
        call sfile%add_column(time%n,'Timestep number')
        call sfile%add_column(time%t,'Time')
        call sfile%add_column(ls_dt,'Particle dt')
        call sfile%add_column(time%cfl,'Maximum CFL')
        call sfile%add_column(ls%np,'Particle number')
        call sfile%add_column(ls%VFmax,'VFmax')
        call sfile%add_column(ls%Umin,'Particle Umin')
        call sfile%add_column(ls%Umax,'Particle Umax')
        call sfile%add_column(ls%Vmin,'Particle Vmin')
        call sfile%add_column(ls%Vmax,'Particle Vmax')
        call sfile%add_column(ls%Wmin,'Particle Wmin')
        call sfile%add_column(ls%Wmax,'Particle Wmax')
        call sfile%add_column(ls%ibmForce(1),'Particle Fx')
        call sfile%add_column(ls%ibmForce(2),'Particle Fy')
        call sfile%add_column(ls%ibmForce(3),'Particle Fz')
        call sfile%write()
      end block create_monitor

    end subroutine simulation_init


    !> Perform an NGA2 simulation
    subroutine simulation_run
      implicit none
      real(WP) :: cfl

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call ls%get_cfl(time%dt,time%cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Advance solid solver
         solid: block
           real(WP) :: dt_done,mydt
           ! Compute divergence of fluid stress
           call fs%get_div_stress(divx=dQdt(:,:,:,1,1),divy=dQdt(:,:,:,2,1),divz=dQdt(:,:,:,3,1))
           ! Sub-iteratore
           call ls%get_cfl(ls_dt,cfl=cfl)
           if (cfl.gt.0.0_WP) ls_dt=min(ls_dt*time%cflmax/cfl,ls_dt_max)
           dt_done=0.0_WP
           do while (dt_done.lt.time%dtmid)
              ! Decide the timestep size
              mydt=min(ls_dt,time%dtmid-dt_done)
              ! Advance particles
              call ls%advance(dt      =mydt,           &
              &               stress_x=dQdt(:,:,:,1,1),&
              &               stress_y=dQdt(:,:,:,2,1),&
              &               stress_z=dQdt(:,:,:,3,1))
              ! Increment
              dt_done=dt_done+mydt
           end do
         end block solid

         ! Remember conserved variables
         fs%Qold=fs%Q

         ! Prepare SGS viscosity models
         call prepare_viscosities()

         ! ! First RK step ====================================================================================
         ! ! Get non-SL RHS and increment
         ! call fs%rhs(dQdt(:,:,:,:,1))
         ! fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         ! ! Apply IBM
         ! call apply_ibm()

         ! ! Second RK step ===================================================================================
         ! ! Get non-SL RHS and increment at midpoint
         ! call fs%rhs(dQdt(:,:,:,:,2))
         ! fs%Q=fs%Qold+time%dt*dQdt(:,:,:,:,2)
         ! ! Apply IBM
         ! call apply_ibm()

         ! First RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,1))
         ! Advance
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,1)
         ! Apply IBM
         call apply_ibm()

         ! Second RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,2))
         ! Advance
         fs%Q=fs%Qold+0.5_WP*time%dt*dQdt(:,:,:,:,2)
         ! Apply IBM
         call apply_ibm()

         ! Third RK step ====================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt=dQdt(:,:,:,:,3))
         ! Advance
         fs%Q=fs%Qold+1.0_WP*time%dt*dQdt(:,:,:,:,3)
         ! Apply IBM
         call apply_ibm()

         ! Fourth RK step ===================================================================================
         ! Get non-SL RHS and increment
         call fs%rhs(dQdt(:,:,:,:,4))
         ! Advance
         fs%Q=fs%Qold+time%dt/6.0_WP*(dQdt(:,:,:,:,1)+2.0_WP*dQdt(:,:,:,:,2)+2.0_WP*dQdt(:,:,:,:,3)+dQdt(:,:,:,:,4))
         ! Apply IBM
         call apply_ibm()

         ! Apply boundary conditions
         call apply_bconds()

         ! Interpolate velocity
         call fs%interp_vel(Ui,Vi,Wi)

         ! Compute local Mach number
         Ma=sqrt(Ui**2+Vi**2+Wi**2)/fs%C

         ! Compute dilatation
         call get_div()

         !> Perform and output monitoring
         call fs%get_info()
         call ls%get_max()
         call mfile%write()
         call cflfile%write()
         call consfile%write()
         call sfile%write()

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
              use lss_class, only: max_bond
              integer :: i,n,nbond
              call ls%update_partmesh(pmesh)
              do i=1,ls%np_
                 nbond=0
                 do n=1,max_bond
                    if (ls%p(i)%ibond(n).gt.0) nbond=nbond+1
                 end do
                 if (ls%p(i)%nbond.gt.0) then
                    pmesh%var(1,i)=1.0_WP-real(nbond,WP)/real(ls%p(i)%nbond,WP)
                 else
                    pmesh%var(1,i)=0.0_WP
                 end if
                 pmesh%var(2,i)  =ls%p(i)%dil
                 pmesh%var(3,i)  =ls%p(i)%flag
                 pmesh%vec(:,1,i)=ls%p(i)%vel
                 pmesh%vec(:,2,i)=ls%p(i)%Abond
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

      end do

 end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(dQdt,Ui,Vi,Wi,Ma,beta,visc,visc_t,div)
      
   end subroutine simulation_final
   
   
end module simulation
