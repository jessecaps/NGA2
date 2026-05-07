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
   type(monitor) :: mfile,cflfile,consfile,sfile,dispfile
   
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

   integer :: target_index
   real(WP), dimension(3) :: target_position


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

   subroutine get_tracked_particle()
      use mpi_f08
      implicit none
      integer :: i, ierr
      real(WP) :: local_pos(3), global_pos(3)

      local_pos = 0.0_WP

      do i=1,ls%np_
         if (ls%p(i)%i.eq.target_index) then
            local_pos = ls%p(i)%pos            
         end if
      end do
      call MPI_ALLREDUCE(local_pos, global_pos, 3, MPI_DOUBLE_PRECISION, MPI_SUM, ls%cfg%comm, ierr)

      target_position = global_pos
   end subroutine

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


   !> Overwrite cosnerved variables using volume-of-solid IBM
   subroutine apply_ibm()
     implicit none
     integer :: i,j,k,ii,jj,kk
     real(WP) :: sum_VF,sum_VFQ1,sum_VFQ2
     do k=cfg%kmin_,cfg%kmax_
        do j=cfg%jmin_,cfg%jmax_
           do i=cfg%imin_,cfg%imax_
              if (ls%VF(i,j,k).eq.0.0_WP) cycle
              ! Neumann: VF-weighted neighbor average for Q(1) and Q(2)
              sum_VF=0.0_WP; sum_VFQ1=0.0_WP; sum_VFQ2=0.0_WP
              do kk=-1,1; do jj=-1,1; do ii=-1,1
                 if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                 sum_VF  =sum_VF  +(1.0_WP-ls%VF(i+ii,j+jj,k+kk))
                 sum_VFQ1=sum_VFQ1+(1.0_WP-ls%VF(i+ii,j+jj,k+kk))*fs%Q(i+ii,j+jj,k+kk,1)
                 sum_VFQ2=sum_VFQ2+(1.0_WP-ls%VF(i+ii,j+jj,k+kk))*fs%Q(i+ii,j+jj,k+kk,2)
              end do; end do; end do
              if (sum_VF.gt.0.0_WP) then
                 fs%Q(i,j,k,1)=(1.0_WP-ls%VF(i,j,k))*fs%Q(i,j,k,1)+ls%VF(i,j,k)*sum_VFQ1/sum_VF
                 fs%Q(i,j,k,2)=(1.0_WP-ls%VF(i,j,k))*fs%Q(i,j,k,2)+ls%VF(i,j,k)*sum_VFQ2/sum_VF
              end if
              ! No-slip now that density is determined
              fs%Q(i,j,k,3)=(1.0_WP-0.5_WP*(ls%VF(i-1,j,k)+ls%VF(i,j,k)))*fs%Q(i,j,k,3)+0.5_WP*(fs%Q(i-1,j,k,1)+fs%Q(i,j,k,1))*ls%VFU(i,j,k)
              fs%Q(i,j,k,4)=(1.0_WP-0.5_WP*(ls%VF(i,j-1,k)+ls%VF(i,j,k)))*fs%Q(i,j,k,4)+0.5_WP*(fs%Q(i,j-1,k,1)+fs%Q(i,j,k,1))*ls%VFV(i,j,k)
              fs%Q(i,j,k,5)=(1.0_WP-0.5_WP*(ls%VF(i,j,k-1)+ls%VF(i,j,k)))*fs%Q(i,j,k,5)+0.5_WP*(fs%Q(i,j,k-1,1)+fs%Q(i,j,k,1))*ls%VFW(i,j,k)
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

      ! Allocate work arrays
       allocate_work_arrays: block
         allocate(dQdt  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:5,1:4))
       end block allocate_work_arrays

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
        time=timetracker(amRoot=cfg%amRoot)
        call param_read('Max timestep size',time%dtmax)
        call param_read('Max cfl number',time%cflmax)
        call param_read('Max time',time%tmax)
        time%dt=time%dtmax
        time%itmax=2
      end block initialize_timetracker

       

      ! ! Initialize Lagrangian solid solver
      ! initialize_lss: block
      !    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_INTEGER
      !    real(WP) :: dx,mu,kk,max_stretch,Lx,Ly,Lz
      !    real(WP) :: xmin,xmax,ymin,ymax,zmin,zmax,ratio,P_load
      !    integer :: np,nt,nx,ny,nz,ierr,global_index
      !    type triangle_type
      !       real(WP), dimension(3) :: norm
      !       real(WP), dimension(3) :: v1
      !       real(WP), dimension(3) :: v2
      !       real(WP), dimension(3) :: v3
      !    end type triangle_type
      !    type(triangle_type), dimension(:), allocatable :: t
         
         
      !    ! Create solver
      !    ls=lss(cfg=cfg,name='solid')
      !    !call fs%initialize(cfg=cfg,name='Compressible NS')
     
         
      !    ! Set material properties
      !    call param_read('Elastic Modulus',ls%elastic_modulus)
      !    call param_read('Poisson Ratio',ls%poisson_ratio)
      !    call param_read('Solid density',ls%rho)
      !    call param_read('Critical Energy Release Rate',ls%crit_energy)

      !    ! Maximum timestep size used for particles
      !    call param_read('Particle timestep size',ls_dt_max,default=huge(1.0_WP))
      !    ls_dt=min(ls_dt_max,time%dtmax)
         
      !    ! Discretization
      !    ! ls%delta=fs%cfg%min_meshsize*1.01
      !    ! Load',P_load)
      !    call param_read('Lx',Lx)
      !    call param_read('Ly',Ly)
      !    call param_read('Lz',Lz)
      !    call param_read('Subdivisions',ny)
      !    nz = ny
      !    nx = NINT(Lx/Lz)*ny
      !    call param_read('Horizon Ratio',ratio)
      !    ls%delta = Ly/real(ny,WP)*ratio
      !    ! Output some info on stretch
      !    mu=ls%elastic_modulus/(2.0_WP+2.0_WP*ls%poisson_ratio)
      !    kk=ls%elastic_modulus/(3.0_WP-6.0_WP*ls%poisson_ratio)
      !    max_stretch=sqrt(ls%crit_energy/((3.0_WP*mu+(kk-5.0_WP*mu/3.0_WP)*0.75_WP**4)*ls%delta))
         
      !    ! Only root process initializes solid particles
      !    if (ls%cfg%amRoot) then
      !       ! Read the STL file and get domain extents and levelset
      !       print*, Lx * Ly * Lz / real(ny*nz*nx,WP)
      !       read_bin: block
              
      !         use messager, only: die
      !         integer :: p,iunit,ierr, wall_np, i, j, k
      !         global_index = 0
      !         target_index = 0
              
      !         ! Read in grid definition
      !         wall_np = ny*nz*(nx+3)
      !       !   call ls%resize(np+wall_np)
      !         call ls%resize(wall_np)
      !         p=0
      !         do i=1,nx+3
      !           do j=1,ny
      !             do k=1,nz
      !               p = p+1
      !               ls%p(p)%pos(1) = (i-1) * (Lx/real(nx,WP))
      !               ls%p(p)%pos(2) = (j) * (Ly/real(ny,WP)) - Ly/2.0_WP
      !               ls%p(p)%pos(3) = (k) * (Lz/real(nz,WP)) - Lz/2.0_WP
      !               ls%p(p)%vol    = Lx * Ly * Lz / real(ny*nz*nx,WP)
      !               ls%p(p)%id=1
      !               if(i.le.3) ls%p(p)%id=-2
                    
      !               ls%p(p)%vel=[0.0_WP,0.0_WP,0.0_WP]
      !               ! Zero out force
      !               ls%p(p)%Abond=0.0_WP
      !               ! Zero out fluid unless end, using this for the load
      !               ls%p(p)%Afluid=0.0_WP
      !               !if(i.eq.nx+3) ls%p(p)%Afluid=[(P_load*Ly*Lz)/(ls%rho*ls%p(p)%vol),0.0_WP,0.0_WP]
      !               ! Locate the particle on the mesh
      !               ls%p(p)%ind=ls%cfg%get_ijk_global(ls%p(p)%pos,[ls%cfg%imin,ls%cfg%jmin,ls%cfg%kmin])
      !               ! Assign a unique integer to particle
      !               ls%p(p)%i=p
      !               ! Activate the particle
      !               ls%p(p)%flag=0
      !               if(i.eq.(nx/2+1).and.j.eq.(ny/2+1).and.k.eq.(nz/2+1)) target_index = p
      !             end do
      !           end do
      !         end do
             
      !       np = wall_np
      !       end block read_bin
      !    end if

      !    ! Allreduce with MPI_MAX ensures the nonzero index propagates to all
      !    call MPI_ALLREDUCE(target_index, global_index, 1, MPI_INTEGER, MPI_MAX, ls%cfg%comm, ierr)

      !    ! Update target_index globally
      !    target_index = global_index
         
      
      !    ! Communicate particles
      !    call ls%sync()

      !    call get_tracked_particle()

      !    ! Get initial volume fraction
      !    ! call ls%update_VF()
         
      !    ! Initalize bonds
      !    call ls%bond_init()

      !    if (ls%cfg%amRoot) then
      !       print*,"===== Solid Setup Description ====="
      !       print*,'Number of particles', np
      !       print*,'Maximum stretching =',max_stretch
      !    end if
         
      ! end block initialize_lss

      initialize_lss: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_INTEGER
         real(WP) :: dx,mu,kk,max_stretch,Lx,Ly,Lz,R,x,y,z,load_rate
         real(WP) :: xmin,xmax,ymin,ymax,zmin,zmax,ratio,dist
         integer :: np,nt,nx,ny,nz,ierr,global_index,N
         type triangle_type
            real(WP), dimension(3) :: norm
            real(WP), dimension(3) :: v1
            real(WP), dimension(3) :: v2
            real(WP), dimension(3) :: v3
         end type triangle_type
         type(triangle_type), dimension(:), allocatable :: t

         
         
        
         ! Create solver
         ls=lss(cfg=cfg,name='solid')
         !call fs%initialize(cfg=cfg,name='Compressible NS')
     
         
         ! Set material properties
         call param_read('Elastic Modulus',ls%elastic_modulus)
         call param_read('Poisson Ratio',ls%poisson_ratio)
         call param_read('Solid density',ls%rho)
         call param_read('Critical Energy Release Rate',ls%crit_energy)
         call param_read('Solid Damping Constant',ls%beta)
         call param_read('Cool down time',ls%cool_down_time)
         call param_read('Continuous damping',ls%continuous_damping)

         ! Maximum timestep size used for particles
         call param_read('Particle timestep size',ls_dt_max,default=huge(1.0_WP))
         ls_dt=min(ls_dt_max,time%dtmax)
         
         ! Discretization
         ! ls%delta=fs%cfg%min_meshsize*1.01
         ! Load',P_load)
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         call param_read('Lz',Lz)
         call param_read('R',R)
         ! call param_read('Solid Spacing',dist)
         call param_read('N Across',N)
         call param_read('Load Rate',load_rate)
         ! Lx = 1.0_WP
         ! Ly = 1.0_WP
         ! dist = 0.01_WP ! Space between particles
         ! Lx = Lx + 6.0_WP * dist
         ! Ly = Ly + 3.0_WP * dist
         dist = Lz/N                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
         
         nz = N
         nx = ceiling(Lx/Lz)*N + 12
         ny = ceiling(Ly/Lz)*N
         call param_read('Horizon Ratio',ratio)
         ls%delta = dist*ratio
         
         ! Output some info on stretch
         mu=ls%elastic_modulus/(2.0_WP+2.0_WP*ls%poisson_ratio)
         kk=ls%elastic_modulus/(3.0_WP-6.0_WP*ls%poisson_ratio)
         max_stretch=sqrt(ls%crit_energy/((3.0_WP*mu+(kk-5.0_WP*mu/3.0_WP)*0.75_WP**4)*ls%delta))
         
         ! Only root process initializes solid particles
         if (ls%cfg%amRoot) then
            read_bin: block
              
              use messager, only: die
              integer :: p,iunit,ierr, wall_np, i, j, k
              real(WP) :: net_vol
              net_vol = 0.0_WP
              global_index = 0
              target_index = 0
              ! Read in grid definition
              wall_np = (ny)*(nz)*(nx)
            !   call ls%resize(np+wall_np)
              call ls%resize(wall_np)
              p=0
              do i=1,nx
                do j=1,ny
                  do k=1,nz
                    x = (i-7) * dist - (Lx/2.0_WP - dist/2.0_WP);
                    y = (j-1) * dist - (Ly/2.0_WP - dist/2.0_WP);
                    z = (k-1) * dist - (Lz/2.0_WP - dist/2.0_WP);
                    !if ((x*x + y*y).gt.R*R) cycle;
                    if (((x)*(x) + y*y).lt.R*R) cycle;
                    p = p+1
                    ls%p(p)%pos(1) = x
                    ls%p(p)%pos(2) = y
                    ls%p(p)%pos(3) = z
                    ls%p(p)%ipos=ls%p(p)%pos
                    ls%p(p)%displacement=0.0_WP
                    ls%p(p)%vol    = dist*dist*dist
                    ls%p(p)%id=1
                    if(i.le.6) ls%p(p)%id=-1
                    if(i.ge.nx-5) ls%p(p)%id=-1
                    ls%p(p)%vel=[0.0_WP,0.0_WP,0.0_WP]
                    if(i.le.6) ls%p(p)%vel=[-load_rate/2.0_WP,0.0_WP,0.0_WP]
                    if(i.ge.nx-5) ls%p(p)%vel=[load_rate/2.0_WP,0.0_WP,0.0_WP]
                    ! if(i.ge.nx-5) net_vol=net_vol+ls%p(p)%vol
                    ! Zero out force
                    ls%p(p)%Abond=0.0_WP
                    ! Zero out fluid unless end, using this for the load
                    ls%p(p)%Afluid=0.0_WP
                    ! Locate the particle on the mesh
                    ls%p(p)%ind=ls%cfg%get_ijk_global(ls%p(p)%pos,[ls%cfg%imin,ls%cfg%jmin,ls%cfg%kmin])
                    ! Assign a unique integer to particle
                    ls%p(p)%i=p
                    ! Activate the particle
                    ls%p(p)%flag=0
                    if(i.eq.(nx/2+2).and.j.eq.(ny/2+1).and.k.eq.(nz/2+1)) target_index = p
                  end do
                end do
              end do
             
            np = wall_np
            print*, "Nx: ", nx
            print*, "Ny: ", ny
            print*, "Nz: ", nz
            print*, "Net Force Volume", net_vol
            print*, "Used Volume", (dist**3 * ny * nz * 6)
            end block read_bin
         end if

         ! Allreduce with MPI_MAX ensures the nonzero index propagates to all
         call MPI_ALLREDUCE(target_index, global_index, 1, MPI_INTEGER, MPI_MAX, ls%cfg%comm, ierr)

         ! Update target_index globally
         target_index = global_index
         
      
         ! Communicate particles
         call ls%sync()

         call get_tracked_particle()

         ! Get initial volume fraction
         ! call ls%update_VF()
         
         ! Initalize bonds
         call ls%bond_init()
         call ls%get_bond_force()
         call ls%sync()

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
         pmesh=partmesh(nvar=5,nvec=3,name='solid')
         pmesh%varname(1)='failfrac'
         pmesh%varname(2)='id'
         pmesh%varname(3)='nbond'
         pmesh%varname(4)='von-Mises'
         pmesh%varname(5)='correcMag'
 

         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='bond_force'
         pmesh%vecname(3)='disp'
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
            pmesh%var(2,i)  =ls%p(i)%id
            pmesh%vec(:,1,i)=ls%p(i)%vel
            pmesh%vec(:,2,i)=ls%p(i)%Abond*ls%rho*ls%p(i)%vol
            pmesh%var(3,i)  =ls%p(i)%nbond
            pmesh%var(4,i)  =ls%p(i)%vonMises
            pmesh%var(5,i)  =ls%p(i)%correcMag
            pmesh%vec(:,3,i)  =ls%p(i)%displacement

         end do
      end block create_pmesh

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='shock')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create monitor files
      create_monitor: block
        real(WP) :: cfl
        ! Prepare some info about fields
        call ls%get_cfl(time%dt,time%cfl)
        call ls%get_max()
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
        dispfile=monitor(ls%cfg%amRoot,'displacement')
        call dispfile%add_column(time%n,'Timestep number')
        call dispfile%add_column(time%t,'Time')
        call dispfile%add_column(ls_dt,'Particle dt')
        call dispfile%add_column(target_position(1),'X')
        call dispfile%add_column(target_position(2),'Y')
        call dispfile%add_column(target_position(3),'Z')
        call dispfile%write()
      end block create_monitor

    end subroutine simulation_init


    !> Perform an NGA2 simulation
    subroutine simulation_run
      implicit none
      real(WP) :: cfl
      logical :: cool_down_time

      cool_down_time = .false.
      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call ls%get_cfl(time%dt,time%cfl)
         ! call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Advance solid solver
         solid: block
           real(WP) :: dt_done,mydt
           ! Sub-iteratore
           call ls%get_cfl(ls_dt,cfl=cfl)
           if (cfl.gt.0.0_WP) ls_dt=min(ls_dt*time%cflmax/cfl,ls_dt_max)
           dt_done=0.0_WP
           do while (dt_done.lt.time%dtmid)
              ! Decide the timestep size
               if (time%t.gt.ls%cool_down_time) cool_down_time = .true.
               mydt=min(ls_dt,time%dtmid-dt_done)
               !  ! Advance particles
                call ls%advance(dt      =mydt,cool_down = cool_down_time, continuous = ls%continuous_damping)
               !  ! Increment
               dt_done=dt_done+mydt



            ! mydt=min(ls_dt,time%dtmid-dt_done)
            !     ! Advance particles
            !    call ls%advance(dt      =mydt)
            !    !  ! Increment
            !    dt_done=dt_done+mydt
              
           end do 
         end block solid

         !> Perform and output monitoring
         call ls%get_max()
         call get_tracked_particle()
         call sfile%write()
         call dispfile%write()
         

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
                 pmesh%var(2,i)  =ls%p(i)%id
                 pmesh%vec(:,1,i)=ls%p(i)%vel
                 pmesh%vec(:,2,i)=ls%p(i)%Abond*ls%rho*ls%p(i)%vol
                 pmesh%var(3,i)  =ls%p(i)%nbond
                 pmesh%var(4,i)  =ls%p(i)%vonMises
                 pmesh%var(5,i)  =ls%p(i)%correcMag
                 pmesh%vec(:,3,i)  =ls%p(i)%displacement


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
      deallocate(dQdt)
   end subroutine simulation_final
   
   
end module simulation
