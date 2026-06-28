!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use lss_class,         only: lss
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, an incompressible flow solver and corresponding time tracker
   type(fft2d),       public :: ps
   type(ddadi),       public :: vs
   type(incomp),      public :: fs
   type(lss),         public :: ls
   type(partmesh),    public :: pmesh
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,sfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: div_x,div_y,div_z
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Uib,Vib,Wib,srcM
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU

   !> Max timestep size for solid solver
   real(WP) :: ls_dt,ls_dt_max
   
   
contains
   
   
   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain
   
   
   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(div_x(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(div_y(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(div_z(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Uib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcM(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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

      ! Initialize Lagrangian solid solver
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
         
         ! Set material properties
         call param_read('Elastic Modulus',ls%elastic_modulus)
         call param_read('Poisson Ratio',ls%poisson_ratio)
         call param_read('Solid density',ls%rho)
         call param_read('Critical Energy Release Rate',ls%crit_energy)

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
         call param_read('Solid Spacing',dist)

         nx = ceiling(2.0_WP*R/dist)
         ny=nx
         nz=  floor(Lz/dist)
         dist = Lz/nz
         print *, "dist: ", dist

         call param_read('Horizon Ratio',ratio)
         ls%delta = dist*ratio

         ! Output some info on stretch
         mu=ls%elastic_modulus/(2.0_WP+2.0_WP*ls%poisson_ratio)
         kk=ls%elastic_modulus/(3.0_WP-6.0_WP*ls%poisson_ratio)
         max_stretch=sqrt(ls%crit_energy/((3.0_WP*mu+(kk-5.0_WP*mu/3.0_WP)*0.75_WP**4)*ls%delta))
         
         ! Only root process initializes solid particles
         if (ls%cfg%amRoot) then
            ! Read the STL file and get domain extents and levelset
            read_bin: block
               use messager, only: die
               integer :: p,iunit,ierr, wall_np, i, j, k
               real(WP) :: net_vol
               net_vol = 0.0_WP
               ! Read in grid definition
               wall_np = (ny)*(nz)*(nx) + ceiling(1.5_WP*R/dist) * nz
               call ls%resize(wall_np)
               p=0
               do i=1,nx
                  do j=1,ny
                     do k=1,nz
                       
                        x = (i-1) * dist - (R)!  - dist/2.0_WP);
                        y = (j-1) * dist - (R)! - dist/2.0_WP);
                        z = (k-1) * dist - Lz/2.0_WP + dist/2.0_WP;
                        ! print *, "looping"
                        ! print *, "dist: ", ((x)*(x) + y*y + z*z) 
                        if (((x)*(x) + y*y).ge.R*R) cycle;    ! Forming a sphere here
                        p = p+1
                        ls%p(p)%pos(1) = x
                        ls%p(p)%pos(2) = y
                        ls%p(p)%pos(3) = z
                        ls%p(p)%ipos=ls%p(p)%pos
                        ls%p(p)%displacement=0.0_WP
                        ls%p(p)%vol    = dist*dist*dist
                        ls%p(p)%id=-1
                        ls%p(p)%vel=[0.0_WP,0.0_WP,0.0_WP]
                        net_vol=net_vol+ls%p(p)%vol
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
                  end do
               end do

               ! ! Add a rectangular flap extending in +x direction behind the cylinder
               ! nx = ceiling(4.0_WP*R/dist)
               ! do i=1,nx
               !    do j=1,3
               !       do k=1,nz
               !          x = (i-1) * dist + (R)
               !          y = (j-2) * dist
               !          z = (k-1) * dist - Lz/2.0_WP
               !          p = p+1
               !          ls%p(p)%pos(1) = x
               !          ls%p(p)%pos(2) = y
               !          ls%p(p)%pos(3) = z
               !          ls%p(p)%ipos        = ls%p(p)%pos
               !          ls%p(p)%displacement= 0.0_WP
               !          ls%p(p)%vol         = dist*dist*dist
               !          ls%p(p)%id          = 1
               !          ls%p(p)%vel         = [0.0_WP,0.0_WP,0.0_WP]
               !          net_vol             = net_vol + ls%p(p)%vol
               !          ls%p(p)%Abond       = 0.0_WP
               !          ls%p(p)%Afluid      = 0.0_WP
               !          ls%p(p)%ind         = ls%cfg%get_ijk_global(ls%p(p)%pos,[ls%cfg%imin,ls%cfg%jmin,ls%cfg%kmin])
               !          ls%p(p)%i           = p
               !          ls%p(p)%flag        = 0
               !       end do
               !    end do
               ! end do
            np = p
            print*, "Nx: ", nx
            print*, "Ny: ", ny
            print*, "Nz: ", nz
            print*, "Net Force Volume", net_vol
            end block read_bin
         end if

         ! Communicate particles
         call ls%sync()





         ! Get initial volume fraction
         call ls%update_VF()
         
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
         pmesh=partmesh(nvar=5,nvec=7,name='solid')
         pmesh%varname(1)='failfrac'
         pmesh%varname(2)='id'
         pmesh%varname(3)='nbond'
         pmesh%varname(4)='von-Mises'
         pmesh%varname(5)='correcMag'
 

         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='bond_force'
         pmesh%vecname(3)='disp'
         pmesh%vecname(4)='t1'
         pmesh%vecname(5)='t2'
         pmesh%vecname(6)='tc'
         pmesh%vecname(7)='Afluid'
        
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
            pmesh%vec(:,4,i)  =ls%p(i)%t1
            pmesh%vec(:,5,i)  =ls%p(i)%t2
            pmesh%vec(:,6,i)  =ls%p(i)%tc
            pmesh%vec(:,7,i)  =ls%p(i)%Afluid
             

         end do
      end block create_pmesh
      
      
      ! Create a flow solver with inflow-outflow
      create_flow_solver: block
         use incomp_class, only: dirichlet,clipped_neumann
         real(WP) :: visc
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define boundary conditions
         call fs%add_bcond(name='inflow', type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )
         ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use random,       only: random_normal
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         real(WP) :: Uin
         ! Read inflow velocity
         call param_read('Inlet velocity',Uin)
         ! IB arrays
         Uib=0.0_WP; Vib=0.0_WP; Wib=0.0_WP; srcM=0.0_WP
         ! Make initial velocity field random to trigger transition
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  fs%U(i,j,k)=0.0_WP
                  fs%V(i,j,k)=0.0_WP
                  fs%W(i,j,k)=0.0_WP
               end do
            end do
         end do
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Set inflow velocity
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=Uin
         end do
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr(src=srcM)
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         resU=srcM/fs%rho           !< Careful, we need to provide
         call fs%get_div(src=resU)  !< a volume source term to div
         
      end block initialize_velocity
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='cylinder')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('velocity_s',Uib,Vib,Wib)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('VFs',ls%VF)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
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
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none
      real(WP) :: cfl
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call ls%get_cfl(time%dt,time%cfl)
         call fs%get_cfl(time%dt,cfl)
         call fs%get_cfl(time%dt,cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Advance solid solver
         solid: block
           real(WP) :: dt_done,mydt
           ! Compute divergence of fluid stress
           call fs%get_div_stress(divx=div_x(:,:,:),divy=div_y(:,:,:),divz=div_z(:,:,:))
           ! Sub-iteratore
           call ls%get_cfl(ls_dt,cfl=cfl)
           if (cfl.gt.0.0_WP) ls_dt=min(ls_dt*time%cflmax/cfl,ls_dt_max)
           dt_done=0.0_WP
           do while (dt_done.lt.time%dtmid)
              ! Decide the timestep size
              mydt=min(ls_dt,time%dtmid-dt_done)
              ! Advance particles
              call ls%advance(dt      =mydt,           &
              &               stress_x=div_x(:,:,:),&
              &               stress_y=div_y(:,:,:),&
              &               stress_z=div_z(:,:,:))
              ! Increment
              dt_done=dt_done+mydt
           end do
         end block solid


         ! Evaluate IB velocity and mass source
         calc_ib_velocity: block
            integer :: i,j,k
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     ! VF based velocity
                     Uib(i,j,k)=ls%VFU(i,j,k)/(sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))+epsilon(1.0_WP))
                     Vib(i,j,k)=ls%VFV(i,j,k)/(sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))+epsilon(1.0_WP))
                     Wib(i,j,k)=ls%VFW(i,j,k)/(sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))+epsilon(1.0_WP))
                  end do
               end do
            end do
            call cfg%sync(Uib)
            call cfg%sync(Vib)
            call cfg%sync(Wib)
            ! Compute IB mass source
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     srcM(i,j,k)=fs%rho*(ls%VF(i,j,k)*sum(fs%divp_x(:,i,j,k)*Uib(i:i+1,j,k))+&
                     &                                sum(fs%divp_y(:,i,j,k)*Vib(i,j:j+1,k))+&
                     &                                sum(fs%divp_z(:,i,j,k)*Wib(i,j,k:k+1)))
                  end do
               end do
            end do
            call cfg%sync(srcM)
         end block calc_ib_velocity

         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply direct IB forcing
            ibforcing: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_; do j=fs%cfg%jmin_,fs%cfg%jmax_; do i=fs%cfg%imin_,fs%cfg%imax_
                  fs%U(i,j,k)=(1.0_WP-sum(fs%itpr_x(:,i,j,k)*ls%VF(i-1:i,j,k)))*fs%U(i,j,k)+ls%VFU(i,j,k)
                  fs%V(i,j,k)=(1.0_WP-sum(fs%itpr_y(:,i,j,k)*ls%VF(i,j-1:j,k)))*fs%V(i,j,k)+ls%VFV(i,j,k)
                  fs%W(i,j,k)=(1.0_WP-sum(fs%itpr_z(:,i,j,k)*ls%VF(i,j,k-1:k)))*fs%W(i,j,k)+ls%VFW(i,j,k)
               end do; end do; end do
               call fs%cfg%sync(fs%U)
               call fs%cfg%sync(fs%V)
               call fs%cfg%sync(fs%W)
            end block ibforcing
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%correct_mfr(src=srcM)
            resU=srcM/fs%rho           !< Careful, we need to provide
            call fs%get_div(src=resU)  !< a volume source term to div
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         resU=srcM/fs%rho           !< Careful, we need to provide
         call fs%get_div(src=resU)  !< a volume source term to div
         
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
                 pmesh%vec(:,4,i)  =ls%p(i)%t1
                 pmesh%vec(:,5,i)  =ls%p(i)%t2
                 pmesh%vec(:,6,i)  =ls%p(i)%tc
                 pmesh%vec(:,7,i)  =ls%p(i)%Afluid
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call sfile%write()
         
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
      deallocate(div_x,div_y,div_z,resU,resV,resW,Ui,Vi,Wi,Uib,Vib,Wib,srcM,gradU)
      
   end subroutine simulation_final
   
   
end module simulation
