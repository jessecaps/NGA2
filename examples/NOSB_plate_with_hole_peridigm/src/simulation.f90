!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP, I8
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use lsspd_class,       only: lss, pd_partition, PDC_MOVES,PDC_INTEGRATES,PDC_BONDS
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
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   type(partmesh),    public :: pmesh
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
      use param, only: param_read,param_exists
      use parallel, only: amRoot
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
         use mathtools, only: Pi
         
         integer(I8), allocatable :: gids(:),rgid(:)
         real(WP), allocatable :: pos(:,:),vel(:,:),voll(:),rpos(:,:),rvel(:,:),rvol(:)
         integer, allocatable :: flags(:),owner(:),rflag(:)
         integer :: i,j,k,n,nn,nr,nx,ny,nz,N_w
         real(WP) :: dx,x0,y0,z0,x1,y1,z1,x2,y2,z2,cx,sx,cy,sy,cz,sz
         real(WP) :: R,L,H,W,dist
         real(WP) :: rho,E,nu,elem,delta,contract,ratio

         ls=lss(cfg=cfg,name='solid')

         call param_read('R',R,default=0.015_WP)
         call param_read('L',L,default=0.1_WP)
         call param_read('H',H,default=0.1_WP)
         call param_read('W',W,default=0.02_WP)
         call param_read('N_w',N_w,default=5)
         elem = W/real(N_w,WP)
         call param_read('Horizon',         delta,default=3.0125_WP*elem)

         call param_read('Material density',rho,default=7850.0_WP)
         call param_read('Elastic modulus', E,default=2.0e11_WP)
         call param_read('Poisson ratio',   nu,default=0.3_WP)
         call param_read('Tau',             ls%tau, default=huge(1.0_WP))

         call param_read('Particle timestep size',ls_dt_max,default=huge(1.0_WP))

         call param_read('Horizon Ratio',ratio)

         call param_read('Unfreeze time', ls%unfreeze_time)
         ls%damping_rate=0.0_WP
         ls_dt=min(ls_dt_max,time%dtmax)
         ! Configure by field assignment (grid-free: no domain, no periodicity)
         
         ls%rho=rho; ls%elastic_modulus=E; ls%poisson_ratio=nu
         ls%delta=delta; ls%dV=elem**3
         ! Root builds the whole lattice; pd_partition routes it (gids are
         ! simply 1..n -- any unique positive keys work)
         
         

         nz = N_w
         elem = W/real(N_w,WP)
         ny = int(L/W*N_w)
         nx = int(H/W*N_w) + 12 ! 6 on each side, where we pull from
         nn=0 
         if(amRoot) then
            
            do k=1,nz; do j=1,ny; do i=1,nx
               
               x0 = (real(i,WP) - 0.5_WP*real(nx+1,WP))*elem
               y0 = (real(j,WP) - 0.5_WP*real(ny+1,WP))*elem
               z0 = (real(k,WP) - 0.5_WP*real(nz+1,WP))*elem
               if (((x0)*(x0) + y0*y0).le.R*R) cycle;
               nn=nn+1
            end do; end do; end do
         end if

         
         allocate(gids(max(nn,1)),pos(3,max(nn,1)),vel(3,max(nn,1)),flags(max(nn,1)),voll(max(nn,1)),owner(max(nn,1)))
         allocate(ls%icell(3,max(nn,1)))
         n=0
         do k=1,nz; do j=1,ny; do i=1,nx
            if (.not.amRoot) exit
            x0 = (real(i,WP) - 0.5_WP*real(nx+1,WP))*elem
            y0 = (real(j,WP) - 0.5_WP*real(ny+1,WP))*elem
            z0 = (real(k,WP) - 0.5_WP*real(nz+1,WP))*elem
            if (((x0)*(x0) + y0*y0).le.R*R) cycle;
            n=n+1
            pos(:,n)=[x0, y0, z0]
            vel(:,n)=[0.0_WP, 0.0_WP, 0.0_WP]
            flags(n)=PDC_MOVES+PDC_INTEGRATES+PDC_BONDS !< IVM, bitwise, this should keep it still?
            gids(n)=int(n,I8)
            voll(n)=elem**3
            if (i.lt.7) then; flags(n)=PDC_MOVES+PDC_BONDS; vel(:,n)=[-1.0e-3_WP, 0.0_WP, 0.0_WP]; end if
            if (i.gt.nx-6) then; flags(n)=PDC_MOVES+PDC_BONDS; vel(:,n)=[1.0e-3_WP, 0.0_WP, 0.0_WP]; end if
            
            ! ls%icell(:,n)=ls%cfg%get_ijk_global(pos(:,n),[ls%cfg%imin,ls%cfg%jmin,ls%cfg%kmin])
         end do; end do; end do
         call pd_partition(nn,gids,pos,vel,flags,voll,owner,nr,rgid,rpos,rvel,rflag,rvol)
         call ls%set_nodes(nr,rgid,rpos,rvel,rflag,rvol)
         
         call ls%detect_families()
         ! call ls%update_VF()

         
      end block initialize_lss



     ! Create partmesh object for visualizing Lagrangian particles
      create_pmesh: block
         use lss_class, only: max_bond
         integer :: i,n,nbond
         pmesh=partmesh(nvar=2,nvec=3,name='solid')
         pmesh%varname(1)='damage'
         pmesh%varname(2)='flag'
         ! pmesh%varname(3)='nbond' ! IVM, seems like we don't currently track this?
         ! mesh%varname(4)='von-Mises'
 

         pmesh%vecname(1)='velocity'
         ! pmesh%vecname(2)='bond_force'
         pmesh%vecname(2)='fluid_force'
         pmesh%vecname(3) = 'displacement'
         call ls%update_partmesh(pmesh)
         
         do i=1,ls%nown ! IVM, probably not the right thing 
            
            pmesh%vec(:,1,i)=ls%v(:,i)
            pmesh%vec(:,2,i)=ls%ff(:,i)
            pmesh%vec(:,3,i)=ls%y(:,i)-ls%x0(:,i)
            pmesh%var(2,i)=ls%flag(i)

            

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
      
      
      ! ! Initialize our velocity field
      ! initialize_velocity: block
      !    use random,       only: random_normal
      !    use incomp_class, only: bcond
      !    type(bcond), pointer :: mybc
      !    integer :: n,i,j,k
      !    real(WP) :: Uin
      !    ! Read inflow velocity
      !    call param_read('Inlet velocity',Uin)
      !    ! IB arrays
      !    Uib=0.0_WP; Vib=0.0_WP; Wib=0.0_WP; srcM=0.0_WP
      !    ! Make initial velocity field random to trigger transition
      !    do k=fs%cfg%kmin_,fs%cfg%kmax_
      !       do j=fs%cfg%jmin_,fs%cfg%jmax_
      !          do i=fs%cfg%imin_,fs%cfg%imax_
      !             fs%U(i,j,k)=0.0_WP
      !             fs%V(i,j,k)=0.0_WP
      !             fs%W(i,j,k)=0.0_WP
      !         end do
      !      end do
      !   end do
      !    call fs%cfg%sync(fs%U)
      !    call fs%cfg%sync(fs%V)
      !    call fs%cfg%sync(fs%W)
      !    ! Set inflow velocity
      !    call fs%get_bcond('inflow',mybc)
      !    do n=1,mybc%itr%no_
      !       i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
      !       fs%U(i,j,k)=Uin
      !    end do
      !    ! Compute MFR through all boundary conditions
      !    call fs%get_mfr()
      !    ! Adjust MFR for global mass balance
      !    call fs%correct_mfr(src=srcM)
      !    ! Compute cell-centered velocity
      !   call fs%interp_vel(Ui,Vi,Wi)
      !    ! Compute divergence
      !    resU=srcM/fs%rho           !< Careful, we need to provide
      !    call fs%get_div(src=resU)  !< a volume source term to div
         
      ! end block initialize_velocity


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
         call ens_out%add_scalar('SRCM',srcM)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      
      
      ! Create monitor files
      create_monitor: block
        real(WP) :: cfl
      !   ! Prepare some info about fields
      !    call fs%get_cfl(time%dt,time%cfl)
      !    call fs%get_max()
      !   ! Create simulation monitor
      !   mfile=monitor(fs%cfg%amRoot,'simulation')
      !   call mfile%add_column(time%n,'Timestep number')
      !   call mfile%add_column(time%t,'Time')
      !   call mfile%add_column(time%dt,'Timestep size')
      !   call mfile%add_column(time%cfl,'Maximum CFL')
      !   call mfile%add_column(fs%Umax,'Umax')
      !   call mfile%add_column(fs%Vmax,'Vmax')
      !   call mfile%add_column(fs%Wmax,'Wmax')
      !    call mfile%add_column(fs%Pmax,'Pmax')
      !    call mfile%add_column(fs%divmax,'Maximum divergence')
      !    call mfile%add_column(fs%psolv%it,'Pressure iteration')
      !    call mfile%add_column(fs%psolv%rerr,'Pressure error')
      !   call mfile%write()
      !   ! Create CFL monitor
      !   cflfile=monitor(fs%cfg%amRoot,'cfl')
      !   call cflfile%add_column(time%n,'Timestep number')
      !   call cflfile%add_column(time%t,'Time')
      !   call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
      !   call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
      !   call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
      !   call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
      !   call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
      !   call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
      !   call cflfile%write()

        ! Create solid monitor
        sfile=monitor(fs%cfg%amRoot,'solid')
        call sfile%add_column(time%n,'Timestep number')
        call sfile%add_column(time%t,'Time')
        call sfile%add_column(ls_dt,'Particle dt')
        call sfile%add_column(time%cfl,'Maximum CFL')
      !   call sfile%add_column(ls%np,'Particle number')
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

      print *, '================== simulation_init COMPLETE =================='

    end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
    subroutine simulation_run
      implicit none
      real(WP) :: cfl
      logical :: freeze_particles

      freeze_particles = .false.

      ! Perform time integration
      do while (.not.time%done())
         ! Increment time
         call ls%get_cfl(time%dt,time%cfl)
         ! call fs%get_cfl(time%dt,cfl); 
         time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()
      
         ! Advance solid solver
         solid: block
           real(WP) :: dt_done,mydt
           ! Compute divergence of fluid stress
         !   call fs%get_div_stress(divx=div_x(:,:,:),divy=div_y(:,:,:),divz=div_z(:,:,:))
           ! Sub-iteratore
           call ls%get_cfl(ls_dt,cfl=cfl)
           if (cfl.gt.0.0_WP) ls_dt=min(ls_dt*time%cflmax/cfl,ls_dt_max)
           dt_done=0.0_WP
           do while (dt_done.lt.time%dtmid)
              ! Decide the timestep size
            
              mydt=min(ls_dt,time%dtmid-dt_done)
              ! Advance particles
              if (time%t.gt.ls%unfreeze_time) freeze_particles=.true.
              call ls%advance(dt      =mydt, & 
              &               unfreeze = freeze_particles,            &
              &               div_stress_x=div_x(:,:,:),&
              &               div_stress_y=div_y(:,:,:),&
              &               div_stress_z=div_z(:,:,:))
              ! Increment
             
              dt_done=dt_done+mydt
           end do
           
         end block solid

         ! ! Evaluate IB velocity and mass source
         ! calc_ib_velocity: block
         !    integer :: i,j,k
         !    do k=fs%cfg%kmin_,fs%cfg%kmax_
         !       do j=fs%cfg%jmin_,fs%cfg%jmax_
         !          do i=fs%cfg%imin_,fs%cfg%imax_
         !             ! VF based velocity
         !             Uib(i,j,k)=0.5_WP*(ls%VFU(i-1,j,k)+ls%VFU(i,j,k))/(sum(fs%itpr_x(:,i,j,k)*ls%VF(i-1:i,j,k))+epsilon(1.0_WP))
         !             Vib(i,j,k)=0.5_WP*(ls%VFV(i,j-1,k)+ls%VFV(i,j,k))/(sum(fs%itpr_y(:,i,j,k)*ls%VF(i,j-1:j,k))+epsilon(1.0_WP))
         !             Wib(i,j,k)=0.5_WP*(ls%VFW(i,j,k-1)+ls%VFW(i,j,k))/(sum(fs%itpr_z(:,i,j,k)*ls%VF(i,j,k-1:k))+epsilon(1.0_WP))
         !          end do
         !       end do
         !    end do
         !    call cfg%sync(Uib)
         !    call cfg%sync(Vib)
         !    call cfg%sync(Wib)
         !    ! Compute IB mass source
         !    do k=fs%cfg%kmin_,fs%cfg%kmax_
         !       do j=fs%cfg%jmin_,fs%cfg%jmax_
         !          do i=fs%cfg%imin_,fs%cfg%imax_
         !             srcM(i,j,k)=fs%rho*(ls%VF(i,j,k)*(sum(fs%divp_x(:,i,j,k)*Uib(i:i+1,j,k))+&
         !             &                                sum(fs%divp_y(:,i,j,k)*Vib(i,j:j+1,k))+&
         !             &                                sum(fs%divp_z(:,i,j,k)*Wib(i,j,k:k+1))))

         !          end do
         !       end do
         !    end do
         !    call cfg%sync(srcM)
         ! end block calc_ib_velocity

         
         ! ! Remember old velocity
         ! fs%Uold=fs%U
         ! fs%Vold=fs%V
         ! fs%Wold=fs%W
         
         ! ! Perform sub-iterations
         ! do while (time%it.le.time%itmax)
            
         !    ! Build mid-time velocity
         !    fs%U=0.5_WP*(fs%U+fs%Uold)
         !    fs%V=0.5_WP*(fs%V+fs%Vold)
         !    fs%W=0.5_WP*(fs%W+fs%Wold)
            
         !    ! Explicit calculation of drho*u/dt from NS
         !    call fs%get_dmomdt(resU,resV,resW)
            
         !    ! Assemble explicit residual
         !    resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dtmid*resU
         !    resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dtmid*resV
         !    resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dtmid*resW
            
         !    ! Form implicit residuals
         !    call fs%solve_implicit(time%dtmid,resU,resV,resW)

         !    ! Apply these residuals
         !    fs%U=2.0_WP*fs%U-fs%Uold+resU
         !    fs%V=2.0_WP*fs%V-fs%Vold+resV
         !    fs%W=2.0_WP*fs%W-fs%Wold+resW
            
         !    ! Apply direct IB forcing
         !    ibforcing: block
         !       integer :: i,j,k
         !       do k=fs%cfg%kmin_,fs%cfg%kmax_; do j=fs%cfg%jmin_,fs%cfg%jmax_; do i=fs%cfg%imin_,fs%cfg%imax_
         !          fs%U(i,j,k)=(1.0_WP-sum(fs%itpr_x(:,i,j,k)*ls%VF(i-1:i,j,k)))*fs%U(i,j,k)+0.5_WP*(ls%VFU(i-1,j,k)+ls%VFU(i,j,k))
         !          fs%V(i,j,k)=(1.0_WP-sum(fs%itpr_y(:,i,j,k)*ls%VF(i,j-1:j,k)))*fs%V(i,j,k)+0.5_WP*(ls%VFV(i,j-1,k)+ls%VFV(i,j,k))
         !          fs%W(i,j,k)=(1.0_WP-sum(fs%itpr_z(:,i,j,k)*ls%VF(i,j,k-1:k)))*fs%W(i,j,k)+0.5_WP*(ls%VFW(i,j,k-1)+ls%VFW(i,j,k))
         !       end do; end do; end do
         !       call fs%cfg%sync(fs%U)
         !       call fs%cfg%sync(fs%V)
         !       call fs%cfg%sync(fs%W)
         !    end block ibforcing
            
         !    ! Apply other boundary conditions
         !    call fs%apply_bcond(time%t,time%dtmid)

         !    ! Solve Poisson equation
         !    call fs%correct_mfr(src=srcM)
         !    resU=srcM/fs%rho           !< Careful, we need to provide
         !    call fs%get_div(src=resU)  !< a volume source term to div
         !    fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dtmid
         !    fs%psolv%sol=0.0_WP
         !    call fs%psolv%solve()
         !    call fs%shift_p(fs%psolv%sol)
            
         !    ! Correct velocity
         !    call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         !    fs%P=fs%P+fs%psolv%sol
         !    fs%U=fs%U-time%dtmid*resU/fs%rho
         !    fs%V=fs%V-time%dtmid*resV/fs%rho
         !    fs%W=fs%W-time%dtmid*resW/fs%rho
            
         !    ! Increment sub-iteration counter
         !    time%it=time%it+1
            
         ! end do
         
         ! ! Recompute interpolated velocity and divergence
         ! call fs%interp_vel(Ui,Vi,Wi)
         ! resU=srcM/fs%rho           !< Careful, we need to provide
         ! call fs%get_div(src=resU)  !< a volume source term to div

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
              use lss_class, only: max_bond
              integer :: i,n,nbond
              call ls%update_partmesh(pmesh)
               do i=1,ls%nown ! IVM, probably not the right thing 
                  pmesh%vec(:,1,i)=ls%v(:,i)
                  pmesh%vec(:,2,i)=ls%ff(:,i)
                  pmesh%vec(:,3,i)=ls%y(:,i)-ls%x0(:,i)
                  pmesh%var(2,i)=ls%flag(i)
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

         ! ! Perform and output monitoring
         ! call fs%get_max()
         ! ! call ls%get_max() ! IVM, need a fix for this guy, I am sure something exists we can pull
         ! call mfile%write()
         ! call cflfile%write()
         ! call sfile%write()
         
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
