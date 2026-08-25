!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use pardata_class,     only: pardata
   implicit none
   private
   
   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type(incomp) :: fs
   type(fft2d)  :: ps
   type(timetracker) :: time
   
   !> Implicit solver
   logical :: use_implicit
   type(ddadi) :: vs
   
   !> Ensight postprocessing
   type(event) :: ens_evt
   type(ensight) :: ens_out
   
   !> Simulation monitor file
   type(monitor) :: mfile,pfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Work arrays
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Gib0,VF0
   real(WP), dimension(3) :: pos,vel,force
   real(WP) :: visc,rhop,dp,mass,CFLp
   
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
   

   !> Compute force on particle
   subroutine get_force()
     use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_IN_PLACE
     use parallel, only: MPI_REAL_WP
     implicit none
     integer :: i,j,k,ii,jj,kk,ierr
     real(WP) :: vol
     real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ

     ! Allocate flux arrays
     allocate(FX(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
     allocate(FY(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
     allocate(FZ(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))

     ! Stress in x
     do kk=fs%cfg%kmin_,fs%cfg%kmax_+1
        do jj=fs%cfg%jmin_,fs%cfg%jmax_+1
           do ii=fs%cfg%imin_,fs%cfg%imax_+1
              ! Fluxes on x-face
              i=ii-1; j=jj-1; k=kk-1
              FX(i,j,k)=+fs%visc(i,j,k)*(sum(fs%grdu_x(:,i,j,k)*fs%U(i:i+1,j,k))+sum(fs%grdu_x(:,i,j,k)*fs%U(i:i+1,j,k)) &
                   &         -2.0_WP/3.0_WP*(sum(fs%divp_x(:,i,j,k)*fs%U(i:i+1,j,k))+sum(fs%divp_y(:,i,j,k)*fs%V(i,j:j+1,k))+sum(fs%divp_z(:,i,j,k)*fs%W(i,j,k:k+1)))) &
                   &         -fs%P(i,j,k)
              ! Fluxes on y-face
              i=ii; j=jj; k=kk
              FY(i,j,k)=+sum(fs%itp_xy(:,:,i,j,k)*fs%visc(i-1:i,j-1:j,k))*(sum(fs%grdu_y(:,i,j,k)*fs%U(i,j-1:j,k))+sum(fs%grdv_x(:,i,j,k)*fs%V(i-1:i,j,k)))
              ! Fluxes on z-face
              i=ii; j=jj; k=kk
              FZ(i,j,k)=+sum(fs%itp_xz(:,:,i,j,k)*fs%visc(i-1:i,j,k-1:k))*(sum(fs%grdu_z(:,i,j,k)*fs%U(i,j,k-1:k))+sum(fs%grdw_x(:,i,j,k)*fs%W(i-1:i,j,k)))
           end do
        end do
     end do
     ! Integrate the stress
     force(1)=0.0_WP
     do k=fs%cfg%kmin_,fs%cfg%kmax_
        do j=fs%cfg%jmin_,fs%cfg%jmax_
           do i=fs%cfg%imin_,fs%cfg%imax_
              vol=sum(fs%itpr_x(:,i,j,k)*(1.0_WP-cfg%VF(i-1:i,j,k))*fs%cfg%vol(i-1:i,j,k))
              force(1)=force(1)+(sum(fs%divu_x(:,i,j,k)*FX(i-1:i,j,k))+&
                   &             sum(fs%divu_y(:,i,j,k)*FY(i,j:j+1,k))+&
                   &             sum(fs%divu_z(:,i,j,k)*FZ(i,j,k:k+1)))*vol
           end do
        end do
     end do

     ! Stress in y
     do kk=fs%cfg%kmin_,fs%cfg%kmax_+1
        do jj=fs%cfg%jmin_,fs%cfg%jmax_+1
           do ii=fs%cfg%imin_,fs%cfg%imax_+1
              ! Fluxes on x-face
              i=ii; j=jj; k=kk
              FX(i,j,k)=+sum(fs%itp_xy(:,:,i,j,k)*fs%visc(i-1:i,j-1:j,k))*(sum(fs%grdv_x(:,i,j,k)*fs%V(i-1:i,j,k))+sum(fs%grdu_y(:,i,j,k)*fs%U(i,j-1:j,k)))
              ! Fluxes on y-face
              i=ii-1; j=jj-1; k=kk-1
              FY(i,j,k)=+fs%visc(i,j,k)*(sum(fs%grdv_y(:,i,j,k)*fs%V(i,j:j+1,k))+sum(fs%grdv_y(:,i,j,k)*fs%V(i,j:j+1,k)) &
                   &         -2.0_WP/3.0_WP*(sum(fs%divp_x(:,i,j,k)*fs%U(i:i+1,j,k))+sum(fs%divp_y(:,i,j,k)*fs%V(i,j:j+1,k))+sum(fs%divp_z(:,i,j,k)*fs%W(i,j,k:k+1)))) &
                   &         -fs%P(i,j,k)
              ! Fluxes on z-face
              i=ii; j=jj; k=kk
              FZ(i,j,k)=+sum(fs%itp_yz(:,:,i,j,k)*fs%visc(i,j-1:j,k-1:k))*(sum(fs%grdv_z(:,i,j,k)*fs%V(i,j,k-1:k))+sum(fs%grdw_y(:,i,j,k)*fs%W(i,j-1:j,k)))
           end do
        end do
     end do
     ! Integrate the stress
     force(2)=0.0_WP
     do k=fs%cfg%kmin_,fs%cfg%kmax_
        do j=fs%cfg%jmin_,fs%cfg%jmax_
           do i=fs%cfg%imin_,fs%cfg%imax_
              vol=sum(fs%itpr_y(:,i,j,k)*(1.0_WP-cfg%VF(i,j-1:j,k))*fs%cfg%vol(i,j-1:j,k))
              force(2)=force(2)+(sum(fs%divv_x(:,i,j,k)*FX(i:i+1,j,k))+&
                   &             sum(fs%divv_y(:,i,j,k)*FY(i,j-1:j,k))+&
                   &             sum(fs%divv_z(:,i,j,k)*FZ(i,j,k:k+1)))*vol
           end do
        end do
     end do

     ! Stress in z
     do kk=fs%cfg%kmin_,fs%cfg%kmax_+1
        do jj=fs%cfg%jmin_,fs%cfg%jmax_+1
           do ii=fs%cfg%imin_,fs%cfg%imax_+1
              ! Fluxes on x-face
              i=ii; j=jj; k=kk
              FX(i,j,k)=+sum(fs%itp_xz(:,:,i,j,k)*fs%visc(i-1:i,j,k-1:k))*(sum(fs%grdw_x(:,i,j,k)*fs%W(i-1:i,j,k))+sum(fs%grdu_z(:,i,j,k)*fs%U(i,j,k-1:k)))
              ! Fluxes on y-face
              i=ii; j=jj; k=kk
              FY(i,j,k)=+sum(fs%itp_yz(:,:,i,j,k)*fs%visc(i,j-1:j,k-1:k))*(sum(fs%grdw_y(:,i,j,k)*fs%W(i,j-1:j,k))+sum(fs%grdv_z(:,i,j,k)*fs%V(i,j,k-1:k)))
              ! Fluxes on z-face
              i=ii-1; j=jj-1; k=kk-1
              FZ(i,j,k)=+fs%visc(i,j,k)*(sum(fs%grdw_z(:,i,j,k)*fs%W(i,j,k:k+1))+sum(fs%grdw_z(:,i,j,k)*fs%W(i,j,k:k+1)) &
                   &         -2.0_WP/3.0_WP*(sum(fs%divp_x(:,i,j,k)*fs%U(i:i+1,j,k))+sum(fs%divp_y(:,i,j,k)*fs%V(i,j:j+1,k))+sum(fs%divp_z(:,i,j,k)*fs%W(i,j,k:k+1)))) &
                   &         -fs%P(i,j,k)
           end do
        end do
     end do
     ! Integrate stress
     force(3)=0.0_WP
     do k=fs%cfg%kmin_,fs%cfg%kmax_
        do j=fs%cfg%jmin_,fs%cfg%jmax_
           do i=fs%cfg%imin_,fs%cfg%imax_
              vol=sum(fs%itpr_z(:,i,j,k)*(1.0_WP-cfg%VF(i,j,k-1:k))*fs%cfg%vol(i,j,k-1:k))
              force(3)=force(3)+(sum(fs%divw_x(:,i,j,k)*FX(i:i+1,j,k))+&
                   &             sum(fs%divw_y(:,i,j,k)*FY(i,j:j+1,k))+&
                   &             sum(fs%divw_z(:,i,j,k)*FZ(i,j,k-1:k)))*vol
           end do
        end do
     end do

     ! Communicate
     call MPI_ALLREDUCE(MPI_IN_PLACE,force,3,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)

     ! Deallocate flux arrays
     deallocate(FX,FY,FZ)

   end subroutine get_force

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
        use incomp_class,   only: dirichlet
        ! Create flow solver
        fs=incomp(cfg=cfg,name='Incompressible NS')
        ! Define boundary conditions
        call fs%add_bcond(name= 'bottom',type=dirichlet,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
        call fs%add_bcond(name= 'top',   type=dirichlet,locator=right_of_domain,face='x',dir=+1,canCorrect=.false.)
        ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Check if implicit velocity solver is used
         call param_read('Use implicit solver',use_implicit)
         if (use_implicit) then
            ! Configure implicit solver
            vs=ddadi(cfg=cfg,name='Velocity',nst=7)
            ! Finish flow solver setup
            call fs%setup(pressure_solver=ps,implicit_solver=vs)
         else
            ! Finish flow solver setup
            call fs%setup(pressure_solver=ps)
         end if
      end block create_flow_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))   
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Gib0(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(VF0(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize our velocity field
      initialize_velocity: block
        ! Zero out velocity
        fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
        ! Compute cell-centered velocity
        call fs%interp_vel(Ui,Vi,Wi)
        ! Compute divergence
        call fs%get_div()
        ! Compute MFR through all boundary conditions
        call fs%get_mfr()
        ! Store initial levelset and VF (absent of particles)
        Gib0=cfg%Gib; VF0=cfg%VF
      end block initialize_velocity


      ! Initialize our particle
      create_particle: block
        use mathtools, only: Pi
        use ibconfig_class, only: bigot,sharp
        integer :: i,j,k
        ! Read in particle properties
        call param_read('Particle density',rhop)
        call param_read('Particle diameter',dp)
        call param_read('Particle position',pos)
        ! Compute mass
        mass=rhop*Pi*dp**3/6.0_WP
        ! Zero-out
        force=0.0_WP; vel=0.0_WP; CFLp=0.0_WP
        ! Overwrite levelset and volume fraction
        do k=cfg%kmino_,cfg%kmaxo_
           do j=cfg%jmino_,cfg%jmaxo_
              do i=cfg%imino_,cfg%imaxo_
                 cfg%Gib(i,j,k)=0.5_WP*dp-sqrt((cfg%xm(i)-pos(1))**2+(cfg%ym(j)-pos(2))**2+(cfg%zm(k)-pos(3))**2)
              end do
           end do
        end do
        ! Get normal vector
        call cfg%calculate_normal()
        ! Get VF field
        call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
      end block create_particle
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='pipe')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('Gwall',Gib0)
         call ens_out%add_scalar('Gpart',cfg%Gib)
         call ens_out%add_scalar('VFwall',VF0)
         call ens_out%add_scalar('VFpart',cfg%VF)
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
          ! Createparticle monitor
         pfile=monitor(fs%cfg%amRoot,'particle')
         call pfile%add_column(time%n,'Timestep number')
         call pfile%add_column(time%t,'Time')
         call pfile%add_column(time%dt,'Timestep size')
         call pfile%add_column(time%cfl,'Maximum CFL')
         call pfile%add_column(pos(1),'X')
         call pfile%add_column(pos(2),'Y')
         call pfile%add_column(pos(3),'Z')
         call pfile%add_column(vel(1),'U')
         call pfile%add_column(vel(2),'V')
         call pfile%add_column(vel(3),'W')
         call pfile%add_column(force(1),'Fx')
         call pfile%add_column(force(2),'Fy')
         call pfile%add_column(force(3),'Fz')
         call pfile%write()
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
         call cflfile%add_column(CFLp,'Particle CFL')
         call cflfile%write()
      end block create_monitor
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         CFLp=abs(maxval(vel))*time%dt/cfg%min_meshsize
         time%cfl=max(time%cfl,CFLp)
         call time%adjust_dt()
         call time%increment()
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Particle solve
            update_particle: block
              use ibconfig_class, only: bigot,sharp
              integer :: i,j,k
              real(WP), dimension(3) :: pos_old,vel_old,acc
              ! Remember the old particle
              pos_old=pos; vel_old=vel
              ! Advance with Euler prediction
              !---------------------------------------------------
              call get_force()
              acc=force/mass+fs%gravity
              pos=pos_old+0.5_WP*time%dtmid*vel
              vel=vel_old+0.5_WP*time%dtmid*acc
              ! Overwrite levelset and volume fraction
              do k=cfg%kmino_,cfg%kmaxo_
                 do j=cfg%jmino_,cfg%jmaxo_
                    do i=cfg%imino_,cfg%imaxo_
                       cfg%Gib(i,j,k)=0.5_WP*dp-sqrt((cfg%xm(i)-pos(1))**2+(cfg%ym(j)-pos(2))**2+(cfg%zm(k)-pos(3))**2)
                    end do
                 end do
              end do
              ! Get normal vector
              call cfg%calculate_normal()
              ! Get VF field
              call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
              ! Correct with midpoint rule
              !---------------------------------------------------
               call get_force()
              acc=force/mass+fs%gravity
              pos=pos_old+time%dtmid*vel
              vel=vel_old+time%dtmid*acc
              ! Overwrite levelset and volume fraction
              do k=cfg%kmino_,cfg%kmaxo_
                 do j=cfg%jmino_,cfg%jmaxo_
                    do i=cfg%imino_,cfg%imaxo_
                       cfg%Gib(i,j,k)=0.5_WP*dp-sqrt((cfg%xm(i)-pos(1))**2+(cfg%ym(j)-pos(2))**2+(cfg%zm(k)-pos(3))**2)
                    end do
                 end do
              end do
              ! Get normal vector
              call cfg%calculate_normal()
              ! Get VF field
              call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
            end block update_particle

            ! Remember old velocity
            fs%Uold=fs%U
            fs%Vold=fs%V
            fs%Wold=fs%W

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)

            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Finish update
            if (use_implicit) then
               ! Form implicit residuals
               call fs%solve_implicit(time%dt,resU,resV,resW)
               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU
               fs%V=2.0_WP*fs%V-fs%Vold+resV
               fs%W=2.0_WP*fs%W-fs%Wold+resW
            else
               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho
               fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho
               fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho
            end if
            
            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
              integer :: i,j,k
              real(WP) :: VF
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        ! No-slip on tube
                        fs%U(i,j,k)=fs%U(i,j,k)*sum(fs%itpr_x(:,i,j,k)*VF0(i-1:i,j,k))
                        fs%V(i,j,k)=fs%V(i,j,k)*sum(fs%itpr_y(:,i,j,k)*VF0(i,j-1:j,k))
                        fs%W(i,j,k)=fs%W(i,j,k)*sum(fs%itpr_z(:,i,j,k)*VF0(i,j,k-1:k))
                        ! No-slip on particle
                        VF=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                        fs%U(i,j,k)=fs%U(i,j,k)*VF+(1.0_WP-VF)*vel(1)
                        VF=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
                        fs%V(i,j,k)=fs%V(i,j,k)*VF+(1.0_WP-VF)*vel(2)
                        VF=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
                        fs%W(i,j,k)=fs%W(i,j,k)*VF+(1.0_WP-VF)*vel(3)
                     end do
                  end do
               end do
               call fs%cfg%sync(fs%U)
               call fs%cfg%sync(fs%V)
               call fs%cfg%sync(fs%W)
            end block ibforcing
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)      
            
            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
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
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call pfile%write()
         call cflfile%write()

         ! End simuilation if particle leaves the domain
         if (pos(1).gt.fs%cfg%x(fs%cfg%imax+1)) time%event_done=.true.
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,gradU,Gib0,VF0)
      
   end subroutine simulation_final
   
end module simulation
