!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use vfs_class,         only: vfs
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Only get a VOF solver and corresponding time tracker
   type(vfs),         public :: vf
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: U,V,W
   
contains
   
   
   !> Function that defines a level set function for Zalesak's notched circle problem
   function levelset_zalesak(xyz,center,radius,width,height) result(G)
      implicit none
      real(WP) :: G
      real(WP), dimension(3), intent(in) :: xyz,center
      real(WP), intent(in) :: radius,width,height
      real(WP) :: c,b,b1,b2,h1,h2
      c=radius-sqrt(sum((xyz-center)**2))
      b1=center(1)-0.5_WP*width; b2=center(1)+0.5_WP*width
      h1=center(2)-radius*cos(asin(0.5_WP*width/radius)); h2=center(2)-radius+height
      if     (c.ge.0.0_WP.and.xyz(1).le.b1.and.xyz(2).le.h2) then
         b=b1-xyz(1); G=min(c,b)
      elseif (c.ge.0.0_WP.and.xyz(1).ge.b2.and.xyz(2).le.h2) then
         b=xyz(1)-b2; G=min(c,b)
      elseif (c.ge.0.0_WP.and.xyz(1).ge.b1.and.xyz(1).le.b2.and.xyz(2).ge.h2) then
         b=xyz(2)-h2; G=min(c,b)
      elseif (c.ge.0.0_WP.and.xyz(1).le.b1.and.xyz(2).ge.h2) then
         b=sqrt(sum((xyz-(/b1,h2,0.0_WP/))**2)); G=min(c,b)
      elseif (c.ge.0.0_WP.and.xyz(1).ge.b2.and.xyz(2).ge.h2) then
         b=sqrt(sum((xyz-(/b2,h2,0.0_WP/))**2)); G=min(c,b)
      elseif (xyz(1).ge.b1.and.xyz(1).le.b2.and.xyz(2).le.h2.and.xyz(2).ge.h1) then
         G=-min(abs(xyz(1)-b1),abs(xyz(1)-b2),abs(xyz(2)-h2))
      elseif (xyz(1).ge.b1.and.xyz(1).le.b2.and.xyz(2).le.h1) then
         G=-min(sqrt(sum((xyz-(/b1,h1,0.0_WP/))**2)),sqrt(sum((xyz-(/b2,h1,0.0_WP/))**2)))
      else
         G=c
      end if
   end function levelset_zalesak
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Initialize our VOF
      initialize_vof: block
         integer :: i,j,k,ii,jj
         real(WP), dimension(3) :: xyz
         integer, parameter :: nsub=20
         ! Create a VOF solver
         vf=vfs(cfg=cfg,name='VOF')
         ! Initialize to Zalesak disk
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  vf%VF(i,j,k)=0.0_WP
                  do jj=1,nsub
                     do ii=1,nsub
                        xyz=[vf%cfg%x(i)+real(ii,WP)*vf%cfg%dx(i)/real(nsub+1,WP),vf%cfg%y(j)+real(jj,WP)*vf%cfg%dy(j)/real(nsub+1,WP),0.0_WP]
                        if (levelset_zalesak(xyz=xyz,center=[0.0_WP,0.25_WP,0.0_WP],radius=0.15_WP,width=0.05_WP,height=0.25_WP).gt.0.0_WP) vf%VF(i,j,k)=vf%VF(i,j,k)+1.0e-4_WP
                     end do
                  end do
               end do
            end do
         end do
      end block initialize_vof
      
      
      ! Prepare our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         integer :: i,j,k
         ! Allocate arrays
         allocate(U(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_))
         allocate(V(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_))
         allocate(W(vf%cfg%imino_:vf%cfg%imaxo_,vf%cfg%jmino_:vf%cfg%jmaxo_,vf%cfg%kmino_:vf%cfg%kmaxo_))
         ! Initialize to solid body rotation
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  U(i,j,k)=-twoPi*vf%cfg%ym(j)
                  V(i,j,k)=+twoPi*vf%cfg%xm(i)
                  W(i,j,k)=0.0_WP
               end do
            end do
         end do
      end block initialize_velocity
      
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=vf%cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         time%dt=time%dtmax
         time%itmax=1
      end block initialize_timetracker
      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='zalesak')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('VOF',vf%VF)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call vf%get_max()
         ! Create simulation monitor
         mfile=monitor(vf%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform explicit Euler time integration
      do while (.not.time%done())
         
         ! Increment time
         call time%increment()
         
         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! ================ VOF SOLVER =======================
            
            ! Time advance
            call vf%advance(dt=time%dt,U=U,V=V,W=W)
            
            ! Reconstruct interface
            
            
            ! ===================================================
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call vf%get_max()
         call mfile%write()
         
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
      deallocate(U,V,W)
      
   end subroutine simulation_final
   
   
   
   
   
end module simulation