!> Various definitions and tools for initializing NGA2 config
module geometry
   use ibconfig_class, only: ibconfig
   use precision,      only: WP
   implicit none
   private
   
   !> Single config
   type(ibconfig), public :: cfg
   
   !> Pipe diameter
   real(WP), public :: D1,D2
   
   public :: geometry_init
   
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx,Dmax
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('L',Lx)
         call param_read('D1',D1)
         call param_read('D2',D2)
         call param_read('ny',ny); allocate(y(ny+1))
         call param_read('nx',nx); allocate(x(nx+1))
         call param_read('nz',nz); allocate(z(nz+1))

         Dmax=max(D1,D2)
         dx=Lx/real(nx,WP)
         no=4
         if (ny.gt.1) then
            Ly=Dmax+real(2*no,WP)*Dmax/real(ny-2*no,WP)
         else
            Ly=dx
         end if
         if (nz.gt.1) then
            Lz=Dmax+real(2*no,WP)*Dmax/real(ny-2*no,WP)
         else
            Lz=dx
         end if
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='tube')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=ibconfig(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Geometry: stepped cylinder, D1 < D2
      !
      !                 r
      !                 ^
      !                 |
      !        D2/2  ---+                    +------------------+
      !                 |                    |                  |
      !                 |                    |                  |
      !        D1/2  ---+--------------------+                  |
      !                 |                                       |
      !                 |                                       |
      !       -D1/2  ---+--------------------+                  |
      !                 |                    |                  |
      !                 |                    |                  |
      !       -D2/2  ---+                    +------------------+
      !                 |
      !                 +--------------------+------------------+---> x
      !                 0                    H                  L
      create_walls: block
        use ibconfig_class, only: bigot,sharp
        use messager,  only: die
        integer :: i,j,k
        real(WP) :: H,r,x,R1,R2,Rlo,Rhi
        real(WP) :: xcA,xcC,rcB,distA,distB,distC,dist,sgn
        call param_read('H',H)
        if (H.gt.cfg%xL) call die('H cannot be larger than Lx!')
        R1=0.5_WP*D1; R2=0.5_WP*D2
        Rlo=min(R1,R2); Rhi=max(R1,R2)
        do k=cfg%kmino_,cfg%kmaxo_
           do j=cfg%jmino_,cfg%jmaxo_
              do i=cfg%imino_,cfg%imaxo_
                 x=cfg%xm(i)
                 r=sqrt(cfg%ym(j)**2+cfg%zm(k)**2)
                 ! Distance to side wall below H: ray {r=R1, x<=H}
                 xcA=min(x,H)
                 distA=sqrt((x-xcA)**2+(r-R1)**2)
                 ! Distance to side wall above H: ray {r=R2, x>=H}
                 xcC=max(x,H)
                 distC=sqrt((x-xcC)**2+(r-R2)**2)
                 ! Distance to flat shoulder at x=H, connecting R1 to R2
                 rcB=min(max(r,Rlo),Rhi)
                 distB=sqrt((x-H)**2+(r-rcB)**2)

                 dist=min(distA,distB,distC)

                 ! Inside/outside test
                 if (x.lt.H) then
                    if (r.lt.R1) then
                       sgn=-1.0_WP
                    else
                       sgn=+1.0_WP
                    end if
                 else
                    if (r.lt.R2) then
                       sgn=-1.0_WP
                    else
                       sgn=+1.0_WP
                    end if
                 end if

                 cfg%Gib(i,j,k)=sgn*dist

              end do
           end do
        end do
        ! Get normal vector
        call cfg%calculate_normal()
        ! Get VF field
        call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
