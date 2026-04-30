!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class,   only: config
   use precision,      only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg

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
         integer :: i,j,k,nx,ny,nz,N
         real(WP) :: Lx,Ly,Lz,dist,dx,R
         real(WP), dimension(:), allocatable :: x,y,z
         
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         call param_read('Lz',Lz)
         call param_read('R',R)
         call param_read('N Across',N)

         ! Lx = 1.0_WP ! beam length
         ! Ly = 1.0_WP ! beam length

         ! dist = 0.01_WP ! Space between particles
         dist = Ly/N
         Lx = Lx + 6.03_WP * dist ! total length of the beam
         ! Ly = Ly + 3.0_WP * dist ! total length of the beam

         dx = 3.015_WP*dist ! grid spacing
         ! print*, "Grid Spacing : ", dx

         nx = ceiling(Lx/dx)+2 ! number of division in x 
         ny = ceiling(Ly/dx)+2 
         nz = ceiling(Lz/dx)+2 

         allocate(x(nx+1))
         allocate(y(ny+1))
         allocate(z(nz+1))

         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-2,WP)*dx - Lx/2.0_WP - 1.5_WP*dist
         end do
         do j=1,ny+1
            y(j)=real(j-2,WP)*dx - Ly/2.0_WP - 1.5_WP*dist
         end do
         do k=1,nz+1
            z(k)=real(k-2,WP)*dx - Lz/2.0_WP - 1.5_WP*dist
         end do
         ! General serial grid object (no=3 needed to support ghost/image point interpolation/extrapolation)
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='box')
         
      end block create_grid

      ! create_grid: block
      !    use sgrid_class, only: cartesian
      !    integer :: i,j,k,nx,ny,nz
      !    real(WP) :: Lx,Ly,Lz,dist
      !    real(WP), dimension(:), allocatable :: x,y,z
         
      !    ! Read in grid definition
      !    call param_read('Lx',Lx); Lx=Lx
      !    call param_read('Ly',Ly); Ly=Ly
      !    call param_read('Lz',Lz); Lz=Lz
      !    call param_read('Subdivisions',ny)
      !    dist = 3.0_WP * Ly / real(ny,WP)
      !    Lx = Lx + 3.0_WP * dist
      !    nx = ceiling(Lx / dist) + 4
      !    ny = ceiling(Ly / dist) + 2
      !    nz = ceiling(Lz / dist) + 2

      !    Lx = real(nx,WP) * dist
      !    Ly = real(ny,WP) * dist
      !    Lz = real(nz,WP) * dist
        

      !    allocate(x(nx))
      !    allocate(y(ny+1))
      !    allocate(z(nz+1))

         
      !    ! Create simple rectilinear grid
      !    do i=1,nx
      !       x(i)=real(i-2,WP)*dist
      !    end do
      !    do j=1,ny+1
      !       y(j)=real(j-1,WP)*dist-0.5_WP*Ly
      !    end do
      !    do k=1,nz+1
      !       z(k)=real(k-1,WP)*dist-0.5_WP*Lz
      !    end do
      

        
         
      !    ! General serial grid object (no=3 needed to support ghost/image point interpolation/extrapolation)
      !    grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='box')
         
      ! end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      
      ! Create walls for this config
      create_walls: block
        cfg%VF=1.0_WP
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
