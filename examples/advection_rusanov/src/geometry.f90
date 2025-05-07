!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
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
         logical :: periodic
         integer :: i, j, k, nx, ny, nz
         real(WP) :: Lx, Ly, Lz
         real(WP), dimension(:), allocatable :: x, y, z

         ! read in grid definition
         call param_read('Lx', Lx)
         call param_read('nx', nx)
         call param_read('Ly', Ly)
         call param_read('ny', ny)
         call param_read('Lz', Lz)
         call param_read('nz', nz)
         call param_read('Periodic', periodic)

         ! allocate
         allocate(x(nx+1), y(ny+1), z(nz+1))

         ! create simple rectilinear grid
         do i = 1, nx+1
            x(i) = real(i-1,WP) / real(nx, WP) * Lx
         end do
         do j = 1, ny+1
            y(j) = real(j-1,WP) / real(ny, WP) * Ly
         end do
         do k = 1, nz+1
            z(k) = real(k-1,WP) / real(nz, WP) * Lz
         end do

         ! generate serial grid object
         grid=sgrid(coord=cartesian, no=1, x=x, y=x, z=x, xper=periodic,      &
           & yper=periodic, zper=periodic)

      end block create_grid

      ! create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition

         ! read in partition
         call param_read('Partition', partition, short='p')

         ! create partitioned grid
         cfg = config(grp=group, decomp=partition, grid=grid)

      end block create_cfg

      ! Create masks for this config
      create_walls: block
         cfg%VF=1.0_WP
      end block create_walls

   end subroutine geometry_init

end module geometry

