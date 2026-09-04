!> Lagrangian solid solver object
!> Attempt at integrating the pdsolver_class without AMR
module lsspd_class
   use precision,      only: WP
   use string,         only: str_medium
   use config_class,   only: config
   use ddadi_class,    only: ddadi
   use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
   use pdsolver_class, only: pdsolver, PDC_IS_DEAD, PDC_BONDS, PDC_INTEGRATES, PDC_MOVES, pd_partition
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: lss, PDC_MOVES, PDC_IS_DEAD, PDC_BONDS, PDC_INTEGRATES, pd_partition
   
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   

   !> I/O chunk size to read at a time
   integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing

   
   !> Lagrangian solid solver object definition
   !> Extends the existing pdsolver_class, incorporating the coupling functions
   type, extends(pdsolver) :: lss
   
      ! This config is used for parallelization and for calculating bond/collision forces
      class(config), pointer :: cfg

      type(ddadi) :: implicit                             !< Implicit solver for filtering
            
      ! Solid volume fraction and momentum
      real(WP), dimension(:,:,:), allocatable :: VF       !< Volume fraction, cell-centered
      real(WP), dimension(:,:,:), allocatable :: VFU      !< Solid velocity, U-face
      real(WP), dimension(:,:,:), allocatable :: VFV      !< Solid velocity, V-face
      real(WP), dimension(:,:,:), allocatable :: VFW      !< Solid velocity, W-face
      
       ! CFL numbers
      real(WP) :: CFLp_x,CFLp_y,CFLp_z,CFLp_a
      
      real(WP) :: VFmax                              !< Volume fraction info
      real(WP), dimension(3) :: ibmForce             !< Total force due to IBM

      ! Filtering operation
      real(WP) :: filter_width                       !< Characteristic filter width
      real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z    !< Divergence operator
      real(WP), dimension(:,:,:,:), allocatable :: grd_x,grd_y,grd_z    !< Gradient operator

      ! Compatibility with the old non-amr version
      integer, dimension(:,:), allocatable :: icell   !< Index of cell containing the particle 
                                                      !< (this might be unnecessary or already exist somewhere, 
                                                      !<  but not in the pdsolver alone I think)

      ! Moving or not (allow flow to setup)
      real(WP) :: unfreeze_time
      
      
   contains
      procedure :: advance                           !< Step forward the particle ODEs
      procedure :: update_partmesh                   !< Update a partmesh object using current particles
      procedure :: update_VF                         !< Compute volume fraction
      procedure :: filter                            !< Apply volume filtering to field
      ! procedure :: get_cfl
   end type lss
   
   
   !> Declare lss constructor
   interface lss
      procedure constructor
   end interface lss
   
contains
   
   
   !> Default constructor for Lagrangian solid solver
   function constructor(cfg,name) result(self)
      implicit none
      type(lss) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! ! Initialize MPI derived datatype for a particle
      ! call prepare_mpi_part() ! IVM, do we need this still? I think that pdsolver handles communication...

      
      ! Allocate VF array on cfg mesh
      allocate(self%VF(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF=0.0_WP
      allocate(self%VFU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFU=0.0_WP
      allocate(self%VFV(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFV=0.0_WP
      allocate(self%VFW(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFW=0.0_WP

      ! Allocate finite volume divergence operators
      allocate(self%div_x(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      allocate(self%div_y(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      allocate(self%div_z(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               self%div_x(:,i,j,k)=self%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< Divergence from [x ,ym,zm]
               self%div_y(:,i,j,k)=self%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,y ,zm]
               self%div_z(:,i,j,k)=self%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,ym,z ]
            end do
         end do
      end do

      ! Allocate finite difference velocity gradient operators
      allocate(self%grd_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
      allocate(self%grd_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
      allocate(self%grd_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
      ! Create gradient coefficients to cell faces
      do k=self%cfg%kmin_,self%cfg%kmax_+1
         do j=self%cfg%jmin_,self%cfg%jmax_+1
            do i=self%cfg%imin_,self%cfg%imax_+1
               self%grd_x(:,i,j,k)=self%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< Gradient in x from [xm,ym,zm] to [x,ym,zm]
               self%grd_y(:,i,j,k)=self%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< Gradient in y from [xm,ym,zm] to [xm,y,zm]
               self%grd_z(:,i,j,k)=self%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< Gradient in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do

      ! Loop over the domain and zero divergence in walls
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) then
                  self%div_x(:,i,j,k)=0.0_WP
                  self%div_y(:,i,j,k)=0.0_WP
                  self%div_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Zero out gradient to wall faces
      do k=self%cfg%kmin_,self%cfg%kmax_+1
         do j=self%cfg%jmin_,self%cfg%jmax_+1
            do i=self%cfg%imin_,self%cfg%imax_+1
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i-1,j,k).eq.0.0_WP) self%grd_x(:,i,j,k)=0.0_WP
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j-1,k).eq.0.0_WP) self%grd_y(:,i,j,k)=0.0_WP
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j,k-1).eq.0.0_WP) self%grd_z(:,i,j,k)=0.0_WP
            end do
         end do
      end do

      ! Adjust metrics to account for lower dimensionality
      if (self%cfg%nx.eq.1) then
         self%div_x=0.0_WP
         self%grd_x=0.0_WP
      end if
      if (self%cfg%ny.eq.1) then
         self%div_y=0.0_WP
         self%grd_y=0.0_WP
      end if
      if (self%cfg%nz.eq.1) then
         self%div_z=0.0_WP
         self%grd_z=0.0_WP
      end if

      ! Create implicit solver object for filtering
      self%implicit=ddadi(cfg=self%cfg,name='Filter',nst=7)
      self%implicit%stc(1,:)=[ 0, 0, 0]
      self%implicit%stc(2,:)=[+1, 0, 0]
      self%implicit%stc(3,:)=[-1, 0, 0]
      self%implicit%stc(4,:)=[ 0,+1, 0]
      self%implicit%stc(5,:)=[ 0,-1, 0]
      self%implicit%stc(6,:)=[ 0, 0,+1]
      self%implicit%stc(7,:)=[ 0, 0,-1]
      call self%implicit%init()

      ! Set default filter width
      self%filter_width=1.0_WP*self%cfg%min_meshsize

      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (self%cfg%amRoot) then
            write(message,'("LSS object [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end function constructor

   

   !> Advance the particle equations by a specified time step dt
   subroutine advance(this,dt,unfreeze,div_stress_x,div_stress_y,div_stress_z)
      implicit none
      class(lss), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div_stress_x  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div_stress_y  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div_stress_z  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n,i
      logical, intent(in) :: unfreeze
      
      
      
      do i=1,this%nown

         if (this%flag(i).eq.PDC_IS_DEAD) cycle
            ! this%ff(:,i)=this%cfg%get_velocity( &  ! we do (div_stress)/rho later to make it acc inside of the pd_advnace routine
            !    pos=this%y(:,i), &
            !    i0=this%icell(1,i), &
            !    j0=this%icell(2,i), &
            !    k0=this%icell(3,i), &
            !    U=div_stress_x,V=div_stress_y,W=div_stress_z) ! interpolates the divergence of stress to the location of the particle (this term is a force density now)
            
            this%ff(:,i)=0.0_WP! we do (div_stress)/rho later to make it acc inside of the pd_advnace routine
            
      end do
      call this%pd_advance(dt) ! use fluid forces and compute bond forces, and update position due to verlet scheme
      
    
      ! do i=1,this%nown
      !    if (this%flag(i).eq.PDC_IS_DEAD) cycle

      !    this%icell(:,i)=this%cfg%get_ijk_global(this%y(:,i),this%icell(:,i)) ! do we need to do this?

      !    ! if(unfreeze) this%flag(i) = PDC_BONDS + PDC_INTEGRATES + PDC_MOVES
      ! end do

      if (unfreeze) then
         do i = 1,this%nown
            this%damping_rate = 0.0005_WP
            if (this%flag(i).eq.(PDC_MOVES+PDC_BONDS)) this%v(:,i)= 0.0_WP
         end do
      end if


      
      ! call this%update_VF() ! now we update the volume fraction

      
      ! Log/screen output (do we need to do this still?)
      ! logging: block
      !    use, intrinsic :: iso_fortran_env, only: output_unit
      !    use param,    only: verbose
      !    use messager, only: log
      !    use string,   only: str_long
      !    character(len=str_long) :: message
      !    if (this%cfg%amRoot) then
      !       write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
      !       if (verbose.gt.1) write(output_unit,'(a)') trim(message)
      !       if (verbose.gt.0) call log(message)
      !    end if
      ! end block logging
      
   end subroutine advance


   !> Update particle volume fraction using our current particles
   subroutine update_VF(this)
      implicit none
      class(lss), intent(inout) :: this
      integer :: i
      ! Reset volume fraction and momentum
      this%VF=0.0_WP; this%VFU=0.0_WP; this%VFV=0.0_WP; this%VFW=0.0_WP
      ! Transfer particle volume
      do i=1,this%nown ! halo included here? 
         ! Skip inactive particle
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         ! Transfer volume to mesh
         call this%cfg%set_scalar(Sp=this%vol(i),           pos=this%y(:,i),i0=this%icell(1,i),j0=this%icell(2,i),k0=this%icell(3,i),S=this%VF ,bc='n')
         call this%cfg%set_scalar(Sp=this%vol(i)*this%v(1,i), pos=this%y(:,i),i0=this%icell(1,i),j0=this%icell(2,i),k0=this%icell(3,i),S=this%VFU,bc='n')
         call this%cfg%set_scalar(Sp=this%vol(i)*this%v(2,i), pos=this%y(:,i),i0=this%icell(1,i),j0=this%icell(2,i),k0=this%icell(3,i),S=this%VFV,bc='n')
         call this%cfg%set_scalar(Sp=this%vol(i)*this%v(3,i), pos=this%y(:,i),i0=this%icell(1,i),j0=this%icell(2,i),k0=this%icell(3,i),S=this%VFW,bc='n')
      end do
      this%VF =this%VF /this%cfg%vol
      this%VFU=this%VFU/this%cfg%vol
      this%VFV=this%VFV/this%cfg%vol
      this%VFW=this%VFW/this%cfg%vol
      ! Sum at boundaries
      call this%cfg%syncsum(this%VF )
      call this%cfg%syncsum(this%VFU)
      call this%cfg%syncsum(this%VFV)
      call this%cfg%syncsum(this%VFW)
      ! Apply volume filter
      call this%filter(this%VF )
      call this%filter(this%VFU)
      call this%filter(this%VFV)
      call this%filter(this%VFW)
      ! Clip
      where (this%VF.lt.0.0_WP) this%VF=0.0_WP
      this%VF=min(this%VF,1.0_WP-epsilon(1.0_WP))

    end subroutine update_VF

   !  subroutine get_cfl(this,dt,cfl)
   !    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
   !    use parallel, only: MPI_REAL_WP
   !    implicit none
   !    class(lss), intent(inout) :: this
   !    real(WP), intent(in)  :: dt
   !    real(WP), intent(out) :: cfl
   !    integer :: i,ierr
   !    real(WP) :: my_CFLp_x,my_CFLp_y,my_CFLp_z,kk,mu,a
      
   !    ! Set the CFLs to zero
   !    my_CFLp_x=0.0_WP; my_CFLp_y=0.0_WP; my_CFLp_z=0.0_WP
   !    do i=1,this%nown
   !       my_CFLp_x=max(my_CFLp_x,abs(this%v(1,i))*this%cfg%dxi(this%icell(1,i)))
   !       my_CFLp_y=max(my_CFLp_y,abs(this%v(2,i))*this%cfg%dyi(this%icell(2,i)))
   !       my_CFLp_z=max(my_CFLp_z,abs(this%v(3,i))*this%cfg%dzi(this%icell(3,i)))
   !    end do
   !    my_CFLp_x=my_CFLp_x*dt; my_CFLp_y=my_CFLp_y*dt; my_CFLp_z=my_CFLp_z*dt
      
   !    ! Get the parallel max
   !    call MPI_ALLREDUCE(my_CFLp_x,this%CFLp_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(my_CFLp_y,this%CFLp_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(my_CFLp_z,this%CFLp_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)

   !    ! CFL based on elastic wave speed in material
   !    kk=this%elastic_modulus/(3.0_WP-6.0_WP*this%poisson_ratio)
   !    mu=this%elastic_modulus/(2.0_WP+2.0_WP*this%poisson_ratio)      
   !    a=sqrt((kk+4.0_WP*mu/3.0_WP)/this%rho)
   !    this%CFLp_a=dt*a/this%delta
      
   !    ! Return the maximum CFL
   !    cfl=max(this%CFLp_x,this%CFLp_y,this%CFLp_z,this%CFLp_a)
      
   ! end subroutine get_cfl
    

    !> Laplacian filtering operation
    subroutine filter(this,A)
      implicit none
      class(lss), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: filter_coeff
      integer :: i,j,k,n,nstep
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ

      ! Recompute filter coefficient
      filter_coeff=max(this%filter_width**2-this%cfg%min_meshsize**2,0.0_WP)/(16.0_WP*log(2.0_WP))
      if (filter_coeff.le.0.0_WP) return

      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      if (.not.this%implicit%setup_done) then
         ! Prepare diffusive operator (only need to do this once)
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%implicit%opr(1,i,j,k)=1.0_WP-(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x(-1,i+1,j,k)+&
                  &                                  this%div_x( 0,i,j,k)*filter_coeff*this%grd_x( 0,i  ,j,k)+&
                  &                                  this%div_y(+1,i,j,k)*filter_coeff*this%grd_y(-1,i,j+1,k)+&
                  &                                  this%div_y( 0,i,j,k)*filter_coeff*this%grd_y( 0,i,j  ,k)+&
                  &                                  this%div_z(+1,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k+1)+&
                  &                                  this%div_z( 0,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k  ))
                  this%implicit%opr(2,i,j,k)=      -(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x( 0,i+1,j,k))
                  this%implicit%opr(3,i,j,k)=      -(this%div_x( 0,i,j,k)*filter_coeff*this%grd_x(-1,i  ,j,k))
                  this%implicit%opr(4,i,j,k)=      -(this%div_y(+1,i,j,k)*filter_coeff*this%grd_y( 0,i,j+1,k))
                  this%implicit%opr(5,i,j,k)=      -(this%div_y( 0,i,j,k)*filter_coeff*this%grd_y(-1,i,j  ,k))
                  this%implicit%opr(6,i,j,k)=      -(this%div_z(+1,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k+1))
                  this%implicit%opr(7,i,j,k)=      -(this%div_z( 0,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k  ))
               end do
            end do
         end do
      end if
      ! Explicit step
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FX(i,j,k)=filter_coeff*sum(this%grd_x(:,i,j,k)*A(i-1:i,j,k))
               FY(i,j,k)=filter_coeff*sum(this%grd_y(:,i,j,k)*A(i,j-1:j,k))
               FZ(i,j,k)=filter_coeff*sum(this%grd_z(:,i,j,k)*A(i,j,k-1:k))
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%rhs(i,j,k)=sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Implicit step
      call this%implicit%setup()
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      A=A+this%implicit%sol
      call this%cfg%sync(A)

      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)

    end subroutine filter
   
   !> Update particle mesh using our current particles
   subroutine update_partmesh(this,pmesh)
      use partmesh_class, only: partmesh
      implicit none
      class(lss), intent(inout) :: this
      class(partmesh), intent(inout) :: pmesh
      integer :: i
      ! Reset particle mesh storage
      call pmesh%reset()
      ! Nothing else to do if no particle is present
      if (this%nown.eq.0) return
      ! Copy particle info
      call pmesh%set_size(this%nown)
      do i=1,this%nown !< IVM, this might not be good, I think we will get duplicates this way
         pmesh%pos(:,i)=this%y(:,i)
      end do
   end subroutine update_partmesh
   
   
   ! !> Creation of the MPI datatype for particle ! IVM, Maybe we dont need this, since comm is handled by pdsolver?
   ! subroutine prepare_mpi_part()
   !    use mpi_f08
   !    use messager, only: die
   !    implicit none
   !    integer(MPI_ADDRESS_KIND), dimension(part_nblock) :: disp
   !    integer(MPI_ADDRESS_KIND) :: lb,extent
   !    type(MPI_Datatype) :: MPI_PART_TMP
   !    integer :: i,mysize,ierr
   !    ! Prepare the displacement array
   !    disp(1)=0
   !    do i=2,part_nblock
   !       call MPI_Type_size(part_tblock(i-1),mysize,ierr)
   !       disp(i)=disp(i-1)+int(mysize,MPI_ADDRESS_KIND)*int(part_lblock(i-1),MPI_ADDRESS_KIND)
   !    end do
   !    ! Create and commit the new type
   !    call MPI_Type_create_struct(part_nblock,part_lblock,disp,part_tblock,MPI_PART_TMP,ierr)
   !    call MPI_Type_get_extent(MPI_PART_TMP,lb,extent,ierr)
   !    call MPI_Type_create_resized(MPI_PART_TMP,lb,extent,MPI_PART,ierr)
   !    call MPI_Type_commit(MPI_PART,ierr)
   !    ! If a problem was encountered, say it
   !    if (ierr.ne.0) call die('[lss prepare_mpi_part] MPI Particle type creation failed')
   !    ! Get the size of this type
   !    call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   ! end subroutine prepare_mpi_part
   
   
end module lsspd_class
