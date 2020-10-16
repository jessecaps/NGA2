!> This is a typical NGA driver.
!> @version 1.0
!> @author O. Desjardins
program nga
   use random,     only: random_init
   use param,      only: param_init,param_final
   use parallel,   only: parallel_init,parallel_final
   use monitor,    only: monitor_init,monitor_final
   use geometry,   only: geometry_init
   use simulation, only: simulation_init,simulation_run
   implicit none
   
   
   ! =======================================
   ! Initialization sequence ===============
   ! =======================================
   ! Initialize parallel environment =======
   call parallel_init
   ! Initialize monitoring capabilities ====
   call monitor_init
   ! Initialize user interaction ===========
   call param_init
   ! Initialize random number generator ====
   call random_init
   ! =======================================
   
   
   
   
   ! =======================================
   ! Setup the problem grid ================
   ! =======================================
   call geometry_init
   ! =======================================
   
   
   
   
   
   
   ! =======================================
   ! Setup the solver ======================
   ! =======================================
   call simulation_init
   ! =======================================
   
   
   
   
   ! =======================================
   ! Run the solver ========================
   ! =======================================
   call simulation_run
   ! =======================================
   
   
   
   
   ! Code termination ======================
   ! Terminate solvers, data, and geometry
   !call simulation_final
   !call data_final
   !call geometry_final
   
   ! Terminate I/O, parser, timing, and monitor
   !call io_final
   
   
   
   
   
   ! =======================================
   ! Termination sequence ==================
   ! =======================================
   ! Clean up user parameters ==============
   call param_final
   ! Clean up monitoring tools =============
   call monitor_final
   ! Clean parallel exit ===================
   call parallel_final
   ! =======================================
   
   
end program nga
