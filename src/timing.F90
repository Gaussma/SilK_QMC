module timing
  use vartypes
  integer(int8) :: time_start, time_rho, time_start_evol
  integer(int8) :: time_start_evol_it, time_start_ext_diag, time_start_debug, &
               time_start_part1, time_start_part2, time_start_part3, &
               time_start_part4
  real(dp)    :: ttimer_evol, ttimer_ext_diag, ttimer_debug, &
               ttimer_evol_part1, ttimer_evol_part2, ttimer_evol_part3, ttimer_evol_part4
  real(dp)    :: timer_evol, timer_ext_diag, timer_debug, &
               timer_evol_part1,  timer_evol_part2,  timer_evol_part3,  timer_evol_part4
  contains

  subroutine checkpoint_timing()
    use checkpointing
    implicit none
    call checkpoint_real_scalar("ttimer_evol",       ttimer_evol)
    call checkpoint_real_scalar("ttimer_ext_diag",   ttimer_ext_diag)
    call checkpoint_real_scalar("ttimer_debug",      ttimer_debug)
    call checkpoint_real_scalar("ttimer_evol_part1", ttimer_evol_part1)
    call checkpoint_real_scalar("ttimer_evol_part2", ttimer_evol_part2)
    call checkpoint_real_scalar("ttimer_evol_part3", ttimer_evol_part3)
    call checkpoint_real_scalar("ttimer_evol_part4", ttimer_evol_part4)
    call checkpoint_real_scalar("timer_evol",        timer_evol)
    call checkpoint_real_scalar("timer_ext_diag",    timer_ext_diag)
    call checkpoint_real_scalar("timer_debug",       timer_debug)
    call checkpoint_real_scalar("timer_evol_part1",  timer_evol_part1)
    call checkpoint_real_scalar("timer_evol_part2",  timer_evol_part2)
    call checkpoint_real_scalar("timer_evol_part3",  timer_evol_part3)
    call checkpoint_real_scalar("timer_evol_part4",  timer_evol_part4)
  end subroutine

  subroutine restart_timing()
    use checkpointing
    implicit none
    call restart_real_scalar("ttimer_evol",       ttimer_evol)
    call restart_real_scalar("ttimer_ext_diag",   ttimer_ext_diag)
    call restart_real_scalar("ttimer_debug",      ttimer_debug)
    call restart_real_scalar("ttimer_evol_part1", ttimer_evol_part1)
    call restart_real_scalar("ttimer_evol_part2", ttimer_evol_part2)
    call restart_real_scalar("ttimer_evol_part3", ttimer_evol_part3)
    call restart_real_scalar("ttimer_evol_part4", ttimer_evol_part4)
    call restart_real_scalar("timer_evol",        timer_evol)
    call restart_real_scalar("timer_ext_diag",    timer_ext_diag)
    call restart_real_scalar("timer_debug",       timer_debug)
    call restart_real_scalar("timer_evol_part1",  timer_evol_part1)
    call restart_real_scalar("timer_evol_part2",  timer_evol_part2)
    call restart_real_scalar("timer_evol_part3",  timer_evol_part3)
    call restart_real_scalar("timer_evol_part4",  timer_evol_part4)
  end subroutine

end module timing

