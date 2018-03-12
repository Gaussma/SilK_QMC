#ifndef KINK_GSL
#define SPRNG
#define SIMPLE_SPRNG
#else
#define SINGLE_STREAM
#endif

module mod_random
  use vartypes

#ifdef SPRNG
#include "sprng_f.h"
  SPRNG_POINTER sprng_p
#else
#define RAND
#endif

#ifndef SPRNG
interface
  function rand_d(i)
    use vartypes
    real(dp) :: rand_d
    integer, intent(in) :: i
  end function
end interface
#endif

  contains

  subroutine init_random(seed_index)
    implicit none
    integer :: seed_index
#ifdef RAND
    integer :: ierr
#ifdef SINGLE_STREAM
    call rand_init(1, ierr)
#else
    call rand_init(1, ierr) ! TODO: use MPI ranks here (ranks+1)
#endif
#else
    integer:: root_seed, seed, k1

    root_seed = 293847589
    sprng_p = init_sprng(1,root_seed,SPRNG_DEFAULT)

    seed=isprng()
    do k1=1,seed_index
      seed=isprng()
    end do
    sprng_p = init_sprng(1,seed,SPRNG_DEFAULT)
#endif
  end subroutine init_random

  function random()
    implicit none
    real(dp) :: random
#ifdef RAND
#ifdef SINGLE_STREAM
    random = rand_d(0)
#else
    random = rand_d(0) ! TODO: use MPI rank here
#endif
#else
    random = sprng()
#endif
  end function

  subroutine checkpoint_random()
    use checkpointing 
    implicit none
#ifdef SPRNG
    character :: buffer(MAX_PACKED_LENGTH)
    integer*8 :: size

    size = pack_sprng(buffer)
    call checkpoint_int_scalar("sprng_size", size)
    call checkpoint_char_array("sprng_state", buffer(1:size), size)
#endif
  end subroutine

  subroutine restart_random()
    use checkpointing 
    implicit none
#ifdef SPRNG
    character :: buffer(MAX_PACKED_LENGTH)
    integer*8 :: size

    call restart_int_scalar("sprng_size", size)
    call restart_char_array("sprng_state", buffer(1:size), size)
#endif
  end subroutine

end module
