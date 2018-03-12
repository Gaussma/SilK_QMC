! We need sparse matrices for different types, but we don't want to copy
! the implementation. Thus, we here use #define to declare a type, let this
! file include itself, now using this type, and doing this as many times as
! we have types.
! Also, note that this type does _not_ handle the question which elements
! should be stored (e.g., the ones that are not zero). It only handles storing
! elements of a matrix in a way that is efficient for sparse matrixes

#ifndef MOD_SPARSE_MATRIX
#define MOD_SPARSE_MATRIX

#define CACHESIZE 100
!#define CDEBUG

#define  TYPE real(dp)
#ifdef TRAD_CPP
#define STYPE(X) X/**/_dp
#else
#define STYPE(X) X ## _dp
#endif
#include "mod_sparse_matrix.F90"

#define  TYPE integer
#ifdef TRAD_CPP
#define STYPE(X) X/**/_int
#else
#define STYPE(X) X ## _int
#endif
#include "mod_sparse_matrix.F90"

#define  TYPE type(ham_element)
#ifdef TRAD_CPP
#define STYPE(X) X/**/_ham
#else
#define STYPE(X) X ## _ham
#endif
#include "mod_sparse_matrix.F90"

#else

module STYPE(mod_sparse_matrix)
  use vartypes

  ! type of one matrix element
  type STYPE(mat_element)
    integer :: i ! index within a row/column
    TYPE    :: e ! element
  end type

  ! Blame Fortran: the only way to get an array of pointers is by using a type
  type pointer_array
    type(STYPE(mat_element)), dimension(:), pointer :: p
  end type pointer_array

  ! type for the whole matrix
  type STYPE(sparse_matrix)
    integer                          :: vlen           ! size of v
    type(pointer_array), allocatable :: v(:)           ! Array of pointers (to arrays)
    integer,             allocatable :: v_size(:)      ! valid portion of each array
    integer,             allocatable :: v_allocated(:) ! allocated size of each array
    integer,             allocatable :: cache_idx(:)   ! cache indexes (cols)
    integer,             allocatable :: cache(:,:)     ! cache column indexes of cols
    integer(int8)                    :: cache_hits
    integer(int8)                    :: cache_misses

    contains
    procedure, pass :: STYPE(init)
    procedure, pass :: STYPE(grow)
    procedure, pass :: STYPE(set_elem)
    procedure, pass :: STYPE(del_elem)
    procedure, pass :: STYPE(get_elem)
    procedure, pass :: STYPE(find_index)
    procedure, pass :: STYPE(update_cache)
    procedure, pass :: STYPE(output_info)
    procedure, pass :: STYPE(get_col_arrays)
    ! only for debugging
#ifdef CDEBUG
    procedure, pass :: STYPE(test_cache)
    procedure, pass :: STYPE(test_sorted)
#endif
  end type STYPE(sparse_matrix)

  contains

  subroutine STYPE(init)(this, llen)
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,                     intent(in)    :: llen
    this%vlen = llen
    allocate(this%v          (llen))
    allocate(this%v_size     (llen))
    allocate(this%v_allocated(llen))
    this%v_size(:)      = 0
    this%v_allocated(:) = 0

    allocate(this%cache_idx(CACHESIZE))
    allocate(this%cache(CACHESIZE, llen))
    this%cache_idx = 0
    this%cache_hits = 0
    this%cache_misses = 0
  end subroutine STYPE(init)

  ! update the caches elements to be from within this list (only these)
  subroutine STYPE(update_cache)(this, ipass, list_len, list)
    use io
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,              intent(in) :: ipass
    integer,              intent(in) :: list_len
    integer, allocatable, intent(in) :: list(:)

    integer :: ci, li, hi, i
    logical :: found, new

    if (ipass == 0) then
      open (unit=8675, file=trim(adjustl(outdir))//'/'//"cache.dat" ,status="unknown")
    else
      open (unit=8675, file=trim(adjustl(outdir))//'/'//"cache.dat" ,status="old",access='append')
    endif
    write(*,*) "CACHE STATS", this%cache_hits, this%cache_misses, &
               real(this%cache_hits,dp)/(this%cache_misses+this%cache_hits)
    write(8675,*) ipass, this%cache_hits, this%cache_misses, &
                  real(this%cache_hits,dp)/(this%cache_misses+this%cache_hits)
    close(8675)
    this%cache_hits   = 0
    this%cache_misses = 0

    ! Mark newly unused cache entries
    do ci = 1, CACHESIZE
      if (this%cache_idx(ci) > 0) then
        found = .false.
        do li = 1, list_len
          if (this%cache_idx(ci) == list(li)) found = .true.
        end do
        if (.not. found) this%cache_idx(ci) = 0
      end if
    end do
#ifdef CDEBUG
    write(*,*) "before update"
    call this%STYPE(test_cache)()
    write(*,*) "ok"
#endif
    ! now make sure that columns that are new are copied into the cache
    ! first find which indexes are new
    do li = 1, list_len
      new = .true.
      ! TODO: use 'while' here
      do ci = 1, CACHESIZE
        if (this%cache_idx(ci) == list(li)) then
          new = .false.
        endif
      end do
      ! now, if this is a new entry, fill cache (otherwise skip)
      if (new) then
        ! for that, we first need to find an empty entry in the cache
        ci = 1
        do while (this%cache_idx(ci) > 0)
          ci = ci+1
          if (ci > CACHESIZE) then
            write(*,*) "Cache size exceeded."
            write(*,*) "list_len", list_len
            stop
          end if
        end do
        hi = 1
        do i = 1, size(this%cache(ci,:))
          if (hi <= this%v_size(list(li))) then
            if (this%v(list(li))%p(hi)%i == i) then
              this%cache(ci, i) = hi
              hi = hi + 1
            else
              this%cache(ci, i) = -hi
            end if
          else
            this%cache(ci, i) = -hi
          end if
        end do
        this%cache_idx(ci) = list(li)
      end if
    end do
#ifdef CDEBUG
    write(*,*) "after update"
    call this%STYPE(test_cache)()
    write(*,*) "ok"
#endif
  end subroutine

  ! output some info about this matrix to a file
  subroutine STYPE(output_info)(this, ipass)
    use io
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,              intent(in) :: ipass
    integer :: i, sum_nonzero
    sum_nonzero = 0
    do i = 1, this%vlen
      sum_nonzero = sum_nonzero + this%v_size(i)
    end do
    if (ipass == 0) then
      open (unit=8676, file=trim(adjustl(outdir))//'/'//"cache_info.dat" ,&
            status="unknown")
      write(8676,*) "# 0:pass 1: sum nonzero 2: relative nonzero"
    else
      open (unit=8676, file=trim(adjustl(outdir))//'/'//"cache_info.dat" ,&
            status="old",access='append')
    endif
    write(8676,*) ipass, sum_nonzero, real(sum_nonzero) / this%vlen
    close(8676)
  end subroutine

#ifdef CDEBUG
  subroutine STYPE(test_cache)(this)
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer :: c, ci, i, ii, hi, val
    return
    i=0
    do ci = 1, CACHESIZE
      c = this%cache_idx(ci)
      if (c > 0) then
        hi = 1
        do i = 1, size(this%cache(ci,:))
          if (hi <= this%v_size(c)) then
            if (this%v(c)%p(hi)%i == i) then
              val = hi
              hi = hi + 1
            else
              val = -hi
            end if
          else
            val = -hi
          end if
          if (this%cache(ci, i) /= val) then
            write(*,*) "Error in cache detected"
            write(*,*) "  c, i", c, i
            write(*,*) "  this%cache(ci, i), val", this%cache(ci, i), val
            hi = 1
            do ii = 1, i
              if (this%v(c)%p(hi)%i == ii) then
                write(*,*) ii, this%cache(ci, ii), hi
                hi = hi + 1
              else
                write(*,*) ii, this%cache(ci, ii)
              end if
            end do
            stop
          end if
        end do
      end if
    end do
  end subroutine
#endif

  ! TODO: should try Ridders' method
  ! 'bisect' this%v(c)%p, and assume roughly linear distribution (assume it to be sorted)
  ! returns index I within it so that this%v(c)%p(I) == i.
  ! If no such element could be found, return the negated index into which
  ! a new element would need to be inserted for the resulting array to remain sorted.
  integer function STYPE(find_index)(this, c, i)
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,intent(in) :: c, i
    integer(int8)      :: k,i1,i2,it,ci
    logical            :: cached

#ifdef CDEBUG
    call this%STYPE(test_cache)()
#endif
    cached = .false.
    STYPE(find_index) = 0
    do ci = 1, CACHESIZE
      if (this%cache_idx(ci) == c) then
        this%cache_hits = this%cache_hits + 1
        STYPE(find_index) = this%cache(ci, i)
        cached = .true.
        return
        exit
      end if
    end do
    this%cache_misses = this%cache_misses + 1
    it = 0
    k  = 1
    i1 = 1
    i2 = this%v_size(c)
    ! See if we have storage at all, and if not return early
    if (i2 == 0) then
      if (cached .and. STYPE(find_index) /= -1) then
        write(*,*) "Error in cache 1"
        stop
      endif
      STYPE(find_index) = -1
      return
    end if
    ! catch cases where the needed index is on the boundary
    if (this%v(c)%p(i1)%i == i)  k = i1
    if (this%v(c)%p(i2)%i == i)  k = i2
    ! catch cases where the needed index is outside the possible range
    if (this%v(c)%p(i1)%i > i) then
      if (cached .and. STYPE(find_index) /= -1) then
        write(*,*) "Error in cache 2"
        stop
      endif
      STYPE(find_index) = -1
      return
    end if
    if (this%v(c)%p(i2)%i < i) then
      if (cached .and. STYPE(find_index) /= -(this%v_size(c)+1)) then
        write(*,*) "Error in cache 3", STYPE(find_index), -(this%v_size(c)+1)
        write(*,*) "Error in cache 3", this%v(c)%p(i2)%i, i
        write(*,*) "BB", this%cache(ci, i), this%v(c)%p(this%cache(ci, i))%i
        stop
      endif
      STYPE(find_index) = -(this%v_size(c)+1)
      return
    end if
    do while (i1 < i2 .and. this%v(c)%p(k)%i .ne. i)
      ! Compute an index in between i1 and i2 depending on this%v(c)%p(i?)%i and i
      ! that is hopefully close to the index we are looking for
      ! (assume linear distribution instead of simple bisection) 
      k = (i1*(this%v(c)%p(i2)%i-i) + &
           i2*(i-this%v(c)%p(i1)%i)) / (this%v(c)%p(i2)%i-this%v(c)%p(i1)%i)
      if (k < i1 .or. k > i2) then
        write(*,*) "error!", i1, i2, k
        write(*,*) "As", this%v(c)%p(i1)%i, this%v(c)%p(i2)%i
        write(*,*) "Looking for ",i
        stop
      end if
      ! make sure we don't end up exactly on a1 or a2 due to round-off
      ! (or otherwise we loop forever)
      if (k == i1) k = i1 + 1
      if (k == i2) k = i2 - 1
      ! If the function value at the new index is smaller than the one we are
      ! looking for, make the new index the new left boundary
      if (this%v(c)%p(k)%i < i) then
        i1 = k
      else
        ! Otherwise, make the new index the new right boundary (if it happens to
        ! be the index we are looking for, the loop will break next iteration
        i2 = k
      end if
      ! Exit loop if interval too small
      if ( i1==i2-1 ) then
        if (i < this%v(c)%p(i1)%i) then
          k = i1
        else
          k = i2
        endif
        exit
      endif
      it = it + 1
      ! end if we found it, or we ended up with an empty interval (didn't find it)
    end do
    if (this%v(c)%p(k)%i == i) then
      STYPE(find_index) =  int(k)
    else
      STYPE(find_index) = -int(k)
    end if
  end function STYPE(find_index)

  ! increase storage for a specific column
  subroutine STYPE(grow)(this, i)
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,                     intent(in)    :: i

    type(STYPE(mat_element)), allocatable :: tmp(:)
    ! if still unallocated - allocate and fill with new element
    if (this%v_allocated(i) == 0) then
      this%v_allocated(i) = 8 ! start with some (8) elements
      allocate(this%v(i)%p(this%v_allocated(i)))
      this%v_size(i) = 0
      return
    end if
    ! allocate temporary storage to hold data
    allocate(tmp(this%v_size(i)))
    ! copy data to temporary storage
    tmp = this%v(i)%p(:this%v_size(i))
    ! increase storage (deallocate, reallocate)
    deallocate(this%v(i)%p)
    this%v_allocated(i) = this%v_allocated(i) * 2
    allocate(this%v(i)%p(this%v_allocated(i)))
    ! copy data back from temporary storage
    this%v(i)%p(:this%v_size(i)) = tmp
    ! deallocate temporary storage
    deallocate(tmp)
  end subroutine STYPE(grow)

#ifdef CDEBUG
  ! a debug routine
  subroutine STYPE(test_sorted)(this, vi, c)
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,                     intent(in)    :: vi
    character,                   intent(in)    :: c
    integer :: i, j, old_i
    return
    old_i = -1
    do i=1, this%v_size(vi)
      if (this%v(vi)%p(i)%i <= old_i) then
        write(*,*) "unsorted!", c, this%v(vi)%p(1:this%v_size(i))
        do j=1, this%v_size(vi)
          write(*,*) j, this%v(vi)%p(j)%i
        end do
        stop
      end if
      old_i = this%v(vi)%p(i)%i
    end do
  end subroutine
#endif

  subroutine STYPE(del_elem)(this, i, j)
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,                     intent(in)    :: i, j
    integer                                    :: this_i, ci, li

#ifdef CDEBUG
    write(*,*) "before del"
    call this%STYPE(test_cache)()
    write(*,*) "ok"
#endif
    ! check that we have storage and at least one element.
    ! If not, return (nothing to delete)
    if (this%v_allocated(i) == 0 .or. this%v_size(i) == 0) then
      return
    endif
    this_i = this%STYPE(find_index)(i, j)
    ! return if the to-be-deleted index does not exist (nothing to delete)
    if (this_i <= 0) then
      return
    endif
#ifdef CDEBUG
    call this%STYPE(test_sorted)(i, 'd')
#endif
    ! TODO: we might need to optimize this is it turns out that the compiler
    ! generates a temporary array on the stack for this overlapping array
    ! assignment
    this%v_size(i) = this%v_size(i) - 1
    ! Handle the special cases of removing the last element (simply return)
    if (this_i == this%v_size(i) + 1) then
#ifdef CDEBUG
      call this%STYPE(test_sorted)(i, 'D')
#endif
      do ci = 1, CACHESIZE
        if (this%cache_idx(ci) == i) then
          this%cache(ci, j:)    = -this%cache(ci, j)
        end if
      end do
#ifdef CDEBUG
      write(*,*) "after early del"
      call this%STYPE(test_cache)()
      write(*,*) "ok"
#endif
      return
    endif
    this%v(i)%p (this_i : this%v_size(i)) = this%v(i)%p (this_i + 1 : this%v_size(i) + 1)
#ifdef CDEBUG
    call this%STYPE(test_sorted)(i, 'D')
#endif
    do ci = 1, CACHESIZE
      if (this%cache_idx(ci) == i) then
        this%cache(ci, j)    = -this%cache(ci, j)
        do li = j+1, size(this%cache(ci, :))
          this%cache(ci, li) = sign(abs(this%cache(ci, li)) - 1, this%cache(ci, li))
        end do
      end if
    end do
#ifdef CDEBUG
    write(*,*) "after del"
    call this%STYPE(test_cache)()
    write(*,*) "ok"
#endif
  end subroutine STYPE(del_elem)

  ! set a non-zero element in the i-th column, index j
  ! if 'j' negative: shortcut to skip index search. In this case the element is
  !                  known to exist, and '-j' is it's index. No error checking
  !                  in that case. If these assumptions don't hold, bad things
  !                  will happen.
  subroutine STYPE(set_elem)(this, i, j, val)
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,                     intent(in)    :: i, j
    TYPE,                        intent(in)    :: val
    integer                                    :: this_i, ci, hi

#ifdef CDEBUG
    write(*,*) "before set"
    call this%STYPE(test_cache)()
    write(*,*) "ok"
#endif
    ! check that we have storage. If not, allocate.
    if (this%v_allocated(i) == 0) then
      call this%STYPE(grow)(i)
    endif
    ! check that we have at least one element, if not, simply assign and return
    if (this%v_size(i) == 0) then
      this%v(i)%p(1)%i = j
      this%v(i)%p(1)%e = val
      this%v_size(i)   = 1
#ifdef CDEBUG
      call this%STYPE(test_sorted)(i, 'c')
#endif
      do ci = 1, CACHESIZE
        if (this%cache_idx(ci) == i) this%cache(ci, j) = 1
      end do
      return
    end if
    ! TODO: some old cruft to be removed. For now it is here to see if this is
    ! still called this way.
    if (j < 0) then
      write(*,*) "remove this (mod_sparse_matrix)"
      stop
    else
      ! find the index of the element, if it already exists
      this_i = this%STYPE(find_index)(i, j)
    end if
    ! directly set element if it already exists, and return
    if (this_i > 0) then
      this%v(i)%p(this_i)%e = val
#ifdef CDEBUG
      call this%STYPE(test_sorted)(i, 's')
#endif
      return
    end if
    this_i = -this_i
    ! insert new element
    ! if new size would be too large, grow the array
    if (this%v_size(i) + 1 > this%v_allocated(i)) then
      call this%STYPE(grow)(i)
    endif
#ifdef CDEBUG
    call this%STYPE(test_sorted)(i, 'p')
#endif
    this%v_size(i) = this%v_size(i) + 1
    ! Handle the special cases of appending an element
    if (this_i == this%v_size(i)) then
      this%v(i)%p(this%v_size(i))%i = j
      this%v(i)%p(this%v_size(i))%e = val
#ifdef CDEBUG
      call this%STYPE(test_sorted)(i, 'a')
#endif
      do ci = 1, CACHESIZE
        if (this%cache_idx(ci) == i) then
          do hi = j, size(this%cache(ci,:))
            this%cache(ci, hi) = sign(abs(this%cache(ci,hi)) + 1, this%cache(ci,hi))
          end do
          this%cache(ci, j) = this_i
        end if
      end do
#ifdef CDEBUG
      write(*,*) "after early set"
      call this%STYPE(test_cache)()
      write(*,*) "ok"
#endif
      return
    endif
    ! TODO: we might need to optimize this is it turns out that the compiler
    ! generates a temporary array on the stack for this overlapping array
    ! assignment
    this%v(i)%p (this_i+1 : this%v_size(i)) = this%v(i)%p (this_i : this%v_size(i)-1)
    this%v(i)%p (this_i)%i = j
    this%v(i)%p (this_i)%e = val
    do ci = 1, CACHESIZE
      if (this%cache_idx(ci) == i) then
        do hi = j, size(this%cache(ci,:))
          this%cache(ci, hi) = sign(abs(this%cache(ci,hi)) + 1, this%cache(ci,hi))
        end do
        this%cache(ci, j) = this_i
      end if
    end do
#ifdef CDEBUG
    call this%STYPE(test_sorted)(i, 'i')
    write(*,*) "after set"
    call this%STYPE(test_cache)()
    write(*,*) "ok"
#endif
  end subroutine STYPE(set_elem)

  ! get element, but return "default" if not found
  ! If 'j' is negative, it is assumed to be the sparse-matrix index,
  ! skipping the index search
  function STYPE(get_elem)(this, i, j, default)
    implicit none
    TYPE                                       :: STYPE(get_elem)
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer,                     intent(in)    :: i, j
    TYPE,                        intent(in)    :: default
    integer                                    :: this_i

    ! check that we have at least one element. If not, return default
    if (this%v_allocated(i) == 0 .or. this%v_size(i) == 0) then
      STYPE(get_elem) = default
    else
      if (j < 0) then
        this_i = -j
      else
        ! find the index of the element, if it already exists
        this_i = this%STYPE(find_index)(i, j)
      endif
      ! return default if it does not exit
      if (this_i <= 0) then
        STYPE(get_elem) = default
      else
        ! finally, return element in question
        STYPE(get_elem) = this%v(i)%p(this_i)%e
      end if
    end if
  end function STYPE(get_elem)

  subroutine STYPE(get_col_arrays)(this, i, col_size, indices, elements)
    implicit none
    class(STYPE(sparse_matrix)), intent(inout) :: this
    integer                    , intent(in)    :: i
    integer                    , intent(out)   :: col_size
    integer, allocatable       , intent(inout) :: indices(:)
    TYPE,    allocatable       , intent(inout) :: elements(:)
    integer :: j
    col_size = this%v_size(i)
    !allocate(indices(this%v_size(i)))
    !allocate(elements(this%v_size(i)))
    do j = 1, this%v_size(i)
      indices(j)  = this%v(i)%p(j)%i
      elements(j) = this%v(i)%p(j)%e
    end do
  end subroutine

end module
#undef TYPE
#undef STYPE
#endif
