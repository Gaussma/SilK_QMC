! This module handles the storage of the Hamiltonian. Any code outside of this
! file does not, and should not know about how the Hamiltonian is actually
! stored

#define INVALID -1.d99
#define NONZERO_SMALL ham_rho_min

#define CACHESIZE 20

module hamiltonian
  use vartypes
  use kink
  use mod_sparse_matrix_int
  use mod_sparse_matrix_ham

  real(dp)                          :: ham_rho_min

  integer,                  private :: ham_len              ! size of Hamiltonian
  real(dp),    allocatable, private :: ham_diag(:)          ! diagonal elements
  
  type(sparse_matrix_ham)       :: ham
  logical                       :: ham_initialized

contains

  function ham_size()
    implicit none
    integer :: i, ham_size

    ham_size = 0
    do i=1, nstate
      ham_size = ham_size + ham%v_size(i)
    end do
  end function ham_size

  subroutine test_new_ham()
    implicit none
    integer  :: i, j, idx
    type(ham_element) :: dummy
    real(dp) :: tmp_new
    if (.not. ham_initialized) return
    ! check
    do i = 1, nstate
      do j = 1, nstate
        idx = ham%find_index_ham(i, j)
        if (idx > 0) then
          dummy = ham%get_elem_ham(i, -idx, dummy)
          tmp_new = dummy%v
        else
          tmp_new = 0.d0
        end if
        if (get_ham(i,j) .ne. tmp_new) then
          write(*,*) "ham inconsistency", tmp_new, get_ham(i,j)
          stop
        end if
      end do
    end do
  end subroutine

  ! handle allocation
  subroutine init_hamiltonian(len_)
    implicit none
    integer, intent(in) :: len_
    integer  :: i, j

    ham_initialized = .false.
    ham_len = len_
    allocate(ham_diag(ham_len))

    write(*,*) " allocation done, now initializing."
    ham_diag    = INVALID
    call ham%init_ham(ham_len)
    write(*,*) "Calculating all Hamiltonian elements"
    !OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(GUIDED)
    !OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(DYNAMIC)
    do i = 1, nstate
      write(*,*) i, nstate
      do j = 1, i
        call set_ham(i, j, gen_simple(i,j))
      end do
    end do
    ham_initialized = .true.
    write(*,*) "done."
  end subroutine

  ! "internal" routine to set an Hamiltonian element. This routine does _NOT_
  ! make sure that the Hamiltonian stays symmetric (do that yourself), but in
  ! return it is thread-safe in 'i'.
  subroutine set_ham_(i, j, val)
    implicit none
    integer, intent(in)   :: i, j
    real(dp),  intent(in) :: val
    type(ham_element)     :: ham_elem
    integer               :: idx
    real(dp)              :: rho

    ! only store actual element if nonzero
    if (i == j) then
      rho = exp(-beta/p*val)
    else
      rho = exp(-beta/p*val) - 1.d0
    endif
    ! only store actual element if nonzero
    if (abs(rho) > NONZERO_SMALL) then
      ham_elem%v    = val
      ham_elem%rho1 = rho
      !if (i /= j) ham_elem%rho1 = ham_elem%rho1 - 1.d0
      if (i == j) ham_diag(i) = val
      call ham%set_elem_ham(i, j, ham_elem)
    else
      if (i == j) ham_diag(i) = 0.d0 ! This should not really happen
      ! If we get a request to set an element to zero that was previously
      ! non-zero, remove it. Otherwise, ignore the request.
      idx = ham%find_index_ham(i, j)
      if (idx > 0) then
        call ham%del_elem_ham(i, j)
      end if
    endif
  end subroutine

  subroutine set_ham(i, j, val)
    implicit none
    integer, intent(in)   :: i, j
    real(dp),  intent(in) :: val
    type(ham_element)     :: ham_elem
    integer               :: idx
    real(dp)              :: rho

    if ( i > 0 .and. j > 0 .and. i <= ham_len .and. j <= ham_len ) then
      if (i == j) then
        rho = exp(-beta/p*val)
      else
        rho = exp(-beta/p*val) - 1.d0
      endif
      ! only store actual element if nonzero
      if (abs(rho) > NONZERO_SMALL) then
        ham_elem%v    = val
        ham_elem%rho1 = rho
        !if (i /= j) ham_elem%rho1 = ham_elem%rho1 - 1.d0
        if (i == j) ham_diag(i) = val
        call ham%set_elem_ham(i, j, ham_elem)
        call ham%set_elem_ham(j, i, ham_elem)
      else
        if (i == j) ham_diag(i) = 0.d0 ! This should not really happen
        ! If we get a request to set an element to zero that was previously
        ! non-zero, remove it. Otherwise, ignore the request.
        idx = ham%find_index_ham(i, j)
        if (idx > 0) then
          call ham%del_elem_ham(i, j)
          call ham%del_elem_ham(j, i)
        end if
      endif
    else
      write(*,*) "Called set_ham with ",i,j,", but size is ", ham_len
      stop
    endif
  end subroutine

function gen_simple(k1,k2)
  use kink
  use gen_data
  implicit none

  integer,intent(in) :: k1, k2
  real(dp)             :: gen_simple

  logical found
  integer k3, k4, ne,itemp
  integer gen_simple_sign, ndiff, i, j
  integer it_list(ne_up + ne_dn), it_diff(ne_up + ne_dn)

  ne = ne_up + ne_dn
  ! First, we need to find, by how much jer(:,k1) and jer(:,k2) differ (ndiff),
  ! but since we also have to place them in maximum coincidence we need to
  ! reshuffle the entries in one of them a bit
  ! See [1], page 70 (section 2.3.3)
  it_list(:)      = jer(:, k2)
  ndiff           = 0
  gen_simple = 0.d0
  gen_simple_sign = 1
  do k3 = 1, ne
    ! Find index k4 in it_list with it_list(k4) == jer(k3,k1), and swap
    ! it with it_list(k3)
    found = .false.
    k4 = 1
    do while(k4 <= ne .and. .not. found)
      if (jer(k3,k1) == it_list(k4)) then
        found = .true.
        if (k4 .ne. k3) then
          gen_simple_sign = -gen_simple_sign
          ! swap it_list(k4) and it_list(k3)
          itemp       = it_list(k4)
          it_list(k4) = it_list(k3)
          it_list(k3) = itemp
        endif
      endif
      k4 = k4+1
    enddo
    ! If there was no such index (see above), save this in it_diff()
    if (.not. found) then
      ndiff = ndiff + 1
      it_diff(ndiff) = k3
      if (ndiff > 2) return
    endif
  enddo

  if (ndiff == 0) then
    gen_simple = repulsion_energy
    do k3 = 1, ne
      gen_simple = gen_simple + one_e_ham(jer(k3,k1), jer(k3,k1))
      do k4 = 1, ne
        gen_simple = gen_simple - two_e_ham(jer(k3,k1),jer(k4, 1),jer(k3,k1),jer(k4, 1)) &
                            + 0.5*two_e_ham(jer(k3,k1),jer(k4,k1),jer(k3,k1),jer(k4,k1))
      end do
    end do
  else if(ndiff == 1) then
    i = it_diff(1)
    j = it_list(i)
    gen_simple = one_e_ham(jer(i,k1),j)
    do k3 = 1, ne
      gen_simple = gen_simple - two_e_ham(jer(i,k1),jer(k3, 1),j,jer(k3, 1)) &
                              + two_e_ham(jer(i,k1),jer(k3,k1),j,jer(k3,k1))
    end do
  else if(ndiff .eq. 2) then
    gen_simple = two_e_ham(jer(it_diff(1), k1), jer(it_diff(2), k1), &
                           it_list(it_diff(1)), it_list(it_diff(2)))
  end if

  gen_simple = gen_simple * gen_simple_sign
end function gen_simple

  ! get element
  function get_ham(i, j)
    implicit none
    integer, intent(in) :: i, j
    real(dp)            :: get_ham
    type(ham_element)   :: result

    result%v = 0.d0
    result = ham%get_elem_ham(i, j, result)
    get_ham = result%v

    if (i == j .and. ham_diag(i) == INVALID) then
      ham_diag(i) = get_ham
    end if
  end function

  ! get rho(element)
  function rho_1(i, j)
    implicit none
    integer, intent(in) :: i, j
    real(dp)            :: rho_1
    type(ham_element)   :: result

    ! set default for when element is not stored
    result%rho1 = 0.d0
    result = ham%get_elem_ham(i, j, result)
    rho_1 = result%rho1
  end function

  ! get multiple elements by reference (but don't use this to change ham!)
  subroutine get_hams(i1, i2, j1, j2, arr, k1, k2, l1, l2)
    implicit none
    integer, intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    real(dp), allocatable :: arr(:,:)
    integer             :: k, l

    if ( i1 > 0 .and. j1 > 0 .and. i1 <= ham_len .and. j1 <= ham_len .and. &
         i2 > 0 .and. j2 > 0 .and. i2 <= ham_len .and. j2 <= ham_len ) then
      !OMP PARALLEL DO PRIVATE(k,l)    DISABLED
      do k = 0, k2-k1
        do l = 0, l2-l1
          arr(k1+k, l1+l) = get_ham(i1+k, j1+l)
        end do
      end do
    else
      write(*,*) "Called get_hams with ",i1,i2,j1,j2,", but size is ", ham_len
      stop
    endif
  end subroutine

  ! return the intersection of two columns of the Hamiltonian
  ! c1, c2: the two columns
  ! C: the resulting array, with the used size stored in nc
  ! These are only indexes to the actual data in the Hamiltonian, and what is
  ! returned is an array which contains index entries for which both input
  ! arrays contain non-zero data.
  recursive subroutine ham_column_intersection(c1, c2, C, nc)
    implicit none
    integer,              intent(in)  :: c1, c2
    integer,              intent(out) :: C(ham_len*2), nc
    integer                           :: i, i2, c1_idx

    ! First, make sure that column c1 is at least as large as column c2
    ! (for later optimizations)
    if (ham%v_size(c2) > ham%v_size(c1)) then
      call ham_column_intersection(c2, c1, C, nc)
      return
    end if

    nc = 0
    ! For about equal sized arrays, simply walk both arrays simultaniously and
    ! save the common elements. For arrays of vastly different size, searching
    ! can be faster, and can be done in parallel. Assume the latter for now.
    ! Go through smaller array
    do i = 1, ham%v_size(c2)
      i2 = ham%v(c2)%p(i)%i
      ! For each element of smaller array, find corresponding element in
      ! larger rarray
      c1_idx = ham%find_index_ham(c1, i2)
      ! if found, save and continue (if not found, do nothing)
      if (c1_idx > 0) then
        nc = nc + 1
        C(nc) = i2
      end if
    end do
  end subroutine ham_column_intersection

  ! Save hamiltonian in matrix form - so far only diagonal as proof of concept
  subroutine checkpoint_ham()
    use checkpointing
    implicit none
    integer(int8)                  :: len8
    integer                        :: col_size, i, j
    integer,           allocatable :: indices(:)
    integer(int8),     allocatable :: indices8(:)
    integer(int8),     allocatable :: col_size8(:)
    real(dp),          allocatable :: ham_col(:)
    type(ham_element), allocatable :: elements(:)
    character(len=100)             :: fieldname

    len8 = ham_len
    allocate(indices(ham_len))
    allocate(elements(ham_len))
    allocate(col_size8(ham_len))
    allocate(indices8(ham_len))
    allocate(ham_col(ham_len))
    call checkpoint_int_scalar("ham_len", len8)
    call checkpoint_real_array("ham_diag",    ham_diag,    len8)
    do i = 1, ham_len
      call ham%get_col_arrays_ham(i, col_size, indices, elements)
      indices8(:col_size) = indices(:col_size)
      do j = 1, col_size
        ham_col(j)   = elements(j)%v
        if (indices(j) <= i) col_size8(i) = j
      end do
      write(fieldname, "(A6,I0.10)") "ham_i_", i
      call checkpoint_int_array(trim(fieldname), indices8, col_size8(i))
      write(fieldname, "(A6,I0.10)") "ham_e_", i
      call checkpoint_real_array(trim(fieldname), ham_col, col_size8(i))
    end do
    call checkpoint_int_array("v_size", col_size8, len8)
    deallocate(ham_col)
    deallocate(indices8)
    deallocate(col_size8)
    deallocate(elements)
    deallocate(indices)
  end subroutine

  ! Save hamiltonian in matrix form - so far only diagonal as proof of concept
  subroutine restart_ham()
    use checkpointing
    implicit none
    integer(int8) :: len8
    integer                        :: i, j
    integer(int8),     allocatable :: col_size8(:)
    integer(int8),     allocatable :: indices8(:)
    real(dp),          allocatable :: ham_col(:)
    character(len=100)             :: fieldname

    call restart_int_scalar("ham_len", len8)
    ham_len = int(len8)
    allocate(ham_diag(ham_len))
    ham_diag = INVALID
    call ham%init_ham(ham_len)

    allocate(col_size8(ham_len))
    allocate(indices8(ham_len))
    allocate(ham_col(ham_len))
    call restart_real_array("ham_diag",    ham_diag,   len8)
    call restart_int_array("v_size",      col_size8,   len8)
    do i = 1, ham_len
      write(fieldname, "(A6,I0.10)") "ham_i_", i
      call restart_int_array(trim(fieldname), indices8, col_size8(i))
      write(fieldname, "(A6,I0.10)") "ham_e_", i
      call restart_real_array(trim(fieldname), ham_col, col_size8(i))
      do j = 1, int(col_size8(i))
        call set_ham(i, int(indices8(j)), ham_col(j))
      end do
    end do
    deallocate(ham_col)
    deallocate(indices8)
    deallocate(col_size8)
  end subroutine

  ! output some info about the hamiltonian elements to a file
  subroutine output_ham_vinfo(ipass, nbins)
    use io
    implicit none
    integer, intent(in) :: ipass, nbins
    integer  :: i, j, bin, col_size
    integer  :: bin_count(0:nbins+1)
    integer  :: bin_count_nondiag(0:nbins+1)
    real(dp) :: bin_bound(0:nbins+1)
    integer,           allocatable :: indices(:)
    type(ham_element), allocatable :: elements(:)

    if (ipass == 0) then
      open (unit=8677, file=trim(adjustl(outdir))//'/'//"ham_vinfo.dat" ,&
            status="unknown")
      write(8677,*) "# 0:pass 1: bin nr 2: lower bin bound 3: bin count 4: bin count nondiag"
    else
      open (unit=8677, file=trim(adjustl(outdir))//'/'//"ham_vinfo.dat" ,&
            status="old",access='append')
    endif
    ! determine bin boundaries from NONZERO_SMALL (and 1)
    bin_bound(0) = -1.d99
    do i = 1, nbins+1
      bin_bound(i) = ((log(10.d0)-log(NONZERO_SMALL)) / nbins) * (i-1) + &
                     log(NONZERO_SMALL)
    end do
    ! get elements from sparse matrix in a "fast" fashion, and count them in bins
    allocate(indices(ham_len)) 
    allocate(elements(ham_len))
    bin_count = 0
    bin_count_nondiag = 0
    do i = 1, ham_len
      call ham%get_col_arrays_ham(i, col_size, indices, elements)
      do j = 1, col_size
        bin = int((log(abs(elements(j)%rho1)) - log(NONZERO_SMALL)) / &
                  (log(10.d0) - log(NONZERO_SMALL)) * nbins) + 1
        if (bin < 0)       bin = 0
        if (bin > nbins+1) bin = nbins+1
        bin_count(bin) = bin_count(bin) + 1
        if (i /= indices(j)) then
          bin_count_nondiag(bin) = bin_count_nondiag(bin) + 1
        end if
      end do
    end do
    write(8677,*) ipass, 0, 0.d0, bin_count(0), bin_count_nondiag(0)
    do i = 1, nbins+1
      write(8677,*) ipass, i, exp(bin_bound(i)), bin_count(i), &
                    bin_count_nondiag(i)
    end do
    write(8677,*) ""
    close(8677)
    deallocate(elements)
    deallocate(indices)
  end subroutine


end module hamiltonian

