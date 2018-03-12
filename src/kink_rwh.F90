#define BETAX 0.405

#define PARAM_VERSION 17

module mpi_data
  integer:: my_id, ierr
  integer:: num_procs
  integer:: root_process

  contains





  ! Determine array bounds of decomposed array of size N; bounds are inclusive
  ! and are using Fortran-array bounds (1..N)
  subroutine decompose(N, low, high, N_local)
    implicit none
    integer, intent(in)  :: N
    integer, intent(out) :: low, high, N_local
  
    integer :: small_size ! size of the small parts
    integer :: N_small    ! number of processes using small_size
    integer :: num_procs_used ! might be smaller than num_procs
  
    num_procs_used = num_procs
    if (N < num_procs) then
      num_procs_used = N
      if (my_id >= N) then
        low     = 0
        high    = -1
        N_local = 0
        return
      endif
    endif
    small_size = N/num_procs_used
    N_small = num_procs_used + small_size*num_procs_used - N
    if (my_id < N_small) then
      low  =  my_id   * small_size + 1
      high = (my_id+1)* small_size
      N_local = small_size
    else
      low  = N_small * small_size + (my_id  -N_small) * (small_size+1) + 1
      high = N_small * small_size + (my_id+1-N_small) * (small_size+1)
      N_local = small_size + 1
    endif
  end subroutine decompose
end module mpi_data

module kink_temp
  use vartypes
  integer,allocatable :: list_of_states(:) 
  real(dp),allocatable:: Ha_initial_a(:,:)
  integer,allocatable::Ha_initial_b(:,:) ,Ha_initial_max(:)
  real(dp),allocatable:: Ha_temp(:,:)
  integer,allocatable:: list_of_find(:)
  real(dp),allocatable::list_choice(:)
  integer:: nappear_find
  integer:: nappear_happen
  integer,allocatable :: list_of_happen(:)
  integer,allocatable :: list_of_happen_sort(:)
end module kink_temp


module sfunc_mod
  use vartypes
  integer, allocatable :: ig(:),itaken(:) ! size n
  real(dp),  allocatable :: wvec(:)         ! size n
  real(dp),  allocatable :: fvec(:)         ! size 0..n-1
end module


Program Kink_QMC
  use vartypes
  use kink
  use kink_temp

  use mpi_data
  use gen_data
  use io
  use timing
  use hamiltonian
  use checkpointing
  use mod_random


  implicit none
  include 'mpif.h'

  interface
   function mytime()
     use vartypes
     integer(int8) :: mytime
   end function
  end interface

  real(dp),allocatable:: xv(:)
  integer,allocatable:: j_old(:)

  integer(int8) :: time_h, time_min, time_sec
  real(dp)  :: time_dbl

  integer :: l,k
  integer :: i,j
  real(dp) :: energy,sfunc,sign_sfunc, maxtmp
  integer:: output_count
  real(dp):: accume_e,accume_sign,accume_kink
  integer:: k1,k2,k4

  real(dp),allocatable :: xvv(:)
  real(dp),allocatable :: rho_add(:,:)
  real(dp),allocatable :: rho_minus(:)
  real(dp),allocatable :: rho_evec(:,:)
  ! Removal,Insertion variables
  integer,allocatable :: jtrial_rem(:),jtrial_insert(:)
  integer,allocatable :: j_insert_add(:)
  integer,allocatable :: j_insert_minus(:)
  real(dp) ::rn3,rn4,rn5
  real(dp) :: xnorm,xnorm_add,xnorm_rem,xnorm_insert,xnorm_rem_add, xnorm_insert_add
  real(dp) :: tmp_sum
  integer :: ipick,ipick_n,ipick_l

  integer:: iplus,iminus
  integer:: seed_index
  real(dp)::sfunc_old
  real(dp):: sum_t
  integer:: jtmp
  integer:: ip1,ip2,ip3,ip4

  real(dp),allocatable::Ha_local(:,:) 
  logical,allocatable::iused(:) 
  integer,allocatable::jused(:)
  real(dp),allocatable::rho_mat(:,:)
  logical:: done
!  logical:: changed_list_of_find
  integer,allocatable:: times_appear(:)
  integer,allocatable:: list_of_states_init(:)
  
  integer:: threshold

  integer::nappear
  integer::nappear_init
  real(dp),allocatable::w(:),work(:)
  integer::INFO
  real(dp):: sf_both

  REAL*8:: ar4,ar5,ar6
  integer :: iar4,iar5,iar6
  real(dp)::sum5,sum3,sum4
  real(dp):: q_est
  integer :: sign_final
  integer :: ver_number
  integer :: calc_rho
  integer::turn_diag
  real(dp):: big_rho, big_rho_insert, big_rho_rem
  real(dp):: big,big_rem,big_insert
  real(dp)::rho_big,rho_big_rem,rho_big_insert 
  real(dp) :: E_tr
  character(len=1000) :: arg
  character(len=4) :: restrict_X ! SD/SDT/none
!!! to calcualte partition function , add on Jan 15
  !!add some new variable 
  real(dp) ::energy_hartree_fock ,Z_total,Z_total_insert
  real(dp),allocatable ::T_add(:),T_minus(:),T_rem_minus(:),T_insert_minus(:)
  real(dp) ::sum_exp_T_add
  integer::ipick_plus,ipick_minus
 
  integer:: k_rem,k_rem_mid,k_rem_s
  integer:: k_insert,k_insert_mid,k_insert_s 
  integer::k3 

  integer*8:: getRSS
#if 0
  !we temporarily do not use those variables.
  real(dp):: magnetization_ave,magn_run
  real(dp):: magnetization_nave


  real(dp):: magnetization4_ave
  real(dp):: magnetization3_ave
  real(dp):: magnetization2_ave
  real(dp):: magnetization0_ave
  real(dp):: magnetization4_nave
  real(dp):: magnetization3_nave
  real(dp):: magnetization2_nave
  real(dp):: magnetization0_nave
  integer::output4_count
  integer::output3_count
  integer::output2_count
  integer::output0_count
#endif

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

  allocate(fact(0:10000))
  ! Calculate the factorial
  fact(0)=1
  do k1=1,10000
    fact(k1)=k1*fact(k1-1)
  end do
  ar4=0.d0
  ar5=0.d0
  ar6=0.d0
  iar4=0
  iar5=0
  iar6=0
  call_size=0 
  call getarg(1, arg)

  open (unit=688, file=trim(adjustl(arg)) ,status="unknown")
  !write(*,*) "please write the version number"
  read (688,*) ver_number

  ! make judgement on version Number
  if (ver_number .eq. PARAM_VERSION) then
     if( my_id.eq.0)then
     write(*,*) "version number is" ,ver_number
     end if
  else
     write(*,*) "wrong version number"
     stop
  end if

  !write(*,*) ' specify the symmetry'
  read(688,*) symmetry
  !check for supported symmetries to avoid uncaught errors because of typos
  if (symmetry /= 'c2v' .and. symmetry /= 'd2h' .and. symmetry /= 'sto') then
    write(*,*) '>',symmetry, '< is not a supported symmetry'
    stop
  endif
  !write(*,*) 'basis set'
  read(688,*) basisset
  if (basisset /= 'ST ' .and. basisset /= 'DZ ' .and.\
      basisset /= 'STH'.and. basisset /= 'DZH') then
    if (my_id == 0) write(*,*) '>',basisset, '< is a unrecognized basis set. There is nothing wrong with that, but it might be a typo.'
  endif
  restrict_SD  = .false.
  restrict_SDT = .false.
  read(688,*) restrict_X ! SD/SDT/none
  if (restrict_X == 'none') then
  else if (restrict_X == 'SD') then
    restrict_SD = .true.
    write(*,*) "Using SD"
  else if (restrict_X == 'SDT') then
    restrict_SD = .true.
    restrict_SDT = .true.
    write(*,*) "Using SDT"
  else
    write(*,*) "Unknown string for restriction, supported: SD/SDT/none, got: ", restrict_X
    stop
  end if
  if (my_id == 0) write(*,*) 'Using symmetry ', symmetry, ' and basis set ', basisset
  !Write(*,*) ' please write the number of rows of the lattice'
  Read(688,*) nrows
  !Write(*,*) ' please write the number of columns of the lattice'
  Read(688,*) ncols
  !Write(*,*) ' please write the number of up e-'
  Read(688,*) ne_up
  !Write(*,*) ' please write the number of down e-'
  Read(688,*) ne_dn
  ! e_up_state = int(fact(nrows*ncols)/fact(nrows*ncols)/fact(ne_up));
  ! e_dn_state = int(fact(nrows*ncols)/fact(nrows*ncols)/fact(ne_dn));
  !Write(*,*) ' please write the number of kinks'
  read(688,*) nkink
  !Write(*,*) ' please write the number of passes'
  read(688,*) npass
  !write(*,*) ' For how many passes should we diagonalize? (npass_diag)'
  read(688,*) npass_diag
  !write(*,*) ' Every how many passes should we diagonalize before npass_diag? (npass_diag)'
  read(688,*) diag_every
  !write(*,*) ' When should we do the first diagonalization (independend of other parameters)?'
  read(688,*) first_diag
  !write(*,*) " please write max number of kinks at any one time "
  read(688,*) nmax_kinks
  !write(*,*) " please write nmax_states to diagonalize"
  read(688,*) nmax_states
  !write(*,*) ' please enter beta '
  read(688,*) beta
  !write(*,*) ' please enter p '
  read(688,*) p
  !write(*,*( ' please enter the minimum value of rho that should be stored in the hamiltonian'
  read(688,*) ham_rho_min
  !write(*,*) ' please enter number of passes to combine'
  read(688,*) ncomb
  !write(*,*) ' please enter the name of the directory output should go to'
  read(688,*) outdir
  !write(*,*) ' please enter the name of the directory containing any input data files'
  read(688,'(A)') indir
  !write(*,*) ' please enter the file name of the rho and log_rho'
  read(688,*) rho_filename
  !write(*,*) ' please indicate whether you want to calculate rho (1) or not (0)'
  read(688,*) calc_rho
  !write(*,*) " diagonalizaiton, 1 for do diagonalization, 0 for skip diagonalization"
  read(688,*)  turn_diag
  !write(*,*) " please input seed index (0 if unsure)"
  read(688,*) seed_index
  read(688,*) x_bar 
  read(688,*) n_occupied 
  read(688,*) switch_frozen 
  read(688,*) x_dimension

  ! if outdir is the special string '$parfile', then outdir will be whatever
  !  the name of the parameter file was, minus the extension .par
  if (trim(adjustl(outdir)) == '$parfile' .and. len_trim(adjustl(arg)) > 4 ) then
    outdir = arg(1:len_trim(adjustl(arg))-4)
    ! Only use the name of the file, not the whole path to the parameter file
    k1 = scan(outdir, '/', .true.)
    if (k1 > 0) then
      outdir = outdir(k1+1:)
    endif
  endif
  call system('mkdir -p '//trim(adjustl(outdir)))
  time_start = mytime()

  close(688)

  ! Read repulsion energy
  open(unit=126,file=trim(adjustl(indir)) // '/repulsion_energy.dat',status='OLD')
  read(126,*) repulsion_energy
  close(126)

!#define PDEB \
!  open (unit=4242+my_id, file=trim(adjustl(outdir))//'/'//"debug_proc"//trim(debugfilename)//".dat" ,status="unknown",access='append'); \
!  write(4242+my_id, *)
!#define PDEBC close(4242+my_id)
!
!  write(debugfilename, '(I5.5)') my_id
!  open (unit=4242+my_id, file=trim(adjustl(outdir))//'/'//"debug_proc"//trim(debugfilename)//".dat" ,status="unknown")
!  write(4242+my_id, *)
   
  e_up_state = int(fact(nrows*ncols-n_occupied)/fact(nrows*ncols-ne_up)/fact(ne_up-n_occupied));
  e_dn_state = int(fact(nrows*ncols-n_occupied)/fact(nrows*ncols-ne_dn)/fact(ne_dn-n_occupied));

  ! Note: this is changed further down
  nstate = e_up_state * e_dn_state;

  if (my_id == 0) then
    write(*,*) ' That makes ',nstate, ' total states (nstate without SY).'
  endif

  call init_random(seed_index)
 
  open (unit=1951, file=trim(adjustl(outdir))//'/'//"screen.dat" ,status="unknown")
  Write(1951,*)    '# screen data from windows computer'


  num_orbital=2*nrows*ncols
  allocate(vect(ne_up)) 
  vect = 0

  ! NOTE: This call changes 'nstate'
  open(unit=1900,file=trim(adjustl(outdir))//'/'//'1900_jer.dat',status='unknown')
  open(unit=1901,file=trim(adjustl(outdir))//'/'//'1901_jer.dat',status='unknown')
  call setup_jer(1,1, my_id)
  beta_x=BETAX

 

  if (my_id == 0) then
    write(*,*) ' Symmetries reduced number of states to ',nstate, '.'
  endif

    open(unit=2037,file=trim(adjustl(outdir))//'/'//'2037.dat',status='unknown')
    open(unit=2038,file=trim(adjustl(outdir))//'/'//'2038.dat',status='unknown')
  
#if 0
  do k1=1,nstate 
   if( mod(k1,1000).eq.0) then 
    write(*,*) k1, "Year "
   end if  
    i=0
    do k2=1,nstate
      if( gen_simple(k1,k2).ne.0) then 
        i=i+1 
       write(2037,*) gen_simple(k1,k2),k2 
      end if 
    end do 
    write(2038,*) i,k1 
   ! if (my_id == 0) write(*,*), i, k1 , "i,k1"
  end do 
#endif 
 



  if (my_id == 0) write(*,*) "ham initialization"
  if (.not. load_ham_init()) then
    call init_hamiltonian(nstate)
    call save_ham_init()
  end if
  !if (my_id == 0) write (*,*) "writing Hamiltonian to disk"
  !call checkpoint()
  ! Just for testing
  !call restart()
  if (my_id == 0) write (*,*) "done"
  
  allocate(Ha_temp(nmax_states,nstate))  
  !write(*,*) "nstate is: ", nstate

  allocate(Ha_local(nmax_states,nstate))
  !Ha_local=0 
  allocate(iused(nstate)) 
  allocate(jused(nstate))
  ! Allocate persistant memory for sfunc_calc()
  allocate(T_add(nstate))
  allocate(T_minus(nstate)) 
  allocate(T_rem_minus(nstate))
  allocate(T_insert_minus(nstate))
  call sfunc_init(nstate)

  allocate(j_old(nstate))
  allocate(xv(nstate))
  allocate(xvv(nstate))
  allocate(jtrial_rem(nstate-1))
  allocate(jtrial_insert(nstate+1))
  allocate(j_insert_add(nstate+2))
  allocate(j_insert_minus(nstate))
  allocate(rho_minus(nmax_kinks))
  j_old = 0
  xv = 0
  xvv = 0
  jtrial_rem = 0
  jtrial_insert = 0
  j_insert_add = 0
  j_insert_minus = 0

  allocate(rho_add(nmax_kinks+5, nstate))
  !rho_add = 0

  allocate(rho_mat(nmax_states,nmax_states))
  rho_mat = 0
  allocate(Ha_diag(nstate))
  allocate(inv_exp_betax_Ha_diag(nstate))
  Ha_diag = 0
  inv_exp_betax_Ha_diag = 1
  !allocate(Ha_temp(nmax_states,nstate))
  !Ha_temp = 0

  !allocate(rho_evec(nstate,nstate))

  allocate(list_of_states(nstate))
  list_of_states = 0
  call ham%update_cache_ham(0, 0, list_of_states)
  call ham%output_info_ham(0)
  call output_ham_vinfo(0, 100)
  allocate(list_of_states_init(nstate))
  list_of_states_init=0
  allocate(list_of_find(nstate))
  list_of_find = 0
  allocate(times_appear(nstate))
  times_appear = 0
  allocate(nnz_low_bound_record(nmax_kinks+1)) 
  allocate(nnz_high_bound_record(nmax_kinks+1)) 
  
 
 !@ allocate(w(nstate),work(nstate*3))
  !keep it for exact diag purpose


 
  allocate(Ha_low(nstate))
  Ha_low = 0
  allocate(list_of_happen(nstate)) ! the list of all states happens in diagonalization
  list_of_happen = 0
  allocate(list_of_happen_sort(nstate))
 
  k_rem=0
  k_rem_mid=0
  k_rem_s=0
  k_insert=0
  k_insert_mid=0
  k_insert_s=0  



!********************* MODULE CSR

  ! Open data files
  if (my_id.eq.0) then
    open (unit=102, file=trim(adjustl(outdir))//'/'//"kink_diag_npass_nkink.dat" ,status="unknown")
    Write(102,*)    '# npass nkink'
    open (unit=103, file=trim(adjustl(outdir))//'/'//"kink_diagonalization.dat" ,status="unknown")
    Write(103,*)    '# diagonalization'
    open (unit=333, file=trim(adjustl(outdir))//'/'//"test_diag.dat" ,status="unknown")
    Write(333,*)    '# state'
    open (unit=1011, file=trim(adjustl(outdir))//'/'//"timing.dat" ,status="unknown")
    Write(1011,*)    '# 1: ipass, 2: total time, 3: total evolution time, time since last output(4:evol, 5: part1, 6: part2, 7: part3, 8: part4, 9:ext_diag, 10:debug), total time (11: evol, 12: part1, 13: part2, 14: part3, 15: part4, 16:ext_diag, 17:debug), 18: memory usage'
  
    open(unit=1934,file=trim(adjustl(outdir))//'/'//'interval_2nd.dat',status='unknown')
    open(unit=1935,file=trim(adjustl(outdir))//'/'//'1935_of_check.dat',status='unknown')
    open(unit=1936,file=trim(adjustl(outdir))//'/'//'nappear_happen_vs_ipass.dat',status='unknown')
    open(unit=1939,file=trim(adjustl(outdir))//'/'//'truncation_Q2_ipass.dat',status='unknown')
    open(unit=1949,file=trim(adjustl(outdir))//'/'//'acceptance_ratio_rem.dat',status='unknown')
      write(1949,*) "# k_rem/10*ipass, k_rem_mid/k_rem,k_rem_s/k_rem,  "
    open(unit=1950,file=trim(adjustl(outdir))//'/'//'acceptance_ratio_insert.dat',status='unknown')
      write(1950,*) "# k_insert/10*ipass, k_insert_mid/k_insert,k_insert_s/k_insert,  "

   ! open(unit=2037,file="/project/xma9/H2O_2R/2037.dat",status='unknown')
   ! open(unit=2038,file="/project/xma9/H2O_2R/2038.dat",status='unknown')



  end if
  ground_state_index = 0
  root_process=0  ! Define Root Process is Zero in MPI.

  num_orbital=2*nrows*ncols
  num_total=int(fact(num_orbital)/fact(num_orbital-2)/fact(2),4)

  call decompose(num_total, two_e_low_bound, two_e_high_bound, two_e_block_size)

  

  ! TODO: Is the first size correct for two_e()?
  allocate(two_e(two_e_block_size,num_orbital,num_orbital))
  two_e = 0
  if (my_id == 0) then
    write(*,*) "Allocating two_e"
    write(*,*) "  dimensions :", two_e_block_size, num_orbital, num_orbital
    write(*,*) "  size(bytes):", two_e_block_size*num_orbital*num_orbital*8
    write(*,*) "  total size :", num_orbital*num_orbital*num_orbital*num_orbital*8
  endif
  allocate(list_choice(0:num_procs))
  list_choice = 0

  !call  get_hamil()
  !call gen(1,1,h_e)  ! we must call this to initialize the one_e_ham and two_e_ham and jer! very important !
  ! write(*,*) h_e, "gen(1,1) "

  !call gen_simple(1,1,h_new)
  !write(*,*) h_e ,h_new, "h_e"

  if (nstate < 150) then
    call save_gen()
  endif

  !call read_2038()
  !call read_2037()
  
  !write(*,*)  Ha_initial_a( 1,1), Ha_initial_a(1,2), "2037"
  !stop 
  
  ! we read in the non zero element of initial Hamiltonian ! 

  ! do some debugging output

 
  time_rho = mytime()
  ! either load rho from file, or compute it
  if (calc_rho .eq. 0) then
    if (my_id == 0) write(*,*) "loading rho from file"
    call load_rho()
    energy_hartree_fock=gen_simple(1,1)
    inv_exp_betax_Ha_diag(1:nstate) = exp(-BETAX*(Ha_diag(1:nstate)-energy_hartree_fock))
    sum_exp_T_add = sum(inv_exp_betax_Ha_diag(1:nstate))
  else
    if (my_id == 0) write(*,*) "calculating rho"
    do i=1,nstate
      Ha_diag(i) = gen_simple(i,i)
    end do    ! this is the right method, rather than using get_hamil()
    energy_hartree_fock=gen_simple(1,1)
    ! precompute these exponential functions (and update later as necessary)
    inv_exp_betax_Ha_diag(1:nstate) = exp(-BETAX*(Ha_diag(1:nstate)-energy_hartree_fock))
    ! we only need to calculate sum_exp_T_add next time when inv_exp_betax_Ha_diag update!  
    sum_exp_T_add = sum(inv_exp_betax_Ha_diag(1:nstate))
    call calc_exact_energy()
  end if ! compute rho
  time_rho = mytime() - time_rho
  if (my_id == 0) write(*,*) "rho done (", time_rho, "ms), preparing evolution"

  ! TODO: comment this re-use of these variables
  ! Change the dimension of rho_evec, from 225 to 15 for example
  ! set namx_kinks is the dimension for sub diagonalization
   !@deallocate(rho_evec)
  !@ deallocate(w,work)
 !!@ allocate(rho_evec(nmax_states,nmax_states))
 !!@ rho_evec = 0

 !!@ allocate(w(nmax_states),work(nmax_states*3))
  !  w = 0
  !  work = 0

  ! WARNING: This only works for nkink=1, otherwise we need to initialize
  ! j_old(2) and j_old(n_kink) as well
  ! TODO: implement for generic nkink
  ! set for the initial condition
  j_old=0
  j_old(1) = 1
  list_of_find=0
  nappear_find=1 
  list_of_find(1)=j_old(1)

  call MPI_BARRIER (MPI_COMM_WORLD,ierr) 

   do k2=1,num_procs

  if ( (my_id+1).eq.k2) then 
    ! call MPI_BARRIER (MPI_COMM_WORLD,ierr)  
 ! else 
    ! call MPI_BARRIER (MPI_COMM_WORLD,ierr)  
  end if 
  call MPI_BARRIER (MPI_COMM_WORLD,ierr)  

  end do  
  do k2= 1, nstate
    Ha_temp(nappear_find,k2) = get_ham(list_of_find(nappear_find), k2)
  end do
  !call get_hams(list_of_find(nappear_find), &
  !              list_of_find(nappear_find), &
  !              1, nstate, &
  !              Ha_temp, &
  !              nappear_find, &
  !              nappear_find, &
  !              1, nstate)
   !call MPI_BARRIER (MPI_COMM_WORLD,ierr)  
  !call mpi_finalize(ierr) 

  ! running QMC
  output_count=0

  ! initialize what is the value of sfunc_old
  ! (calculate S: one part of the weight)
  call sfunc_calc(j_old,nstate,nkink,p,sfunc_old,sign_sfunc)

  ! set initial value of IP (counters for the number of states)
  ip1=0
  ip2=0
  ip3=0
  ip4=0

  nappear=0
  nappear_init=0

  accume_e=0.d0
  accume_sign=0.d0
  accume_kink=0.d0
  ipass=0

  nappear_happen=0
  list_of_happen=0

  ! nappear_find=0
  ! list_of_find=0 
  
  timer_evol       = 0.d0
  timer_evol_part1 = 0.d0
  timer_evol_part2 = 0.d0
  timer_evol_part3 = 0.d0
  timer_evol_part4 = 0.d0
  timer_ext_diag   = 0.d0
  timer_debug      = 0.d0
  ttimer_evol       = 0.d0
  ttimer_evol_part1 = 0.d0
  ttimer_evol_part2 = 0.d0
  ttimer_evol_part3 = 0.d0
  ttimer_evol_part4 = 0.d0
  ttimer_ext_diag   = 0.d0
  ttimer_debug      = 0.d0

  if (my_id == 0) write(*,*) "starting evolution"
  time_start_evol = mytime()
  do while(ipass .le. npass)  
    ipass=ipass+1
    time_start_evol_it = mytime() 

    ! Timing
    if (.false. .and. modulo(ipass*10, npass) .eq. 0) then
      time_dbl = (mytime() - time_start) * (dble(npass)/ipass - 1)
      time_h   = int(time_dbl/60/60)
      time_min = int((time_dbl-time_h*60*60) / 60)
      time_sec = int( time_dbl-time_h*60*60-time_min*60)
      write(*,*)  (mytime() - time_start), "time",ipass
    end if

    time_start_part1 = mytime() 
    do imulti=1,10
      do i= 1,nkink
        iplus  = 1+modulo(i-1 +1,nkink)
        iminus = 1+modulo(i-1 -1,nkink)

        ! l: random state
        l = min(nstate, 1+int(random()*nstate))

        ! first criteria is l is not equal to these three value
        if (j_old(iplus)  .ne. l .and. &
            j_old(iminus) .ne. l .and. &
            j_old(i)      .ne. l) then
          ! second criteria is if l is ok, the next configuration not zero
          if (nkink .eq. 1 .or.  (rho_1(j_old(iminus),l).ne.0.0 .and. &
                                  rho_1(j_old(iplus ),l).ne.0.0)) then
            ! this is to calculate the energy difference first part, log,
            energy = 0.d0
            if (nkink > 1) then
              ! calculate the energy difference in log
              energy = log_rho_1(j_old(iminus), l       ) + log_rho_1(j_old(iplus ), l       ) - &
                       log_rho_1(j_old(iminus), j_old(i)) - log_rho_1(j_old(iplus ), j_old(i))

            endif
            ! We call sfunc_calc with a slightly changed j_old, so we need to save
            ! what we changed
            jtmp     = j_old(i)
            j_old(i) = l
            call sfunc_calc(j_old,nstate,nkink,p,sfunc,sign_sfunc)
            ! restore change
            j_old(i) = jtmp

            energy = energy + sfunc - sfunc_old

            ! do QMC update only if energy got higher than random number 
            if (energy .ge. log(random())) then
              j_old(i)  = l
         
              sfunc_old = sfunc
              ! update list_of_find with new states
              ! list_of_find is list that contains all states that appeared between
              !   two diagonalizations
              ! first place to update list_of_find, guarantee that no number repeated in
              !   list_of_find
              if (.not. any(list_of_find(1:nappear_find) == l) .and. &
                  nappear_find < nmax_states) then
                nappear_find = nappear_find+1 
                list_of_find(nappear_find) = l
                do k2 = 1, nstate
                  Ha_temp(nappear_find,k2) = get_ham(list_of_find(nappear_find),k2)
                end do
!  call get_hams(list_of_find(nappear_find), &
!                list_of_find(nappear_find), &
!                1, nstate, &
!                Ha_temp, &
!                nappear_find, &
!                nappear_find, &
!                1, nstate )
              end if 
            end if ! update
          end if ! second criteria
        end if ! first criteria
      end do ! nkink
    enddo ! imulti
    timer_evol_part1 = timer_evol_part1 + (mytime() - time_start_part1)/1000.d0
       
    !check off diagonal element in rho_1 matrix 

    ! ----------------------- Second Part ----------------------
    time_start_part2 = mytime()
    beta_x=0.405
    if (ipass.ge.10000000) then  
      energy_hartree_fock = gen_simple(1,1) 
      T_add(1:nstate) = -beta_x * (Ha_diag(1:nstate) - energy_hartree_fock)

      ! from there , rewrite Second Part using new T(c'|c) !
      do imulti = 1,10
        if (nkink > 1) then
          do k1 = 1,nkink 
            T_minus(k1) = beta_x*(Ha_diag(j_old(k1)) - energy_hartree_fock)
          end do 
          big = maxval(T_minus(1:nkink))
        else
          T_minus(1) = -8800.00
          big        = -8800.00
        end if
        big = max(big, maxval(T_add(1:nstate)))

        Z_total = sum_exp_T_add / exp(big)
        xnorm_add = (nkink+1) * Z_total 
        xnorm = xnorm_add + sum(exp(T_minus(1:nkink) - big))
        call weight_kink(nstate,j_old,rho_big)

        rn3 = random()
        ! Do we remove or insert a random kink - here: remove?
        ! TODO: figure out if we can elininate calculating xnorm
        !       if rn3 > xnorm_add/xnorm .and. nkink == 1 (which might be often)
        if (rn3 > xnorm_add/xnorm .and. nkink > 1) then
          k_rem=k_rem+1 
          ! Try to remove the kink - find out which to remove based on rn3
          ! we must include the case nkink.gt.1, if nkink=1, we can not remove the kink 
          tmp_sum = xnorm_add/xnorm
          done = .false.
          ipick = 0
          do while (.not.done .and. ipick < nkink-1)
            ipick = ipick + 1
            tmp_sum = tmp_sum + T_minus(ipick)/xnorm
            if (tmp_sum > rn3) done = .true.
          end do
          if (.not.done) ipick = nkink
          ! in this step, we make judgement whether ipick is acceptable or not !  
          ipick_plus  = modulo(ipick-1 +1, nkink)+1
          ipick_minus = modulo(ipick-1 -1, nkink)+1

          ! note: keep this as two separate 'ifs' to prevent execution of the
          ! second expression when the first is not true
          if (j_old(ipick_minus) .ne. j_old(ipick_plus)) then
            if (get_ham(j_old(ipick_minus), j_old(ipick_plus)) .ne. 0) then 
              k_rem_mid = k_rem_mid + 1
              ! remove the chosen kink (ipick)
              jtrial_rem(    1:ipick-1) = j_old(      1:ipick-1)
              jtrial_rem(ipick:nkink-1) = j_old(ipick+1:nkink  )
              nkink = nkink-1 
              big_rem = big
              ! nkink has now been changed, so we need to recalculate stuff
              if (nkink == 1) then 
                ! No removal weight for 1 state case!  
                xnorm_rem = (nkink+1) * Z_total  
              else
                do k1=1,nkink 
                  T_rem_minus(k1) = beta_x*(Ha_diag(j_old(k1)) - energy_hartree_fock)
                end do 
              end if
              xnorm_rem = (nkink+1)*Z_total + sum(exp(T_rem_minus(1:nkink)-big_rem))
              call weight_kink(nstate,jtrial_rem,rho_big_rem)
              rn4 = random()
              ! Do we we accept this removal?
              if (log(abs(xnorm/xnorm_rem)) + big - big_rem - &
                2*beta_x*(Ha_diag(j_old(ipick)) - energy_hartree_fock) + rho_big_rem-rho_big > log(rn4)) then
                k_rem_s=k_rem_s+1 
                ! must satisfy two criteria
                ! we can move the IF statement but can not Remove to do correct sampling! 
                ! So in this step, we must let the removal satisfy criteria.  
                ! set ipick_plus, ipick_minus first , this can be done in first step !   
                ! accept removal
                j_old(1:nkink) = jtrial_rem(1:nkink)
              else
                ! revert removal instead
                nkink = nkink+1
              end if
            end if 
          end if
        ! insert kink ?
        else if (rn3 < xnorm_add/xnorm .and. nkink < nmax_kinks) then 
          ! need to locate ipick_n and ipick_l  
          k_insert=k_insert+1 
          done = .false.
          ipick_l = 1
          ipick_n = 1
          tmp_sum = 0.d0
          do while (.not.done)
            tmp_sum = tmp_sum + abs(T_add(ipick_l))/xnorm
            if(tmp_sum >= rn3) then
              done = .true.
            else
              ipick_n = ipick_n+1
              if(ipick_n > nkink+1) then
                ipick_n = 1
                ipick_l = ipick_l+1
              end if
              if(ipick_l > nstate) then
                ipick_n = 0
                ipick_l = 0
                done = .true.
              endif
            endif
          end do
          if(.not. done) then
            ipick_n = nkink+1
            ipick_l = nstate
          end if
          ! is ipick acceptable?
          ipick_plus  = modulo(ipick_n-1 +1, nkink)+1
          ipick_minus = modulo(ipick_n-1 -1, nkink)+1
          if (ipick_l .ne. j_old(ipick_minus) .and. &
              ipick_l .ne. j_old(ipick_plus )) then   
            if (get_ham(ipick_l, j_old(ipick_minus)) .ne. 0 .and. &
                get_ham(ipick_l, j_old(ipick_plus )) .ne. 0) then 
              k_insert_mid = k_insert_mid + 1
              nkink = nkink + 1
              jtrial_insert(        1:ipick_n-1) = j_old(      1:  ipick_n-1)
              jtrial_insert(      ipick_n      ) = ipick_l
              jtrial_insert(ipick_n+1:    nkink) = j_old(ipick_n:    nkink-1)

              big_insert = big
              do k1=1,nkink 
                T_insert_minus(k1) = beta_x*(Ha_diag(jtrial_insert(k1)) - energy_hartree_fock)
                if (big_insert < T_insert_minus(k1)) then 
                  big_insert=T_insert_minus(k1) 
                end if 
              end do 
              Z_total_insert = sum_exp_T_add / exp(big_insert)
              xnorm_insert = sum(exp(T_insert_minus(1:nkink)-big_insert)) + &
                             (nkink+1) * Z_total_insert
              call weight_kink(nstate,jtrial_insert,rho_big_insert) 
              rn5 = random()
              if (log(abs(xnorm/xnorm_insert)) + big - big_insert + &
                  2*beta_x*(Ha_diag(ipick_l) - energy_hartree_fock) + rho_big_insert - rho_big > log(rn5) .and.&
                  nkink <= nmax_kinks) then
                k_insert_s=k_insert_s+1
                ! must include the two criteria !  
                j_old(1:nkink) = jtrial_insert(1:nkink)
              else
                nkink=nkink-1 
              end if 
            end if  
          end if  !criteria 2
        end if   !criteria 1 
      end do !imulti
    else ! ipass.lt.10000  

    ! dynamically change the number of kinks.

    ! removal or insertion kinks begin
    do imulti=1,10
      ! refresh rho
      ! rho_add(:nkink+1, :nstate) = 0.d0
      ! rho_minus(:nkink)          = 0.0

      ! Add a kink temporarily to see what happens

      ! First, get largest energy from kink removal
      call removal_kink (nstate, j_old, rho_minus)
      big_rho = maxval(rho_minus(1:nkink))

      ! Now see if some kink insertion has an even higher energy,
      ! and if so, save that instead
      ! This function modifies big_rho if necessary, and uses it to return
      ! an already "normalized" xnorm_add
      call insertion_kink(nstate,j_old,rho_add,big_rho,xnorm_add)
      ! use the potentially new big_rho to get xnorm from rho_minus too
      xnorm = xnorm_add + sum(exp(rho_minus(1:nkink)-big_rho))
      ! Now xnorm contains a sum of possible energy changes by kink changes,
      ! xnorm_add is the sum of changes by inserting a kink
      if (xnorm .ne. 0) then
        ! Do we remove or insert a random kink?
        rn3 = random()
        if (rn3 > xnorm_add/xnorm .and. nkink > 1) then
          ! Try to remove the kink - find out which to remove based on rn3
          tmp_sum = xnorm_add/xnorm
          done = .false.
          ipick = 0
          do while (.not.done .and. ipick < nkink-1)
            ipick = ipick + 1
            tmp_sum = tmp_sum + abs(exp(rho_minus(ipick)-big_rho))/xnorm
            if(tmp_sum > rn3) then
              done  = .true.
            end if
          end do
          if (.not.done) ipick = nkink

          ! remove the chosen kink (ipick)
          jtrial_rem(    1:ipick-1) = j_old(      1:ipick-1)
          jtrial_rem(ipick:nkink-1) = j_old(ipick+1:nkink  )
          nkink = nkink-1
          ! calcuate the energies associated with removing and/or adding another
          ! kink from new configuration
          call removal_kink(nstate-1, jtrial_rem, rho_minus)
          big_rho_rem = maxval(rho_minus(1:nkink))

          call insertion_kink(nstate-1, jtrial_rem, rho_add, big_rho_rem, xnorm_rem_add)
          xnorm_rem = xnorm_rem_add + sum(exp(rho_minus(1:nkink) - big_rho_rem))

          ! Do we we accept this removal - based on a random number and the
          ! energy difference of the new configuration and possible other new
          ! configurations
          rn4 = random()
          if (log(abs(xnorm/xnorm_rem)) + big_rho - big_rho_rem > log(rn4)) then
            ! accept removal
            j_old(1:nkink) = jtrial_rem(1:nkink)
          else
            ! revert removal
            nkink = nkink+1
          end if

        else if (rn3 < xnorm_add/xnorm .and. nkink < nmax_kinks) then 
          ! Try to insert a kink
          ! Now pick one kink to insert, based on a random number (rn3)
          ! ipick_n : position at which kink is to be inserted
          ! ipick_l : which state should the new kink "go to"
          done = .false.
          ipick_l = 1
          ipick_n = 1
          tmp_sum = 0.d0
!define OLD_INS
#ifdef OLD_INS
          do while (.not. done .and. ipick_l <= nstate)
            tmp_sum = tmp_sum + abs(exp(rho_add(ipick_n,ipick_l)-big_rho))/xnorm
            if(tmp_sum >= rn3) then
              done = .true.
            else
              ipick_n = ipick_n + 1
              if(ipick_n > nkink + 1) then
                ipick_n = 1
                ipick_l = ipick_l + 1
              end if
            endif
          end do
#else
          do while (.not. done .and. ipick_n <= nkink + 1)
            tmp_sum = tmp_sum + abs(exp(rho_add(ipick_n,ipick_l)-big_rho))/xnorm
            if(tmp_sum >= rn3) then
              done = .true.
            else
              ipick_l = ipick_l + 1 
              if( ipick_l > nstate) then 
                 ipick_l = 1
                 ipick_n = ipick_n + 1
              end if    
            end if
          end do
#endif
          ! prevent roundoff error
          if(.not. done) then
            ipick_n = nkink+1
            ipick_l = nstate
          end if
       
          ! Insert kink
          nkink = nkink + 1
          jtrial_insert(        1:ipick_n-1) = j_old(      1:  ipick_n-1)
          jtrial_insert(      ipick_n      ) = ipick_l
          jtrial_insert(ipick_n+1:    nkink) = j_old(ipick_n:    nkink-1)

          ! Now get energies for removal/insertion of kinks from that new configuration
          call removal_kink(nstate+1, jtrial_insert, rho_minus)
          big_rho_insert = maxval(rho_minus(1:nkink))

          call insertion_kink(nstate+1,jtrial_insert,rho_add,big_rho_insert,xnorm_insert_add)
          xnorm_insert = xnorm_insert_add + sum(exp(rho_minus(1:nkink)-big_rho_insert))

          ! Do we accept this insertion - based on a random number and the
          ! energy differnce of the new configuration  and possible other new
          ! configurations
          rn5 = random()
          if (log(abs((xnorm/xnorm_insert))) + big_rho - big_rho_insert > log(rn5) .and. &
              nkink <= nmax_kinks) then
            ! accept insertion
            j_old(1:nkink) = jtrial_insert(1:nkink)

            if (.not. any(list_of_find(1:nappear_find) == ipick_l) .and. &
              nappear_find < nmax_states) then 
              nappear_find = nappear_find+1 
              list_of_find(nappear_find)    = ipick_l
#if 1
              do k2= 1, nstate
                Ha_temp(nappear_find,k2)=get_ham(list_of_find(nappear_find),k2)
              end do
#else
              call get_hams(list_of_find(nappear_find), &
                            list_of_find(nappear_find), &
                            1, nstate, &
                            Ha_temp, &
                            nappear_find, &
                            nappear_find, &
                            1, nstate)
#endif
            end if
          else
            ! revert insertion
            nkink = nkink-1
          end if
        end if ! for initial criteria (rn3.gt.xnom_add/xnorm)
      end if ! for xnorm.ne.0 case
    enddo ! imulti
    end if 
    ! removal or insertion kinks end
    if (nkink.eq.1) then
      ip1=ip1+1
    else if (nkink.eq.2) then
      ip2=ip2+1
    else if (nkink.eq.3) then
      ip3=ip3+1
    else if (nkink.eq.4) then
      ip4=ip4+1
    end if
    timer_evol_part2 = timer_evol_part2 + (mytime() - time_start_part2)/1000.d0

    !------------------------End Second Part ---------------------
    ! The third part of code begin

    time_start_part3 = mytime()
    ! diagonalize the selected states, update rho
    if (turn_diag == 1) then
      do j=1, max(1,nkink)
        times_appear(j_old(j)) = times_appear(j_old(j))+1
        ! if there wasn't something in list_of_states_init, add it to the end
        if (.not. any(list_of_states_init(1:nappear_init) == j_old(j)) .and. &
            nappear_init < nmax_states) then
          nappear_init = nappear_init+1
          list_of_states_init(nappear_init) = j_old(j)
        endif
      enddo

      if( ipass == first_diag .or. &
         (ipass <= npass_diag .and. modulo(ipass, diag_every) == 0) .or. &
          (nkink .ge. nmax_kinks)) then
        write (*,*) "KN nappear_init 1", nappear_init
        if(ipass >= npass_diag.and.nkink.ge.nmax_kinks) then
          npass=ipass+1000
          write(102,*) ' new max passes ',nkink, npass
        endif

        if (my_id.eq.0) then
        !  write(102,*) ipass, nappear_init,nkink ,"npass,nappear_init,nkink"
       
          do k1=1,nkink
         !   write (333,*) j_old(k1)
          end do
       !   write(333, *) "one group"
        end if
      

        if(ipass == 100) then
          threshold=1
        else if(ipass == 500) then
          threshold=1
        else if(ipass == 5000) then
          threshold=1
        else
          threshold=1
        endif
  
        nappear=0
        list_of_states=0
        ! Filter from list_of_states_init() to list_if_states() 
        write (*,*) "KN nappear", nappear, nappear_init, nmax_kinks
        do k1 = 1, nappear_init
          if (nappear < nmax_states .and. times_appear(list_of_states_init(k1)) >= threshold) then
            nappear = nappear+1
            list_of_states(nappear) = list_of_states_init(k1)
          endif
        enddo
        call ham%update_cache_ham(ipass, nappear, list_of_states)
        call ham%output_info_ham(ipass)
        call output_ham_vinfo(ipass, 100)

        do k1=1,nstate
            iused(k1)=.false.
        enddo
        write (*,*) "KN nappear_init 2", nappear_init

!     list_of_states span all the dimension from 1 to nstate, the same dimension with iused, jused 

       do k2=1,nappear
        iused(list_of_states(k2)) = .true.
        jused(list_of_states(k2))=k2
       enddo

        ! List_of_happen() is a list of all the states that happened in diagonalizations 
        do j=1, nappear
          k=0
          done=.false.
          do while (.not.done.and.k.lt.nappear_happen)
            k=k+1
            if ( list_of_states(j).eq.list_of_happen(k)) then
              done=.true.
            end if
          end do
          if(.not.done) then
            nappear_happen = nappear_happen+1
            list_of_happen(nappear_happen)=list_of_states(j)
          end if
        end do
        ! nappear_happen is increasing
        list_of_happen_sort(1:nappear_happen)= list_of_happen( 1:nappear_happen) 
        ! sorting the list_of_happen_sort, for future use
        call PIKSRT ( nappear_happen, list_of_happen_sort)  
        if( my_id.eq.0) then 
          write(102,*) ipass,nappear,nkink,nappear_happen,"npass,Real nappear,nkink,nappear_happen"

           write(103,*) ipass,nappear,nappear_happen,"ipass,sub_size, total"
          write(1936,*) nappear_happen, ipass,"nappear_happen, ipass"
         ! write(1936,*) call_size, ipass, "call size ipass=100" 
        end if 


        

        ! set up matrix for diagonalization, rho_mat is sub_mat to diagonalize.
        do k1=1,nappear
          do k2=1,nappear
            rho_mat(k1,k2)=get_ham(list_of_states(k1),list_of_states(k2))
          enddo
        enddo

        allocate(rho_evec(nappear,nappear))
        allocate(w(nmax_states),work(nmax_states*3))
         w = 0
         work = 0
        ! copy rho_mat to rho_evec to diagonalize, so that the eigenvector will fill the
        !   matrix rho_evec, to be the coefficient to the next step.
        rho_evec(1:nappear, 1:nappear) = rho_mat(1:nappear,1:nappear)
 
        ! if( ipass.ge.1) then 
       ! write(*,*) gen_new(1902,39965), gen_new(39965,1902), "1902-39965"
        !  do  k2=1,nappear
       !       do k4=1, nappear 
          !  write(1934,*) rho_evec(k2,k4), k2,k4, "rho_evec Before",nappear
        !      end do 
        !    end do 
        !  write(1934,*)  list_of_states(1:nappear), nappear , "nappear" ,ipass
        !  write(1934,*) list_of_happen(1:nappear_happen),nappear_happen, "nappear_happen",ipass 
        !  end if 

        time_start_ext_diag = mytime()

      
       call DSYEV( 'V', 'U', nappear, rho_evec,nappear, W, WORK,nappear*3, INFO )
        timer_ext_diag = timer_ext_diag + (mytime() - time_start_ext_diag)/1000.d0
        !!extremely important !  the eigenvector will overwrite the rho_mat to rho_evec

       
        ! Renormalization of the rho_evec
        do k1=1,nappear
          sum3 = sum(rho_evec(:nappear, k1)**2)
          rho_evec(1:nappear,k1) = rho_evec(1:nappear,k1)/dsqrt(sum3)
        enddo
      !if( my_id.eq.0) then 
      !  write(*,*) sum3, nappear, nmax_states,"normalization",ipass
      !  write(*,*) INFO, "INFO"
      !end if 


        ! test when ipass=75: 
        if (my_id.eq.0.and.ipass.eq.75) then  
          do k1=1,nappear 
            do k2=1,nappear   
              sum5=0 
              do k3=1,nappear 
                sum4=0  
                do k4=1,nappear 
                  sum4=sum4+ rho_mat(k3,k4) * rho_evec(k4,k2)
                end do
                sum5=sum5+rho_evec(k3,k1)*sum4
              end do   
              if( k2==k1) then
                write(*,*)  sum5,w(k1), k1,"w(k1)",ipass 
              else 
                write(*,*)  sum5,"should be zero Year",k1,k2,ipass
              end if 
            end do 
          end do   
        end if

        do k1=1, nappear 
          do k2=1,nstate 
            sum5=0.d0
            if (.not.iused(k2)) then 
              do k4=1,nappear 
                sum5=sum5+rho_evec(k4,k1)*get_ham(list_of_states(k4),k2) 
              end do
            else if ( k2.eq.list_of_states(k1)) then 
              sum5=w(k1)   
            else 
              sum5=0.d0 
            end if
            ! use Ha_local to temporary save all the changed  elements in this diagonalization.  
            Ha_local(k1,k2)=sum5
          end do 
        end do  

        do k1 = 1 , nappear
          do k2 = 1 , nstate
            call set_ham(list_of_states(k1), k2, Ha_local(k1, k2))
          end do
        end do  

        do k1= 1, nappear 
          Ha_diag(list_of_states(k1))=w(k1) ! update eigenvalues! awesome
          !this step only update limit number of vector ! 
        end do
       
        ! the rest of diagonal elements do not change
        times_appear=0

        ! Disable optimization in rho_1
        ground_state_index = 0
  
        nkink    = 1
        j_old(1) = 1
        maxtmp = rho_1(j_old(1),j_old(1))

       
        do k2=2, nstate
          if(rho_1(k2,k2) > maxtmp) then
            j_old(1) = k2
            maxtmp = rho_1(k2, k2)
          endif
        enddo
        ! in here , j_old(1) get updated after diagonalization.  
        ground_state_index = j_old(1) 
        ! we update ground_state_index after every diagonalization

       ! write(1935,*) ground_state_index , "ground state index is ",ipass 
        
        ! 2014-04-08: This is the bulk of the computation of part III
        time_start_debug = mytime()
        do k1= 1, nstate
          Ha_low(k1) = get_ham(j_old(1), k1)
        end do 

        timer_debug = timer_debug + (mytime() - time_start_debug)/1000.d0
        ! the ground state and excitation is dense, that's why it is time consuming.  

        rho_mat(:nappear, :nappear)  = 0.d0
        rho_evec(:nappear, :nappear) = 0.d0
        list_of_states(1:nappear)    = 0
        nappear                      = 0
        list_of_states_init(1:nappear_init) = 0 
   
        !list_of_states               = 0 
        write (*,*) "KN nappear_init 3", nappear_init
        nappear_init                 = 0   

        list_of_find(1:nappear_find) = 0 
        ! after reset, we need list_of_find save ground state
        nappear_find = 1 
        list_of_find(1) = j_old(1) 

      
        Ha_temp(1,1:nstate) = Ha_low(1:nstate)
        
        deallocate(rho_evec) 
        deallocate(w)
        deallocate(work) 
     
      end if
    end if ! turn_diagonalization
    ! calculate the energy

    ! track j_old 
#if 0
    if( ipass.gt.100000) then 
      if (my_id.eq.0) then 
  
     
        do k1=1,nkink
         ! write(1935,*)  j_old(k1), "j_old",k1,nkink,ipass
         ! write(1935,*) get_ham( j_old(k1),j_old(k1)) , "Diagonal element"

          do k2=1,nstate
            if( abs(get_ham(j_old(k1),k2)).gt.1.d-14) then 
            !  write(1935,*) get_ham(j_old(k1),k2), "off_diagonal", j_old(k1),k2 
            end if 
          end do  
        end do 

      end if   

    end if  

#endif
    call sfunc_calc(j_old,nstate,nkink,p,sfunc_old,sign_sfunc)
    timer_evol_part3 = timer_evol_part3 + (mytime() - time_start_part3)/1000.d0
 
    ! ------------------- Third Part END, Update RHO Finish ----------------

    time_start_part4 = mytime()
    !part4 main calculating energy! 
    q_est=0                   ! Q_estimator initial value
    if (nkink.eq.1) then
      q_est=0.0
      q_est=rho_1(j_old(1),j_old(1))**p
      if (q_est.gt.0.0) then
        sign_final=1
      else
        sign_final=-1
      end if
      sf_both=rho_dev_1(j_old(1),j_old(1))/abs(rho_1(j_old(1),j_old(1)))
      accume_e    = accume_e    + sf_both
      accume_sign = accume_sign + sign_final
      accume_kink = accume_kink + nkink
!print*, sf_both , ipass 
    else ! nkink>2
      do i=1,nkink
        xv(i) = rho_1    (j_old(i),j_old(i))
        xvv(i)= rho_dev_1(j_old(i),j_old(i))
      end do
  
      call efunc_calc_rwh(xv,xvv,nkink,p,sfunc,sign_sfunc,sf_both)
      sum_t=0.0
      q_est=1.0
      do i=1,nkink
        ! calcuate a function f
        iplus=1+modulo(i,nkink)
        sum_t=sum_t+rho_dev_1(j_old(i),j_old(iplus))/ rho_1(j_old(i),j_old(iplus))/p
        sign_sfunc=sign_sfunc*sign(1.d0,rho_1(j_old(i),j_old(iplus)))
        q_est=q_est*rho_1(j_old(i),j_old(iplus))
   
      end do
    !    if(ipass.eq.228) then 
!print*,sum_t , ipass 
     !   end if 
      sign_sfunc=dfloat(nint(sign_sfunc))
      sfunc=exp(sfunc)

      q_est=q_est*sfunc*sign_sfunc

      accume_e    = accume_e    + sf_both+sum_t
      accume_sign = accume_sign + sign_sfunc
      accume_kink = accume_kink + nkink

      sign_final=nint(sign_sfunc)
 !     do  i=1,nkink
 !       iplus=1+modulo(i,nkink)
 !       if(rho_1(j_old(i),j_old(iplus)) .lt. 0) sign_final=-sign_final
 !     enddo
      sf_both=0
    endif ! nkink==1
 
    ! diagonalize
    if (modulo(ipass, diag_every) == 0 .and. nkink == 1) then 
      call truncation_Q2(nstate,j_old, E_tr)  
      if (modulo(ipass-1,ncomb) == 0 .and. my_id == 0) then 
      !  write(1939,*) E_tr, nkink,ipass , "E_tr,ipass"
      end if
    end if 
    timer_evol_part4 = timer_evol_part4 + (mytime() - time_start_part4)/1000.d0
  
    ! Some output from time to time
    if (modulo(ipass-1,ncomb) == 0 .and. my_id == 0) then
      if(ipass .eq. 1) then
        open(unit=40,file=trim(adjustl(outdir))//'/'//'kink.out',status='unknown')
      else
        open(unit=40,file=trim(adjustl(outdir))//'/'//'kink.out',status='old',access='append')
      endif
      write(40,5000) accume_e/ncomb,accume_sign/ncomb,accume_kink/ncomb,ipass
      write(*,*) ipass, accume_e/ncomb, nkink, ham_size(), ham_size()/(1.d-2*nstate*nstate)
5000  format(3e20.10,2i10,e20.10)
      close(40)
      accume_e=0.d0
      accume_sign=0.d0
      accume_kink=0.d0
     ! write(1950,*) real(k_insert)/real(10*ipass), real(k_insert_mid)/real(k_insert+1), &
     !               real(k_insert_s)/real(k_insert+1)
     ! write(1949,*) real(k_rem)/real(10*ipass), real(k_rem_mid)/real(k_rem+1), &
     !               real(k_rem_s)/real(k_rem+1)
    endif
    timer_evol = timer_evol + (mytime() - time_start_evol_it)/1000.d0
    ! timing output
    if (modulo(ipass, diag_every) == 0 .and. my_id == 0) then
      ttimer_evol       = ttimer_evol       + timer_evol
      ttimer_evol_part1 = ttimer_evol_part1 + timer_evol_part1
      ttimer_evol_part2 = ttimer_evol_part2 + timer_evol_part2
      ttimer_evol_part3 = ttimer_evol_part3 + timer_evol_part3
      ttimer_evol_part4 = ttimer_evol_part4 + timer_evol_part4
      ttimer_ext_diag   = ttimer_ext_diag   + timer_ext_diag
      ttimer_debug      = ttimer_debug      + timer_debug
      write(1011,*) ipass, mytime()/1000.d0-time_start/1000.d0, &
                           mytime()/1000.d0-time_start_evol/1000.d0, &
                    timer_evol, timer_evol_part1,  timer_evol_part2,  timer_evol_part3,  timer_evol_part4, timer_ext_diag, timer_debug, &
                    ttimer_evol, ttimer_evol_part1, ttimer_evol_part2, ttimer_evol_part3, ttimer_evol_part4, ttimer_ext_diag, ttimer_debug, getRSS()
      timer_evol       = 0.d0
      timer_evol_part1 = 0.d0
      timer_evol_part2 = 0.d0
      timer_evol_part3 = 0.d0
      timer_evol_part4 = 0.d0
      timer_ext_diag   = 0.d0
      timer_debug      = 0.d0
      FLUSH(1011)
    end if

  end do ! ipass


  
  

#if 0
! calculate Q_0,Q_2,Q_3
! calculate Q_0_bar, Q_2_bar,Q_3_bar,
!previously we use this part to check whether QMC is working correctly ! QMC and theoretical result should match!

      output_count=0
      magn_run=0.0
      magnetization_ave=0.0
      magnetization_nave=0.0
      big_rho=0.d0
      do nkink=1,nstate
        big_rho=max(big_rho,p*log(abs(rho_1(nkink,nkink))))
      enddo
      do nkink=1,3

!  when nkink=4 case
      output4_count=0
      magn_run=0.0
      magnetization4_ave=0.0
      magnetization4_nave=0.0



       if( nkink.eq.4) then

      do i=1,nstate
       do j=1,nstate
        do k=1,nstate
         do l=1,nstate
          if(i.ne.j.and.j.ne.k .and.k.ne.l.and.l.ne.i) then
           if(rho_1(i,j) .ne. 0.d0 .and. rho_1(j,k) .ne. 0.d0  .and. &
               rho_1(k,l) .ne. 0.d0 .and. rho_1(l,i) .ne. 0.d0) then


            output4_count=output4_count+1

            j_old(1)=i
            j_old(2)=j
            j_old(3)=k
            j_old(4)=l
            call sfunc_calc(j_old,nstate,nkink,p,sfunc,sign_sfunc)


            magnetization4_ave=magnetization4_ave+ &
                nint(sign_sfunc)*rho_1(i,j)*rho_1(j,k)*rho_1(k,l)*rho_1(l,i)*exp(sfunc-big_rho)*p/4

            magnetization4_nave=magnetization4_nave+&
              abs(rho_1(i,j)*rho_1(j,k)*rho_1(k,l)*rho_1(l,i))*exp(sfunc-big_rho)*p/4




            if( sign_sfunc.lt.1) then
             write(*,*) "negative value", sign_sfunc
            end if


           end if
          end if
         end do
        end do
       end do
      end do

       end if

! when nknik=4 finish

! when nkink=3 start

    if (nkink.eq.3) then
      output3_count=0
      magn_run=0.0
      magnetization3_ave=0.0
      magnetization3_nave=0.0
      j_old(1)=0
            j_old(2)=0
            j_old(3)=0
            j_old(4)=0



      do i=1,nstate
       do j=1,nstate
        do k=1,nstate
!         do l=1,nstate
          if(i.ne.j.and.j.ne.k .and.k.ne.i) then
           if(rho_1(i,j) .ne. 0.d0 .and. rho_1(j,k) .ne. 0.d0  .and. &
               rho_1(k,i) .ne. 0.d0) then
            output3_count=output3_count+1
            j_old(1)=i
            j_old(2)=j
            j_old(3)=k
!            j_old(4)=rho_1(l,l)
            call sfunc_calc(j_old,nstate,nkink,p,sfunc,sign_sfunc)

            magnetization3_ave=magnetization3_ave+&
               nint(sign_sfunc)*rho_1(i,j)*rho_1(j,k)*rho_1(k,i)*exp(sfunc-big_rho)*p/3
            magnetization3_nave=magnetization3_nave+ &
               abs(rho_1(i,j)*rho_1(j,k)*rho_1(k,i))*exp(sfunc-big_rho)*p/3

            if( sign_sfunc.lt.1) then
             write(*,*) "negative value", sign_sfunc
            end if
           end if
          end if
         end do
        end do
       end do
!      end do
    end if

! When nkink=3 finish
! when nkink=2 start

   if (nkink.eq.2) then
      output2_count=0
      magn_run=0.0
      magnetization2_ave=0.0
      magnetization2_nave=0.0

      j_old(1)=0
      j_old(2)=0
      j_old(3)=0
      j_old(4)=0

      do i=1,nstate
       do j=1,nstate
!        do k=1,nstate
!         do l=1,nstate
          if(i.ne.j.and.j.ne.i) then
           if(rho_1(i,j).ne. 0.d0 .and. rho_1(j,i) .ne. 0.d0) then

         !   write(*,*) rho_1(i,j) , i,j,"rho_1(i,j) nkink=2"
            output2_count=output2_count+1

            j_old(1)=i
            j_old(2)=j

!            j_old(3)=rho_1(k,k)
!            j_old(4)=rho_1(l,l)
            call sfunc_calc(j_old,nstate,nkink,p,sfunc,sign_sfunc)


            magnetization2_ave=magnetization2_ave+ &
                nint(sign_sfunc)*rho_1(i,j)*rho_1(j,i)*exp(sfunc-big_rho)*p/2

            magnetization2_nave=magnetization2_nave+ &
               abs(rho_1(i,j)*rho_1(j,i))*exp(sfunc-big_rho)*p/2

            if( sign_sfunc.lt.1) then
             write(*,*) "negative value", sign_sfunc
            end if
           end if
          end if
         end do
        end do
!       end do
!      end do
     end if

! when nkink=2 finish
! when nkink=1 start

    if ( nkink.eq.1) then
      output0_count=0
      magn_run=0.0
      magnetization0_ave=0.0
      magnetization0_nave=0.0

      j_old(1)=0
      j_old(2)=0
      j_old(3)=0
      j_old(4)=0

      do i=1,nstate
        if (rho_1(i,i).ne.0.0) then
          output0_count=output0_count+1
          magnetization0_ave=magnetization0_ave+exp(p*log(abs(rho_1(i,i)))-big_rho)
          magnetization0_nave=magnetization0_nave+exp(p*log(abs(rho_1(i,i)))-big_rho)
        end if
      end do
    end if

! when nkink=1 finish
! calculate total
  end do
  magnetization_ave=magnetization4_ave+ &
  magnetization3_ave+magnetization2_ave+ &
  magnetization0_ave

  output_count=output4_count+ output3_count+ output2_count+output0_count

  write(*,*) magnetization4_ave,magnetization3_ave,magnetization2_ave,&
             magnetization0_ave ,"Q4-Q0"
  write(*,*) magnetization4_nave,magnetization3_nave,magnetization2_nave,&
             magnetization0_nave , "|Qi| "
  write(*,*)output_count,magnetization_ave/magnetization_nave
  write(*,*)output4_count,output3_count,output2_count,output0_count
#endif

  close(102)
  close(333)
  !deallocate(rho_evec)
  !deallocate(Ha)
  deallocate(j_old)

  deallocate(xv)
  deallocate(xvv)

  deallocate(jtrial_rem)
  deallocate(jtrial_insert)

  deallocate(j_insert_add)
  deallocate(j_insert_minus)
  if (calc_rho .ne. 0) then
!    call save_rho()
  endif

  deallocate(rho_add)
  deallocate(rho_minus)


  call MPI_finalize(ierr)
CONTAINS

! Initiate checkpointing of all data
subroutine checkpoint()
  call checkpoint_start()
  call checkpoint_ham()
  call checkpoint_timing()
  call checkpoint_random()
  call checkpoint_end()
end subroutine

! Initiate restart using checkpointed data
subroutine restart()
  use checkpointing
  if (.not. restart_start()) then
    write(*,*) "Error restarting"
    stop
  end if
  call restart_timing()
  call restart_ham()
  call restart_random()
  call restart_end()
end subroutine

! Save only the Hamiltonian to a file
subroutine save_ham_init()
  call checkpoint_start_to_file(trim(adjustl(outdir))//'/ham_init.h5')
  call checkpoint_ham()
  call checkpoint_end()
end subroutine

! Load only the Hamiltonian from a file
function load_ham_init()
  use checkpointing
  implicit none
  logical load_ham_init
  integer herr
  load_ham_init = .false.
  ! disable hdf5 error printing (we are probing here)
  call h5eset_auto_f(0, herr)
  if (.not. restart_start_from_file(trim(adjustl(indir))//'/ham_init.h5')) then
    write(*,*) "Ham not found saved as "//trim(adjustl(indir))//'/ham_init.h5'
    if (.not. restart_start_from_file(trim(adjustl(outdir))//'/ham_init.h5')) then
      write(*,*) "Ham not found saved as "//trim(adjustl(outdir))//'/ham_init.h5'
      call h5eset_auto_f(1, herr)
      return
    else
      write(*,*) "Loading ham from "//trim(adjustl(outdir))//'/ham_init.h5'
      load_ham_init = .true.
    end if
  else
    write(*,*) "Loading ham from "//trim(adjustl(indir))//'/ham_init.h5'
    load_ham_init = .true.
  end if
  ! enable hdf5 error printing
  call h5eset_auto_f(1, herr)
  call restart_ham()
  call restart_end()
end function

! use module for sfunc_calc to avoid memory allocation every time sfunc_calc
! is called (which is very often)

subroutine test_exchange ( k_1,k_2) 
 implicit none 
  integer, intent(in) :: k_1,k_2 
  integer:: k1,k2,tmp 

  k1=k_1
  k2=k_2 
  
  if( k1>k2) then 
  
   tmp=k2 
   k2=k1 
   k1=tmp 
  write(*,*) k_1,k_2, "suc exchange",k1,k2 
  else 
  write(*,*) k_1,k_2,"notexchange",k1,k2 
  end if 
end subroutine  


subroutine sfunc_init(n)
  use sfunc_mod
  implicit none
  integer, intent(in) :: n
  ! initialize cache of sfunc
  allocate(ig(n+1))
  allocate(itaken(n+1))
  allocate(wvec(n+1))
  allocate(fvec(0:n))
end subroutine

! Calculate S: one part of the weight w
! Has side-effects: fect, wvec
subroutine sfunc_calc(iarr,iarr_size,n,part_n,sfunc,sign_sfunc)
  use kink
  use sfunc_mod
  implicit none

  integer, intent(in) :: iarr_size
  integer, intent(in) :: iarr(iarr_size)
  integer, intent(in) :: n
  integer(int8), intent(in) :: part_n

  real(dp), intent(out) :: sfunc, sign_sfunc

  integer k1, k2, k3, k4, m
  real(dp) tol, big_x, gsum, sf, sum3
  real(dp) :: rho(n)
  parameter (tol=1.d-12)

  do k1=1,n
    rho(k1) = rho_1(iarr(k1), iarr(k1))
  end do
  if (n.eq.1) then
    sfunc=part_n*log(rho(1))
    sign_sfunc=1
  else
    m=0
    itaken(:n) = 0    ! equivalent to use a do loop to default itaken(n)
    big_x  = 0.d0
    do k1=1,n
      if(itaken(k1) .eq. 0) then
        m=m+1
        ig(m)=1
        wvec(m)=rho(k1)   !in here, we get fresh wvec by rho_1 based on nkink(input)
        big_x=max(big_x,abs(wvec(m)))
        itaken(k1)=1
        do k2=k1+1,n
          if(itaken(k2) .eq. 0) then
            if(abs(rho(k2) - wvec(m)) .lt. tol) then
              itaken(k2)=1
              ig(m)=ig(m)+1
            endif
          endif
        enddo
      endif
    enddo
    sf = 0.d0
    do k1=1,m
      fvec(0)=(wvec(k1)/big_x)**(part_n-1)
      do k2=1,m
        if (k1 .ne. k2) then
          fvec(0)=fvec(0)/(wvec(k1)-wvec(k2))**ig(k2)
        endif
      enddo

      do k2=1,ig(k1)-1
        sum3=0.d0
        do k3=0,k2-1
          gsum=0.d0
          do k4=1,m
            if(k1 .ne. k4) then
              gsum = gsum+ig(k4)/(wvec(k1)-wvec(k4))**(k3+1)
            endif
          enddo ! k4
          sum3 = sum3 + fact(k2-1) / fact(k2-1-k3) * fvec(k2-1-k3) * (-1)**k3 * &
                      ((part_n-1)/wvec(k1)**(k3+1)-gsum)
        enddo ! k3
        fvec(k2)=sum3
      enddo ! k2
      sf=sf+fvec(ig(k1)-1)/fact(ig(k1)-1)

    enddo ! k1
    sfunc=(part_n-1)*log(big_x)+log(abs(sf))

    if(sf .gt. 0) then
      sign_sfunc=1
    else
      sign_sfunc=-1
    endif
  end if
return
end subroutine

! Comment, this is an adaptive method to calculate S
! Please see equation 9-13 in Chemical Physics Letters 362 (2002) Page 550

! in here, we empty Xiaoyao Method to calculate energy

! Read in 'temp.dat': the Hailtonian and write into 'rho_evec(i,j)'

! Dr Hall's code to calcuate energy, it is using in the main code
subroutine efunc_calc_rwh(xvec,dxvec,n, &
           part_n,sfunc,sign_sfunc,efunc)
  use kink


!**********Xiaoyao begin
  implicit none

  real(dp),intent(in)::xvec(n),dxvec(n)
  integer,intent(in)::n
   integer(int8),intent(in)::part_n
  !@integer, intent(in):: part_n
  real(dp), intent(out)::sfunc,sign_sfunc,efunc

  integer m,k1,k2,k3,k4,k5
  real(dp)   tol, big_x,gsum,dgsum,sum3,tmp_sum,sf,dsf
  parameter (tol=1.d-12)
  integer  ig(n),itaken(n)
  real(dp)    wvec(n),dwvec(n)
  real(dp)    fvec(0:n-1),evec(0:n-1)

!**********Xiaoyao achieve
  if (n.eq.1) then
    sfunc = part_n * log(xvec(1))
    efunc = dxvec(1) / xvec(1)
  else
    m=0          ! index of unique kinks
    itaken = 0   ! array of kinks 'we have already seen'
    big_x = 0.d0 ! maximum of xvec
    do k1=1,n
      if(itaken(k1) .eq. 0) then
        m = m + 1
        ig(m) = 1  ! count how many kinks there are for a specify xvec
        wvec(m)  = xvec(k1)  ! copy data from xvec to wvec
        dwvec(m) = dxvec(k1)
        big_x = max(big_x, abs(xvec(k1))) ! remember maximum
        itaken(k1) = 1 ! mark it 'as seen'
        ! loop over the remaining array to find similar entries
        do k2=k1+1,n
          ! maybe not necessary but maybe faster
          if (itaken(k2) .eq. 0) then
            ! check if there is another, similar entry
            if (abs(xvec(k2)-xvec(k1)) .lt. tol) then
              ! if so, count them and mark them as seen
              itaken(k2) = 1
              ig(m) = ig(m)+1
            endif
          endif
        enddo
      endif
    enddo
    ! eq 10: F_l^(p) := d^p/(dx_l^p) * x_l^(_P-1) / prod{k!=l}(x_l-x_k)^g_k
    ! eq 14: E...
    ! note: p!=P  :-)
    sf=0.d0
    dsf=0.d0
    ! loop over distinct values (m, not n)
    do k1=1,m
      fvec(0) = (wvec(k1)/big_x)**(part_n-1)
      evec(0) = (part_n-1.d0)*dwvec(k1)/wvec(k1)
      do k2=1,m
        if(k1 .ne. k2) then
          fvec(0) = fvec(0) / (wvec(k1)-wvec(k2))**ig(k2)
          evec(0) = evec(0)-ig(k2)*(dwvec(k1)-dwvec(k2))/ &
                    (wvec(k1)-wvec(k2))
        endif
      enddo
      evec(0)=evec(0)*fvec(0)/part_n

     ! based on notes, eq(36),
     ! fvec(0)=F_l^(0)= x_l^(_P-1) / prod{k!=l}(x_l-x_k)^g_k
     ! evec(0)=E_l^(0)=F_l^(0)*[(p-1)*x_l^(1)/x_l-sum{k!=l}(g_k*[x_l^(1)-x_k^(1)]/[x_l-x_k]]
      do k3=1,ig(k1)-1
        sum3=0.d0
        tmp_sum=0.d0
        do k4=0,k3-1
          gsum=0.d0
          dgsum=0.d0
          do k5=1,m
            if(k1 .ne. k5) then
              gsum  = gsum+ig(k5)/(wvec(k1)-wvec(k5))**(k4+1)
              dgsum = dgsum+ig(k5)*(k4+1)*(dwvec(k1)-dwvec(k5))/ &
                      (wvec(k1)-wvec(k5))**(k4+2)
            endif
          end do
          sum3 = sum3+fact(k3-1)/fact(k4)/fact(k3-1-k4)* &
                      fvec(k3-1-k4)* &
           (-1)**k4*fact(k4)*((part_n-1)/wvec(k1)**(k4+1)- &
           gsum)
           tmp_sum=tmp_sum+fact(k3-1)/fact(k4)/fact(k3-1-k4)* &
           (fvec(k3-1-k4)* &
           (-1)**k4*fact(k4)/part_n* &
           (-(part_n-1.d0)*dwvec(k1)/wvec(k1)**(k4+2)+dgsum)+ &
           evec(k3-1-k4)* &
           (-1)**k4*fact(k4)*((part_n-1)/wvec(k1)**(k4+1)- &
           gsum))

         ! comment, sum and tmp_sum is the equation (12) and deriation of (12)
         ! over beta.
          !         ===
         ! F_l^(n)=\\ (n-1) * G_l^(m)*F_l^(n-1-m)
         !         //  (m)
         !          ===
        end do
        fvec(k3)=sum3
        evec(k3)=tmp_sum

        !using the adaptive method to get the value fvec and evec
      enddo
      sf=sf+fvec(ig(k1)-1)/fact(ig(k1)-1)
      dsf=dsf+evec(ig(k1)-1)/fact(ig(k1)-1)
      ! sf is S in equation 11 in the paper
      ! dsf is the second term(have not divide S yet) in equation 14 in the paper, the first term
      ! is calculated in the main program
      ! S= sum{l=1,m} F_l^(g_l-1)/(g_l-1)!
      ! dsf=sum{l=0,m}*{[(g_l-1)!]^(-1)}*&
      ! sum{j=0,g_l-2}*combinator(gl-2,j)*&
      !{D_l^(j)*F_l^(g_l-2-j)}+G_l^(j)*E_l^(gl-2-j)}
    enddo
    efunc=dsf/sf
    sfunc=(part_n-1)*log(big_x)+log(abs(sf))

    ! sfunc=log(abs(S)), sign_sfunc is the sign of S, we calculate log(abs(S)) is
    ! for judgement in QMC with log (random number)
    ! However, efunc keep the same, even for big_x term, dsf/sf, the ratio will cancel big_x effect
    if(sf .gt. 0) then
      sign_sfunc=1
    else
      sign_sfunc=-1
    endif
  !       print*,dsf,sf,efunc,sfunc,' values '
  endif
return
end subroutine
! comment, calculate energy in an adaptive method, see eq(14) in
! chemical Physics Letter 362 (2002) page 549-553

! calculate the estimated energy for only zero and two kinks
! In the two kink case, one of kinks is the ground state, the other is
! any other that can form a kink with the ground state.
subroutine truncation_Q2(j_size, j_old, E_tr)
  use kink
  use mpi_data
  use timing
  implicit none

  include 'mpif.h'

  integer, intent(in) :: j_size
  integer, intent(in) :: j_old(j_size)
  real(dp):: big
  real(dp) :: sum_weight,sum_weight_y

  real(dp) :: rho_add_x(nmax_kinks+5,nstate),rho_add_y(nmax_kinks+5,nstate)
  integer:: j_new_add(j_size+1)
  integer:: k1,iplus,iminus,i1,j1,i2
  real(dp) :: energy, sfunc, sign_sfunc
  real(dp) :: sum_node,sum_node_y, big_node
  real(dp) :: e_zero, q_zero, E_tr
  integer :: iplus_2

  big_node = 0
  do k1 = 1, 1
    iplus  = modulo(k1-1 +1, nkink)+1
    iminus = modulo(k1-1 -1, nkink)+1

    ! copy all values of j_old into j_new_add, but leave the k1'th free
    j_new_add(    :k1-1)    = j_old(  :k1-1)
    j_new_add(k1+1:nkink+1) = j_old(k1:nkink)

    ! TODO: DISABLED
    !OMP PARALLEL DO PRIVATE(j1,i1,i2,energy,sfunc,sign_sfunc) FIRSTPRIVATE(j_new_add)
    do j1= 1, nstate
      j_new_add(k1) = j1
      if (j_new_add(iminus) .ne. j1 .and. &
          j_new_add(iplus)  .ne. j1 .and. &
          (nkink+1) <= nmax_kinks) then
        if (abs(rho_1(j_new_add(iminus),j1)*rho_1(j1,j_new_add(iminus))) >= 1.d-20 ) then
        ! from here we calculate rho
          energy = log(p/dfloat(nkink+1))
          do i1 = 1,nkink+1
            i2 = 1 + modulo(i1,nkink+1)
            energy = energy + log_rho_1(j_new_add(i1), j_new_add(i2))
          end do
          ! in the system, nkink=1, we use this to calculate Q2, so we have nkink+1
          !$OMP CRITICAL
          do i=1,nkink+1
            xv(i) = rho_1    (j_new_add(i),j_new_add(i))
            xvv(i)= rho_dev_1(j_new_add(i),j_new_add(i))
          end do
          call efunc_calc_rwh(xv,xvv,(nkink+1),p,sfunc,sign_sfunc,sf_both)
          sum_t=0.0
          do i=1,nkink+1
            ! calcuate a function f
            iplus_2 = 1 + modulo(i, nkink+1)
            sum_t = sum_t + rho_dev_1(j_new_add(i), j_new_add(iplus_2)) / &
                                rho_1(j_new_add(i), j_new_add(iplus_2)) / p
          end do
          !$OMP END CRITICAL
          energy = energy + sfunc
          rho_add_x(k1,j1) = energy        ! rho_add here is the log(rho_add)
          rho_add_y(k1,j1) = sum_t + sf_both
        else
          rho_add_x(k1,j1) = -8800.d0
          rho_add_y(k1,j1) =     0.d0
        end if
      else
        rho_add_x(k1,j1) = -8800.d0
        rho_add_y(k1,j1) =     0.d0  
        ! consider all the possible nkink+1 insertions, each insertion and nstate choice.
      end if
      !$OMP CRITICAL
      if ( big_node < rho_add_x(k1,j1) ) then
        big_node = rho_add_x(k1,j1)
      end if
      !$OMP END CRITICAL
    end do ! j1 (nstate)
  end do ! k1 (nkink+1)
  ! so in here, we have big_node in each process, we need big now,
  ! which is the largest element
  call  MPI_reduce(big_node, big, 1, &
                   MPI_DOUBLE_PRECISION,MPI_MAX,root_process,MPI_COMM_WORLD,ierr)
  ! this step we get the big from max in each process,  next step is bcast
  ! this value to all process
  call MPI_BCAST(big, 1, &
                 MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr)
  ! next step is to calculate the Summation of weight in each process,
  ! exp( rho_add_x(k1,j1)-big)
  sum_node=0.d0
  do k1=1,1
    do j1= 1, nstate
      sum_node=sum_node+exp(rho_add_x(k1,j1)-big)
    end do
  end do
  sum_weight=0.d0
  Call MPI_reduce(sum_node, sum_weight,1, &
                  MPI_DOUBLE_PRECISION,MPI_SUM,root_process,MPI_COMM_WORLD,ierr)
  ! sum_weight in here is xnrom_add.
  ! now we need to have a flag to know whether we need this subroutine to figure out ::
  ! if try to insert a kink, which kink and which state in insert.
  ! for convenience , we use MPI_gather,to get list_sum for all procoess, for next step
  ! judgement to know which process it will be picked up.
  !call mpi_gather(sum_node, 1, MPI_DOUBLE_PRECISION, list_choice, 1, &
  !               MPI_DOUBLE_PRECISION,root_process,mpi_comm_world,ierr)
 
  sum_node_y=0.d0
  do k1=1, 1
    do j1= 1, nstate
      sum_node_y=sum_node_y+exp(rho_add_x(k1,j1)-big)*rho_add_y(k1,j1)    
  !write(*,*) sum_node_y,j1, "sum_node_y"
  ! write(*,*) 
    end do
  end do
  sum_weight_y=0
  Call MPI_reduce(sum_node_y, sum_weight_y,1, &
                  MPI_DOUBLE_PRECISION,MPI_SUM,root_process,MPI_COMM_WORLD,ierr)

  
 ! calculate Q0 and E0,  Q0=exp(-beta*E0) 
 
  e_zero=rho_dev_1(j_old(1),j_old(1))/abs(rho_1(j_old(1),j_old(1)))
  !q_zero=rho_1(j_old(1),j_old(1))**p  
  q_zero= exp(p*log(rho_1(j_old(1),j_old(1)))-big)
  !q_zero in here is times exp(-big) as a factor 

 !  write(*,*) sum_weight, sum_weight_y, e_zero, q_zero,rho_1(j_old(1),j_old(1)),big,exp(p*log(rho_1(j_old(1),j_old(1)))-big),"Q_90"  
  E_tr= (e_zero*q_zero+sum_weight_y)/( q_zero+sum_weight )
  
  
end subroutine Truncation_Q2

! Try to remove a kink at every possible location, and return in rho_minus the
! energy for each case

subroutine weight_kink(j_size,j_old,rho)
  use kink
  implicit none

  integer, intent(in) :: j_size
  integer, intent(in) :: j_old(j_size)
  real(dp), intent(out) :: rho

 ! integer:: j_new_minus(nkink-1)
  integer i1
  real(dp)  energy, sfunc, sign_sfunc

 
   if (nkink==1) then 
      rho = p*log(rho_1(j_old(1),j_old(1)))
   else
     energy = log(p/dble(nkink))
      do i1=1,nkink 
         energy=energy+log_rho_1(j_old(i1),j_old(1+modulo(i1,nkink))) 
      end do 
      call sfunc_calc(j_old,j_size,nkink,p,sfunc,sign_sfunc)
       rho=energy+sfunc 
   endif   

  ! rho=exp(rho)

         
end subroutine weight_kink

subroutine removal_kink(j_size, j_old, rho_minus)
  use kink
  implicit none

  integer, intent(in) :: j_size
  integer, intent(in) :: j_old(j_size)
  real(dp), intent(out) :: rho_minus(j_size)

  integer  :: j_new_minus(nkink-1)
  integer  :: ikink, iplus, iminus, i
  real(dp) :: energy, sfunc, sign_sfunc

  ! shortcut for just to kinks
  if (nkink == 2) then
    rho_minus(1) = p*log(rho_1(j_old(2),j_old(2)))
    rho_minus(2) = p*log(rho_1(j_old(1),j_old(1)))
  else
    ! go through each possible kink location
    do ikink = 1, nkink
      iplus =  modulo(ikink-1 +1, nkink)+1
      iminus = modulo(ikink-1 -1, nkink)+1
      ! If removal here would create one fewer kink (and not two), and
      ! the transision between the previous and next state is possible (Ham>0)
      ! then calculate the energy (otherwise set something very small)
      !if (j_old(iplus) .ne. j_old(iminus) .and. &
      !    rho_1(j_old(iplus),j_old(iminus)).ne.0) then
      if (j_old(iplus) .ne. j_old(iminus) .and. &
          abs(rho_1(j_old(iplus),j_old(iminus))).ge.1.d-20) then

        j_new_minus(    1:ikink-1) = j_old(      1:ikink-1)
        j_new_minus(ikink:nkink-1) = j_old(ikink+1:nkink)

        energy = log(p/dble(nkink-1))
        do i = 1, nkink-1
          energy = energy + log_rho_1(j_new_minus(i), j_new_minus(1+modulo(i,nkink-1)))
        enddo
        call sfunc_calc(j_new_minus, j_size-1, nkink-1, p, sfunc, sign_sfunc)
        rho_minus(ikink) = energy + sfunc
      else
        rho_minus(ikink) = -8800.0
      end if
    end do ! ikink
  end if
end subroutine removal_kink

subroutine insertion_kink(j_size, j_old, rho_add_x, max_rho, sum_weight)
  use kink
  implicit none

  integer, intent(in)   :: j_size
  integer, intent(in)   :: j_old(j_size)
  real(dp), intent(out) :: rho_add_x(nmax_kinks+5,nstate)
  real(dp),intent(inout):: max_rho
  real(dp),intent(out)  :: sum_weight

  integer  :: j_new_add(j_size+1)
  integer  :: k1,iplus,iminus,i1,j1,i2,n1
  real(dp) :: energy, sfunc, sign_sfunc
  integer, allocatable :: rho1_nonzero_plusminus(:,:)

  allocate(rho1_nonzero_plusminus(nkink+1, 0:nstate))
  ! Try to add a kink in any of the positions (between/before/after existing kinks)
  do k1=1, nkink+1
    iplus  = modulo(k1-1 +1, nkink+1)+1
    iminus = modulo(k1-1 -1, nkink+1)+1

    ! copy all values of j_old into j_new_add, but leave the k1'th free
    j_new_add(    :k1-1)    = j_old(  :k1-1)
    j_new_add(k1+1:nkink+1) = j_old(k1:nkink)

    rho_add_x(k1,:) = -8800.0
    call ham_column_intersection(j_new_add(iminus), j_new_add(iplus), &
                                 rho1_nonzero_plusminus(k1,1:),      &
                                 rho1_nonzero_plusminus(k1,0))
    ! TODO: DISABLED
    !OMP PARALLEL DO PRIVATE(n1,j1,i1,i2,energy,sfunc,sign_sfunc) FIRSTPRIVATE(j_new_add)
    do n1 = 1, rho1_nonzero_plusminus(k1,0)
      j1 = rho1_nonzero_plusminus(k1,n1)
      j_new_add(k1)=j1

      if (j_new_add(iminus) .ne. j1 .and. &
          j_new_add(iplus)  .ne. j1 .and. &
          (nkink+1) <= nmax_kinks) then
        energy = log(p/dfloat(nkink+1))
        do i1 = 1, nkink+1
          i2 = 1 + modulo(i1,nkink+1)
          energy = energy + log_rho_1(j_new_add(i1), j_new_add(i2))
        end do
        !$OMP CRITICAL
        ! This is expensive: called way too often
        call sfunc_calc(j_new_add,j_size+1,nkink+1,p,sfunc,sign_sfunc)
        !$OMP END CRITICAL
        energy = energy + sfunc
        rho_add_x(k1,j1) = energy        ! rho_add here is the log(rho_add)
      end if
      !$OMP CRITICAL
      if ( max_rho < rho_add_x(k1,j1) ) then
        max_rho = rho_add_x(k1,j1)
      end if
      !$OMP END CRITICAL
    end do ! j1 (nstate)
  end do ! k1 (nkink+1)
  ! next step is to calculate the Summation of weight in each process,
  ! exp( rho_add_x(k1,j1)-max_rho)
  sum_weight=0.d0
  do k1=1, nkink+1
    do n1 = 1, rho1_nonzero_plusminus(k1,0)
      j1 = rho1_nonzero_plusminus(k1,n1)
      sum_weight=sum_weight+exp(rho_add_x(k1,j1)-max_rho)
    end do
  end do
  ! sum_weight in here is xnorm_add.
  list_choice = sum_weight
  deallocate(rho1_nonzero_plusminus)
end subroutine insertion_kink


!!!Input file
! TODO: intent properly
! TODO: make parameters general (istate, m?, allocation statements)

! returns some incrementing time measurement. Don't use this for an absolute
! time - only for time _differences_.

subroutine save_gen()
  use kink
  use mpi
  implicit none
  integer i,j
  if (my_id > 0) then
    return
  endif
  write(*,*) 'Dumping gen_simple(:,:)'
  call system('mkdir -p '//trim(adjustl(outdir)))
  open(unit=555, file=trim(adjustl(outdir))//'/gen.dat')
  do i=1,nstate
    do j=1,nstate
      write(555,*) gen_simple(i,j)
    end do
  end do
  close(555)
end subroutine

subroutine save_rho()
  use kink
  implicit none
  integer i,j
  call system('mkdir -p '//trim(adjustl(outdir)))
  open(unit=555, file=trim(adjustl(outdir))//'/'//trim(adjustl(rho_filename)))
  do i=1,nstate
    do j=1,nstate
      write(555,*) rho_1(i,j)
     ! write(555,*) rho_dev_1(i,j)
    end do
  end do
  close(555)
  open(unit=555, file=trim(adjustl(outdir))//'/input')
  write (555,'(i23,a)') PARAM_VERSION, ' # version'
  write (555,'(a,a)')   symmetry, '# symmetry'
  write (555,'(a,a)')   basisset, '# basisset'
  write (555,'(i23,a)') nrows, ' # nrows'
  write (555,'(i23,a)') ncols, ' # ncols'
  write (555,'(i23,a)') npass_diag, ' # npass_diag'
  write (555,'(i23,a)') diag_every, ' # diag_every'
  write (555,'(i23,a)') first_diag, ' # first_diag'
  write (555,'(i23,a)') ne_up, ' # ne_up'
  write (555,'(i23,a)') ne_dn, ' # ne_dn'
  write (555,'(i23,a)') nkink, ' # nkink'
  write (555,'(i23,a)') npass, ' # npass'
  write (555,'(i23,a)') nmax_kinks, ' # nmax_kinks'
  write (555,'(ES23.16,a)') beta, ' # beta'
  write (555,'(i23,a)') p, ' # p'
  write (555,'(i23,a)') ncomb, ' # ncomb'
  write (555,'(a,a)') '.', ' # outdir'
  write (555,'(a,a)') trim(adjustl(rho_filename)), ' # rho_filename'
  write (555,'(i23,a)') 0, ' # calc_rho'


  close(555)
end subroutine

subroutine load_rho()
  use kink
  implicit none
  integer i,j
  open(unit=555, file=rho_filename)
  do i=1,nstate
    do j=1,nstate
    end do
  end do
  close(555)
end subroutine

! The following function log_rho_1() is really just:
!real(dp) function log_rho_1(i,j)
!  integer,intent(in)::i,j
!  log_rho_1 = log(abs(rho_1(i,j)))
!end function log_rho_1
! but it is left as full function to optimization later

real(dp) function log_rho_1(i,j)
  use kink
  use mpi_data
  implicit none

  integer,intent(in)::i,j
  integer :: k
  logical :: found

  if (i == j) then
    log_rho_1 = -beta/p*Ha_diag(i)
  else if (i == ground_state_index) then 
    log_rho_1 = log(abs(exp(-beta/p*Ha_low(j))-1.d0))
  else
    found = .false.
    k = 1
    do while (.not. found .and. k < nappear_find)
      if (list_of_find(k) == i) then
        found = .true.
      else
        k = k+1
      end if
    end do
    if (found) then
      log_rho_1 = log(abs(exp(-beta/p*Ha_temp(k,j))-1.d0))
    else
      log_rho_1 = log(abs(exp(-beta/p*get_ham(i,j))-1.d0))
    end if
  end if
end function log_rho_1

! we do not use this function extensively
real(dp) function rho_dev_1(i,j)
  use kink
  use mpi_data
  implicit none

  integer,intent(in)::i,j
  real(dp):: h_new

  if (i == j) then
    rho_dev_1 = exp(-beta/p*Ha_diag(i))*Ha_diag(i)
  else
    h_new = get_ham(i,j)
    rho_dev_1 = exp(-beta/p*h_new)*h_new
  end if
end function rho_dev_1

!!!!!!!!NEW GROUP


End Program Kink_QMC

! Get time as accurate as possible, in milliseconds
integer(int8) function mytime()
  use vartypes
  implicit none

  integer(int8) :: count, count_rate, count_max
  call system_clock(count, count_rate, count_max)

  if (count_rate == 0) then
    mytime = 0
  else
    mytime = (1000*count)/count_rate
  endif
end function

#if 0
      subroutine get_hamil()
       ! this structrue is for Hubbard Model, to construct Hamiltonian  ! , please do not delete it !
      use kink
      implicit none
      real(dp), allocatable :: one_e_ham(:,:),store(:,:),two_e_ham(:,:,:,:)
      real(dp), allocatable :: w(:),work(:)
      real(dp) gen_simple,prod_up,prod_dn,u,elow
      integer k1,k2,k3,k4,k5,k6,icount,ne,jcount,ipick,it,itemp
      integer ix,iy,nboxes,iprod_up,iprod_dn
      integer info,ia,nstart,nend,m1,ndiff,ii1,j1,ii2,j2
      integer it_list(ne_up+ne_dn),it_diff(ne_up+ne_dn)
      integer, allocatable :: ker_up(:,:),ker_dn(:,:),jer(:,:)
      logical done,done2
      allocate(one_e_ham(nrows*ncols,nrows*ncols))
      allocate(store(nrows*ncols,nrows*ncols))
      allocate(two_e_ham(nrows*ncols,nrows*ncols,nrows*ncols,nrows*ncols))
      allocate(w(nrows*ncols))
      allocate(work(nrows*ncols*3))
      one_e_ham = 0
      store = 0
      two_e_ham = 0
      w = 0
      work = 0
! u should really be read in from an input file
      u=4.d0
      nboxes=nrows*ncols
      prod_up=1.d0
      do k1=1,ne_up
        prod_up=prod_up*dfloat(nboxes+1-k1)/k1
      enddo
      iprod_up=nint(prod_up)
      prod_dn=1.d0
      do k1=1,ne_dn
        prod_dn=prod_dn*dfloat(nboxes+1-k1)/k1
      enddo
      iprod_dn=nint(prod_dn)
      allocate(ker_up(ne_up,iprod_up+1))
      allocate(ker_dn(ne_dn,iprod_dn+1))
      ne=ne_up+ne_dn
      allocate(jer(ne,iprod_up*iprod_dn))
      ker_up = ker_dn = jer = 0

      one_e_ham = 0.d0

      icount=0
      do k2=1,nrows
       do k1=1,ncols
        icount=icount+1
        if(k1-1 .ne. 0) then
         one_e_ham(icount,icount-1)=-1.d0
         one_e_ham(icount-1,icount)=-1.d0
 !add
        else
 !        one_e_ham(icount,icount+(ncols-1))=-1.d0
 !        one_e_ham(icount+(ncols-1),icount)=-1.d0
        endif
        if(k1+1 .le. ncols) then
         one_e_ham(icount,icount+1)=-1.d0
         one_e_ham(icount+1,icount)=-1.d0
         else
 !        one_e_ham(icount,icount-(ncols-1))=-1.d0
 !        one_e_ham(icount-(ncols-1),icount)=-1.d0
        endif
        if(k2-1 .ne. 0) then
         one_e_ham(icount,icount-ncols)=-1.d0
         one_e_ham(icount-ncols,icount)=-1.d0
        else
 !        one_e_ham(icount,icount+(nrows-1)*ncols)=-1.d0
 !        one_e_ham(icount+(nrows-1)*ncols,icount)=-1.d0
        endif
        if(k2+1 .le. nrows) then
         one_e_ham(icount,icount+ncols)=-1.d0
         one_e_ham(icount+ncols,icount)=-1.d0
        else
   !      one_e_ham(icount,icount-(nrows-1)*ncols)=-1.d0
   !      one_e_ham(icount-(nrows-1)*ncols,icount)=-1.d0
        endif
       enddo
      enddo
      store=one_e_ham
      call DSYEV( 'V', 'U', nrows*ncols, store,nrows*ncols, w, work, nrows*ncols*3, info )
!      do k1=1,icount
!        print*,k1,w(k1)
!      enddo

      do k1=1,ne_up
        ker_up(k1,1)=k1
      enddo
      done=.false.
      jcount=1
      ipick=ne_up
      do while(.not. done)
        done2=.false.
        do while(.not. done2)
          it=ker_up(ipick,jcount)+1
          if(it .gt. icount-(ne_up-ipick)) then
           ipick=ipick-1
          else
           done2=.true.
          endif
          if(ipick .eq. 0) done2=.true.
        enddo
        if(ipick .eq. 0) done=.true.
        if(.not. done) then
          jcount=jcount+1
          do k1=1,ne_up
            ker_up(k1,jcount)=ker_up(k1,jcount-1)
          enddo
          ker_up(ipick,jcount)=ker_up(ipick,jcount)+1
          do k1=ipick+1,ne_up
            ker_up(k1,jcount)=ker_up(k1-1,jcount)+1
          enddo
        ipick=ne_up
        endif
      enddo

      do k1=1,ne_dn
        ker_dn(k1,1)=k1
      enddo
      done=.false.
      jcount=1
      ipick=ne_dn
      do while(.not. done)
        done2=.false.
        do while(.not. done2)
          it=ker_dn(ipick,jcount)+1
          if(it .gt. icount-(ne_dn-ipick)) then
           ipick=ipick-1
          else
           done2=.true.
          endif
          if(ipick .eq. 0) done2=.true.
        enddo
        if(ipick .eq. 0) done=.true.
        if(.not. done) then
          jcount=jcount+1
          do k1=1,ne_dn
            ker_dn(k1,jcount)=ker_dn(k1,jcount-1)
          enddo
          ker_dn(ipick,jcount)=ker_dn(ipick,jcount)+1
          do k1=ipick+1,ne_dn
            ker_dn(k1,jcount)=ker_dn(k1-1,jcount)+1
          enddo
        ipick=ne_dn
        endif
      enddo
      jcount=0
      do k1=1,iprod_up
        do k2=1,iprod_dn
          jcount=jcount+1
          do k3=1,ne_up
            jer(k3,jcount)=ker_up(k3,k1)
          enddo
          do k3=1,ne_dn
            jer(k3+ne_up,jcount)=ker_dn(k3,k2)
          enddo
        enddo
      enddo

      do k1=1,icount
        do k2=1,icount
          do k3=1,icount
            do k4=1,icount
              gen_simple=0.d0
              do k5=1,icount
!                do k6=1,icount
!                  do k7=1,icount
!                    do k8=1,icount
!                     if(k5 .eq. k7 .and. k6 .eq. k8 .and. k5 .eq. k6)
                gen_simple=gen_simple+store(k5,k1)*store(k5,k2)*u*store(k5,k3)*store(k5,k4)
              enddo
            two_e_ham(k1,k2,k3,k4)=sumcheck
            enddo
          enddo
        enddo
      enddo
    do k3=1, 9
      do k4=1,9
  ! write(*,*) two_e_ham(1,2,k3,k4), two_e_ham(1,2,k4,k3) ,k3,k4
    write(*,*) one_e_ham(k3,k4), one_e_ham(k4,k3),k3,k4, "LOOP"
      end do
     end do
   !stop
      elow=0.d0
      do k1=1,jcount
        do k2=k1,jcount
        sumcheck=0.d0
        ndiff=0
        do k5=1,ne_up+ne_dn
         it_list(k5)=jer(k5,k2)
        enddo

        m1=1
        do ia=1,2
          if(ia .eq. 1) then
            nstart=1
            nend=ne_up
          else
            nstart=ne_up+1
            nend=ne
          endif
          do k5=nstart,nend
            k3=1
            k4=nstart
            do while(k4 .le. nend .and. k3 .eq. 1)
              if(jer(k5,k1) .eq. it_list(k4)) then
                k3=0
                itemp=it_list(k4)
                it_list(k4)=it_list(k5)
                it_list(k5)=itemp
                if(k4 .ne. k5) m1=-m1           !Notes, a1,a2,...,an,,,,am,,,,
              endif                             ! we exchange an, am,  2*(m-n) +1 order, always odd,so change sign
              k4=k4+1
            enddo
            ndiff=ndiff+k3
            if(k3 .eq. 1) then
              it_diff(ndiff)=k5
            endif
          enddo
        enddo




        if(ndiff .eq. 0) then
          do k3=1,ne
            sumcheck=sumcheck+w(jer(k3,k1))        ! the kinetic energy each electrons (ne=ne_up+ne_dn)
          enddo
!print*,sumcheck
          do k3=1,ne_up
            do k4=1,ne_dn
              sumcheck=sumcheck+two_e_ham(jer(k3,k1),jer(ne_up+k4,k1),jer(k3,k2),jer(ne_up+k4,k2))
!print*,sumcheck
           enddo
         enddo
        else if(ndiff .eq. 1) then
          ii1=it_diff(1)
          j1=it_list(ii1)

          do k3=1,ne
            if(k3 .ne. ii1) then
              sumcheck=sumcheck+two_e_ham(jer(k3,k1),jer(ii1,k1),it_list(k3),j1)
              if((k3 .le. ne_up .and. ii1 .le. ne_up) .or. (k3 .gt. ne_up .and. ii1 .gt. ne_up))&
                sumcheck=sumcheck-two_e_ham(jer(k3,k1),jer(ii1,k1),j1,it_list(k3))
            endif
          enddo          ! we k3 and ii1 on the same side, we need to consider slater determinant,|12>=|12>-|21>
        else if(ndiff .eq. 2) then
          ii1=(it_diff(1))
          ii2=(it_diff(2))
          j1=it_list(ii1)
          j2=it_list(ii2)
          sumcheck=sumcheck+two_e_ham(jer(ii1,k1),jer(ii2,k1),j1,j2)
          if((ii1 .le. ne_up .and. ii2 .le. ne_up) .or. (ii1 .gt. ne_up .and. ii2 .gt. ne_up)) &
            sumcheck=sumcheck-two_e_ham(jer(ii1,k1),jer(ii2,k1),j2,j1)
        endif
!elow=min(elow,sumcheck)
!print*,k1,k2,sumcheck,elow,ndiff
    !    ha(k1,k2)=sumcheck*m1
    !    ha(k2,k1)=sumcheck*m1
      enddo
      enddo
      deallocate(store,one_e_ham,w,work,ker_up,ker_dn,jer,two_e_ham)
      return
      end
#endif


subroutine read_2038()
 
 use kink 
 use kink_temp
 use gen_data
 use io
 implicit none
! integer,dimension(20000):: r1
 integer:: j
 integer:: x1
 
 integer:: k_count ! the variable to count how many freeze oribtals 
 k_count=0
 open(unit=2038,file=trim(adjustl(indir)) // '/2038.dat',status='OLD')
   j = 0
   do
     read(2038,*,end=30) x1
     j = j+1
     Ha_initial_max(j) = x1
    ! if( j<10) then 
   !  write(*,*)   Ha_initial_max(j),j
    ! end if 
   end do
30 Close(2038)
   

end subroutine 


subroutine read_2037()
 
 use kink 
 use kink_temp
 use gen_data
 use io
 implicit none
! integer,dimension(20000):: r1
 integer:: j
 integer:: x1

 real(dp) :: x2 
 integer:: k_count ,j_count
 k_count=1
 open(unit=2037,file=trim(adjustl(indir)) // '/2037.dat',status='OLD')
   j = 0
   j_count=0
   do
     read(2037,*,end=30) x2,x1
     j_count=j_count+1 
     j = j+1
     Ha_initial_a(k_count,j)=x2
     Ha_initial_b(k_count,j)=x1 
     
    
     if( Ha_initial_max(k_count) ==j ) then 
        k_count=k_count+1
        j=0
        if( mod(k_count,10000).eq.0) then 
          write(*,*) "cheers",k_count 
        end if
     end if 
   end do
30 Close(2037)
   
! now we can successfully read 2037.dat and 2038.dat,
! the non zero element from inputfile, rather than calcuation again! 
end subroutine 




subroutine read_orbital_energy()
 use kink 
 use gen_data
 use io
 implicit none
 real(dp),dimension(20000):: r1
 integer:: j,k1
 real(dp):: x1
 integer:: x3
 integer:: k_count ! the variable to count how many freeze oribtals 
 k_count=0
   open(unit=139,file=trim(adjustl(indir)) // '/orbital_energy.dat',status='OLD')
       j = 0
     do
       read(139,*,end=30) x3,x1
         j = j+1
         r1(j) = x1
         if( switch_frozen.eq.1) then 
           if(x1.lt.x_bar.and.j.le.ne_up) then 
             k_count=k_count+1
             vect(k_count)=j
           end if 
         end if 
     end do
30 Close(139)

   e_orb=0.d0
   do k1=1, j
      e_orb(k1)=r1(k1)
   end do
   n_occupied=k_count

   if (n_occupied > 0) then
    !  write(*,*) n_occupied, vect(n_occupied) , "kk"
   endif

end subroutine

subroutine read_casscf()
  use kink
  use mpi_data
  use gen_data
  use io
  implicit none
  include 'mpif.h'

  integer:: j,j_max
 ! integer:: vc(ne)
  integer:: a1
  real(dp):: a2
  integer:: ker_up_order(ne_up), ker_dn_order(ne_dn)
  integer:: tran(num_orbital/2)
  integer:: k1,k2,k4
  integer:: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
  integer:: ne

  ne=ne_up+ne_dn

  open(unit=128,file=trim(adjustl(indir)) // '/casscf_important_states.dat',status='OLD')
  j = 0
  !two_e=0.d0
  !two_e_ham=0.d0
  casscf_vector=0

  do
  
   read(128,*,end=30) a1,a2, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
      j = j+1

     casscf_vector(1,j)=x1
     casscf_vector(2,j)=x2
     casscf_vector(3,j)=x3
     casscf_vector(4,j)=x4
     casscf_vector(5,j)=x5
     casscf_vector(6,j)=x6
     casscf_vector(7,j)=x7
     casscf_vector(8,j)=x8
     casscf_vector(9,j)=x9
     casscf_vector(10,j)=x10
     casscf_coeff(j)= a2 

  end do

30  close(128)
 ! write(*,*) "hello"
   j_max=j
  ! Transformation from CASSCF to HF, the data from NWCHEM
   tran(1)=1
   tran(2)=2
   tran(9)=3
   tran(3)=4
   tran(13)=5
   tran(4)=6
   tran(10)=7
   tran(14)=8
   tran(5)=9
   tran(11)=10
   tran(6)=11
   tran(12)=12
   tran(7)=13
   tran(8)=14
  ! begin the tranformation from CASSCF to HF
  ! write(*,*) "hello"

  do k2=1,j_max
    do k1=1,10
  ! write(*,*) k1,k2
      casscf_vector( k1, k2) = tran( casscf_vector (k1,k2))

   ! write(*,*) casscf_vector(k1,k2)
    end do
 !  write(*,*) casscf_vector(:,k2) , "full", k2 
 end do

 !stop
! after transformation, we need do the sorting up electron and dn electron




  do k4=1, j_max
     ker_up_order(1:ne_up)=casscf_vector (1:ne_up,k4)
     call PIKSRT ( ne_up, ker_up_order(:))
     casscf_vector( 1:ne_up, k4) = ker_up_order(1:ne_up)
  end do

  do k4=1, j_max
     ker_dn_order(1:ne_dn)=casscf_vector (ne_up+1:ne,k4)
     call PIKSRT ( ne_dn, ker_dn_order(:))
     casscf_vector(ne_up+1:ne, k4) = ker_dn_order(1:ne_dn)
  end do

  !write(*,*)  casscf_vector ( 1:10, 2)

!stop


  ! now up and down electron position are ascending order !

  casscf_vector(ne_up+1:ne, 1:j_max )=casscf_vector(ne_up+1: ne, 1:j_max) +ne_up

 ! write(*,*)  casscf_vector ( :, 2)

  do k1=1,j_max
     do k2=1, ne_up
        if ( casscf_vector(k2,k1)>ne_up ) then
           casscf_vector(k2,k1)=casscf_vector(k2,k1) +ne_up
        end if
     end do
  end do

   do k1=1,j_max
     do k2=ne_up+1,ne
        if ( casscf_vector(k2,k1)>(ne_up+ne_dn) ) then
           casscf_vector(k2,k1)=casscf_vector(k2,k1)+num_orbital/2 -ne_up
        end if
     end do
   end do

 ! write(*,*) casscf_vector( :,3) , "mark"
 ! stop

  !allocate(jer(ne,j_max))
  !jer(1:10, 1:j_max) = casscf_vector( 1:10, 1:j_max)
  !nstate=j_max


end subroutine read_casscf
! symmetry_reduction_c2v(a,b) returns the symmetry of the product of a and b,
! according to the C_{2v} group
! (see, e.g., http://www.webqc.org/symmetrypointgroup-c2v.html)
!
! a, b and the returned value are integers, denoting a symmetry.
! The symmetry integers come from definitions in nwchem (hardcoded),
!  see (src/symmetry/sym_char_tab.F).
!  These are:   symmetry  int-code
!                  A1        0
!                  A2        1
!                  B1        2
!                  B2        3
! The product table this function implements is:
!       | A1  A2  B1  B2
!    ---+---------------
!    A1 | A1  A2  B1  B2
!    A2 | A2  A1  B2  B1
!    B1 | B1  B2  A1  A2
!    B2 | B2  B1  A2  A1
function symmetry_reduction_c2v(a,b)
  implicit none
  integer,intent(in) :: a,b
  integer::symmetry_reduction_c2v
  integer, dimension(4,4)::matrix
  matrix=reshape((/ 0,1,2,3, 1,0,3,2, 2,3,0,1, 3,2,1,0 /), shape(matrix))
  symmetry_reduction_c2v = matrix(a+1, b+1)
end function symmetry_reduction_c2v

! symmetry_reduction_d2h(a,b) returns the symmetry of the product of a and b,
! according to the D_{2h} group
! (see, e.g., http://www.webqc.org/symmetrypointgroup-d2h.html)
!
! a, b and the returned value are integers, denoting a symmetry.
! The symmetry integers come from definitions in nwchem (hardcoded),
!  see (src/symmetry/sym_char_tab.F).
!  These are:   symmetry  int-code
!                  Ag        0
!                  Au        1
!                  B1g       2
!                  B1u       3
!                  B2g       4
!                  B2u       5
!                  B3g       6
!                  B3u       7
! The product table this function implements is:
!      | Ag  Au  B1g B1u B2g B2u B3g B3u
! -----+-----------------------------------
!  Ag  | Ag  Au  B1g B1u B2g B2u B3g B3u
!  Au  | Au  Ag  B1u B1g B2u B2g B3u B3g
!  B1g | B1g B1u Ag  Au  B3g B3u B2g B2u
!  B1u | B1u B1g Au  Ag  B3u B3g B2u B2g
!  B2g | B2g B2u B3g B3u Ag  Au  B1g B1u
!  B2u | B2u B2g B3u B3g Au  Ag  B1u B1g
!  B3g | B3g B3u B2g B2u B1g B1u Ag  Au
!  B3u | B3u B3g B2u B2g B1u B1g Au  Ag
function symmetry_reduction_d2h(a,b)
  implicit none
  integer,intent(in) :: a,b
  integer::symmetry_reduction_d2h
  integer, dimension(8,8)::matrix
  matrix=reshape((/ 0,1,2,3,4,5,6,7, 1,0,3,2,5,4,7,6, 2,3,0,1,6,7,4,5, 3,2,1,0,7,6,5,4, &
                    4,5,6,7,0,1,2,3, 5,4,7,6,1,0,3,2, 6,7,4,5,2,3,0,1, 7,6,5,4,3,2,1,0 /), &
                 shape(matrix))
  symmetry_reduction_d2h = matrix(a+1, b+1)
end function symmetry_reduction_d2h

! wrapper function for all symmetry functions. Other code only needs to call
! symmetry_reduction(), and this function then figures out which symmetry to
! use by looking at the 'symmetry' variable within the kink module.
! We might potentially optimize the string comparisons if this turns out to be
! a bottleneck.
function symmetry_reduction(a,b)
  use kink
  implicit none
  integer,intent(in) :: a,b
  integer::symmetry_reduction
  integer::symmetry_reduction_c2v
  integer::symmetry_reduction_d2h

  if (symmetry == 'c2v') then
    symmetry_reduction = symmetry_reduction_c2v(a,b)
  elseif (symmetry == 'd2h') then
    symmetry_reduction = symmetry_reduction_d2h(a,b)
  else
    write(*,*) "unknown symmetry: ", symmetry
    stop
  endif
end function symmetry_reduction


!!! add a subroutine

subroutine compute_electron(icount,ne_up_x,jcount_size,ker_up_x)
     use kink
     implicit none
     integer:: k1,ipick,it
     integer,intent(in):: icount, ne_up_x
     integer,intent(in):: jcount_size
     integer::jcount
     logical:: done,done2

   ! integer,allocatable,intent(out):: ker_up_x(:,:)



   ! allocate(ker_up_x(ne_up_x,jcount_size))

     !allocate(ker_up_x(2,7))
    integer :: ker_up_x(ne_up_x,jcount_size)
    ! write(*,*) jcount_size, ne_up_x,"jcount_size, ne_up_x in side"
      do k1=1,ne_up_x
        ker_up_x(k1,1)=k1
      enddo

     ! write(*,*) ker_up_x(1,1), ker_up_x(2,1) , "too"
      done=.false.
      jcount=1
      ipick=ne_up_x
      do while(.not. done)
        done2=.false.
        do while(.not. done2)
          it=ker_up_x(ipick,jcount)+1
          if(it .gt. icount-(ne_up_x-ipick)) then
           ipick=ipick-1
          else
           done2=.true.
          endif
          if(ipick .eq. 0) done2=.true.
        enddo
        if(ipick .eq. 0) done=.true.
        if(.not. done) then
          jcount=jcount+1
          do k1=1,ne_up_x
            ker_up_x(k1,jcount)=ker_up_x(k1,jcount-1)
          enddo
          ker_up_x(ipick,jcount)=ker_up_x(ipick,jcount)+1
          do k1=ipick+1,ne_up_x
            ker_up_x(k1,jcount)=ker_up_x(k1-1,jcount)+1
          enddo
        ipick=ne_up_x
        endif
      enddo
  end subroutine compute_electron

! Sort array ARR of length N
SUBROUTINE PIKSRT(N,ARR)
 implicit none
  integer:: N, i,j,a
  integer ARR(N)
  do j=2, N
    a=ARR(j)
    do i=j-1,1,-1
      if (ARR(i)<=a) goto 10
      ARR(i+1)=ARR(i)
    end do
        i=0
10  ARR(i+1)=a
  end do
  return
END subroutine PIKSRT

! Sort array ARR of length N, and ARR_NNZ of the same size in the same way as ARR
SUBROUTINE PIKSRT_two(N,ARR,ARR_NNZ)
  use vartypes
  implicit none
  integer:: N, i,j,a
  real(dp)::b 
  integer ARR(N)
  real(dp):: ARR_NNZ(N) 

  do j=2, N
    a=ARR(j)
    b=ARR_NNZ(j) 
    do i=j-1,1,-1
      if (ARR(i)<=a) goto 10
 
      ARR(i+1)=ARR(i)
      ARR_NNZ(i+1)=ARR_NNZ(i) 
  ! very important!
    end do
        i=0
10  ARR(i+1)=a ; ARR_NNZ(i+1)=b     
  end do
  return
END subroutine PIKSRT_two 


subroutine array_union (A,na,B,nb,C,nc) 
  implicit none 
  integer,intent(in)  :: A(na), na, B(nb), nb
  integer,intent(out) :: C(na+nb), nc
  integer:: ia, ib

  ia = 1 
  ib = 1 
  nc = 0
  ! loop as long as we are within at least A or B
  do while (ia <= na .or. ib <= nb)
    nc = nc+1
    ! If we are done with A, save the remainder of B
    if (ia > na) then
   ! write(*,*) "hi"
      C(nc:nc+nb-ib) = B(ib:nb)
      nc = nc + nb-ib
      ib = nb + 1
    ! If we are done with B, save the remainder of A
    elseif (ib > nb) then
  !  write(*,*) "hello"
      C(nc:nc+na-ia) = A(ia:na)
      nc = nc + na-ia
      ia = na + 1
    ! If we we are still within both A and B, compare the elements and save
    ! the smaller and advance that array - and in case of equality save only
    ! once and advance both A and B
    elseif (A(ia) == B(ib)) then
      C(nc) = A(ia)
      ia = ia + 1
      ib = ib + 1
    elseif (A(ia) < B(ib)) then
      C(nc) = A(ia)
      ia = ia + 1
    else
      C(nc) = B(ib)
      ib = ib + 1
    end if
  enddo
end subroutine array_union

! TODO: should try Ridders' method
! 'bisect' A, and assume roughly linear distribution,
integer function get_linear_subsection_index(A,na, B,nb,ib)
  use vartypes
  use mod_sparse_matrix_int
  implicit none
  integer,intent(in)    :: na,nb,ib
  type(mat_element_int) :: A(na),B(nb)
  integer(int8)         :: k_b,ia1,ia2,it

  get_linear_subsection_index = 0
  it  = 0
  k_b = 1
  ia1 = 1
  ia2 = na
  ! catch cases where the needed index is on the boundary
  if (A(ia1)%e == B(ib)%e) k_b = ia1
  if (A(ia2)%e == B(ib)%e) k_b = ia2
  ! catch cases where the needed index is outside the possible range
  if (A(ia1)%e >  B(ib)%e .or. A(ia2)%e <  B(ib)%e) return
  do while (ia1 < ia2 .and. A(k_b)%e .ne. B(ib)%e)
    ! Compute an index in between ia1 and ia2 depending on A(ia?) and B(ib)
    ! that is hopefully close to the index we are looking for
    ! (assume linear distribution instead of simple bisection) 
    k_b = (ia1*(A(ia2)%e-B(ib)%e) + ia2*(B(ib)%e-A(ia1)%e)) / (A(ia2)%e-A(ia1)%e)
    if (k_b < ia1 .or. k_b > ia2) then
      write(*,*) "error!", ia1, ia2, k_b
      write(*,*) "As", A(ia1)%e, A(ia2)%e
      write(*,*) "Looking for ",B(ib)%e, "in ", A(1:na)
      stop
    end if
    ! make sure we don't end up exactly on a1 or a2 due to round-off
    if (k_b == ia1) k_b = ia1 + 1
    if (k_b == ia2) k_b = ia2 - 1
    ! If the function value at the new index is smaller than the one we are
    ! looking for, make the new index the new left boundary
    if (A(k_b)%e < B(ib)%e) then
      ia1 = k_b
    else
    ! Otherwise, make the new index the new right boundary
      ia2 = k_b
    end if
    if (ia1 == ia2-1 .and. A(ia1)%e < B(ib)%e .and. A(ia2)%e > B(ib)%e) exit
    it = it + 1
    ! end if we found it, or we ended up with an empty interval (didn't find it)
  end do
  if (A(k_b)%e == B(ib)%e) then
    get_linear_subsection_index = int(k_b)
  end if
end function

! return the intersection of two sorted arrays (elements are in both arrays)
! A, B: the two input arrays
! C: the resulting array, with the used size stored in nc
recursive subroutine array_intersection(A,na, B,nb, C,nc)
  use mod_sparse_matrix_int
  implicit none 
  integer,              intent(in)  :: na,nb
  type(mat_element_int),intent(in)  :: A(na),B(nb) 
  integer,              intent(out) :: nc, C(na+nb)  
  integer                           :: ia,ib,k_b
  integer :: get_linear_subsection_index

  ! First, make sure that A is at least as large as B (for later optimizations)
  if (nb > na) then
    call array_intersection(B,nb, A, na, C, nc)
    return
  end if

  ia = 1 
  ib = 1 
  nc = 0 
  ! For about equal sized arrays, simply walk both arrays simultaniously and
  ! save the common elements. For arrays of vastly different size, searching
  ! can be faster, and can be done in parallel.
  ! try to find each element of B in A
  do ib = 1, nb
    k_b = get_linear_subsection_index(A,na,B,nb,ib)
    if (k_b > 0) then
      nc = nc + 1
      C(nc) = B(ib)%e
    end if
  end do
end subroutine array_intersection

subroutine array_intersection_indexB(A,na,B,nb,C,nc,D) 
 ! two array, we find the element that A and B both have
 ! D store index of B(which element of B) contribute 

 implicit none 
 integer,intent(in):: na,nb,A(na),B(nb) 
 integer,intent(out)::nc,C(na+nb),D(na+nb) 
 integer:: ia,ib,k_b

 ia=1 
 ib=1 
 k_b=1
 nc=0 
 do while (ia.le.na)    
   if( A(ia).eq.B(ib)) then 
      nc=nc+1 
      C(nc)=B(ib)
      D(nc)=ib
     if (ia.eq.na) then 
        ia=ia+1
     else if (ib.eq.nb) then 
        ia=ia+1 
        ib=k_b        
     else 
        ia=ia+1
        ib=ib+1 
        k_b=ib  
     end if  
   else 
     if ( ia.eq.na) then 
        ia=ia+1 
     else if ( ib.eq.nb) then 
        ia=ia+1 
        ib=k_b        
     else 
        ib=ib+1 
     end if  
   end if 
 end do     
end subroutine array_intersection_indexB

subroutine calc_exact_energy()
  use vartypes
  use hamiltonian
  use kink
  use io
  implicit none
  real(dp),allocatable::exact_w(:),exact_work(:)
  real(dp):: e_0, qsum, esum
  real(dp),allocatable :: rho_exact(:,:)
  integer::INFO
  integer:: k1,k2

  ! This only works for small systems. Large systems either run out of memory or
  ! compute themselves to death. :)
  ! The number for nstate is kind of arbitrary. It should be large enough so
  ! that the exact energy is calculated for small systems where this is
  ! feasable, and it is not for large systems where this would take too long
  ! or not work at all due to memory limits
  if (nstate < 10000) then
    write(*,*) "Starting calculation of exact energy, stand by..."
    allocate(rho_exact(nstate, nstate))
    allocate(exact_w(nstate))
    allocate(exact_work(nstate*3))
    do k1=1,nstate
      do k2=k1,nstate
        rho_exact(k1,k2) = gen_simple(k1,k2)
        rho_exact(k2,k1) = rho_exact(k1,k2)
      end do
    end do
    ! diagonalization, to calcuate ground state energy
    call DSYEV( 'N', 'U', nstate, rho_exact, nstate, exact_w, exact_work, nstate*3, INFO )
    e_0 = exact_w(1)
    qsum = sum(exp(-beta*(exact_w-e_0)))
    esum = sum(exp(-beta*(exact_w-e_0))*exact_w)
    ! output
    write(*,*) '  Q, <E>_Q, <E>', qsum, esum/qsum, e_0,nstate 
    open(unit=4321,file=trim(adjustl(outdir))//'/'//'energy_exact.dat',status='unknown')
    write(4321,*) '# 1: Q, 2: <E>_Q, 3: <E>'
    write(4321,*) qsum, esum/qsum, e_0
    close(4321)
    ! cleanup
    deallocate(rho_exact)
    deallocate(exact_w)
    deallocate(exact_work)
    write(*,*) "...done", gen_simple(1,1)

  else
    write(*,*) "Skipping calculation of exact energy because the number of states is too large."
  endif
end subroutine


! [1] "Modern Quantum Chemistry, Introduction to Advanced Electronic Structure",
!     by Attila Szabo and Neil S. Ostlund

