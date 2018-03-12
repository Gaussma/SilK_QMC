module gen_data
  use vartypes
  integer,allocatable, protected :: sym_orb(:)
  real(dp), allocatable, protected :: one_e_ham(:,:)
  real(dp), allocatable, protected :: two_e_ham(:,:,:,:)
  integer:: casscf_vector(10,41)
  real(dp),allocatable:: two_e(:,:,:)
  integer:: two_e_low_bound, two_e_high_bound
  integer:: nstate_low_bound, nstate_high_bound
  integer:: nnz_low_bound,nnz_high_bound 
  integer,allocatable::nnz_low_bound_record(:),nnz_high_bound_record(:) 

  integer, allocatable, protected :: jer(:,:)
  integer,allocatable:: jer_new(:,:) 
  integer::la(2,92)
  real(dp),allocatable :: e_orb(:)
  real(dp)::x_bar
  integer::n_occupied 
  integer,allocatable::vect(:)
  integer::list_initial(36) !"LUMO and HOMO"
  integer,allocatable:: list_ccsd(:), list_ccsdt(:)  
  integer::lh_initial 

  contains

  




  subroutine read_one_integral(num_orbital)
    use io
    implicit none
    integer, intent(in) :: num_orbital
    integer :: r1, r2
    real(dp)  :: c
    allocate(one_e_ham(num_orbital,num_orbital))
    one_e_ham = 0.d0
    open(unit=126,file=trim(adjustl(indir)) // '/one_electron_integral.dat',status='OLD')
    do
      read(126,*, end=30) c, r1, r2
      one_e_ham(r1, r2) = c
    end do
30  close(126)
  end subroutine

  subroutine read_two_integral(num_orbital)
!    use mpi_data
    use io
    implicit none
    include 'mpif.h'
    integer, intent(in) :: num_orbital
    integer :: i1, i2, i3, i4, IOstatus
    real(dp)  :: e
    allocate(two_e_ham(num_orbital,num_orbital,num_orbital,num_orbital))
    two_e_ham = 0.d0
    open(unit=127,file=trim(adjustl(indir)) // '/two_electron_integral.dat',status='OLD')
    IOstatus = 0
    do while (IOstatus == 0)
      read(127, *, iostat=IOstatus) e, i1, i2, i3, i4
#if 1
      two_e_ham(i1, i2, i3, i4) = e
#else
      if (i1 < i2 .and. abs(e) > 0.000001) then
      block
        integer :: i
        i = tridiag_index(i1, i2)
        if( two_e_low_bound <= i .and. i <= two_e_high_bound) then
          two_e(i-two_e_low_bound+1, i3, i4) = e
        endif
      end block
      endif
#endif
    end do
    close(127)
  end subroutine

  subroutine read_orbital_symmetry(num_orbital)
    use io
    implicit none
    integer, intent(in) :: num_orbital
    integer             :: i, sym, IOstatus
    allocate(sym_orb(num_orbital))
    sym_orb = 0
    open (unit=109, file=trim(adjustl(indir)) // '/symmetry.dat', status='OLD')
    IOstatus = 0
    do while (IOstatus == 0)
      read(109, *, iostat=IOstatus) i, sym
      sym_orb(i) = sym
    end do
    close(109)
end subroutine  

subroutine setup_jer(k_1,k_2, my_id)
!  use kink
!  use mpi_data
  use io 
  use kink
  implicit none
  integer,intent(in) :: k_1
  integer,intent(in) :: k_2
  integer,intent(in) :: my_id

  real(dp) :: HF_E
  integer symmetry_reduction
  real(dp)  sumcheck,prod_up,prod_dn,elow
  integer k1,k2,k3,k4,k5,ne,jcount,ipick,itemp,jcount_up,jcount_dn
  integer half_num_orbital,iprod_up,iprod_dn
  integer ia,nstart,nend,m1,ndiff,ii1,j1,ii2,j2
  integer it_list(ne_up+ne_dn),it_diff(ne_up+ne_dn)
  integer, allocatable :: ker_up(:,:),ker_dn(:,:)

  real(dp):: sum_effect(ne_up+ne_dn), sum_ndiff1
  logical done
  integer:: jer_order(ne_up+ne_dn)
  integer::k_pick 
  integer:: sym_count,sym_final
  integer:: ne_up_free, ne_dn_free, half_num_orbital_free
  integer,allocatable::jer_1(:,:)

  integer,allocatable:: ker_up_order(:), ker_dn_order(:)
 
  integer:: ip_casscf

  integer:: k_count_up, k_count_dn,k_count
  integer:: up_s,dn_s,up_d,dn_d,up_t,dn_t
  integer,allocatable::list_up_s(:)
  integer,allocatable::list_dn_s(:)
  integer,allocatable::list_up_d(:)
  integer,allocatable::list_dn_d(:)
  integer,allocatable::list_up_t(:)
  integer,allocatable::list_dn_t(:)
 

  allocate(ker_up_order(ne_up))
  allocate(ker_dn_order(ne_dn))
  allocate (e_orb(num_orbital))
  ker_up_order = 0
  ker_dn_order = 0
  e_orb = 0

  call read_orbital_symmetry(num_orbital)
  call read_orbital_energy()

  call read_one_integral(num_orbital)
  call read_two_integral(num_orbital)
 

 
  ip_casscf=0
  if (ip_casscf.eq.1) then 
    call read_casscf()
    sumcheck=0.d0
    ne=ne_up+ne_dn
    nstate=41
    allocate(jer(10,41))
    do k1=1,41
      jer(1:10,k1)= casscf_vector(1:10, k1)
      write(*,*) jer(:,k1), "show"
    end do
  else

  ! since we allocate two_e_ham in gen subroutine, so it is better to call read_two_integral in here 
  half_num_orbital=nrows*ncols

  ! here , we do the changes   
  ne_up_free = ne_up-n_occupied
  ne_dn_free = ne_dn-n_occupied
  half_num_orbital_free = half_num_orbital-n_occupied 

  ! This calculates half_num_orbita_free! / ne_up_free! / (half_num_orbital_free-ne_up_free)!
  prod_up=1.d0
  do k1=1,ne_up_free
    prod_up=prod_up*(half_num_orbital_free+1.d0-k1)/k1
  enddo
  iprod_up=nint(prod_up)
  prod_dn=1.d0
  do k1=1,ne_dn_free
    prod_dn = prod_dn*(half_num_orbital_free+1.d0-k1)/k1
  enddo
  iprod_dn=nint(prod_dn)

  allocate(ker_up(ne_up,iprod_up+1))
  allocate(ker_dn(ne_dn,iprod_dn+1))
  ne = ne_up + ne_dn
  if (restrict_SD) then
    allocate(jer(ne, 197814))
    allocate(jer_new(ne,197814))
  else
    allocate(jer(ne,iprod_up*iprod_dn))
    allocate(jer_new(ne,iprod_up*iprod_dn)) 
  end if
  jer = 0
  jer_new = 0

  if (restrict_SD) then
    allocate(list_up_s(iprod_up))
    allocate(list_dn_s(iprod_dn))
    allocate(list_up_d(iprod_up))
    allocate(list_dn_d(iprod_dn))
    allocate(list_up_t(iprod_up))
    allocate(list_dn_t(iprod_dn))
  endif



 
  sumcheck=0.d0

  ! up electron inital setup
  ker_up=0
  ! Start with a simple sequence of integers
  ker_up(1:ne_up_free, 1) = (/ (k1, k1 = 1, ne_up_free) /)

  jcount_up=1
  ipick = ne_up_free
  do while(ipick > 0)
    ! pick an index from which on later overwriting values. Search from the
    ! end of ker_up for the first value that is smaller than what the maximum
    ! of orbitals allows (using that ker_up is ordered).
    ipick = ne_up_free
    do while(ker_up(ipick,jcount_up) > half_num_orbital_free-(ne_up_free-ipick)-1)
      ipick = ipick - 1
      if (ipick .eq. 0) exit ! Don't move this into the 'while' above
    enddo
    ! Now create a new 'jcount_up', copy the old, but reset all values starting
    ! with ipick by a sequence, starting with one more than what was at ipick
    if (ipick > 0) then
      jcount_up = jcount_up + 1
      ker_up(1:ipick, jcount_up) = ker_up(1:ipick, jcount_up-1)
      k2 = ker_up(ipick, jcount_up) + 1
      ker_up(ipick:ne_up_free,jcount_up) = (/ (k1, k1=k2, k2+ne_up_free-ipick) /)
    endif
  enddo

  !down electron initial setup        
  ker_dn=0 
  ! Start with a simple sequence of integers
  ker_dn(1:ne_dn_free, 1) = (/ (k1, k1 = 1, ne_dn_free) /)

  jcount_dn=1
  ipick = ne_dn_free
  do while(ipick > 0)
    ! pick an index from which on later overwriting values. Search from the
    ! end of ker_rn for the first value that is smaller than what the maximum
    ! of orbitals allows (using that ker_dn is ordered).
    ipick = ne_dn_free
    do while(ker_dn(ipick,jcount_dn) > half_num_orbital_free-(ne_dn_free-ipick)-1)
       ipick=ipick - 1
      if(ipick .eq. 0) exit ! Don't move this into the 'while' above
    enddo
    ! Now create a new 'jcount_dn', copy the old, but reset all values starting
    ! with ipick by a sequence, starting with one more than what was at ipick
    if (ipick > 0) then
      jcount_dn = jcount_dn + 1
      ker_dn(1:ipick, jcount_dn) = ker_dn(1:ipick, jcount_dn-1)
      k2 = ker_dn(ipick, jcount_dn) + 1
      ker_dn(ipick:ne_dn_free,jcount_dn) = (/ (k1, k1=k2, k2+ne_dn_free-ipick) /)
    endif
  enddo

  ! Step one, vect(1) need position.     
  ! frozen case, from inputfile 
  if (switch_frozen.eq.1) then 
    do k4=1,iprod_up
     ! k_pick=1
      do k1=1,n_occupied 
        done = .false.
        k2=1
        do while (.not. done .and. k2.le.ne_up_free)
          if( ker_up(k2,k4).ge.vect(k1)) then 
            ker_up(k2:ne_up_free,k4) = ker_up(k2:ne_up_free,k4)+1 
            done = .true.
           ! k_pick=k2+1 
          end if 
          k2=k2+1 
       end do 
      end do   
    end do

    do k4=1,iprod_up
      ker_up(ne_up_free+1:ne_up,k4) =vect(1:n_occupied) 
      ker_up_order(:)=ker_up(:,k4) 
      call PIKSRT((ne_up),ker_up_order(:))
      ker_up(:,k4)=ker_up_order(:)
    end do 
 
    do k4=1,iprod_dn
      k_pick=1
      do k1=1,n_occupied 
        done = .false.
        k2=k_pick
        do while (.not. done .and. k2.le.ne_dn_free)
          if( ker_dn(k2,k4).ge.vect(k1)) then 
            ker_dn(k2:ne_dn_free,k4) = ker_dn(k2:ne_dn_free,k4)+1 
            done = .true.
            k_pick=k2+1 
          end if 
          k2=k2+1 
       end do 
      end do   
    end do
    do k4=1,iprod_dn
      ker_dn(ne_dn_free+1:ne_dn,k4) =vect(1:n_occupied) 
      ker_dn_order(:)=ker_dn(:,k4) 
      call PIKSRT((ne_dn),ker_dn_order(:))
      ker_dn(:,k4)=ker_dn_order(:)
    end do 
  end if  ! for switch_frozen 

  ker_dn = ker_dn + ne_up

  where (ker_up(1:ne_up, 1:jcount_up) > ne_up)
    ker_up = ker_up + ne_dn
  end where
  where (ker_dn(1:ne_dn, 1:jcount_dn) > ne_up + ne_dn)
    ker_dn = ker_dn + num_orbital/2-ne_up
  end where
  if (restrict_SD) then
    k_count_up=0 
    k_count_dn=0 
    up_s=0  !the value of single excitation from ker_up
    dn_s=0  !the value of signle excitation from ker_dn 
    up_d=0  !the value of double excitation from ker_up
    dn_d=0  ! the value of double excitation from ker_dn
    list_up_s=0
    list_dn_s=0
    list_up_d=0
    list_dn_d=0
  !calculate list_up_s and list_dn_s , the number of single excitation from ker_up and ker_dn !
    !calculate single excitation from ker_up 
    k3=0 
    do k1=1,iprod_up 
      k_count_up=0 
      do k2=1, ne_up 
        if (1.le.ker_up(k2,k1).and.ker_up(k2,k1).le.ne_up) then 
         k_count_up=k_count_up+1 
        end if 
      end do 
      if ( k_count_up.eq.(ne_up-1)) then 
       k3=k3+1
       list_up_s(k3)=k1
      end if 
    end do 
      up_s=k3  ! get up_s value  !
   !calcuate single excitation from ker_dn !
      k3=0 
    do k1=1,iprod_dn 
      k_count_dn=0 
      do k2=1, ne_dn
        if ((1+ne_up).le.ker_dn(k2,k1).and.ker_dn(k2,k1).le.ne) then 
         k_count_dn=k_count_dn+1 
        end if 
      end do 
      if ( k_count_dn.eq.(ne_dn-1)) then 
       k3=k3+1
       list_dn_s(k3)=k1
      end if 
    end do 
    dn_s=k3   ! get dn_s value !
  !calculate list_up_d and list_dn_d , the number of double excitation from ker_up and ker_dn !
    !calculate double excitation from ker_up 
    k3=0  
    do k1=1,iprod_up 
      k_count_up=0 
      do k2=1, ne_up 
        if (1.le.ker_up(k2,k1).and.ker_up(k2,k1).le.ne_up) then 
         k_count_up=k_count_up+1 
        end if 
      end do 
      if ( k_count_up.eq.(ne_up-2)) then 
       k3=k3+1
       list_up_d(k3)=k1
      end if 
    end do 
      up_d=k3 ! get up_d value, double excitation from ker_up 
    !calcualte double excitation from ker_dn   
     k3=0  
      do k1=1,iprod_dn 
      k_count_dn=0 
      do k2=1, ne_dn 
        if ((1+ne_up).le.ker_dn(k2,k1).and.ker_dn(k2,k1).le.ne) then 
         k_count_dn=k_count_dn+1 
        end if 
      end do 
      if ( k_count_dn.eq.(ne_dn-2)) then 
       k3=k3+1
       list_dn_d(k3)=k1
      end if 
    end do 
     dn_d=k3  ! get dn_d value, double excitation from ker_dn
     k3=0  
  
  
    !write(*,*) up_s,dn_s,up_d,dn_d, "check"
   
    jcount=0
  !add HF determinant first
  
    jer(1:ne_up,1)=ker_up(1:ne_up,1)
    jer(1+ne_up:ne,1)=ker_dn(1:ne_dn,1)
    jcount=jcount+1
  !add single excitations for jer 
   !1)  1 up, 0 down 
     do k1=1, up_s 
     jcount=jcount+1
     jer(1:ne_up,jcount)=ker_up(1:ne_up,list_up_s(k1)) 
     jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,1)
     end do 
   !2) 0 up, 1 down 
      do k1=1,dn_s
      jcount=jcount+1
       jer(1:ne_up,jcount)=ker_up(1:ne_up,1)
        jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,list_dn_s(k1))
      end do  
  !add double ecitation for Jer 
    !1) 2 Up, 0 down 
    do k1=1, up_d
     jcount=jcount+1
     jer(1:ne_up,jcount)=ker_up(1:ne_up,list_up_d(k1)) 
     jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,1)
     end do 
    !2) 0 Up, 2 down  
     do k1=1,dn_d
      jcount=jcount+1
       jer(1:ne_up,jcount)=ker_up(1:ne_up,1)
        jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,list_dn_d(k1))
      end do  
  
    !3) 1 Up, 1 down 
     do k1=1,up_s
       do k2=1,dn_s
        jcount=jcount+1 
       jer(1:ne_up,jcount)=ker_up(1:ne_up,list_up_s(k1))
       jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,list_dn_s(k2))
      end do 
     end do   
    write(*,*) jcount, "SD"
  end if ! SD

  if (restrict_SDT) then
    up_t=0
    dn_t=0 
    list_up_t=0
    list_dn_t=0 
  !calculate list_up_t and list_dn_t , the number of triple excitation from ker_up and ker_dn !
    !calcualte triple excitation from ker_up 
    k3=0  
    do k1=1,iprod_up 
      k_count_up=0 
      do k2=1, ne_up 
        if (1.le.ker_up(k2,k1).and.ker_up(k2,k1).le.ne_up) then 
         k_count_up=k_count_up+1 
        end if 
      end do 
      if ( k_count_up.eq.(ne_up-3)) then 
       k3=k3+1
       list_up_t(k3)=k1
      end if 
    end do 
      up_t=k3 ! get up_d value, double excitation from ker_up 
    !calculate triple excitation from ker_dn    
     k3=0  
      do k1=1,iprod_dn 
      k_count_dn=0 
      do k2=1, ne_dn
        if ((1+ne_up).le.ker_dn(k2,k1).and.ker_dn(k2,k1).le.ne) then 
         k_count_dn=k_count_dn+1 
        end if 
      end do 
      if ( k_count_dn.eq.(ne_dn-3)) then 
       k3=k3+1
       list_dn_t(k3)=k1
      end if 
    end do 
     dn_t=k3  ! get dn_d value, double excitation from ker_dn
     k3=0 
  
     write(*,*) up_t,dn_t ,"up_t,dn_t"
  !calculate triple excitation for Jer 
  
   !1) 3 UP , 0 Down 
   
    do k1=1, up_t
     jcount=jcount+1
     jer(1:ne_up,jcount)=ker_up(1:ne_up,list_up_t(k1)) 
     jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,1)
     end do 
    !2) 0 Up, 2 down  
     do k1=1,dn_t
         jcount=jcount+1
        jer(1:ne_up,jcount)=ker_up(1:ne_up,1)
        jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,list_dn_t(k1))
     end do  
     !3) 1 Up, 2 down 
  
     do k1=1,up_s 
       do k2=1,dn_d 
        jcount=jcount+1
        jer(1:ne_up,jcount)=ker_up(1:ne_up,list_up_s(k1))
        jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,list_dn_d(k2))
       end do 
     end do 
    !4) 2 Up, 1 Down  
    
     do k1=1,up_d
       do k2=1,dn_s
        jcount=jcount+1
        jer(1:ne_up,jcount)=ker_up(1:ne_up,list_up_d(k1))
        jer(1+ne_up:ne,jcount)=ker_dn(1:ne_dn,list_dn_s(k2))
       end do 
     end do 
     write(*,*) jcount, "SDT"
   end if ! SDT
      
  if (.not. restrict_SD .and. .not. restrict_SDT) then
    jcount=0
    do k1=1,iprod_up
      do k2=1,iprod_dn
        jcount=jcount+1
        jer(1      :ne_up      ,jcount)=ker_up(1:ne_up,k1)
        jer(1+ne_up:ne_dn+ne_up,jcount)=ker_dn(1:ne_dn,k2)
      enddo
    enddo
  endif ! .not. restrict_SD .and. .not. restrict_SDT
  end if ! ip_casscf
  k1 = k_1
  k2 = k_2

  elow=0.d0
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
          if(k4 .ne. k5) m1=-m1
        endif
        k4=k4+1
      enddo
      ndiff=ndiff+k3
      if(k3 .eq. 1) then
        it_diff(ndiff)=k5
      endif
    enddo
  enddo

  ! worth to mention is that it_diff and jer is different!!!
  ! write(*,*) ndiff, it_diff(1),it_diff(2)
  ! write(*,*) it_list(it_diff(1)), it_list(it_diff(2))

  if (ndiff .eq. 0) then
    sum_effect=0.d0
    do k3=1, ne
      sum_effect(k3)=one_e_ham(jer(k3,k1), jer(k3,k1))
      do k5=1,ne
        sum_effect(k3)=sum_effect(k3) - two_e_ham(jer(k3,k1),jer(k5,k1),jer(k3,k1),jer(k5,k1))
      end do
    end do
    do k3=1, ne
      sumcheck=sumcheck+ sum_effect(k3)
    enddo
    if (my_id == 0) print*,sumcheck, "one_electron"

    do k3=1,ne_up+ne_dn
      do k4=1,ne_up+ne_dn
        sumcheck=sumcheck+0.5*two_e_ham(jer(k3,k1),jer(k4,k1),jer(k3,k2),jer(k4,k2))
      enddo
    enddo
    sumcheck=sumcheck+ repulsion_energy
  else if(ndiff .eq. 1) then
    sum_ndiff1=0.d0
    ii1=it_diff(1)
    ! new jer(ii1,k2)  match one to one  difference. jer(ii1,k2) and it_list is different.
    j1=it_list(ii1)

    sum_ndiff1=one_e_ham(jer(ii1,k1),j1)
    sumcheck=sumcheck+ sum_ndiff1
  else if(ndiff .eq. 2) then
    ii1=(it_diff(1))
    ii2=(it_diff(2))
    j1=it_list(ii1)
    ! remember it_list(ii1) and jer(ii1,k2) is different
    j2=it_list(ii2)
    sumcheck=sumcheck+two_e_ham(jer(ii1,k1),jer(ii2,k1),j1,j2)
  end if

  HF_E=sumcheck*m1

#if 1
  if (symmetry == 'c2v' .or. symmetry == 'd2h') then
    !using symmetry to reduce the size of nstate
    sym_count=0
    do  k1=1,jcount
      k3=1 
      sym_final=sym_orb(jer(1,k1))
      do while (k3.lt.ne) 
        sym_final = symmetry_reduction(sym_final,sym_orb(jer(k3+1,k1)))
        ! write(*,*) sym_final, k3,jcount , "see"
        k3=k3+1 
      end do 
      if (sym_final.eq.0) then 
        sym_count=sym_count+1
        jer_new(:,sym_count)=jer(:,k1) 
      end if 
    end do 
  
    deallocate(jer) 
    allocate(jer(ne,sym_count)) 
    jer = 0
    jer(1:ne,1:sym_count)=jer_new(1:ne,1:sym_count) 
    deallocate(jer_new) 
    nstate=sym_count 
    if (my_id == 0) write(*,*) nstate, "SYS"
  end if 
#endif

#if 0
  ! to get LUMO-HOMO state from DZ basis set ! do not delete it !
  if (basisset == 'DZH') then
    k3=0
    list_initial=0
    do k1=1,nstate
     if( jer(1,k1).eq.1.and.jer(2,k1).eq.2.and.jer(3,k1).eq.3.and.jer(6,k1).eq.6.and.jer(7,k1).eq.7.and.jer(8,k1).eq.8) then
       if(jer(4,k1).eq.4.or.jer(4,k1).eq.5.or.jer(4,k1).eq.11.or.jer(4,k1).eq.12) then
          if( jer(5,k1).eq.5.or.jer(5,k1).eq.11.or.jer(5,k1).eq.12) then
            if(jer(9,k1).eq.9.or.jer(9,k1).eq.10.or.jer(9,k1).eq.20.or.jer(9,k1).eq.21) then
              if(jer(10,k1).eq.10.or.jer(10,k1).eq.20.or.jer(10,k1).eq.21) then
                 k3=k3+1
                 list_initial(k3)=k1
              end if
            end if
          end if
       end if
     end if
    end do
   write(*,*) list_initial(1:k3), "lomo"
   lh_initial=k3 
   !stop  
  end if
#endif    
  write(*,*) "OK"
  write(*,*) ne_up,ne_dn, "e_up,e_dn "

#ifdef CCSD
  allocate(list_ccsd(nstate))
  allocate(list_ccsdt(nstate))
 ! filter all ccsd slater determinant from FULL space DZ  basis, DONOT DELETE IT , thanks !
 ! consider single excitation first  5C1*9C1*2=90 
 k_count=0
 k3=0 
  list_ccsd(1)=1 !include HF   
  k3=1 
 
 
 !single excitation from up electron   
  k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up-1).and.k_count_dn.eq.(ne_dn)) then 
       k3=k3+1
       list_ccsd(k3)=k1 
    end if 
  end do   

 !single excitation from down electron  
  k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up).and.k_count_dn.eq.(ne_dn-1)) then 
       k3=k3+1
       list_ccsd(k3)=k1 
    end if 
  end do   


 write(*,*) "single excitation", k3 

 !Begin Double excitation    

 ! 2 up, 0 down 


   k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up-2).and.k_count_dn.eq.(ne_dn)) then 
       k3=k3+1
       list_ccsd(k3)=k1 
    end if 
  end do   
 ! 0 UP, 2 DOWN 
    k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up).and.k_count_dn.eq.(ne_dn-2)) then 
       k3=k3+1
       list_ccsd(k3)=k1 
    end if 
  end do   




 ! third option, one excitation up, one excitation down.  

  k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up-1).and.k_count_dn.eq.(ne_dn-1)) then 
       k3=k3+1
       list_ccsd(k3)=k1 
    end if 
  end do   




 !do k2=1,k3 
 ! do k1=1,ne 
! write(1901,*) jer(k1,list_ccsd(k3)), k2  
!  write(*,*) jer(k1,list_ccsd(k3)), k2  
 ! end do 

 !end do 

!write(1901,*) k3, "ccsd "

#endif 

! write(*,*) k3, "CCSD"
 
 

#ifdef CCSDT
   ! add the filter for CCSDT , single, double, triple excitation. 
   ! acutally we only need to write the triple excitation here, since single and double is picked out by ccsd
   ! when running ccsdt, we need to turn on the ccsd first !  
    ! we have four cateogries,  
   !   excitations     e_UP    e_Down 
   !       3            3        0    ---> (5C3)*(9C3)= 
   !       3            2        1    ---->(5c2)*(9C2)*(5C1)*(9C1)
   !       3            1        2 
   !       3            0        3    
  ! part 1: 
   list_ccsdt(1:k3)=list_ccsd(1:k3) 
 !part4  
  k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up).and.k_count_dn.eq.(ne_dn-3)) then 
       k3=k3+1
       list_ccsdt(k3)=k1 
    end if 
  end do   
  ! Part1  



   k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up-3).and.k_count_dn.eq.(ne_dn)) then 
       k3=k3+1
       list_ccsdt(k3)=k1 
    end if 
  end do   

  ! part 2 ( 2 up , 1 down excitation)  

  k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up-2).and.k_count_dn.eq.(ne_dn-1)) then 
       k3=k3+1
       list_ccsdt(k3)=k1 
    end if 
  end do   

  ! part3 ( 1 up, 2 down excitation)  

  k_count_up=0
  k_count_dn=0 

  do k1=1, nstate 
    k_count_up=0
    do k2=1,ne_up 
      if(1.le.jer(k2,k1).and.jer(k2,k1).le.ne_up) then 
        k_count_up=k_count_up+1
      end if
    end do 
    k_count_dn=0
     do k2=ne_up+1,ne 
      if((ne_up+1).le.jer(k2,k1).and.jer(k2,k1).le.ne) then 
        k_count_dn=k_count_dn+1
      end if 
     end do  

    if( k_count_up.eq.(ne_up-1).and.k_count_dn.eq.(ne_dn-2)) then 
       k3=k3+1
       list_ccsdt(k3)=k1 
    end if 
  end do  
 write(*,*) "ccsdt",k3
#endif 

#ifdef CCSD
   !for CCSD in DZ basis 

   allocate(jer_1(ne,k3))
   do k1=1,k3 
    jer_1(1:ne,k1)=jer(1:ne,list_ccsdt(k1))
   end do 
 
    deallocate(jer) 
    allocate(jer(ne,k3))
 
    jer=0 
   write(*,*) k3, "k3"
   do k1=1,k3
   jer(1:ne, k1)=jer_1(1:ne,k1)
   end do 
   nstate=k3
 !do k2=1,k3 
!  do k1=1,ne 
 !write(1900,*) jer(k1,k2), k2  
!  end do 

! end do 
#endif
! write(*,*) k3, "CCSD"
 


#ifdef LUMOHOMO
   !for LUMO-HOMO in DZ basis 
   do k1=1,10 
    jer_1(1:ne,k1)=jer(1:ne,list_initial(k1))
   end do 
 
    deallocate(jer) 
    allocate(jer(10,36))
 
    jer=0 
   write(*,*) k3, "k3"
   do k1=1,10
   jer(1:ne, k1)=jer_1(1:ne,k1)
   end do 
   nstate=10
#endif
  if (basisset == 'STH') then
    !for STO-3g basis  
    k3=0
    list_initial=0
    do k1=1,nstate
     if( jer(1,k1).eq.1.and.jer(2,k1).eq.2.and.jer(3,k1).eq.3.and.jer(6,k1).eq.6.and.jer(7,k1).eq.7.and.jer(8,k1).eq.8) then
       if(jer(4,k1).eq.4.or.jer(4,k1).eq.5.or.jer(4,k1).eq.11.or.jer(4,k1).eq.12) then
          if( jer(5,k1).eq.5.or.jer(5,k1).eq.11.or.jer(5,k1).eq.12) then
            if(jer(9,k1).eq.9.or.jer(9,k1).eq.10.or.jer(9,k1).eq.13.or.jer(9,k1).eq.14) then
              if(jer(10,k1).eq.10.or.jer(10,k1).eq.13.or.jer(10,k1).eq.14) then
                 k3=k3+1
                 list_initial(k3)=k1
              end if
            end if
          end if
       end if
     end if
    end do
   !write(*,*) list_initial(1:k3), "lomo",k3
    do k1=1,k3
      write(1935,*) jer(:, list_initial(k1)), "jer -list", list_initial (k1) 
    end do 
  end if

 ! we sorting jer in bottom, since gen(1,1) keep the same! 
  do k3=1, nstate
    jer_order(:) = jer(:,k3)
    call PIKSRT(ne,jer_order)
    jer(:,k3)=jer_order(:)
  end do

  if (my_id == 0) write(*,*) nstate, "nstate after SYMMETRY"
 !stop
  return
end subroutine setup_jer

end module gen_data
