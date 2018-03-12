
module kink
  use vartypes

!*********************************************
  real(dp),allocatable::Ha(:,:)
  real(dp), allocatable:: Ha_diag(:), inv_exp_betax_Ha_diag(:)
!*********************************************
  ! setting a factorial loop up table
  real(dp),allocatable :: fact(:)
  character(len=3) :: symmetry
  logical :: restrict_SD, restrict_SDT
  character(len=3) :: basisset
  integer:: nrows,ncols, npass_diag, first_diag, diag_every
  integer:: ne_up, ne_dn   ! the number of up_electron and dn_electron
  integer :: e_up_state, e_dn_state, nstate
  integer:: ncomb
  integer:: nkink,npass
  integer::nmax_kinks,nmax_states
  integer(int8) :: p
  !@integer:: p
  integer:: ipass, imulti

  real(dp):: beta, beta_x
!add
  integer:: num_orbital
  integer:: num_total
  integer:: two_e_block_size, nstate_block_size,nnz_block_size
  real(dp):: repulsion_energy
  integer:: switch_frozen, x_dimension
  integer:: ground_state_index
  real(dp),allocatable::Ha_low(:)
  integer:: call_size 
  real(dp):: casscf_coeff(41)
 ! integer:: casscf_vector(10,41) 
end module
