module vartypes
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: int8 = selected_int_kind(18)

  type ham_element
    real(dp) :: v    ! value of ham
    real(dp) :: rho1 ! cache for exp(-beta/p*ham)
  end type
end module
