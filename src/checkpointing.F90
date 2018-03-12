! This module handles all low-level checkpointing

module checkpointing
  use vartypes
  use hdf5
  INTEGER        :: hdferr
  INTEGER(HID_T) :: hdfile
  character(1000):: hdfilename = "checkpoint.h5"
  logical        :: gzip_available = .false.

contains

  ! Initiate hdf5 writing, to be called before writing anything
  subroutine checkpoint_start_to_file(filename)
    implicit none
    character(*), intent(in) :: filename
    logical :: avail
    integer :: filter_info, filter_info_both

    ! open interface
    call h5open_f(hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5open_f"
    ! enable compression if available
    call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, avail, hdferr)
    if (avail) then
      call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, hdferr)
      filter_info_both = ior(H5Z_FILTER_ENCODE_ENABLED_F, H5Z_FILTER_DECODE_ENABLED_F)
      if (filter_info == filter_info_both) then
        gzip_available = .true.
      end if
    end if
    ! create new file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdfile, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5fcreate_f"
  end subroutine

  subroutine enable_gzip(dcpl, len)
    implicit none
    INTEGER(HID_T), intent(in) :: dcpl ! Handle
    integer(int8) , intent(in) :: len
    INTEGER(HSIZE_T)           :: chunk(1)
    integer hdferr

    ! Do not compress really small arrays: this creates much more overhead than
    ! compression saves
    if (gzip_available .and. len > 512) then
      call h5pset_deflate_f(dcpl, 1, hdferr)
      if (hdferr .ne. 0) then
        write(*,*) "error enabling compression in hdf5"
      end if
      chunk(1) = len
      call h5pset_chunk_f(dcpl, 1, chunk, hdferr)
      if (hdferr .ne. 0) then
        write(*,*) "error setting chunk size in hdf5"
      end if
    end if
  end subroutine

  subroutine checkpoint_start()
    use io
    implicit none
    call checkpoint_start_to_file(trim(adjustl(outdir))//'/'//trim(adjustl(hdfilename)))
  end subroutine

  ! Initiate hdf5 writing, to be called before writing anything
  function restart_start_from_file(filename)
    implicit none
    logical                  :: restart_start_from_file
    character(*), intent(in) :: filename

    restart_start_from_file = .true.
    ! open interface
    call h5open_f(hdferr)
    if (hdferr .ne. 0) then
      write(*,*) "error in h5open_f"
      restart_start_from_file = .false.
    end if
    ! create new file
    call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDONLY_F, hdfile, hdferr)
    if (hdferr .ne. 0) then
      write(*,*) "error in h5fcreate_f"
      restart_start_from_file = .false.
    end if
  end function

  function restart_start()
    use io
    implicit none
    logical :: restart_start
    restart_start = restart_start_from_file(trim(adjustl(outdir))//'/'//trim(adjustl(hdfilename)))
  end function


  ! Finalize hdf5 writing, to be called once everything is written and done
  subroutine checkpoint_end()
    call h5fclose_f(hdfile , hdferr)
  end subroutine

  ! Finalize hdf5 reading, to be called once everything is read
  subroutine restart_end()
    call h5fclose_f(hdfile , hdferr)
  end subroutine


  ! save a one-dimensional real array
  subroutine checkpoint_real_array(name, array, len)
    implicit none

    character(len=*), intent(in) :: name
    integer(int8),    intent(in) :: len
    real(dp),         intent(in) :: array(len)

    INTEGER(HID_T)  :: fspace, mspace, dset, dcpl ! Handles
    INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, start, count, stride, block

    dims  = len

    ! create file dataspace
    call h5screate_simple_f(1, dims, fspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create memory dataspace
    call h5screate_simple_f(1, dims, mspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create dataset propery list
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
    call enable_gzip(dcpl, len)
    if (hdferr .ne. 0) write(*,*) "error in h5pcreate_f"
    ! create dataset
    call h5dcreate_f(hdfile, trim(adjustl(name)), H5T_IEEE_F64LE, fspace, dset, hdferr, dcpl)
    if (hdferr .ne. 0) write(*,*) "error in h5dcreate_f"
    if (hdferr .ne. 0) stop

    count  = len
    start  = 0
    stride = 1
    block  = 1
    call h5sselect_hyperslab_f (mspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f diag"
    call h5sselect_hyperslab_f (fspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f fspace"

    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
      array, dims, hdferr, file_space_id=fspace, mem_space_id=mspace)
    if (hdferr .ne. 0) write(*,*) "error in h5dwrite_f"

    ! close and release resources
    call h5pclose_f(dcpl , hdferr)
    call h5dclose_f(dset , hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5sclose_f(fspace, hdferr)
  end subroutine


  ! write a single real(dp) to disk; use checkpoint_real_array()
  subroutine checkpoint_real_scalar(name, scalar)
    implicit none

    character(len=*), intent(in) :: name
    real(dp),         intent(in) :: scalar

    real(dp)      :: scalar_array(1)
    integer(int8) :: scalar_len

    scalar_array = scalar
    scalar_len   = 1

    call checkpoint_real_array(name, scalar_array, scalar_len)
  end subroutine


  ! read a one-dimensional real array
  subroutine restart_real_array(name, array, len)
    implicit none

    character(len=*), intent(in)  :: name
    integer(int8),    intent(in)  :: len
    real(dp),         intent(out) :: array(len)

    INTEGER(HID_T)  :: fspace, mspace, dset ! Handles
    INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, start, count, stride, block

    dims  = len

    ! create file dataspace
    call h5screate_simple_f(1, dims, fspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create memory dataspace
    call h5screate_simple_f(1, dims, mspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create dataset propery list
    call h5dopen_f(hdfile, trim(adjustl(name)), dset, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5dopen_f"
    ! create dataset
    if (hdferr .ne. 0) write(*,*) "error in h5dcreate_f"
    if (hdferr .ne. 0) stop

    count  = len
    start  = 0
    stride = 1
    block  = 1
    call h5sselect_hyperslab_f (mspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f diag"
    call h5sselect_hyperslab_f (fspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f fspace"

    call h5dread_f(dset, H5T_NATIVE_DOUBLE, &
      array, dims, hdferr, file_space_id=fspace, mem_space_id=mspace)
    if (hdferr .ne. 0) write(*,*) "error in h5dread_f"

    ! close and release resources
    call h5dclose_f(dset , hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5sclose_f(fspace, hdferr)
  end subroutine


  ! write a single real(dp) to disk; use checkpoint_real_array()
  subroutine restart_real_scalar(name, scalar)
    implicit none

    character(len=*), intent(in)  :: name
    real(dp),         intent(out) :: scalar

    real(dp)    :: scalar_array(1)
    integer(int8) :: scalar_len

    scalar_len   = 1

    call restart_real_array(name, scalar_array, scalar_len)

    scalar = scalar_array(1)
  end subroutine


  ! save a one-dimensional integer array
  subroutine checkpoint_int_array(name, array, len)
    use iso_c_binding
    implicit none

    character(len=*),      intent(in) :: name
    integer(int8),         intent(in) :: len
    integer(int8), target, intent(in) :: array(len)

    INTEGER(HID_T)  :: fspace, mspace, dset, dcpl ! Handles
    INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, start, count, stride, block

    dims  = len

    ! create file dataspace
    call h5screate_simple_f(1, dims, fspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create memory dataspace
    call h5screate_simple_f(1, dims, mspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create dataset propery list
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
    call enable_gzip(dcpl, len)
    if (hdferr .ne. 0) write(*,*) "error in h5pcreate_f"
    ! create dataset
    call h5dcreate_f(hdfile, trim(adjustl(name)), h5kind_to_type(int8, H5_INTEGER_KIND), &
                     fspace, dset, hdferr, dcpl)
    if (hdferr .ne. 0) write(*,*) "error in h5dcreate_f"
    if (hdferr .ne. 0) stop

    count  = len
    start  = 0
    stride = 1
    block  = 1
    call h5sselect_hyperslab_f (mspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f diag"
    call h5sselect_hyperslab_f (fspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f fspace"

    call h5dwrite_f(dset, h5kind_to_type(int8, H5_INTEGER_KIND), C_LOC(array(1)), hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5dwrite_f"

    ! close and release resources
    call h5pclose_f(dcpl , hdferr)
    call h5dclose_f(dset , hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5sclose_f(fspace, hdferr)
  end subroutine


  ! write a single integer(int8) to disk; use checkpoint_int_array()
  subroutine checkpoint_int_scalar(name, scalar)
    implicit none

    character(len=*), intent(in) :: name
    integer(int8),    intent(in) :: scalar

    integer(int8) :: scalar_array(1)
    integer(int8) :: scalar_len

    scalar_array = scalar
    scalar_len   = 1

    call checkpoint_int_array(name, scalar_array, scalar_len)
  end subroutine


  ! read a one-dimensional int array
  subroutine restart_int_array(name, array, len)
    use iso_c_binding
    implicit none

    character(len=*),      intent(in)  :: name
    integer(int8),         intent(in)  :: len
    integer(int8), target, intent(out) :: array(len)

    INTEGER(HID_T)  :: fspace, mspace, dset ! Handles
    INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, start, count, stride, block
    TYPE(C_PTR) :: f_ptr

    dims  = len

    ! create file dataspace
    call h5screate_simple_f(1, dims, fspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create memory dataspace
    call h5screate_simple_f(1, dims, mspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create dataset propery list
    call h5dopen_f(hdfile, trim(adjustl(name)), dset, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5dopen_f"
    ! create dataset
    if (hdferr .ne. 0) write(*,*) "error in h5dcreate_f"
    if (hdferr .ne. 0) stop

    count  = len
    start  = 0
    stride = 1
    block  = 1
    call h5sselect_hyperslab_f (mspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f diag"
    call h5sselect_hyperslab_f (fspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f fspace"

    f_ptr = C_LOC(array(1))
    call h5dread_f(dset, h5kind_to_type(int8, H5_INTEGER_KIND), f_ptr, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5dread_f"

    ! close and release resources
    call h5dclose_f(dset , hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5sclose_f(fspace, hdferr)
  end subroutine


  ! read a single integer(int8) to disk; uses restart_int_array()
  subroutine restart_int_scalar(name, scalar)
    implicit none

    character(len=*), intent(in)  :: name
    integer(int8),    intent(out) :: scalar

    integer(int8) :: scalar_array(1)
    integer(int8) :: scalar_len

    scalar_len   = 1

    call restart_int_array(name, scalar_array, scalar_len)

    scalar = scalar_array(1)
  end subroutine


  ! save a one-dimensional character array
  subroutine checkpoint_char_array(name, array, len)
    use iso_c_binding
    implicit none

    character(len=*),      intent(in) :: name
    integer(int8),         intent(in) :: len
    character, target,     intent(in) :: array(len)

    INTEGER(HID_T)  :: fspace, mspace, dset, dcpl ! Handles
    INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, start, count, stride, block

    dims  = len

    ! create file dataspace
    call h5screate_simple_f(1, dims, fspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create memory dataspace
    call h5screate_simple_f(1, dims, mspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create dataset propery list
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
    call enable_gzip(dcpl, len)
    if (hdferr .ne. 0) write(*,*) "error in h5pcreate_f"
    ! create dataset
    call h5dcreate_f(hdfile, trim(adjustl(name)), H5T_NATIVE_CHARACTER, fspace, dset, hdferr, dcpl)
    if (hdferr .ne. 0) write(*,*) "error in h5dcreate_f"
    if (hdferr .ne. 0) stop

    count  = len
    start  = 0
    stride = 1
    block  = 1
    call h5sselect_hyperslab_f (mspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f diag"
    call h5sselect_hyperslab_f (fspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f fspace"

    call h5dwrite_f(dset, H5T_NATIVE_CHARACTER, C_LOC(array(1)), hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5dwrite_f"

    ! close and release resources
    call h5pclose_f(dcpl , hdferr)
    call h5dclose_f(dset , hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5sclose_f(fspace, hdferr)
  end subroutine


  ! read a one-dimensional character array
  subroutine restart_char_array(name, array, len)
    use iso_c_binding
    implicit none

    character(len=*),  intent(in)  :: name
    integer(int8),     intent(in)  :: len
    character, target, intent(out) :: array(len)

    INTEGER(HID_T)  :: fspace, mspace, dset ! Handles
    INTEGER(HSIZE_T), DIMENSION(1:1) :: dims, start, count, stride, block
    TYPE(C_PTR) :: f_ptr

    dims  = len

    ! create file dataspace
    call h5screate_simple_f(1, dims, fspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create memory dataspace
    call h5screate_simple_f(1, dims, mspace, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5screate_simple_f"
    ! create dataset propery list
    call h5dopen_f(hdfile, trim(adjustl(name)), dset, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5dopen_f"
    ! create dataset
    if (hdferr .ne. 0) write(*,*) "error in h5dcreate_f"
    if (hdferr .ne. 0) stop

    count  = len
    start  = 0
    stride = 1
    block  = 1
    call h5sselect_hyperslab_f (mspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f diag"
    call h5sselect_hyperslab_f (fspace, H5S_SELECT_SET_F, start, count, hdferr, stride, block)
    if (hdferr .ne. 0) write(*,*) "error in h5sselect_hyperslab_f fspace"

    f_ptr = C_LOC(array(1))
    call h5dread_f(dset, H5T_NATIVE_CHARACTER, f_ptr, hdferr)
    if (hdferr .ne. 0) write(*,*) "error in h5dread_f"

    ! close and release resources
    call h5dclose_f(dset , hdferr)
    call h5sclose_f(mspace, hdferr)
    call h5sclose_f(fspace, hdferr)
  end subroutine

end module checkpointing

