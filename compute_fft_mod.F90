!!$Link to FFTW to compute the Fourier transforms.
!!$Irritatingly Intel fortran (at least v 12.1) does not seem to like 
!!$.f03 files, nor the standards, so link to the old style interface.

module compute_fft_mod

  use real_type_mod

#ifdef __INTEL_COMPILER

#else

  use, intrinsic :: iso_c_binding
  
#endif

  implicit none

#ifdef __INTEL_COMPILER

  include 'fftw3.f'

  real(kind=wp), dimension(:), allocatable :: in
  complex(kind=wp), dimension(:), allocatable :: out
  integer*8 :: plan

#else

  include 'fftw3.f03'

  type(C_PTR) :: plan
  real (C_DOUBLE), dimension(:), allocatable :: in
  complex(C_DOUBLE_COMPLEX), dimension(:), allocatable :: out

#endif

contains

  subroutine fft_setup(np)

    integer, intent(in) :: np
    integer :: istat

    allocate(in (np), &
             out(np/2 + 1), &
             STAT=istat)
    if (istat .ne. 0) then
      write(*,*) 'Allocation error in fft_setup'
      stop
    end if

    call dfftw_plan_dft_r2c_1d(plan, np, in, out, &
         FFTW_MEASURE)

  end subroutine fft_setup

  subroutine fft_shutdown

    deallocate(in, out)
    call dfftw_destroy_plan(plan)

  end subroutine fft_shutdown

  subroutine compute_fft(np, q, fft_out)

    integer, intent(in) :: np
    real(kind=wp), dimension(np), intent(in) :: q
    complex(kind=wp), dimension(np/2 + 1), intent(out) :: fft_out

    in = q
    call dfftw_execute_dft_r2c(plan, in, out)
    fft_out = out

  end subroutine compute_fft

end module compute_fft_mod
