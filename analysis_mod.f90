module analysis_mod

  use real_type_mod
  use variables_mod
  use output_mod
  use parameters_mod
  use compute_fft_mod
  use master_function_mod

  implicit none

contains

  subroutine analysis_create_all(out_dir)

    character(len=200), intent(in) :: out_dir

    real(kind=wp), dimension(:), allocatable :: freq
    integer :: istat, i
    real(kind=wp) :: pi

    pi = 4.0_wp * atan(1.0_wp)

    allocate(freq(np/2 + 1), &
             STAT = istat)
    if (istat .ne. 0) then
      write(*,*) 'Allocation error in analysis_create_all'
      stop
    end if
    do i = 1, np/2 + 1
      freq(i) = 2.0_wp * pi * dble(i-1) / (xmax - xmin)
    end do

    call create_output_file(out_dir, 'n_up',    N_UP_FILE)
    call create_output_file(out_dir, 'n_down',  N_DOWN_FILE)
    call create_output_file(out_dir, 'mu_up',   MU_UP_FILE)
    call create_output_file(out_dir, 'mu_down', MU_DOWN_FILE)
    call create_output_file(out_dir, 'n_sqd',   N_SQD_FILE)
    call create_output_file(out_dir, 'v_up',    V_UP_FILE)
!!$    call create_output_file(out_dir, 'c_cc',    C_CC_FILE)
!!$    call create_output_file(out_dir, 'c_ent',   C_ENT_FILE)
!!$    call create_output_file(out_dir, 'cs2',     CS2_FILE)
    call create_output_file(out_dir, 'constraints', CONSTRAINTS_FILE)

    if (do_fft) then
      if (vtk_output) then
        call create_output_timefreq_vtk_file(out_dir, 'fft', FFT_FILE, &
             np, freq, fft_out_total, fft_out_dt)
      else
        call create_output_file(out_dir, 'fft',    FFT_FILE)
      end if
    end if

    deallocate(freq)

  end subroutine analysis_create_all

  subroutine analysis_close_all()

    close(N_UP_FILE   )
    close(N_DOWN_FILE )
    close(MU_UP_FILE  )
    close(MU_DOWN_FILE)
    close(N_SQD_FILE  )
    close(V_UP_FILE   )
!!$    close(C_CC_FILE   )
!!$    close(C_ENT_FILE  )
!!$    close(CS2_FILE    )
    close(CONSTRAINTS_FILE)

    if (do_fft) then
      close(FFT_FILE)
    end if
    
  end subroutine analysis_close_all

  subroutine analysis_output_all(x, n_up, n_down, mu_up, mu_down, mu_down_p, &
       n_sqd, t, dt, dx)

    real(kind=wp), dimension(1-gz:np+gz), intent(in) :: x
    real(kind=wp), dimension(NF, NV, 1-gz:np+gz), intent(in) :: &
         n_up, n_down, mu_up, mu_down, mu_down_p
    real(kind=wp), dimension(NFF, 1-gz:np+gz), intent(in) :: &
         n_sqd
    real(kind=wp), intent(in) :: t, dt, dx

    real(kind=wp), dimension(:, :), allocatable :: v_up, cs2, constraints
    real(kind=wp), dimension(:), allocatable :: c_cc, c_ent
    real(kind=wp), dimension(NFF) :: dLdn
    real(kind=wp) :: d2Ldn2
    real(kind=wp), dimension(NF) :: d2Ldn2_2

    integer :: istat, f, f1, f2, i

    allocate(v_up(NFF, 1-gz:np+gz), &
             cs2 (2, 1-gz:np+gz), &
             c_cc(     1-gz:np+gz), &
             c_ent(    1-gz:np+gz), &
             constraints(NF, 1-gz:np+gz), &
             STAT=istat)

    call output2(N_UP_FILE   , x, n_up   , NF, NV, np, gz, t)
    call output2(N_DOWN_FILE , x, n_down , NF, NV, np, gz, t)
    call output2(MU_UP_FILE  , x, mu_up  , NF, NV, np, gz, t)
    call output2(MU_DOWN_FILE, x, mu_down, NF, NV, np, gz, t)
    call output (N_SQD_FILE  , x, n_sqd  , NFF   , np, gz, t)

    do f = 1, NF
      v_up(f, :) = n_up(f, 2, :) / n_up(f, 1, :)
    end do
    f = NF
    do f1 = 1, NF
      do f2 = f1 + 1, NF
        f = f + 1
        v_up(f, :) = sqrt(1.0_wp - &
             n_sqd(f1, :) * n_sqd(f2, :) / n_sqd(f, :)**2)
      end do
    end do

    call output (V_UP_FILE   , x, v_up   , NFF   , np, gz, t)

!!$    These quantities are from Lars et al's two stream instability paper
!!$    but in the end did not turn out to be particularly useful
!!$    do i = 1-gz, np+gz
!!$      call MasterFunction_Derivs(n_sqd, dLdn)
!!$      call MasterFunction_SecondDerivs(n_sqd(:, i), d2Ldn2, d2Ldn2_2)
!!$      c_cc(i) = 1.0_wp / ( dLdn(1) * dLdn(2) ) * &
!!$           (2.0_wp * sqrt(n_sqd(1, i) * n_sqd(2, i)) * d2Ldn2)**2
!!$      c_ent(i) = dLdn(3)**2 / ( dLdn(1) * dLdn(2) )
!!$      do f1 = 1, NF
!!$        cs2(f1, i) = sqrt(n_sqd(f1, i)) / dLdn(f1) * d2Ldn2_2(f1) + 1
!!$      end do
!!$    end do
!!$
!!$    call output (C_CC_FILE   , x, c_cc   , 1     , np, gz, t)
!!$    call output (C_ENT_FILE  , x, c_ent  , 1     , np, gz, t)
!!$    call output (CS2_FILE    , x, cs2    , 2     , np, gz, t)

!!$    Compute the constraint
!!$    d_t mu^X_x = d_x mu^X_t
!!$    using first order differencing in time and 
!!$    central second order differencing in space

!!$    This can only be done after the first timestep
    if (t .ge. 0.5_wp * dt) then
      constraints = 0.0_wp
      do f = 1, NF
        do i = 1, np
          constraints(f, i) = (mu_down(f, 2, i) - mu_down_p(f, 2, i)) / dt
          constraints(f, i) = constraints(f, i) - 0.5_wp * &
               (mu_down(f, 1, i+1) - mu_down(f, 1, i-1)) / dx
        end do
      end do

      call output (CONSTRAINTS_FILE    , x, constraints, NF, np, gz, t)
    end if
    
    deallocate(v_up, c_cc, c_ent, cs2, constraints)

  end subroutine analysis_output_all

  subroutine analysis_fft(n_sqd, n_up, t)

    real(kind=wp), dimension(NF, NV, 1-gz:np+gz), intent(in) :: n_up
    real(kind=wp), dimension(NFF, 1-gz:np+gz), intent(in) :: n_sqd
    real(kind=wp), intent(in) :: t

    real(kind=wp), dimension(:), allocatable :: freq
    complex(kind=wp), dimension(:, :), allocatable :: fft
    integer :: istat, i, ff, vv
    real(kind=wp) :: pi

    pi = 4.0_wp * atan(1.0_wp)

    allocate(freq(np/2 + 1), &
             fft(NF*NV, np/2 + 1), &
             STAT = istat)
    if (istat .ne. 0) then
      write(*,*) 'Allocation error in analysis_fft'
      stop
    end if

    do ff = 1, NF
      do vv = 1, NV
        call compute_fft(np, n_up(ff, vv, 1:np), fft(NF*(ff-1)+vv, :))
      end do
    end do
    
    if (vtk_output) then
      call output_timefreq_vtk(FFT_FILE, fft, np)
    else
      do i = 1, np/2 + 1
        freq(i) = 2.0_wp * pi * dble(i-1) / (xmax - xmin)
      end do
      call output3(FFT_FILE, freq, fft, np, NV, NF, t)
    end if
    
    deallocate(freq, fft)

  end subroutine analysis_fft

end module analysis_mod
