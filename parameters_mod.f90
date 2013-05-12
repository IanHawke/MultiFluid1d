module parameters_mod

  use real_type_mod
  use variables_mod

  implicit none

!!$  Parameters to be set (later)

  real(kind=wp) :: &
       xmin = 0.0_wp, xmax = 1.0_wp, &
       courant = 0.2_wp, &
       t_end = 0.4_wp, &
       out_every_dt = -1.0_wp, &
       fft_out_every_dt = -1.0_wp

  real(kind=wp), dimension(2) :: &
       initial_n = 1.0_wp, &
       initial_v = 0.0_wp, &
       initial_ptbn_n = 0.0_wp, &
       initial_ptbn_v = 0.0_wp, &
       initial_period_n = 1.0_wp, &
       initial_period_v = 1.0_wp, &
       initial_phase_n = 0.0_wp, &
       initial_phase_v = 0.0_wp

  real(kind=wp), dimension(2) :: &
       eos_m = 1.0, &
       eos_kappa = 0.12, &
       eos_gamma = 1.5, &
       eos_gamma_delta = 1.0, &
       eos_sigma = 1.1_wp
  real(kind=wp) :: &
       eos_kappa_delta = 0.0, &
       eos_kappa_np = 0.0, &
       eos_kappa_ns = 1.0_wp, &
       eos_m_ns = 1.0_wp, &
       eos_alpha = 1.0_wp, &
       eos_beta = 1.33333333333333333333333333_wp, &
       eos_tau = 0.1_wp
  real(kind=wp) :: &
       eos_lambda_0 = 1.0, &
       eos_lambda_1 = 0.1

  real(kind=wp) :: &
       eps_dissipation = 0.0

  integer :: &
       np = 100, &
       gz = 3, &
       rk_order = 3, &
       differencing_order = 4, &
       dissipation_order = 0, &
       out_every_n = 1, &
       fft_out_every_n = 1

  character(len=200) :: &
       out_dir_base = 'results', &
       perturbation_type = 'trig', &
       master_function_type = 'multi_polytrope'

  logical :: &
       out_dir_date = .true., &
       do_fft = .false., &
       vtk_output = .false., &
       vtk_output_binary = .false.

!!$  Derived parameters

  real(kind=wp) :: dx, dt, I2DX, I12DX, I16DX, I64DX, fft_out_dt
  integer :: out_every, fft_out_every, fft_out_total, master_function_i

contains

  subroutine set_parameters(parameter_file)

    character(len=200), intent(in) :: parameter_file

    namelist /parameters/ &
         xmin, xmax, &
         courant, &
         t_end, &
         initial_n, initial_v, &
         initial_ptbn_n, initial_ptbn_v, &
         initial_period_n, initial_period_v, &
         initial_phase_n, initial_phase_v, &
         eos_m, eos_kappa, eos_gamma, eos_gamma_delta, &
         eos_kappa_delta, eos_kappa_np, &
         eos_alpha, eos_tau, eos_sigma, eos_beta, eos_kappa_ns, eos_m_ns, &
         out_every_dt, out_every_n, &
         fft_out_every_dt, fft_out_every_n, &
         np, gz, &
         rk_order, differencing_order, dissipation_order, &
         eps_dissipation, &
         out_dir_base, out_dir_date, &
         perturbation_type, &
         master_function_type, &
         do_fft, &
         vtk_output, vtk_output_binary

    open(3, file=trim(parameter_file), STATUS='old')
    read(3, nml=parameters)
    
    dx = (xmax - xmin) / dble(np)
    dt = courant * dx
    I2DX  = 1.0_wp / ( 2.0_wp * dx)
    I12DX = 1.0_wp / (12.0_wp * dx)
    I16DX = 1.0_wp / (16.0_wp * dx)
    I64DX = 1.0_wp / (64.0_wp * dx)

    if ( (differencing_order > 2 * gz    ).or.&
         (dissipation_order  > 2 * gz - 1) ) then
      write(*,*) 'Number of ghostzones (gz) too small'
      write(*,'(a10,i2)') 'Number is ', gz
      write(*,'(a29,i2,a24,i2)') 'To match differencing order (', &
           differencing_order, &
           ') gz should be at least ', differencing_order / 2
      write(*,'(a28,i2,a24,i2)') 'To match dissipation order (', &
           dissipation_order, &
           ') gz should be at least ', (dissipation_order + 1) / 2
      stop
    end if

    if (out_every_dt > 0) then
      out_every = floor(out_every_dt / dt)
    else if (out_every_n > 0) then
      out_every = out_every_n
    else
      out_every = 0
    end if
    if (fft_out_every_dt > 0) then
      fft_out_every = floor(fft_out_every_dt / dt)
    else if (fft_out_every_n > 0) then
      fft_out_every = fft_out_every_n
    else
      fft_out_every = 0
    end if
    fft_out_total = ceiling(t_end / dt) / fft_out_every + 1
    fft_out_dt = dble(fft_out_every) * dt
    write(*, *) 'Output every ', out_every, ' steps.'
    write(*, *) 'FFT Output every ', fft_out_every, ' steps.'

    select case (master_function_type)

    case ('multi_polytrope') !!$ Prix, Comer et al

      master_function_i = M_FN_MULTI_POLY
      write(*,*) 'Case now absorbed in ''generic_paper''.'
      stop

    case ('cosmological_poly') !!$ Comer, Peter and Andersson

      master_function_i = M_FN_COSMOLOGICAL_POLY
      write(*,*) 'Case now absorbed in ''generic_paper''.'
      stop

    case ('generic_paper') !!$ Hawke, Comer and Andersson eq (20)
!!$      Note that the parameters listed in the paper include
!!$      kappa_{12} which here is eos_kappa_np

      master_function_i = M_FN_GENERIC_PAPER

    case default

      write(*,*) 'Value of master_function_type parameter not recognized: ', &
           master_function_type
      stop

    end select
  
  end subroutine set_parameters
  
end module parameters_mod
