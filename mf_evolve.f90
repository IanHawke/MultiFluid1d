program mf_evolve

  use real_type_mod
  use output_mod
  use variables_mod
  use parameters_mod
  use evolve_mod
  use conversion_mod
  use initial_mod
  use analysis_mod
  use system_mod

  implicit none

!!$  Variables

  real(kind=wp), dimension(:), allocatable :: x
  real(kind=wp), dimension(:, :, :), allocatable :: n_up, mu_down, &
       n_down, mu_up, mu_down_p
  real(kind=wp), dimension(:, :), allocatable :: n_sqd
  
  integer :: i, istat, nt
  character(len=200) :: out_dir, parameter_file
  real(kind=wp) :: t

!!$  Read the command line and find out the input file

  if (iargc() .ne. 1) then
    write(*,*) 'Useage: mf_evolve <parameter file>'
    stop
  end if
  call getarg(1, parameter_file)
  
!!$  Set the parameters

  call set_parameters(parameter_file)

!!$  Open output files, backup parameter file and code

  out_dir = make_out_dir(out_dir_base)
  call execute_command(trim('cp '//trim(parameter_file)//' '//trim(out_dir)), istat)
  if (istat .ne. 0) then
    write(*,*) 'Copying the parameter file failed, ', istat
    stop
  end if
  call execute_command(trim('tar zcf '//trim(out_dir)//'/code.tgz Makefile *.f90'), &
       istat)
  if (istat .ne. 0) then
    write(*,*) 'Backing up code failed, ', istat
    stop
  end if
  call analysis_create_all(out_dir)
  if (do_fft) then
    call fft_setup(np)
  end if

!!$  Set up arrays

  allocate(n_up   (NF, NV, 1-gz:np+gz), &
           n_down (NF, NV, 1-gz:np+gz), &
           mu_up  (NF, NV, 1-gz:np+gz), &
           mu_down(NF, NV, 1-gz:np+gz), &
           mu_down_p(NF, NV, 1-gz:np+gz), &
           n_sqd  (NFF   , 1-gz:np+gz), &
           STAT=istat)

  allocate(x(1-gz:np+gz), &
           STAT=istat)

!!$  Set up initial data

  call initial(x, n_up)
  call all_from_n_up(n_up, n_down, mu_up, mu_down, n_sqd)
  mu_down_p = mu_down

  t = 0.0_wp
  
!!$  Output initial data

  call analysis_output_all(x, n_up, n_down, mu_up, mu_down, mu_down_p, n_sqd, &
       t, dt, dx)
  if (do_fft) then
    call analysis_fft(n_sqd, n_up, t)
  end if
  
!!$  Evolve

  nt = 0

  do while (t < t_end)

    mu_down_p = mu_down

    call evolve(n_up, n_down, mu_up, mu_down, n_sqd)
    nt = nt + 1
    t = t + dt

    if ( mod(nt, out_every) == 0 ) then
      write(*, *) 'Step ', nt, ' Time ', t
      call analysis_output_all(x, n_up, n_down, mu_up, mu_down, mu_down_p, &
           n_sqd, t, dt, dx)
    end if
!!$    Do the FFT at every step if requested
    if (do_fft) then
      if ( mod(nt, fft_out_every) == 0 ) then
        call analysis_fft(n_sqd, n_up, t)
      end if
    end if

  end do

!!$  End of program - close all files and deallocate

  call analysis_close_all()
  if (do_fft) then
    call fft_shutdown()
  end if

  deallocate(n_up, n_down, mu_up, mu_down, mu_down_p, n_sqd)
  
end program mf_evolve
