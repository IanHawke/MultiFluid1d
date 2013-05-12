module evolve_mod

  use real_type_mod
  use parameters_mod
  use variables_mod
  use conversion_mod
  use master_function_mod

  implicit none

contains

  subroutine evolve(n_up, n_down, mu_up, mu_down, n_sqd)

    real(kind=wp), dimension(NF, NV, 1-gz:np+gz), intent(inout) :: &
         n_up, n_down, mu_up, mu_down
    real(kind=wp), dimension(NFF, 1-gz:np+gz), intent(inout) :: n_sqd

    real(kind=wp), dimension(:, :, :), allocatable :: q_new
    real(kind=wp), dimension(:, :), allocatable :: rhs, auxl, q

    integer :: i, j, nt, istat

    allocate(q_new(rk_order+1, NEVOLVE, 1-gz:np+gz), &
             q(NEVOLVE, 1-gz:np+gz), &
             auxl(NF, 1-gz:np+gz), &
             rhs (NEVOLVE, 1-gz:np+gz), &
             STAT=istat)
    if (istat .ne. 0) then
      write(*,*) 'Allocation error in evolve'
      stop
    end if

    q_new = 0.0_wp
    rhs = 0.0_wp
    
    do i = 1, NF
      q_new(1, i, :) = n_up(i, 1, :)
      do j = 1, ND
        q_new(1, j*NF + i, :) = mu_down(i, 1+j, :)
      end do
    end do
    q = q_new(1, :, :)

    do nt = 1, rk_order

!!$    Auxilliary variables are the timelike components of mu, for the moment
      auxl = mu_down(:, 1, :)

      call compute_rhs(q, n_up, mu_down, rhs)      


      call RK_step(q, q_new, rhs, nt, rk_order)

      q = q_new(nt+1, :, :)

      call boundaries(q, NEVOLVE)
      do i = 1-gz, np+gz

        call all_from_evolved_pt(q(:, i), &
             n_up(:, :, i), n_down(:, :, i), &
             mu_up(:, :, i), mu_down(:, :, i), &
             n_sqd(:, i))
      end do
    end do
    
    deallocate(q_new, q, auxl, rhs)

  end subroutine evolve

  subroutine compute_rhs(q, n_up, mu_down, rhs)

    real(kind=wp), dimension(NEVOLVE, 1-gz:np+gz), intent(in) :: q
    real(kind=wp), dimension(NF, NV, 1-gz:np+gz), intent(in) :: n_up, mu_down
    real(kind=wp), dimension(NEVOLVE, 1-gz:np+gz), intent(out) :: rhs

    real(kind=wp), dimension(:, :), allocatable :: deriv_n, deriv_mu, &
         source, deriv_dissip_q

    integer :: i, istat, f

    rhs = 0.0_wp

    allocate(deriv_n (NF, 1-gz:np+gz), &
             deriv_mu(NF, 1-gz:np+gz), &
             source(NEVOLVE, 1-gz:np+gz), &
             deriv_dissip_q(NEVOLVE, 1-gz:np+gz), &
             STAT=istat)
    if (istat .ne. 0) then
      write(*,*) 'Allocation error in compute_rhs'
      stop
    end if

    deriv_n  = 0.0_wp
    deriv_mu = 0.0_wp
    select case (differencing_order)

      case (2)

        do i = 1, np
          do f = 1, NF
!!$        Second order central differencing
            deriv_n (f, i) = I2DX * (n_up   (f, 2, i+1) - n_up   (f, 2, i-1))
            deriv_mu(f, i) = I2DX * (mu_down(f, 1, i+1) - mu_down(f, 1, i-1))
          end do
        end do
        
      case (4)

        do i = 1, np
          do f = 1, NF
!!$        Fourth order central differencing
            deriv_n (f, i) = I12DX * &
                 (-n_up(f, 2, i+2) + n_up   (f, 2, i-2) + &
         8.0_wp * (n_up(f, 2, i+1) - n_up   (f, 2, i-1)) )
            deriv_mu(f, i) = I12DX * &
                 (-mu_down(f, 1, i+2) + mu_down(f, 1, i-2) + &
         8.0_wp * (mu_down(f, 1, i+1) - mu_down(f, 1, i-1)) )
          end do
        end do
        
      case default

        write(*,*) 'Unknown value of ''differencing order''', &
             differencing_order
        stop

    end select

    deriv_dissip_q = 0.0_wp        
    select case (dissipation_order)

      case (0)

        deriv_dissip_q = 0.0_wp

      case (3)

        do i = 1, np
!!$      Third order dissipation
          deriv_dissip_q (:, i) = I16DX * &
               (          q(:, i+2) + q(:, i-2)  - &
                4.0_wp * (q(:, i+1) + q(:, i-1)) + &
                6.0_wp *  q(:, i)                  )
        end do
        
      case (5)

        do i = 1, np
!!$      Fifth order dissipation
          deriv_dissip_q (:, i) = I64DX * &
               (          q(:, i+3) + q(:, i-3)  - &
                6.0_wp * (q(:, i+2) + q(:, i-2)) + &
               15.0_wp * (q(:, i+1) + q(:, i-1)) - &
               20.0_wp *  q(:, i)                  )
        end do
        
      case default

        write(*,*) 'Unknown value of ''dissipation order''', &
             dissipation_order
        stop

    end select
    
    source = 0.0_wp
!!$      Source terms due to inter-fluid or external forces
    do i = 1, np
      call force_terms(q(:, i), n_up(:, :, i), mu_down(:, :, i), source(:, i))
    end do

!!$    Note that we should add a source term here for the external forces
!!$    Need to consider the mu term correctly - divide by n
    do i = 1, np
      do f = 1, NF
        rhs(f     , i) = -deriv_n (f, i) + source(     f, i)
        rhs(NF + f, i) =  deriv_mu(f, i) + source(NF + f, i)
      end do
      rhs(:, i) = rhs(:, i) + eps_dissipation * deriv_dissip_q(:, i)
    end do

    deallocate(deriv_n, deriv_mu, source, deriv_dissip_q)
    
  end subroutine compute_rhs
  
  subroutine rk_step(q, q_new, rhs, n, order)

    integer, intent(in) :: n, order
    real(kind=wp), dimension(NEVOLVE, 1-gz:np+gz), intent(in) :: q, rhs
    real(kind=wp), dimension(order + 1, NEVOLVE, 1-gz:np+gz), intent(out) :: &
         q_new

    integer :: i

    select case (order)

      case (1)

        if (n .ne. 1) then
          write(*,*) 'In RK_step, order 1, n not 1!? ', n
          stop
        end if

        q_new(2, :, :) = q_new(1, :, :) + dt * rhs

      case (2)

        select case (n)

          case (1)

            q_new(2, :, :) = q_new(1, :, :) + dt * rhs

          case (2)

            q_new(3, :, :) = 0.5_wp * &
                 (q_new(1, :, :) + q_new(2, :, :) + dt * rhs)

          case default

            write(*,*) 'In RK_step, order 2, n not 1, 2!? ', n
            stop
            
        end select
          
      case (3)

        select case (n)

          case (1)

            q_new(2, :, :) = q_new(1, :, :) + dt * rhs

          case (2)

            q_new(3, :, :) = 0.25_wp * &
                 (3.0_wp * q_new(1, :, :) + q_new(2, :, :) + dt * rhs)

          case (3)

            q_new(4, :, :) = 1.0_wp / 3.0_wp * &
                 (q_new(1, :, :) + 2.0_wp * q_new(3, :, :) + &
                 2.0_wp * dt * rhs)

          case default

            write(*,*) 'In RK_step, order 3, n not 1-3!? ', n
            stop
            
        end select

      case (4)

        select case (n)

          case (1)

            q_new(2, :, :) = q_new(1, :, :) + 0.5_wp * dt * rhs

          case (2)

            q_new(3, :, :) = q_new(1, :, :) + 0.5_wp * dt * rhs

          case (3)

            q_new(4, :, :) = q_new(1, :, :) + dt * rhs

          case (4)

            q_new(5, :, :) = 1.0_wp / 6.0_wp * &
                 ( -2.0_wp * q_new(1, :, :) + &
                    2.0_wp * q_new(2, :, :) + &
                    4.0_wp * q_new(3, :, :) + &
                    2.0_wp * q_new(4, :, :) + &
                    dt * rhs )

          case default

            write(*,*) 'In RK_step, order 4, n not 1-4!? ', n
            stop
            
        end select

      case default

        write(*,*) 'Expected the RK order to be in 1-4; found ', order
        stop

    end select
    
  end subroutine rk_step
  
  subroutine boundaries(q, nv)

    integer, intent(in) :: nv
    real(kind=wp), dimension(nv, 1-gz:np+gz), intent(inout) :: q

    integer :: i

!!$    Periodic

    do i = 1, gz
      q(:, 1-i ) = q(:, np-i+1)
      q(:, np+i) = q(:, i)
    end do
    
  end subroutine boundaries
  
end module evolve_mod
