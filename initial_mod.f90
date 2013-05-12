module initial_mod

  use real_type_mod
  use parameters_mod
  use variables_mod
  use conversion_mod

  implicit none

contains

  subroutine initial(x, n_up)

    real(kind=wp), dimension(1-gz:np+gz), intent(out) :: x
    real(kind=wp), dimension(NF, NV, 1-gz:np+gz), intent(out) :: n_up

    real(kind=wp) :: n, v, w
    real(kind=wp) :: pi = 4.d0 * atan(1.d0)

    integer :: i, f

    do i = 1-gz, np+gz
      x(i) = (dble(i) - 0.5_wp) * dx + xmin

!!$      Primitive variables; simple trig functions for now
!!$      This is essentially random - something more useful should be found

      do f = 1, NF

        select case (perturbation_type)

          case ('trig')

            n = initial_n(f) * (1.0_wp + &
                 initial_ptbn_n(f) * sin( pi * initial_period_n(f) * &
                                          (x(i) - initial_phase_n(f)) ) )
            v = initial_v(f) + &
                 initial_ptbn_v(f) * sin( pi * initial_period_v(f) * &
                                          (x(i) - initial_phase_v(f)) ) 

          case ('gaussian')

            n = initial_n(f) * (1.0_wp + &
                 initial_ptbn_n(f) * &
                 exp(-(x(i) - 0.5_wp * (xmin + xmax))**2 / 0.05_wp**2) )
            v = initial_v(f) + initial_ptbn_v(f) * &
                 exp(-(x(i) - 0.5_wp * (xmin + xmax))**2 / 0.15_wp**2) 

          case ('doublegaussian')

            n = initial_n(f) * (1.0_wp + &
                 initial_ptbn_n(f) * &
                 exp(-(x(i) - &
                     1.0_wp / 3.0_wp * (xmin + xmax))**2 / 0.075_wp**2) + &
                 initial_ptbn_n(f) * &
                 exp(-(x(i) -  &
                     2.0_wp / 3.0_wp * (xmin + xmax))**2 / 0.075_wp**2) )
            v = initial_v(f) + initial_ptbn_v(f) * &
                 exp(-(x(i) - &
                     1.0_wp / 3.0_wp * (xmin + xmax))**2 / 0.075_wp**2) + &
                initial_ptbn_v(f) * &
                 exp(-(x(i) -  &
                     2.0_wp / 3.0_wp * (xmin + xmax))**2 / 0.075_wp**2) 

          case default

            write(*,*) &
                 'Parameter ''perturbation_type'' takes an unknown value ', &
                 perturbation_type
            stop

        end select

        w = 1.0_wp / sqrt(1.0_wp - v**2)

        n_up(f, 1, i) = n * w
        n_up(f, 2, i) = n * w * v
      
      end do
      
    end do

  end subroutine initial
  
end module initial_mod
