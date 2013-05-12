module master_function_mod

  use real_type_mod
  use variables_mod
  use parameters_mod

  implicit none

contains

  subroutine MasterFunction_Derivs(n_sqd, dLdn)

    real(kind=wp), dimension(NFF), intent(in) :: n_sqd
    real(kind=wp), dimension(NFF), intent(out) :: dLdn

    real(kind=wp), dimension(NF) :: n
    real(kind=wp) :: nxy, delta
    
!!$    Again, this should be generalized
    n = sqrt(n_sqd(1:NF))
!!$    Note that nxy = n_{12}^2 in the code, but we want n_{12}
    nxy = sqrt(n_sqd(NFF))
    
    delta = sqrt(1.0_wp - (n(1) * n(2) / nxy)**2)

    if (master_function_i .eq. M_FN_GENERIC_PAPER) then
      
      dLdn(1) = -(1.d0 / nxy**2 / n(1) * &
           (eos_m(1) * nxy**2 + &
            eos_kappa(1) * n(1)**(eos_gamma(1) - 1.d0) * eos_gamma(1) * nxy**2 + &
            eos_kappa_np * n(1)**(eos_sigma(1) - 1.d0) * eos_sigma(1) * n(2)**eos_sigma(2) * nxy**2 - &
            2.d0 * eos_kappa_delta * n(1) * n(2)**2)) / 2.d0

      dLdn(2) = -(1.d0 / nxy**2 / n(2) * &
           (eos_m(2) * nxy**2 + &
            eos_kappa(2) * n(2)**(eos_gamma(2) - 1.d0) * eos_gamma(2) * nxy**2 + &
            eos_kappa_np * n(1)**eos_sigma(1) * n(2)**(eos_sigma(2) - 1.d0) * eos_sigma(2) * nxy**2 - &
            2.d0 * eos_kappa_delta * n(1)**2 * n(2))) / 2.d0

      dLdn(3) = -eos_kappa_delta * n(1)**2 * n(2)**2 / nxy**4

    else

      write(*,*) 'Have not coded this master function case!'
      stop

    end if

  end subroutine MasterFunction_Derivs

  subroutine MasterFunction_SecondDerivs(n_sqd, d2Ldn2, d2Ldn2_2)

    real(kind=wp), dimension(NFF), intent(in) :: n_sqd
    real(kind=wp), intent(out) :: d2Ldn2
    real(kind=wp), dimension(NF), intent(out) :: d2Ldn2_2

    real(kind=wp), dimension(NF) :: n
    real(kind=wp) :: nxy, delta

    integer :: f
    
!!$    Again, this should be generalized
    n = sqrt(n_sqd(1:NF))
    nxy = sqrt(n_sqd(NFF))
    
    delta = sqrt(1.0_wp - (n(1) * n(2) / nxy)**2)

!!$    This is the mixed derivative of Lambda wrt n1n2, 
!!$    also called dB^x / dn_y^2 in the relativistic two stream paper

    if (master_function_i .eq. M_FN_GENERIC_PAPER) then

      d2Ldn2 = -(1.d0 / nxy**2 / n(1) / n(2) * &
           (eos_kappa_np * n(1)**(eos_sigma(1) - 1.d0) * &
            eos_sigma(1) * n(2)**(eos_sigma(2) - 1.d0) * eos_sigma(2) * nxy**2 - &
            4.d0 * eos_kappa_delta * n(1) * n(2))) / 4.d0

      d2Ldn2_2(1) = -eos_kappa_delta * n(2)**2 / nxy**4
      d2Ldn2_2(2) = -eos_kappa_delta * n(1)**2 / nxy**4

    else

      write(*,*) 'Have not coded this master function case! (Second derivs)'
      stop

    end if

  end subroutine MasterFunction_SecondDerivs

  subroutine force_terms(q, n_up, mu_down, force)

    real(kind=wp), dimension(NEVOLVE), intent(in) :: q
    real(kind=wp), dimension(NF, NV), intent(in) :: n_up, mu_down
    real(kind=wp), dimension(NEVOLVE), intent(out) :: force

    force = 0.0_wp

  end subroutine force_terms

end module master_function_mod
