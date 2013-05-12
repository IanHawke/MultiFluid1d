module conversion_mod

  use real_type_mod
  use variables_mod
  use parameters_mod
  use master_function_mod

  implicit none

contains

  subroutine compute_scalars_pt(n_up, n_down, n_sqd)

    real(kind=wp), dimension(NF, NV), intent(in) :: n_up, n_down
    real(kind=wp), dimension(NFF), intent(out) :: n_sqd

    integer :: f1, f2, f, k

!!$    The scalars are -n^2. The first NF are the n_x^2; the rest are in
!!$    standard lexicographic order

    do f = 1, NF
      n_sqd(f) = -n_up(f, 1) * n_down(f, 1) 
      do k = 2, NV
        n_sqd(f) = n_sqd(f) - n_up(f, k) * n_down(f, k)
      end do
    end do
    f = NF
    do f1 = 1, NF
      do f2 = f1 + 1, NF
        f = f + 1
        n_sqd(f) = -n_up(f1, 1) * n_down(f2, 1) 
        do k = 2, NV
          n_sqd(f) = n_sqd(f) - n_up(f1, k) * n_down(f2, k)
        end do
      end do
    end do
    
  end subroutine compute_scalars_pt

  subroutine compute_scalars(n_up, n_down, n_sqd)

    real(kind=wp), dimension(NF, NV, 1-gz:np+gz), intent(in) :: n_up, n_down
    real(kind=wp), dimension(NFF, 1-gz:np+gz), intent(out) :: n_sqd

    integer :: i

!!$    Use local routine
    do i = 1-gz, np+gz
      call compute_scalars_pt(n_up(:, :, i), n_down(:, :, i), n_sqd(:, i))
    end do
    
  end subroutine compute_scalars

  subroutine compute_mu_down_pt(n_down, n_sqd, mu_down)

    real(kind=wp), dimension(NF, NV), intent(in) :: n_down
    real(kind=wp), dimension(NFF), intent(in) :: n_sqd
    real(kind=wp), dimension(NF, NV), intent(out) :: mu_down

    real(kind=wp), dimension(NFF) :: dLdn

    integer :: f, f1, f2

!!$    Compute the derivatives of the master function
    call MasterFunction_Derivs(n_sqd, dLdn)

!!$    Compute conjugate momenta
    mu_down = 0_wp
!!$    Bulk terms (B^x)
    do f = 1, NF
      mu_down(f, :) = -2.0_wp * dLdn(f) * n_down(f, :)
    end do
!!$    Other terms (A^{xy})
    f = NF
    do f1 = 1, NF
      do f2 = f1 + 1, NF
        f = f + 1
        mu_down(f1, :) = mu_down(f1, :) - dLdn(f) * n_down(f2, :)
        mu_down(f2, :) = mu_down(f2, :) - dLdn(f) * n_down(f1, :)
      end do
    end do

  end subroutine compute_mu_down_pt

  subroutine all_from_n_up_pt(n_up, n_down, mu_up, mu_down, n_sqd)

    real(kind=wp), dimension(NF, NV), intent(in) :: n_up
    real(kind=wp), dimension(NF, NV), intent(out) :: &
         n_down, mu_up, mu_down
    real(kind=wp), dimension(NFF), intent(out) :: n_sqd

    integer :: i, f, f1, f2

!!$    Use the metric to compute n_down
    n_down(:, 1)    = -n_up(:, 1)
    n_down(:, 2:NV) =  n_up(:, 2:NV)

!!$    Compute the scalars (-n^2)
    call compute_scalars_pt(n_up, n_down, n_sqd)

!!$    Compute the conjugate momenta.
    call compute_mu_down_pt(n_down, n_sqd, mu_down)

!!$    Compute mu_up for completeness
    mu_up(:, 1)    = -mu_down(:, 1)
    mu_up(:, 2:NV) =  mu_down(:, 2:NV)

  end subroutine all_from_n_up_pt

  subroutine all_from_n_up(n_up, n_down, mu_up, mu_down, n_sqd)

    real(kind=wp), dimension(NF, NV, 1-gz:np+gz), intent(in) :: n_up
    real(kind=wp), dimension(NF, NV, 1-gz:np+gz), intent(out) :: &
         n_down, mu_up, mu_down
    real(kind=wp), dimension(NFF, 1-gz:np+gz), intent(out) :: n_sqd

    integer :: i

!!$    Use the local routine
    do i = 1-gz, np+gz
      call all_from_n_up_pt(n_up(:, :, i), n_down(:, :, i), &
           mu_up(:, :, i), mu_down(:, :, i), n_sqd(:, i))
    end do

  end subroutine all_from_n_up

  subroutine all_from_evolved_pt(q, n_up, n_down, mu_up, mu_down, n_sqd)

    use local_cons_mod

    real(kind=wp), dimension(NEVOLVE), intent(in) :: q
    real(kind=wp), dimension(NF, NV), intent(out) :: n_up, n_down, &
         mu_up, mu_down
    real(kind=wp), dimension(NFF), intent(inout) :: n_sqd

    integer :: f, f1, f2, k

    real(kind=wp) :: tol = 1.d-10
    integer :: info
    real(kind=wp), dimension(LWA) :: wa
    
    integer, dimension(NF) :: ipiv
    real(kind=wp), dimension(NFF) :: dLdn
    real(kind=wp), dimension(NF, NF) :: A
    real(kind=wp), dimension(NF) :: b
    real(kind=wp), dimension(NFF) :: f_err
    real(kind=wp), dimension(NFF, NFF) :: f_jac

!!$    Set the local_cons variable for the external routines
    do f = 1, NF
      local_cons(f, 1) = q(f)
      local_cons(f, 2) = q(NF + f)
    end do
    
!!$    The first NF evolved variables are timelike components of n_up
    n_up(:, 1) = q(1:NF)

!!$    Initial guesses for the scalars; that they are the squares of
!!$    the timelike components of n
    do f = 1, NF
!!$      Do not let the guess jump too much
      if (n_up(f, 1)**2 < 2.0_wp * n_sqd(f)) then
        n_sqd(f) = n_up(f, 1)**2 - mu_down(f, 2)**2
      end if
    end do
    f = NF
    do f1 = 1, NF
      do f2 = f1 + 1, NF
        f = f + 1
        if (n_up(f1, 1) * n_up(f2, 1) < 2.0_wp * n_sqd(f)) then
          n_sqd(f) = n_up(f1, 1) * n_up(f2, 1) - mu_down(f1, 2) * mu_down(f2, 2)
        end if
      end do
    end do

!!$  We have to try and find the correct values of the scalars. 
!!$  This requires a multi-dimensional root find.
    call hybrd1(Cons2PrimFn, NFF, n_sqd, f_err, tol, info, wa, LWA)
!!$    The minpack procedure complains if the tolerance requested is
!!$    too tight (info=3), and says if it thinks the solution is fine
!!$    and within the tolerance (info=1); all others should be errors
    if ((info .ne. 1).and.(info .ne. 3)) then
      write(*,*) 'hybrd1 call failed'
      write(*,*) 'Error in ''all_from_evolved_pt'''
      write(*,*) 'Error is ', info
      write(*,*) n_sqd 
      write(*,*) n_up 
      write(*,*) mu_down
      write(*,*) local_cons
      stop
    end if

!!$  We now know the values of the scalars - 
!!$  find the terms depending on the master function
    call MasterFunction_Derivs(n_sqd, dLdn)

!!$  Now define the matrix equation A x = b and invert; 
!!$  solution is in b
    do f = 1, NF
      A(f, f) = -2.0_wp * dLdn(f)
      b(f) = local_cons(f, 2)
    end do
    f = NF
    do f1 = 1, NF
      do f2 = f1 + 1, NF
        f = f + 1
        A(f1, f2) = - dLdn(f)
        A(f2, f1) = - dLdn(f)
      end do
    end do
    call dgesv(NF, 1, A, NF, ipiv, b, NF, info)
    if (info .ne. 0) then
      write(*,*) 'dgesv call failed'
      write(*,*) 'Error in ''all_from_n_up'''
      write(*,*) 'Error is ', info
      stop
    end if
    n_up(:, 2) = b

!!$    We now have the complete n_up vector; find all the rest

    call all_from_n_up_pt(n_up, n_down, mu_up, mu_down, n_sqd)

  end subroutine all_from_evolved_pt

  subroutine Cons2PrimFn(NF_NF, n_sqd, f_err, iflag)

    use local_cons_mod

!!$  This argument is required by the minpack routines. However, as we
!!$  already have this through the module, we never use it
    integer, intent(in) :: NF_NF

    integer :: iflag
    real(kind=wp), dimension(NFF) :: n_sqd, new_n_sqd, f_err

    real(kind=wp), dimension(NFF) :: dLdn
    real(kind=wp), dimension(NF, NF) :: A
    real(kind=wp), dimension(NF) :: b

    integer :: info
    integer, dimension(NF) :: ipiv

    integer :: f, f1, f2

!!$  Given the scalars (and through the external module the conserved
!!$  variables; in particular, the spatial components of the conserved
!!$  variables, which are the spatial components of the conjugate
!!$  momenta) we can compute the master function, and in particular its
!!$  derivatives (dLdn). From that we can write the conjugate momenta
!!$  as a matrix (A) times the spatial components of n. Solving this linear
!!$  system gives n. From n we compute the scalars, and hence an error in
!!$  the initial guess for the scalars

!!$    Check if any guess for n_sqd is unphysical (negative); if so, exit
    if (any(n_sqd < 0.0_wp)) then

      f_err = 1.d10

    else

      call MasterFunction_Derivs(n_sqd, dLdn)

!!$  Now define the matrix equation A x = b and invert; 
!!$  solution is in b
      do f = 1, NF
        A(f, f) = -2.0_wp * dLdn(f)
        b(f) = local_cons(f, 2)
      end do
      f = NF
      do f1 = 1, NF
        do f2 = f1 + 1, NF
          f = f + 1
          A(f1, f2) = - dLdn(f)
          A(f2, f1) = - dLdn(f)
        end do
      end do
      call dgesv(NF, 1, A, NF, ipiv, b, NF, info)
      if (info .ne. 0) then
        write(*,*) 'dgesv call failed'
        write(*,*) 'Error in ''Con2PrimFn'''
        write(*,*) 'Error is ', info
        write(*,*) 'A', A
        write(*,*) 'b', b
        write(*,*) 'local_cons', local_cons
        write(*,*) 'n_sqd', n_sqd
        stop
      end if
    
!!$  The newly computed n is placed in b; find the new scalars
!!$  Note that the timelike component of n is given by the timelike
!!$  component of the conserved variables
    
!!$    Should write this out generally

      new_n_sqd(1) = local_cons(1, 1) * local_cons(1, 1) - b(1) * b(1)
      new_n_sqd(2) = local_cons(2, 1) * local_cons(2, 1) - b(2) * b(2)
      new_n_sqd(3) = local_cons(1, 1) * local_cons(2, 1) - b(1) * b(2)
      
      f_err = n_sqd - new_n_sqd

    end if
    
  end subroutine Cons2PrimFn

  subroutine Cons2PrimFnJac(NF_NF, n_sqd, f_err, f_jac, NF_NF_NF, iflag)

    use local_cons_mod

!!$  This argument is required by the minpack routines. However, as we
!!$  already have this through the module, we never use it
    integer, intent(in) :: NF_NF, NF_NF_NF

    integer :: iflag
    real(kind=wp), dimension(NFF) :: n_sqd, new_n_sqd, f_err
    real(kind=wp), dimension(NFF, NFF) :: f_jac

    real(kind=wp), dimension(NFF) :: dLdn
    real(kind=wp), dimension(NF, NF) :: A
    real(kind=wp), dimension(NF) :: b

    integer :: info
    integer, dimension(NF) :: ipiv

    integer :: f, f1, f2

!!$  Given the scalars (and through the external module the conserved
!!$  variables; in particular, the spatial components of the conserved
!!$  variables, which are the spatial components of the conjugate
!!$  momenta) we can compute the master function, and in particular its
!!$  derivatives (dLdn). From that we can write the conjugate momenta
!!$  as a matrix (A) times the spatial components of n. Solving this linear
!!$  system gives n. From n we compute the scalars, and hence an error in
!!$  the initial guess for the scalars


    if (iflag .eq. 1) then ! Compute the function, not the Jacobian

!!$    Check if any guess for n_sqd is unphysical (negative); if so, exit
      if (any(n_sqd < 0.0_wp)) then
        
        f_err = 1.d10
        
      else
        
        call MasterFunction_Derivs(n_sqd, dLdn)
        
!!$  Now define the matrix equation A x = b and invert; 
!!$  solution is in b
        do f = 1, NF
          A(f, f) = -2.0_wp * dLdn(f)
          b(f) = local_cons(f, 2)
        end do
        f = NF
        do f1 = 1, NF
          do f2 = f1 + 1, NF
            f = f + 1
            A(f1, f2) = - dLdn(f)
            A(f2, f1) = - dLdn(f)
          end do
        end do
        call dgesv(NF, 1, A, NF, ipiv, b, NF, info)
        if (info .ne. 0) then
          write(*,*) 'dgesv call failed'
          write(*,*) 'Error in ''Con2PrimFn'''
          write(*,*) 'Error is ', info
          write(*,*) 'A', A
          write(*,*) 'b', b
          write(*,*) 'local_cons', local_cons
          write(*,*) 'n_sqd', n_sqd
          stop
        end if
        
!!$  The newly computed n is placed in b; find the new scalars
!!$  Note that the timelike component of n is given by the timelike
!!$  component of the conserved variables
    
!!$    Should write this out generally

        new_n_sqd(1) = local_cons(1, 1) * local_cons(1, 1) - b(1) * b(1)
        new_n_sqd(2) = local_cons(2, 1) * local_cons(2, 1) - b(2) * b(2)
        new_n_sqd(3) = local_cons(1, 1) * local_cons(2, 1) - b(1) * b(2)
      
        f_err = n_sqd - new_n_sqd

      end if
    
    else if (iflag .eq. 2) then ! Compute Jacobian

      f_jac = 0_wp
      write(*, *) 'Do not know how to compute Jacobian yet'
      stop

    else
      
      write(*,*) 'Do not understand iflag value in Cons2PrimFnJac', iflag

    end if

  end subroutine Cons2PrimFnJac

end module conversion_mod
