module system_mod

  implicit none

contains

  subroutine execute_command(command, status)

    character(len=*), intent(in)  :: command
    integer,          intent(out) :: status

#ifdef __INTEL_COMPILER

    integer, external :: system

#endif

    status = 0

#ifdef __INTEL_COMPILER

    status = system(trim(command))
#else
    call execute_command_line(trim(command), EXITSTAT=status)
#endif

  end subroutine execute_command

end module system_mod
