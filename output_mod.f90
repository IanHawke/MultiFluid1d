module output_mod

  use real_type_mod
  use parameters_mod
  use system_mod

  implicit none

contains

  function make_out_dir(out_dir_base) result(out_dir)

    character(len=*), intent(in) :: out_dir_base

    character(len=200) :: out_dir
    character(len=8) :: date_string
    character(len=10) :: time_string

    integer :: istat

    call date_and_time(date_string, time_string)
    if (out_dir_date) then
      out_dir = trim(trim(out_dir_base)//'_'//trim(date_string)//'_'//trim(time_string(1:6)))
    else
      out_dir = trim(trim(out_dir_base))
    end if
    call execute_command(trim('mkdir -p '//out_dir), istat)
    if (istat .ne. 0) then
      write(*,*) 'System mkdir call failed, error ', istat
      stop
    end if

  end function make_out_dir

  subroutine create_output_file(out_dir, name, file_number)

    character(len=*), intent(in) :: out_dir, name
    integer, intent(in) :: file_number

    character(len=200) :: filestring
    integer :: istat
    
    filestring = trim(trim(out_dir)//'/'//trim(name)//'.dat')
    open(file_number, FILE=filestring, STATUS='replace', IOSTAT=istat)

    if (istat .ne. 0) then
      write(*,*) 'Error opening file ', istat
      stop
    end if
    
  end subroutine create_output_file

  subroutine create_output_timefreq_vtk_file(out_dir, name, file_number, &
       np, freq, out_total, out_dt)

    character(len=*), intent(in) :: out_dir, name
    integer, intent(in) :: file_number, np, out_total
    real(kind=wp), dimension(np/2 + 1), intent(in) :: freq
    real(kind=wp), intent(in) :: out_dt

    character(len=200) :: filestring
    character(1), parameter :: end_rec = char(10)
    character(len=255) :: out_string

    integer :: istat, i

    logical :: vtk_output_binary = .false.
    
    filestring = trim(trim(out_dir)//'/'//trim(name)//'.vtk')
    if (vtk_output_binary) then
      open(file_number, FILE=filestring, &
           FORM='UNFORMATTED', ACCESS = 'SEQUENTIAL', ACTION = 'WRITE', &
           CONVERT='BIG_ENDIAN', &
           STATUS='replace', IOSTAT=istat)
    else
      open(file_number, FILE=filestring, STATUS='replace', IOSTAT=istat)
    end if
    if (istat .ne. 0) then
      write(*,*) 'Error opening file ', istat
      stop
    end if

    if (vtk_output_binary) then
      out_string = "# vtk DataFile Version 3.0"//end_rec
      write(file_number) trim(out_string)
      out_string = "Data produced by mf_evolve"//end_rec
      write(file_number) trim(out_string)
      out_string = "BINARY"//end_rec
      write(file_number) trim(out_string)
      out_string = "DATASET RECTILINEAR_GRID"//end_rec//end_rec
      write(file_number) trim(out_string)
      write(out_string, '(a10,3i6)') "DIMENSIONS", &
           np/2 + 1, out_total, 1
      write(file_number) trim(trim(out_string)//end_rec)
      write(out_string, '(a13,i6,a6)')    "X_COORDINATES", np/2+1, &
           " FLOAT"
      write(file_number) trim(trim(out_string)//end_rec)
      write(file_number) (freq(i), i=1,np/2+1)
      write(file_number) end_rec//end_rec
      write(out_string, '(a13,i6,a6)')    "Y_COORDINATES", out_total, &
           " FLOAT"
      write(file_number) trim(trim(out_string)//end_rec)
      write(file_number) (out_dt * dble(i-1), i=1,out_total)
      write(file_number) end_rec//end_rec
      write(out_string, '(a13,i6,a6)')    "Z_COORDINATES", 1, &
           " FLOAT"
      write(file_number) trim(trim(out_string)//end_rec)
      write(file_number) 0_wp
      write(file_number) end_rec//end_rec
      write(out_string, '(a10,i12)') "POINT_DATA", (np/2+1)*(out_total)
      write(file_number) trim(trim(out_string)//end_rec)
      out_string = "SCALARS Power FLOAT"//end_rec
      write(file_number) trim(out_string)
      out_string = "LOOKUP_TABLE default"//end_rec
      write(file_number) trim(out_string)
      write(file_number) end_rec//end_rec
    else
      write(file_number, '(a26)') "# vtk DataFile Version 3.0"
      write(file_number, '(a26)') "Data produced by mf_evolve"
      write(file_number, '(a5)')  "ASCII"
      write(file_number, '(a24)') "DATASET RECTILINEAR_GRID"
      write(file_number, '(a10,3i6)') "DIMENSIONS", &
           np/2 + 1, out_total, 1
      write(file_number, '(a13,i6,a6)')    "X_COORDINATES", np/2+1, " FLOAT"
      do i = 1, np/2+1
        write(file_number, '(es21.12E3, a1)', advance = 'NO') freq(i), ' '
      end do
      write(file_number, *)
      write(file_number, '(a13,i6,a6)')    "Y_COORDINATES", out_total, " FLOAT"
      do i = 1, out_total
        write(file_number, '(es21.12E3, a1)', advance = 'NO') &
             out_dt * dble(i-1), ' '
      end do
      write(file_number, *)
      write(file_number, '(a13,i6,a6)')    "Z_COORDINATES", 1     , " FLOAT"
      write(file_number, '(i1)') 0
      write(file_number, '(a10,i12)') "POINT_DATA", (np/2+1)*(out_total)
      write(file_number, '(a19)') "SCALARS Power FLOAT"
      write(file_number, '(a20)') "LOOKUP_TABLE default"
    end if

    call flush(file_number)

  end subroutine create_output_timefreq_vtk_file
  
  subroutine output(file_number, x, q, nv, np, gz, t)

    integer, intent(in) :: file_number, nv, np, gz
    real(kind=wp), intent(in) :: t
    real(kind=wp), dimension(1-gz:np+gz), intent(in) :: x
    real(kind=wp), dimension(nv, 1-gz:np+gz), intent(in) :: q

    integer :: i, k

    write(file_number, '(a8, es21.12E3)') '#Time = ', t
    do i = 1, np
      write(file_number, '(es21.12E3, a1)', advance='NO') x(i), ' '
      do k = 1, nv
        write(file_number, '(es21.12E3,a1)',advance='NO') q(k, i), ' '
      end do
      write(file_number, *)
    end do
    write(file_number, *)
    write(file_number, *)

    call flush(file_number)

  end subroutine output
  
  subroutine output2(file_number, x, q, nf, nv, np, gz, t)

    integer, intent(in) :: file_number, nf, nv, np, gz
    real(kind=wp), intent(in) :: t
    real(kind=wp), dimension(1-gz:np+gz), intent(in) :: x
    real(kind=wp), dimension(nf, nv, 1-gz:np+gz), intent(in) :: q

    integer :: f, i, k

    write(file_number, '(a8, es21.12E3)') '#Time = ', t
    do i = 1, np
      write(file_number, '(es21.12E3, a1)', advance='NO') x(i), ' '
      do f = 1, nf
        do k = 1, nv
          write(file_number, '(es21.12E3,a1)',advance='NO') q(f, k, i), ' '
        end do
      end do
      write(file_number, *)
    end do
    write(file_number, *)
    write(file_number, *)

    call flush(file_number)

  end subroutine output2
  
  subroutine output3(file_number, freq, q, np, NV, NF, t)

    integer, intent(in) :: file_number, np, NV, NF
    real(kind=8), intent(in) :: t
    real(kind=wp), dimension(np/2 + 1), intent(in) :: freq
    complex(kind=wp), dimension(NF*NV, np/2 + 1), intent(in) :: q

    integer :: f, i, k, ffvv

    write(file_number, '(a8, es21.12E3)') '#Time = ', t
    do i = 1, np/2 + 1
!!$      write(file_number, '(3es21.12E3)') t, dble(i-1), abs(q(i))
      write(file_number, '(2es21.12E3)', advance = 'NO') t, dble(i-1)
      do ffvv = 1, NF*NV
        write(file_number, '(es21.12E3, a1, es21.12E3, a1)', advance = 'NO') &
             real(q(ffvv, i)), ' ', aimag(q(ffvv, i)), ' '
      end do
      write(file_number, *)
    end do
    write(file_number, *)

    call flush(file_number)

  end subroutine output3

  subroutine output_timefreq_vtk(file_number, q, np)

    integer, intent(in) :: file_number, np
    complex(kind=wp), dimension(np/2 + 1), intent(in) :: q

    integer :: i

    logical :: vtk_output_binary = .true.

    if (vtk_output_binary) then
      write(file_number) (abs(q(i)),i=1,np/2+1)
    else
      do i = 1, np/2 + 1
        write(file_number, '(es21.12E3)') abs(q(i))
      end do
    end if

  end subroutine output_timefreq_vtk
  
end module output_mod
