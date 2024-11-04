subroutine CountLines(filename, num_lines)
    
    implicit none
    
    character(len=*), intent(in) :: filename
    integer, intent(out) :: num_lines
    integer :: unit_num, io_status
    character(len=256) :: line

    num_lines = 0
    open(newunit=unit_num, file=filename, status="old", action="read")

    do
        read(unit_num, '(A)', iostat=io_status) line
        if (io_status /= 0) exit
        num_lines = num_lines + 1
    end do

    close(unit_num)
end subroutine CountLines
