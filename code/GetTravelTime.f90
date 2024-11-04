program GetTravelTime
    implicit none
    integer :: n
    real, allocatable :: Xi(:), Yi(:)
    character(len = 4), allocatable :: STA_TG_CODES(:)
    real :: a, b, sdv, R, std_a, std_b

    character(len = 100) :: line
    integer :: i, unit_num

    real :: xsec
    integer iy, im, id, ih, mm

    call CountLines("../data/ppfile.txt", n)

    n = n - 1
    allocate(Xi(n), Yi(n), STA_TG_CODES(n))


    open(newunit=unit_num, file = "../data/ppfile.txt", status = "old", action = "read")

    ! Read the datas and skipping header------
    rewind(unit_num)
    read(unit_num,'(1x, i4, 4i2, f6.2)') iy, im, id, ih, mm, xsec
    xsec = mm * 60.0 + xsec

    ! Read data and fill into Xi and Yi
    do i = 1, n
        read(unit_num,'(A4, f6.1, 9x, i3, f6.2)') STA_TG_CODES(i), Xi(i), mm, Yi(i)
        Yi(i) = mm * 60.0 + Yi(i) - xsec
    end do
    close(unit_num)










    deallocate(Xi, Yi, STA_TG_CODES)

end program GetTravelTime