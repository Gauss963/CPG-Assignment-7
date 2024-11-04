program GetTravelTime
    implicit none
    integer :: n
    real, allocatable :: Xi(:), Yi(:)
    real :: a, b, sdv, R, std_a, std_b

    character(len = 100) :: line
    integer :: i, unit_num

    real :: xsec
    integer iy, im, id, ih, mm

    call CountLines("../data/ppfile.txt", n)
    allocate(Xi(n), Yi(n))


    open(newunit=unit_num, file = "../data/ppfile.txt", status = "old", action = "read")
    n = n - 1


    ! Read the datas and skipping header------
    rewind(unit_num)
    read(unit_num,'(1x, i4, 4i2, f6.2)') iy, im, id, ih, mm, xsec
    xsec = mm * 60.0 + xsec


    ! Read data and fill into Xi and Yi
    do i = 1, n
        read(unit_num,'(5x, f6.1, 9x, i3, f6.2)') Xi(i), mm, Yi(i)
        Yi(i) = mm * 60.0 + Yi(i) - xsec
    end do
    ! Read the datas and skipping header------

    close(unit_num)

    print *, Yi


    deallocate(Xi, Yi)

end program GetTravelTime