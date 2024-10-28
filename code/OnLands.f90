program OnLands
    implicit none

    character(len = 100) :: line
    
    integer :: n, m, i
    integer :: io_status, unit_num

    real, allocatable :: CONTOUR_X(:), CONTOUR_Y(:)
    real :: x, y, z, w
    real, allocatable :: EQ_X(:), EQ_Y(:), EQ_Z(:), EQ_M(:)
    real :: ax_x_max, ax_x_min, ax_y_max, ax_y_min, margen

    ! Used for locpt subroutine
    integer :: l_locpt, m_locpt
    real :: x0, y0

    ! Get contour data -------------------------------------------------------------------
    open(newunit = unit_num, file = "../data/Taiwan.txt", status = "old", action = "read")
    n = 0
    do
        read(unit_num, '(A)', iostat=io_status) line
        if (io_status /= 0) exit
        n = n + 1
    end do
    close(unit_num)
    allocate(CONTOUR_X(n), CONTOUR_Y(n))
    open(newunit=unit_num, file = "../data/Taiwan.txt", status = "old", action = "read")
    do i = 1, n
        read(unit_num, *, iostat=io_status) CONTOUR_X(i), CONTOUR_Y(i)
        if (io_status /= 0) exit
    end do
    close(unit_num)


    ! Get earthquake data ----------------------------------------------------------------
    open(newunit = unit_num, file = "../data/1999.lis", status = "old", action = "read")
    m = 0
    do
        read(unit_num, '(A)', iostat=io_status) line
        if (io_status /= 0) exit
        m = m + 1
    end do
    close(unit_num)
    allocate(EQ_X(m), EQ_Y(m), EQ_Z(m), EQ_M(m))

    open(newunit=unit_num, file = "../data/1999.lis", status = "old", action = "read")
    do i = 1, m
        read(unit_num, '(18X, F2.0, F5.2, F3.0, F5.2, F6.2, F5.2)', iostat=io_status) x, y, z, w, EQ_Z(i), EQ_M(i)
        EQ_X(i) = (z + w / 60)
        EQ_Y(i) = (x + y / 60)


        if (io_status /= 0) exit
    end do
    close(unit_num)

    do i = 1, 3
        print *, EQ_X(i), EQ_Y(i)
    end do

    do i = 1, 3
        x0 = EQ_X(i)
        y0 = EQ_Y(i)
        call locpt(x0, y0, CONTOUR_X, CONTOUR_Y, n, l_locpt, m_locpt)
        if (l_locpt == 1) then
            print *, 'Earthquake ', i, ' is inside the contour.'
        else if (l_locpt == 0) then
            print *, 'Earthquake ', i, ' is on the contour boundary.'
        else
            print *, 'Earthquake ', i, ' is outside the contour.'
        end if
    end do








    ! Set x, y lims ----------------------------------------------------------------------
    margen = 0.25  ! 0.25 degree larger
    ax_x_min = minval(EQ_X) - margen
    ax_x_max = maxval(EQ_X) + margen
    ax_y_min = minval(EQ_Y) - margen * 2
    ax_y_max = maxval(EQ_Y) + margen * 2

    call pgopen('1999_event_distribution.ps/VCPS')
    call pgsci(1)
    call pgenv(ax_x_min, ax_x_max, ax_y_min, ax_y_max, 0, 1)
    call pgscf(1)
    call pglab('Longitude (E)', 'Latitude (N)', '1999 Event Distribution')


    do i = 1, m

        x0 = EQ_X(i)
        y0 = EQ_Y(i)
        call locpt(x0, y0, CONTOUR_X, CONTOUR_Y, n, l_locpt, m_locpt)
        if (l_locpt == 1) then
            call pgscr(42, 0.121, 0.466, 0.706)
            call pgsci(42)
            call pgsch(1.0)
            call pgpt1(EQ_X(i), EQ_Y(i), 3)

        else if (l_locpt == 0) then
            call pgscr(41, 0.0862, 0.659, 0.173)
            call pgsci(41)
            call pgsch(1.0)
            call pgpt1(EQ_X(i), EQ_Y(i), 3)
            print *, 'Earthquake ', i, ' is on the contour boundary.'
        else
            call pgscr(43, 1.00, 0.498, 0.055)
            call pgsci(43)
            call pgsch(1.0)
            call pgpt1(EQ_X(i), EQ_Y(i), 3)
        end if


    end do
    
    call pgsci(1)
    call pgline(n, CONTOUR_X, CONTOUR_Y)
    
    call pgclos()


    deallocate(CONTOUR_X, CONTOUR_Y)
    deallocate(EQ_X, EQ_Y, EQ_Z, EQ_M)

end program OnLands