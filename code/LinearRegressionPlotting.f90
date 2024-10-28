program LinearRegressionPlotting
    implicit none
    integer :: n
    integer :: i, unit_num
    real, allocatable :: Xi(:), Yi(:)
    real :: a, b, sdv, R, std_a, std_b
    real :: xsec
    real :: ax_x_max, ax_x_min, ax_y_max, ax_y_min, ax_scale_ratio
    real REGRESSION_X(2), REGRESSION_Y(2), reg_scale_ratio

    character(len = 100) :: line
    character(len=100) :: equation
    integer iy, im, id, ih, mm

    real :: x_legend, y_legend
    real :: legend_char_height, symbol_size, line_length, line_spacing, text_offset

    ! Set the number of points
    n = 30


    ! Allocate Xi and Yi
    allocate(Xi(n), Yi(n))


    open(newunit=unit_num, file = "../data/ppfile.txt", status = "old", action = "read")
    
    
    ! Get the number of lines-----------------
    ! Skipping header
    read(unit_num, '(A)') ! Do nothing
    ! Get the number of lines
    n = 0
    do while (.true.)
        read(unit_num, '(A)', iostat=i) line
        if (i /= 0) exit
        n = n + 1
    end do
    ! Get the number of lines-----------------


    ! Read the datas and skipping header------
    rewind(unit_num)
    read(unit_num,'(1x, i4, 4i2, f6.2)') iy, im, id, ih, mm, xsec
    xsec = mm * 60.0 + xsec


    ! Read data and fill into Xi and Yi
    do i = 1, n
        ! read(unit_num, '(6X, F5.1, F4.0)') Xi(i), Yi(i)
        read(unit_num,'(5x, f6.1, 9x, i3, f6.2)') Xi(i), mm, Yi(i)
        Yi(i) = mm * 60.0 + Yi(i) - xsec
    end do
    ! Read the datas and skipping header------

    close(unit_num)
    
    call LinearRegressionFitting(Xi, Yi, n, a, b, sdv, R, std_a, std_b)

    print *, 'Constant a: ', a
    print *, 'Constant b: ', b
    print *, 'Standard deviation of fit (sdv): ', sdv
    print *, 'Linear correlation coefficient (R): ', R
    print *, 'Standard deviation of a: ', std_a
    print *, 'Standard deviation of b: ', std_b

    ! Get equation
    write(equation, '(A, F6.2, A, F6.2)') 'y = ', a, 'x + ', b


    ! Scatter plot 
    !! Set ax size
    ax_scale_ratio = 0.25
    reg_scale_ratio= 0.15
    ax_x_min = minval(Xi) - abs(minval(Xi)) * ax_scale_ratio
    ax_x_max = maxval(Xi) + abs(maxval(Xi)) * ax_scale_ratio
    ax_y_min = minval(Yi) - abs(minval(Yi)) * ax_scale_ratio
    ax_y_max = maxval(Yi) + abs(maxval(Yi)) * ax_scale_ratio

    REGRESSION_X(1) = minval(Xi) - abs(minval(Xi)) * reg_scale_ratio
    REGRESSION_X(2) = maxval(Xi) + abs(maxval(Xi)) * reg_scale_ratio
    REGRESSION_Y(1) = b + a * REGRESSION_X(1)
    REGRESSION_Y(2) = b + a * REGRESSION_X(2)

    x_legend = ax_x_min + 0.1 * (ax_x_max - ax_x_min)
    y_legend = ax_y_max - 0.1 * (ax_y_max - ax_y_min)

    legend_char_height = 1.0
    symbol_size = 1.0
    line_length = (ax_x_max - ax_x_min) * 0.05
    line_spacing = (ax_y_max - ax_y_min) * 0.05
    text_offset = (ax_x_max - ax_x_min) * 0.02

    !! Plotting
    call pgopen('EpicentralDistance-Time.ps/VCPS')
    call pgsci(1)
    call pgenv(ax_x_min, ax_x_max, ax_y_min, ax_y_max, 0, 1)
    call pgscf(1)
    call pglab('Epicentral Distance (km)', 'Time (s)', 'Epicentral Distance - Time')


    call pgsci(2)
    call pgpt(n, Xi, Yi, 5)
    call pgsci(4)
    call pgslw(2)
    call pgline(2, REGRESSION_X, REGRESSION_Y)

    call pgsls(3)
    call pgline(2, REGRESSION_X, REGRESSION_Y + 2 * std_b)
    call pgline(2, REGRESSION_X, REGRESSION_Y - 2 * std_b)

    call pgsci(1)
    call pgtext(50.0, 8.0, equation)



    call pgsch(legend_char_height)
    call pgsci(2)
    call pgpt1(x_legend, y_legend + 0.1, 5)
    call pgsci(1)
    call pgptxt(x_legend + text_offset, y_legend, 0.0, 0.0, 'Data Points')
    call pgsci(4)
    
    call pgsls(1)
    call pgmove(x_legend, y_legend - line_spacing + 0.1)
    call pgdraw(x_legend + line_length - 5, y_legend - line_spacing + 0.1)
    call pgsci(1)
    call pgptxt(x_legend + text_offset, y_legend - line_spacing, 0.0, 0.0, 'Linear Regression')

    call pgsls(3)
    call pgsci(4)
    call pgmove(x_legend, y_legend - 2 * line_spacing + 0.1)
    call pgdraw(x_legend + line_length, y_legend - 2 * line_spacing + 0.1)
    call pgsci(1)
    call pgptxt(x_legend + text_offset + 5, y_legend - 2 * line_spacing, 0.0, 0.0, '95% C.I.')

    
    call pgclos()

    ! y_err = x.std() * np.sqrt(1/len(x) + (x - x.mean())**2 / np.sum((x - x.mean())**2))

    deallocate(Xi, Yi)
end program LinearRegressionPlotting