program GPSLocating
    
    use, intrinsic :: iso_fortran_env
    
    implicit none

    integer, parameter :: dp = real64
    integer :: N, i, j, k, l, m, INFO, o, tg_num
    
    
    real(dp), allocatable :: v_QC1(:), v_QC2(:), v_QC3(:), v_d(:), v_d_FIXED(:)
    real(dp) :: v_m(3), v_location(3)
    real(dp), allocatable :: m_Q(:,:), m_QT(:,:), m_Q_FIXED(:,:)
    real(dp) :: m_QTQ(3, 3), m_QTQ_INVERSE(3, 3)
    real(dp) :: velocity


    ! real(dp), allocatable :: v_QC1(:), v_QC2(:), v_QC3(:), v_d(:), v_d_FIXED(:)
    ! real(dp) :: v_m(4), v_location(4)
    ! real(dp), allocatable :: m_Q(:,:), m_QT(:,:), m_Q_FIXED(:,:)
    ! real(dp) :: m_QTQ(4, 4), m_QTQ_INVERSE(4, 4)


    character(len = 100) :: line, elev_str, lat_str, lon_str
    character(len = 100) :: data_filename
    real, allocatable :: STATION_X(:), STATION_Y(:), STATION_Z(:)
    real :: lat_deg, lat_min, lon_deg, lon_min
    integer :: io_status, unit_num


    real, allocatable :: TRAVEL_DIST(:), TRAVEL_TIMES(:)
    character(len = 4), allocatable :: STA_TG_CODES(:)
    real :: a, b, sdv, R, std_a, std_b

    real :: xsec
    integer iy, im, id, ih, mm
    
    
    velocity = 7 ! km/s


    ! Count the number of lines.
    call CountLines("../data/nsta.dat", o)
    allocate(STATION_X(o), STATION_Y(o), STATION_Z(o))

    ! Read data (lat, lon, elev)
    open(newunit=unit_num, file = "../data/nsta.dat", status = "old", action = "read")
    do i = 1, m
        read(unit_num, '(A)', iostat=io_status) line
        if (io_status /= 0) exit

        lat_str = line(5:13)
        read(lat_str(1:2), *) lat_deg
        read(lat_str(3:7), *) lat_min

        STATION_Y(i) = lat_deg + lat_min / 60.0

        lon_str = line(13:22)
        read(lon_str(1:3), *) lon_deg
        read(lon_str(4:8), *) lon_min

        STATION_X(i) = lon_deg + lon_min / 60.0

        elev_str = line(22:28)
        if (trim(elev_str) == '') then
            STATION_Z(i) = 0.0
        else
            read(elev_str, *) STATION_Z(i)
        end if
    end do
    close(unit_num)




    call CountLines("../data/ppfile.txt", tg_num)
    
    tg_num = tg_num - 1

    allocate(TRAVEL_DIST(tg_num), TRAVEL_TIMES(tg_num), STA_TG_CODES(tg_num))

    open(newunit=unit_num, file = "../data/ppfile.txt", status = "old", action = "read")
    rewind(unit_num)

    read(unit_num,'(1x, i4, 4i2, f6.2)') iy, im, id, ih, mm, xsec
    xsec = mm * 60.0 + xsec
    do i = 1, tg_num
        read(unit_num,'(A4, f6.1, 9x, i3, f6.2)') STA_TG_CODES(i), TRAVEL_DIST(i), mm, TRAVEL_TIMES(i)
        TRAVEL_TIMES(i) = mm * 60.0 + TRAVEL_TIMES(i) - xsec
    end do
    close(unit_num)
    print *, STA_TG_CODES










    v_location = [20, 20, 20]

    ! Prompt user for the number of stations
    print *, 'Enter the number of stations (N >= 4):'
    read *, N
    if (N < 4) then
        print *, 'At least 4 stations are required.'
        stop
    end if
    allocate(v_QC1(N), v_QC2(N), v_QC3(N), v_d(N))


    ! Input station coordinates and dstance
    print *, 'Enter the coordinates (v_QC1, v_QC2, v_QC3) and dstance v_d for each station:'
    do i = 1, N
        print *, 'Station ', i, ':'
        read *, v_QC1(i), v_QC2(i), v_QC3(i), v_d(i)
    end do


    v_d_FIXED = v_d
    ! print *, 'd_FIXED is', v_d_FIXED


    allocate(m_Q(N, 3), m_QT(3, N), m_Q_FIXED(N, 3))
    do j = 1, N
        m_Q(j, 1) = v_QC1(j) ! (column, row) => (j, 1) -> v_QC1(j)
        m_Q(j, 2) = v_QC2(j) ! Since FORTRAN is column-major
        m_Q(j, 3) = v_QC3(j) ! See `https://en.wikipedia.org/wiki/Row-_and_column-major_order`

        m_Q_FIXED(j, 1) = v_QC1(j)
        m_Q_FIXED(j, 2) = v_QC2(j)
        m_Q_FIXED(j, 3) = v_QC3(j)
    end do
    ! Q FIXED initialized here. Print the matrix to check
    print *, 'Matrix Q FIXED:'
    do k = 1, N
        print *, m_Q_FIXED(k, 1), m_Q_FIXED(k, 2), m_Q_FIXED(k, 3)
    end do



    do i = 1, 1000
        
        ! Update matrix Q
        do j = 1, N
            do k = 1, 3
                m_Q(j, k) = v_location(k) - m_Q_FIXED(j, k)
            end do
        end do

        ! Update vector d
        do m = 1, N
            v_d(m) = v_d_FIXED(m)**2 - ( m_Q_FIXED(m, 1) - v_location(1) )**2 &
                                     - ( m_Q_FIXED(m, 2) - v_location(2) )**2 &
                                     - ( m_Q_FIXED(m, 3) - v_location(3) )**2
            
            v_d(m) = v_d(m) * 0.5
        end do
        ! print *, 'vector d = ', v_d
        ! print *, '---------------------------------------------------------------'

        ! print *, 'Matrix Q:'
        ! do k = 1, N
        !     print *, m_Q(k, 1), m_Q(k, 2), m_Q(k, 3)
        ! end do



        ! Solve linear system
        m_QT = transpose(m_Q)
        m_QTQ = matmul(m_QT, m_Q)
        
        call matrix_inverse(3, m_QTQ, m_QTQ_INVERSE, INFO)
        v_m = matmul(matmul(m_QTQ_INVERSE, m_QT), v_d)
        
        
        ! call MATRIXINV(m_QTQ, 3)
        ! v_m = matmul(matmul(m_QTQ, m_QT), v_d)
        
        v_location = v_location + v_m

        ! print *, 'vector m = ', v_m




        ! Check if \delta x, \delta y, \delta z \leq 10e-6
        if (abs(v_m(1)) < 10e-6 .AND. & 
            abs(v_m(2)) < 10e-6 .AND. &
            abs(v_m(3)) < 10e-6) then
            print *, 'Number of iteration is', i
            goto 42
        end if


    end do


    print *, 'Exceed iteration limit (1000)'
    print *, ''

    ! Output the result
42  print *, 'The estimated coordinates of the target point are:'
    print *, 'X = ', v_location(1)
    print *, 'Y = ', v_location(2)
    print *, 'Z = ', v_location(3)

    print *, 'The errors are:'
    print *, 'delta x = ', v_m(1)
    print *, 'delta y = ', v_m(2)
    print *, 'delta z = ', v_m(3)


    deallocate(v_QC1, v_QC2, v_QC3, v_d)
    deallocate(m_Q)

end program GPSLocating