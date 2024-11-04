program GPSLocating
    
    use, intrinsic :: iso_fortran_env
    
    implicit none

    integer, parameter :: dp = real64
    integer :: N, i, j, k, l, m, INFO, o, tg_num, y


    real, allocatable :: v_QC1(:), v_QC2(:), v_QC3(:), v_d(:), v_d_FIXED(:)
    real :: v_m(4), location(4)
    real, allocatable :: m_Q(:,:), m_QT(:,:), m_Q_FIXED(:,:)
    real :: m_QTQ(4, 4), m_QTQ_INVERSE(4, 4)


    character(len = 100) :: line, elev_str, lat_str, lon_str, code_str
    character(len = 100) :: data_filename
    real, allocatable :: STATION_X(:), STATION_Y(:), STATION_Z(:)
    real, allocatable :: STATION_TG_X(:), STATION_TG_Y(:), STATION_TG_Z(:)
    real :: lat_deg, lat_min, lon_deg, lon_min
    integer :: io_status, unit_num


    real, allocatable :: TRAVEL_DIST(:), ARRIVAL_TIMES(:)
    character(len = 4), allocatable :: STA_CODES(:), STA_TG_CODES(:)
    real :: a, b, sdv, R, std_a, std_b

    real :: xsec
    integer iy, im, id, ih, mm
    real :: velocity
    real :: kx, ky, kd
    real :: rx, ry, rz, rt
    real :: kx_guess, ky_guess, kz_guess
    
    
    velocity = 7 ! km/s







    ! Count the number of lines.
    call CountLines("../data/nsta.dat", o)
    allocate(STATION_X(o), STATION_Y(o), STATION_Z(o), STA_CODES(o))

    ! Read data (lat, lon, elev)
    open(newunit=unit_num, file = "../data/nsta.dat", status = "old", action = "read")
    do i = 1, m
        read(unit_num, '(A)', iostat=io_status) line
        if (io_status /= 0) exit

        read(line(1:4), *) STA_CODES(i)

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

    allocate(TRAVEL_DIST(tg_num), ARRIVAL_TIMES(tg_num), STA_TG_CODES(tg_num))
    allocate(STATION_TG_X(tg_num), STATION_TG_Y(tg_num), STATION_TG_Z(tg_num))
    allocate(v_d(tg_num))

    open(newunit=unit_num, file = "../data/ppfile.txt", status = "old", action = "read")
    rewind(unit_num)

    read(unit_num,'(1x, i4, 4i2, f6.2)') iy, im, id, ih, mm, xsec
    xsec = mm * 60.0 + xsec
    do i = 1, tg_num
        read(unit_num,'(A5, f6.1, 9x, i3, f6.2)') STA_TG_CODES(i), TRAVEL_DIST(i), mm, ARRIVAL_TIMES(i)
        ! ARRIVAL_TIMES(i) = mm * 60.0 + ARRIVAL_TIMES(i)
    end do
    close(unit_num)


    do i = 1, size(STA_TG_CODES)

        do j = 1, size(STA_CODES)

            if (ADJUSTL(TRIM(STA_TG_CODES(i))) == ADJUSTL(TRIM(STA_CODES(j)))) then
                STATION_TG_X(i) = STATION_X(j)
                STATION_TG_Y(i) = STATION_Y(j)
                STATION_TG_Z(i) = STATION_Z(j)
                goto 43
            end if
        end do
43      print *, i, "is done, next station"
    end do






    allocate(m_Q(tg_num, 4), m_QT(4, tg_num), m_Q_FIXED(tg_num, 4))

    ! location = [121, 25, -10, 0]
    location = [0, 0, 10, 0] ! 121, 24

    print *, ARRIVAL_TIMES

    do y = 1, 1000

        do i = 1, tg_num
            ! call delaz(24, 121, location(2), location(1), kx_guess, ky_guess, kz_guess)
            call delaz(24.0, 121.0, STATION_TG_Y(i), STATION_TG_X(i), kx, ky, kd)

            ! Update vector d
            rx = (kx - location(1))
            ry = (ky - location(2))
            rz = (-STATION_TG_Z(i) / 1000 - location(3))
            rt = ( ARRIVAL_TIMES(i) - location(4))
            v_d(i) = 0.5 * (rx**2 + ry**2 + rz**2 - velocity**2 * rt**2)
            ! Update matrix Q
            m_Q(i, 1) = rx
            m_Q(i, 2) = ry
            m_Q(i, 3) = rz
            m_Q(i, 4) = -velocity**2 * rt
        end do
        ! Solve linear system
        m_QT = transpose(m_Q)
        m_QTQ = matmul(m_QT, m_Q)
        call MATRIXINV(m_QTQ, 4)
        v_m = matmul(matmul(m_QTQ, m_QT), v_d)
        location = location + v_m

        print *, v_m

        ! Check if \delta x, \delta y, \delta z \leq 10e-6
        ! if (abs(v_m(1)) < 10e-6 .AND. & 
        !     abs(v_m(2)) < 10e-6 .AND. &
        !     abs(v_m(3)) < 10e-6) then
        !     print *, 'Number of iteration is', i
        !     goto 42
        ! end if


    end do


    print *, 'Exceed iteration limit (1000)'
    print *, ''

    ! Output the result
42  print *, 'The estimated coordinates of the target point are:'
    print *, 'X = ', location(1)
    print *, 'Y = ', location(2)
    print *, 'Z = ', location(3)
    print *, 'T = ', location(4)

    print *, 'The errors are:'
    print *, 'delta x = ', v_m(1)
    print *, 'delta y = ', v_m(2)
    print *, 'delta z = ', v_m(3)


    deallocate(v_d)
    deallocate(m_Q)

end program GPSLocating