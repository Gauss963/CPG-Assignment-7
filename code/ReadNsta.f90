 program ReadNsta
    implicit none
    character(len = 100) :: line, elev_str, lat_str, lon_str
    character(len = 100) :: data_filename
    
    integer :: n, m, i
    integer :: io_status, unit_num

    real, allocatable :: STATION_X(:), STATION_Y(:), STATION_Z(:)
    real, allocatable :: TRAVEL_TIME(:)
    
    real :: lat_deg, lat_min, lon_deg, lon_min


    data_filename = "../data/nsta.dat"


    ! Count the number of lines.
    call CountLines(data_filename, m)
    allocate(STATION_X(m), STATION_Y(m), STATION_Z(m))


    ! Read data (lat, lon, elev)
    open(newunit=unit_num, file = data_filename, status = "old", action = "read")
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

    



    deallocate(STATION_X, STATION_Y, STATION_Z)
end program ReadNsta