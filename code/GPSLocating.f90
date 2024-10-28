program GPSLocating
    use, intrinsic :: iso_fortran_env
    implicit none

    integer, parameter :: dp = real64
    integer :: N, i, j, k, INFO
    real(dp), allocatable :: v_QC1(:), v_QC2(:), v_QC3(:), v_d(:)
    real(dp) :: v_m(3)
    real(dp), allocatable :: m_Q(:,:), m_QT(:,:)
    real(dp) :: m_QTQ(3, 3), m_QTQ_INVERSE(3, 3)

    ! Prompt user for the number of stations
    print *, 'Enter the number of stations (N >= 4):'
    read *, N
    if (N < 4) then
        print *, 'At least 4 stations are required.'
        stop
    end if

    ! Allocate arrays
    allocate(v_QC1(N), v_QC2(N), v_QC3(N), v_d(N))

    ! Input station coorv_dnates and v_dstances
    print *, 'Enter the coordinates (v_QC1, v_QC2, v_QC3) and stance v_d for each station:'
    do i = 1, N
        print *, 'Station ', i, ':'
        read *, v_QC1(i), v_QC2(i), v_QC3(i), v_d(i)
    end do

    allocate(m_Q(N, 3), m_QT(3, N))
    do j = 1, N
        m_Q(j, 1) = v_QC1(j)
        m_Q(j, 2) = v_QC2(j)
        m_Q(j, 3) = v_QC3(j)
    end do

    ! Print the matrix to check
    print *, 'Matrix Q:'
    do k = 1, N
        print *, m_Q(k, 1), m_Q(k, 2), m_Q(k, 3)
    end do
    

    m_QT = transpose(m_Q)
    m_QTQ = matmul(m_QT, m_Q)
    call matrix_inverse(3, m_QTQ, m_QTQ_INVERSE, INFO)
    v_m = matmul(m_QTQ_INVERSE, matmul(m_QT, v_d))


    ! Output the result
    print *, 'The estimated coordinates of the target point are:'
    print *, 'X = ', v_m(1)
    print *, 'Y = ', v_m(2)
    print *, 'Z = ', v_m(3)



    deallocate(v_QC1, v_QC2, v_QC3, v_d)
    deallocate(m_Q)

end program GPSLocating