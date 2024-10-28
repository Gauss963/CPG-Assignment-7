subroutine LinearRegressionFitting(Xi, Yi, n, a, b, sdv, R, std_a, std_b)
    implicit none
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: Xi, Yi
    real, intent(out) :: a, b, sdv, R, std_a, std_b
    real :: sum_X, sum_Y, sum_X2, sum_Y2, sum_XY
    real :: denom
    integer :: i

    ! Initialize
    sum_X = 0.0
    sum_Y = 0.0
    sum_X2 = 0.0
    sum_Y2 = 0.0
    sum_XY = 0.0

    ! Calculate the necessity
    do i = 1, n
        sum_X = sum_X + Xi(i)
        sum_Y = sum_Y + Yi(i)
        sum_X2 = sum_X2 + Xi(i)**2
        sum_Y2 = sum_Y2 + Yi(i)**2
        sum_XY = sum_XY + Xi(i) * Yi(i)
    end do

    ! Calculate constant $a$ and $b$
    denom = n * sum_X2 - sum_X**2
    a = (n * sum_XY - sum_X * sum_Y) / denom
    b = (sum_Y * sum_X2 - sum_X * sum_XY) / denom

    ! Calculate sdv
    sdv = 0.0
    do i = 1, n
        sdv = sdv + (Yi(i) - (a * Xi(i) + b))**2
    end do
    sdv = sqrt(sdv / (n - 2))
    

    ! Calculate R
    R = (n * sum_XY - sum_X * sum_Y) / sqrt((n * sum_X2 - sum_X**2) * (n * sum_Y2 - sum_Y**2))

    ! Calculate sdv of a and b
    std_a = sdv * sqrt(real(n)) / sqrt(denom)
    std_b = sdv * sqrt(sum_X2 / denom)

end subroutine LinearRegressionFitting