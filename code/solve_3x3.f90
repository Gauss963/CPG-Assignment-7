subroutine solve_3x3(A, B, X)
    use, intrinsic :: iso_fortran_env
    implicit none
    
    integer, parameter :: dp = real64

    real(dp), intent(in) :: A(3, 3), B(3)
    real(dp), intent(out) :: X(3)
    real(dp) :: detA, detA1, detA2, detA3
    real(dp), dimension(3, 3) :: A1, A2, A3
    

    ! Compute determinant of A
    detA = determinant_3x3(A)
    if (abs(detA) < 1.0e-12_dp) then
        print *, 'Matrix is singular or ill-conditioned.'
        stop
    end if

    ! Replace columns of A with B to compute determinants
    A1 = A
    A1(:,1) = B
    detA1 = determinant_3x3(A1)

    A2 = A
    A2(:,2) = B
    detA2 = determinant_3x3(A2)

    A3 = A
    A3(:,3) = B
    detA3 = determinant_3x3(A3)

    ! Compute solutions
    X(1) = detA1 / detA
    X(2) = detA2 / detA
    X(3) = detA3 / detA


    contains

    function determinant_3x3(M) result(det)
    real(dp), intent(in) :: M(3,3)
    real(dp) :: det
    det = M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) &
        - M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) &
        + M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1))
    end function determinant_3x3

end subroutine solve_3x3