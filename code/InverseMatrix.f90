subroutine matrix_inverse(N, A, AINV, INFO)
    implicit none
    integer, intent(in) :: N
    double precision, intent(in) :: A(N, N)
    double precision, intent(out) :: AINV(N, N)
    integer, intent(out) :: INFO

    integer, allocatable :: IPIV(:)
    double precision, allocatable :: WORK(:)
    integer :: LWORK

    AINV = A

    allocate(IPIV(N))

    ! LU decomposition
    call dgetrf(N, N, AINV, N, IPIV, INFO)
    if (INFO /= 0) then
        print *, 'DGETRF error, INFO = ', INFO
        deallocate(IPIV)
        return
    end if

    LWORK = -1
    allocate(WORK(1))
    call dgetri(N, AINV, N, IPIV, WORK, LWORK, INFO)
    if (INFO /= 0) then
        print *, 'DGETRI error, INFO = ', INFO
        deallocate(IPIV, WORK)
        return
    end if


    LWORK = int(WORK(1))
    deallocate(WORK)
    allocate(WORK(LWORK))


    call dgetri(N, AINV, N, IPIV, WORK, LWORK, INFO)
    if (INFO /= 0) then
        print *, 'DGETRI error, INFO = ', INFO
        deallocate(IPIV, WORK)
        return
    end if


    deallocate(IPIV, WORK)

end subroutine matrix_inverse
