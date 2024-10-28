program InverseTest
    implicit none
    integer :: N, INFO
    double precision, allocatable :: A(:,:), AINV(:,:)
    integer :: i, j, k

    print *, 'Enter the size of the matrix:'
    read *, N

    allocate(A(N, N), AINV(N, N))

    print *, 'Enter the row elements:'
    do i = 1, N
        do j = 1, N 
            print *, 'Row ', i, ', Column ', j, ':'
            read *, A(j, i)
        end do
    end do

    print *, A

    call matrix_inverse(N, A, AINV, INFO)

    if (INFO == 0) then
        print *, 'Original Matrix, A:'
        do i = 1, N
            print '(3F10.4)', (A(j, i), j = 1, N)
        end do

        print *, 'Inverse Matrix, AINV:'
        do i = 1, N
            print '(3F10.4)', (AINV(j, i), j = 1, N)
        end do
    else
        print *, 'Error, INFO = ', INFO
    end if

    deallocate(A, AINV)

end program InverseTest
