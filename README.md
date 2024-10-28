# Assignment 7

## Task: The `Location program`

The linear system can be represented as:

$$
\overrightarrow{d} = \hat Q \cdot \overrightarrow{m}
$$

where:

$$
q_{x_i} m_x + q_{y_i} m_y + q_{z_i} m_z = d_i
$$

for $i = 1, 2, \cdots, n$. This can be written in matrix form as:

```math
\begin{pmatrix}
q_{x_1} & q_{y_1} & q_{z_1} \\\
\vdots & \vdots & \vdots \\\
q_{x_n} & q_{y_n} & q_{z_n}
\end{pmatrix}
\begin{pmatrix}
m_x \\\
m_y \\\
m_z
\end{pmatrix}
=
\begin{pmatrix}
d_1 \\\
\vdots \\\
d_n
\end{pmatrix}
```

Here, $\hat Q$ is an $n \times 3$ matrix, $\overrightarrow{m}$ is a $3 \times 1$ vector, and $\overrightarrow{d}$ is an $n \times 1$ vector.

To solve for \(M\), we use the equation:

$$ 
\overrightarrow{m} =  \left( \hat Q^{T} \cdot \hat Q \right)^{-1} \cdot \hat Q^{T} \cdot \overrightarrow{d} 
$$


1. **Input**
   - Input N stations $(x_i,y_i,z_i)$ and distance $D_i$ to target point.

2. **Output**
   - Calculating the $(X,Y,Z)$ of target point. 
   - Just like GPS location.
3. **Output**
   - I use my custom subroutine, and call `LAPAC` to do the work.
        ```FORTRAN
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
        ```