# Assignment 6

## Task: Plot the `On-land Earthquakes`.

1. **Coastal line**
   - From `1999.lis`.
2. **SUBROUTINE**
   - `locpt (x0, y0, x, y, n, l, m)`
      ```FORTRAN
      SUBROUTINE locpt (x0, y0, x, y, n, l, m)
      !-----------------------------------------------------------------------
      ! GIVEN A POLYGONAL LINE CONNECTING THE VERTICES (X(I),Y(I)) (I = 1,...,N)
      ! TAKEN IN THIS ORDER.  IT IS ASSUMED THAT THE POLYGONAL PATH IS A LOOP,
      ! WHERE (X(N),Y(N)) = (X(1),Y(1)) OR THERE IS AN ARC FROM (X(N),Y(N)) TO
      ! (X(1),Y(1)).  N.B. The polygon may cross itself any number of times.

      ! (X0,Y0) IS AN ARBITRARY POINT AND L AND M ARE VARIABLES.
      ! On output, L AND M ARE ASSIGNED THE FOLLOWING VALUES ...

      !    L = -1   IF (X0,Y0) IS OUTSIDE THE POLYGONAL PATH
      !    L =  0   IF (X0,Y0) LIES ON THE POLYGONAL PATH
      !    L =  1   IF (X0,Y0) IS INSIDE THE POLYGONAL PATH

      ! M = 0 IF (X0,Y0) IS ON OR OUTSIDE THE PATH.  IF (X0,Y0) IS INSIDE THE
      ! PATH THEN M IS THE WINDING NUMBER OF THE PATH AROUND THE POINT (X0,Y0).

      ! Fortran 66 version by A.H. Morris
      ! Converted to ELF90 compatibility by Alan Miller, 15 February 1997

      !-----------------------

      IMPLICIT NONE
      INTEGER:: n

      REAL :: x0, y0, x(n), y(n)

      INTEGER :: l, m

      !     Local variables
      INTEGER :: i, n0
      REAL    :: angle, eps, pi, pi2, sum, theta, theta1, thetai, tol, u, v

      !     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
      !            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0

      eps = EPSILON(1.0)

      !-----------------------------------------------------------------------
      n0 = n
      IF (x(1) == x(n) .AND. y(1) == y(n)) n0 = n - 1
      pi = ATAN2(0.0, -1.0)
      pi2 = 2.0*pi
      tol = 4.0*eps*pi
      l = -1
      m = 0

      u = x(1) - x0
      v = y(1) - y0
      IF (u == 0.0 .AND. v == 0.0) GO TO 20
      IF (n0 < 2) RETURN
      theta1 = ATAN2(v, u)

      sum = 0.0
      theta = theta1
      DO i = 2, n0
      u = x(i) - x0
      v = y(i) - y0
      IF (u == 0.0 .AND. v == 0.0) GO TO 20
      thetai = ATAN2(v, u)
      
      angle = ABS(thetai - theta)
      IF (ABS(angle - pi) < tol) GO TO 20
      IF (angle > pi) angle = angle - pi2
      IF (theta > thetai) angle = -angle
      sum = sum + angle
      theta = thetai
      END DO

      angle = ABS(theta1 - theta)
      IF (ABS(angle - pi) < tol) GO TO 20
      IF (angle > pi) angle = angle - pi2
      IF (theta > theta1) angle = -angle
      sum = sum + angle

      !     SUM = 2*PI*M WHERE M IS THE WINDING NUMBER

      m = ABS(sum)/pi2 + 0.2
      IF (m == 0) RETURN
      l = 1
      IF (sum < 0.0) m = -m
      RETURN

      !     (X0, Y0) IS ON THE BOUNDARY OF THE PATH

      20 l = 0
      RETURN
      END SUBROUTINE locpt

      ```


## Task: Linear Regression `Plots`

1. **Read data**: 
   - Epicentral distance as $X_i$ 
   - P travel time as $Y_i$
   - Data is provided in `ppfile.txt`

2. **Fit a linear model**: 
   - The equation for the line: $Y_i = a \cdot X_i + b$
   - Fit the data to this linear equation.

3. **Find the following values**:
   - Constant $a$
   - Constant $b$
   - **Standard deviation ($s.d.v.$)** of the fit
   - **Linear correlation coefficient ($R$)**
   - **Standard deviation of $a$ and $b$**

5. **SUBROUTINE**
   - `xyfit(x, y, ndata, cept, slop, rms, r_cor, sdv)`
      ```FORTRAN

      !
      !       fitting line using X,Y
      !
            subroutine xyfit(x,y,ndata,cept,slop,rms,r_cor,sdv)
            real*4 x(ndata),y(ndata)

            xsum  =0.0
            ysum  =0.0
            xxsum =0.0
            yysum =0.0
            xysum =0.0
            sum   =0.0
            errsum=0.0

            do i=1,ndata
               xp  =x(i)
               yp  =y(i)
               sum =sum+1.0
               xsum=xsum+xp
               ysum=ysum+yp
               xxsum=xxsum+xp*xp
               yysum=yysum+yp*yp
               xysum=xysum+xp*yp
            enddo


            r_cor=(sum*xysum-xsum*ysum)/sqrt( (sum*xxsum-xsum*xsum)*
         +        (sum*yysum-ysum*ysum) )


            xavg=xsum/sum
            yavg=ysum/sum
            su2 =xxsum-xsum*xavg
            sv2 =yysum-ysum*yavg
            suv =xysum-xsum*yavg
            sdx =sqrt( su2/(sum-1.0) )
            sdy =sqrt( sv2/(sum-1.0) )
            r2  =suv/sqrt(su2*sv2)
            suv2=sv2-su2
            b2  =( suv2+sqrt( suv2*suv2+4.0*suv*suv ) )/(2.0*suv)
            b1  =yavg-b2*xavg
            sdb2=b2*sqrt( (1.-r2*r2)/sum )/r2
            part1=(sdy-sdx*b2)**2.0
            part2=2.*sdx*sdy+(b2*(1.+r2)*xavg**2.0)/r2**2.
            sdb1=sqrt( ( part1+(1.-r2)*b2*part2 )/sum )

            cept=b1 !-  a     Y=aX+b
            slop=b2 !-  b

            do i=1,ndata
               err =b2*x(i)+b1-y(i)
               err2=err*err
               errsum=errsum+err2
            enddo

            sdv=sqrt(errsum/float(ndata-1))

            rms=sqrt(abs(errsum))/float(ndata)

            return
            end
      ```

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