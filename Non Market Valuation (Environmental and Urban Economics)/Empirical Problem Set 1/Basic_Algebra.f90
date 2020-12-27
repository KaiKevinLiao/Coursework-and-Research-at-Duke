! #########
! This module containing routine performing linear algebra
! #########
    
Module Basic_Algebra

    Contains
    
    
function variance(A,n,flag)

! flag == 1 -> A is population
! flag == 0 -> A is sample

    USE Global_Data
    
    implicit none
    
    integer, intent(in) :: n,flag
    real(KIND=DOUBLE), intent(in) :: A(n)
    real(KIND=DOUBLE) variance, sum, sum_sqr
    integer i
    
    sum = 0
    sum_sqr = 0
    do i = 1,n
        sum = sum + A(i)
        sum_sqr = sum_sqr + A(i)**2
    end do
    
    if (flag == 0) then
        variance = (sum_sqr - sum**2/n)/(n-1)
    else
        variance = (sum_sqr - sum**2/n)/(n)
    end if
end function variance


function inverse(A,n)
    
    USE Global_Data
    
    implicit none
    
    integer, intent(in) :: n
    real(KIND=DOUBLE), intent(in) :: A(n,n)
    real(KIND=DOUBLE) inverse(n,n), work(n), Ainv(n,n)
    integer info1, info2
    
    integer ipiv(n)
    
    Ainv = A
    Call dgetrf ( n , n , Ainv , n , ipiv, info1 )
    
    Call dgetri ( n , Ainv , n , ipiv, work, n, info2 )
    
    if ( ( info1 .ne. 0 ) .or. ( info2 .ne. 0 ) ) then
        write(*,*) 'Problem Computing Inverse'
    end if
    
    inverse = Ainv
    
end function inverse
    
! Matrix Multiplication
! Calls MKL routine dgemm
function matrix_mult(A,B,m,n,k)
    
    USE Global_Data

    implicit none

    ! integer, parameter :: DOUBLE     = SELECTED_REAL_KIND(p=10)

    integer, intent(in) :: m,n,k
    real(KIND=DOUBLE), intent(in) :: A(m,n), B(n,k)

    real(KIND=DOUBLE) alpha0, beta0, matrix_mult(m,k)

    alpha0 = 1.0
    beta0  = 0.0

    call dgemm('N', 'N', m, k, n, alpha0, A, m, B, n, beta0, matrix_mult, m)

end function matrix_mult

subroutine interp1( xData, yData, xVal, yVal )
! Inputs: xData = a vector of the x-values of the data to be interpolated
!         yData = a vector of the y-values of the data to be interpolated
!         xVal  = a vector of the x-values where interpolation should be performed
! Output: yVal  = a vector of the resulting interpolated values

USE Global_Data    

  implicit none

  real(KIND=DOUBLE), DIMENSION(2000), intent(in)             :: xData, yData
  real(KIND=DOUBLE), DIMENSION(numBUYLA), intent(in)         :: xVal
  real(KIND=DOUBLE), DIMENSION(numBUYLA), intent(out)        :: yVal
  integer :: inputIndex, dataIndex
  real(KIND=DOUBLE) :: minXdata, minYdata, xRange, weight,maxXData

  ! Possible checks on inputs could go here
  ! Things you may want to check:
  !   monotonically increasing xData
  !   size(xData) == size(yData)
  !   size(xVal) == size(yVal)

  minXData = xData(1)
  maxXData = xData(2000)
  xRange = maxXData - minXData

  do inputIndex = 1, size(xVal)
      ! possible checks for out of range xVal could go here

      ! this will work if x is uniformly spaced, otherwise increment
      ! dataIndex until xData(dataIndex+1)>xVal(inputIndex)
      dataIndex = min(floor(xVal(inputIndex)),1999);

      weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex));
      yVal(inputIndex) = yData(dataIndex) +  weight*( yData(dataIndex+1)-yData(dataIndex) );
  end do
end subroutine

end module Basic_Algebra