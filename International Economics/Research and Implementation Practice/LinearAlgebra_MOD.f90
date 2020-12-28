! Module containing routine performing linear algebra operations
    
Module LinearAlgebra_MOD
    
Contains

! Obtains eigenvalues and eigenvectors of matrix A    
! Calls MKL routine DGEEV
subroutine eigen(A,N,eigenval,eigenvec,info)

    USE Global_Data
    
    implicit none
    
    real(KIND=DOUBLE), intent(in) :: A(:,:)
    integer, intent(in) :: N
    real(KIND=DOUBLE), intent(out) :: eigenval(:,:)
    real(KIND=DOUBLE), intent(out) :: eigenvec(:,:)
    integer, intent(out) :: info

    INTEGER LWORK

    real(KIND=DOUBLE) VL( N , N ), VR( N , N ), &
                      WR( N ), WI( N ), WORK( 1000 )
    
    integer LDA, LDVL, LDVR, LWMAX
    
    LDA = N
    LDVL = N 
    LDVR = N 
    LWMAX = 1000

    LWORK = -1
    CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, &
                 VR, LDVR, WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

    CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL, &
                 VR, LDVR, WORK, LWORK, INFO )
    
    eigenval(:,1) = WR(:)
    eigenvec      = VR
    
end subroutine eigen

! Computes the inverse of matrix A
! Calls MKL routines dgetrf and dgetri
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

! Creates matrix of dimensions sizevec by ncol 
! replicating ncol columns given by vec
function repmat_col(vec,sizevec,ncol)

    USE Global_Data
    
    implicit none
    
    real(KIND=DOUBLE), intent(in) :: vec(:)
    
    integer, intent(in) :: sizevec, ncol
    
    real(KIND=DOUBLE) repmat_col(sizevec,ncol)
    
    integer j
    
    do j = 1, ncol
        repmat_col(1:sizevec,j) = vec(1:sizevec)
    end do
    
end function repmat_col

! Creates matrix of dimensions nrow by sizevec
! replicating nrow rows given by vec
function repmat_row(vec,sizevec,nrow)

    USE Global_Data
    
    implicit none
    
    real(KIND=DOUBLE), intent(in) :: vec(:)
    
    integer, intent(in) :: sizevec, nrow
    
    real(KIND=DOUBLE) repmat_row(nrow,sizevec), repmat_row_tr(sizevec,nrow)
    
    integer j
    
    do j = 1, nrow
        repmat_row_tr(1:sizevec,j) = vec(1:sizevec)
    end do
    
    repmat_row = transpose(repmat_row_tr)
    
end function repmat_row

end module LinearAlgebra_MOD