#define _HOUSEID 1
#define _PRICE 2
#define _COUNTYID 3
#define _BUILTYEAR 4
#define _SQUAREFEET 5
#define _BATHROOMS 6
#define _BEDROOMS 7
#define _ROOMS 8
#define _STORIES 9
#define _VIOLENT 10
#define _PROPERTY 11
#define _SALEYEAR 12

#define _BUYERID 1
#define _VIOLENT_BUYER 3
#define _PROPERTY_BUYER 4
#define _RACE 5
#define _INCOME 6
#define _CITY 7
    
Module Bootstrapped_Multimarket_Rosen_Estimates

    Contains
    
subroutine Bootstrapped_Rosen_Estimates_Results(data_buyers, numBS)
    
    USE Global_Data
    USE Basic_Algebra
    
    integer i,j,k
    
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB), intent(in)   :: data_buyers
    integer, intent(in)                                       :: numBS
    real(KIND=DOUBLE), DIMENSION(2)                           :: grad_la, grad_sf
    real(KIND=DOUBLE), DIMENSION(12)                          :: f
    
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB)               :: boot_data_buyers
    real(KIND=DOUBLE), DIMENSION(numBS, 6)                    :: beta
    real(KIND=DOUBLE), DIMENSION(6)                           :: beta_temp
    real(KIND=DOUBLE), DIMENSION(6)                           :: beta_data
    ! for variance calculation
    real(KIND=DOUBLE), DIMENSION(numBS)                       :: A
    ! variance
    real(KIND=DOUBLE), DIMENSION(6)                           :: var
    
    ! ###########
    ! Read Hedoic Price regression results
    ! ###########
    
    open(unit = 130 , file = trim(dir_mom)//trim(filename_hedonic_price_function_regression)//'_LA'//'.csv')
    open(unit = 131 , file = trim(dir_mom)//trim(filename_hedonic_price_function_regression)//'_SF'//'.csv')
    
    read(130,*)
    read(130,*) f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11),f(12),grad_la(1),grad_la(2)
    read(131,*)
    read(131,*) f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11),f(12),grad_sf(1),grad_sf(2)
    print*, grad_la(1),grad_la(2),grad_sf(1),grad_sf(2)
    
    !!!!!!!!!!!!!!!!
    ! BIG BIG CHANGE
    !!!!!!!!!!!!!!!!!!!!!!
    grad_la(2) = 2 * grad_la(2)
    grad_sf(2) = 2 * grad_sf(2)
    
    open(unit = 140 , file = trim(dir_mom)//trim(filename_rosen_estimate)//'.csv')
    open(unit = 141 , file = trim(dir_mom)//trim(filename_bootstrapped_rosen_estimate)//'.csv')
    
    write(140, 900) 'Constant, Violent Crime Rate (cases per 100000 residents), Income, Asian/Pacific Islander, Black, Hispanic'
    write(141, 900) 'Constant, Violent Crime Rate (cases per 100000 residents), Income, Asian/Pacific Islander, Black, Hispanic'
    
    ! ###########
    ! Compute point estimate
    ! ###########
    beta_data = 0
    Call Rosen_Estimates_Results(data_buyers, beta_data, grad_la, grad_sf)
    
    ! ###########
    ! Write point estimate
    ! ###########
    
    do j = 1,5
        write(140, 612, advance = 'no') beta_data(j),','
    end do
    write(140, 613) beta_data(6)
    
    ! ###########
    ! Compute bootstrapped estimators
    ! ###########
    beta = 0
    do i = 1,numBS
        beta_temp = 0
        boot_data_buyers = 0
        Call Bootstrapped_Buyers(data_buyers, boot_data_buyers)
        Call Rosen_Estimates_Results(boot_data_buyers, beta_temp, grad_la, grad_sf)
        do j = 1,6
            beta(i,j) = beta_temp(j)
        end do
        ! ###########
        ! Write bootstrapped estimators
        ! ###########
        
        do j = 1,5
            write(141, 612, advance = 'no') beta(i,j),','
        end do
        write(141, 613) beta(i,6)
        
        write(*,901, advance = 'no') 'rosen interation: '
        write(*,*) i
    end do
    
    ! ###########
    ! Compute estimators variances
    ! ###########
    
    do j = 1,6
        A = beta(:,j)
        var(j) = variance(A,numBS,0)
    end do
    
    ! ###########
    ! Write estimators variances
    ! ###########
    do j = 1,5
        write(140, 612, advance = 'no') var(j),','
    end do
    write(140, 613) var(6)
    
    
901        format(a20)
900        format(a107)
612        format(f50.8,a1)
613        format(f50.8) 
end subroutine Bootstrapped_Rosen_Estimates_Results



subroutine Rosen_Estimates_Results(data_buyers, beta, grad_la, grad_sf)

    USE Global_Data
    
    integer i, j, k
    
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB), intent(in)   :: data_buyers
    real(KIND=DOUBLE), DIMENSION(6), intent(out)              :: beta
    
    ! Hedoic Price regression results
    real(KIND=DOUBLE), DIMENSION(2), intent(in)               :: grad_la, grad_sf
    ! Buyer Matrix
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB)               :: A
    ! MWTP Unobservables
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: nu                
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: b
    
    
    ! Parameter for calling regression function
    integer                                                     m, n, nrhs, lda, ldb, rank, lwork, info
    integer, DIMENSION(50)                                   :: jpvt
    real(KIND=DOUBLE), DIMENSION(1000000)                    :: work
    real(KIND=DOUBLE)                                        :: rcond
    


    
    ! ###########    
    ! Generate independent variables matrix, with race dummy (omite white)
    ! ###########
    A = 0
    A(:,1) = 1.
    A(:,2) = data_buyers(:,_VIOLENT_BUYER)
    A(:,3) = data_buyers(:,_INCOME)
    do i = 1,numBUY
        if (int(data_buyers(i,_RACE)) == 2) then
            A(i,4) = 1
        else if (int(data_buyers(i,_RACE)) == 3) then
            A(i,5) = 1
        else if (int(data_buyers(i,_RACE)) == 4) then
            A(i,6) = 1
        end if
    end do

    
    
    ! #########
    ! Generate dependent variable: Implicit Price
    ! #########
    b = 0
    do i = 1,numBUY
        ! LA
        if (data_buyers(i,_CITY) == 1) then
            b(i) = data_buyers(i,_VIOLENT_BUYER) * grad_la(2) + grad_la(1)
        else
            b(i) = data_buyers(i,_VIOLENT_BUYER) * grad_sf(2) + grad_sf(1)
        end if
    end do
    
    print*, b(10), b(555555)
    
    ! #########
    ! Do regression
    ! #########
    
    m = numBUY
    n = 6
    nrhs = 1
    lda = numBUY
    ldb = numBUY
    lwork = 1000000
    call dgelsy(m, n, nrhs, A, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
    
    do i = 1,6
        beta(i) = b(i)
    end do
       
end subroutine Rosen_Estimates_Results



subroutine Bootstrapped_Buyers(data_buyers, boot_data_buyers)
    
    USE Global_Data
    
    integer i, j, k
    
    real(KIND=DOUBLE) :: u
    
    ! Data
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB), intent(in)   :: data_buyers
    
    ! Bootstrapped output
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB), intent(out)  :: boot_data_buyers
    
    ! Bootstrapped index
    integer   r_buyers(numBUY)
    
    
    ! Generate random index array
    do i = 1,numBUY
        call random_number(u)
        k = 1 + FLOOR(numBUY*u)  
        r_buyers(i) = k
    end do
    
    
    ! Construct bootstrapped dataset
    do i = 1, numBUY
        boot_data_buyers(i,:) = data_buyers(r_buyers(i),:)
    end do
    

end subroutine Bootstrapped_Buyers
    
    
end Module Bootstrapped_Multimarket_Rosen_Estimates