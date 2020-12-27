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
    
    
Module Non_Parametric_MWTP_Estimates

    Contains

subroutine Non_Parametric_MWTP_Estimates_Results(data_buyers_la)

    USE Global_Data
    USE Basic_Algebra
    
    integer i, j

    ! Data
    real(KIND=DOUBLE), DIMENSION(numBUYLA, numVBLA), intent(in)         :: data_buyers_la
    real(KIND=DOUBLE), DIMENSION(numBUYLA)                              :: data_VC
    real(KIND=DOUBLE), DIMENSION(numBUYLA)                              :: beta
    
    ! working matrix for data read
    real(KIND=DOUBLE), DIMENSION(2000,9)                                :: f
    ! Variables read from non-parametric estimation of hedonic price gradient with h = 1
    real(KIND=DOUBLE), DIMENSION(2000)                                  :: beta_src
    real(KIND=DOUBLE), DIMENSION(2000)                                  :: VC_src
    
    ! independent variables matrix
    real(KIND=DOUBLE), DIMENSION(numBUYLA, 5)                           :: A
    
    ! Parameter for calling regression function
    integer                                                     m, n, nrhs, lda, ldb, rank, lwork, info
    integer, DIMENSION(50)                                   :: jpvt
    real(KIND=DOUBLE), DIMENSION(1000000)                    :: work
    real(KIND=DOUBLE)                                        :: rcond
    
    ! Regression results
    real(KIND=DOUBLE), DIMENSION(5)                          :: beta_results
    real(KIND=DOUBLE), DIMENSION(5)                          :: beta_variance
    
    ! ###########
    ! Read non-parametric estimation of hedonic price gradient
    ! ###########
    open(unit = 170 , file = trim(dir_mom)//trim(filename_non_paramatric_hedonic_price_gradient)//'.csv')
    read(170,*)
    read(170,*)
    do i = 1,2000
        read(170,*) f(i,:)
    end do 
    VC_src = f(:,1)
    beta_src = f(:,3)
    
    ! ################
    ! Interpolate beta
    ! ################
    
    data_VC(:) = data_buyers_la(:,_VIOLENT_BUYER)
    call interp1( VC_src, beta_src, data_VC, beta )
    
    
    ! ###########
    ! Generate independent variables matrix
    ! ###########
    
    A(:,1) = 1
    A(:,2) = data_buyers_la(:,_INCOME)
    do i = 1,numBUYLA
        if (int(data_buyers_la(i,_RACE)) == 2) then
            A(i,3) = 1
        else if (int(data_buyers_la(i,_RACE)) == 3) then
            A(i,4) = 1
        else if (int(data_buyers_la(i,_RACE)) == 4) then
            A(i,5) = 1
        end if
    end do
    
    ! ##########
    ! Regression
    ! ##########
    
    m = numBUYLA
    n = 5
    nrhs = 1
    lda = numBUYLA
    ldb = numBUYLA
    lwork = 1000000
    call dgelsy(m, n, nrhs, A, lda, beta, ldb, jpvt, rcond, rank, work, lwork, info)
    
    do i = 1,5
        beta_results(i) = beta(i)
    end do
    
    ! ##########
    ! Write regression results
    ! ##########

    open(unit = 171 , file = trim(dir_mom)//trim(filename_non_paramatric_MWTP)//'.csv')
    write(171,870) 'Constant,Income, Asian/Pacific, Black, Hispanic'
    write(171,871) beta_results(1), ',',beta_results(2), ',',beta_results(3), ',',beta_results(4), ',',beta_results(5)
    
870 format(a47)
871 format(4(f50.8,a1), f50.8)    
    
end subroutine Non_Parametric_MWTP_Estimates_Results


end Module Non_Parametric_MWTP_Estimates