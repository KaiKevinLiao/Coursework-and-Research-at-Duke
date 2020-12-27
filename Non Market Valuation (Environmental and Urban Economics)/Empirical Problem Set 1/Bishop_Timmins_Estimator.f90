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
    
#define _buyerid 1
#define _violent_buyer 3
#define _property_buyer 4
#define _race 5
#define _income 6
#define _city 7


Module Bishop_Timmins_Estimator

    Contains
    
subroutine Bootstrapped_MLE(data_buyers, numBS)
    
    USE Global_Data
    USE Bootstrapped_Hedonic_Price_Function
    
    implicit none

    integer i,j,k
    
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB), intent(in)   :: data_buyers
    integer, intent(in)                                       :: numBS
    real(KIND=DOUBLE), DIMENSION(6)                           :: grad_la, grad_sf
    real(KIND=DOUBLE), DIMENSION(12)                          :: f
    
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB)               :: boot_data_buyers
    real(KIND=DOUBLE), DIMENSION(numBS, 6)                    :: beta
    real(KIND=DOUBLE), DIMENSION(6)                           :: beta_data
    ! for variance calculation
    real(KIND=DOUBLE), DIMENSION(numBS)                       :: A
    ! variance
    real(KIND=DOUBLE), DIMENSION(6)                           :: var
    
    real(KIND=DOUBLE), DIMENSION(11)                          :: para
    
    real(KIND=DOUBLE)                                         :: alpha2, alpha2_LB, alpha2_UB
    real(KIND=DOUBLE), DIMENSION(10)                          :: alpha
    
    
    open(unit = 530 , file = trim(dir_mom)//trim(filename_bootstrapped_six_order)//'_LA'//'.csv')
    open(unit = 531 , file = trim(dir_mom)//trim(filename_bootstrapped_six_order)//'_SF'//'.csv')
    
    read(530,*)
    read(531,*)
    
    open(unit = 140 , file = trim(dir_mom)//trim(filename_timmis_estimate)//'.csv')
    open(unit = 141 , file = trim(dir_mom)//trim(filename_bootstrapped_timmis_estimate)//'.csv')
    
    write(140, 900) 'Violent Crime Rate (cases per 100000 residents), Constant,  Income, Asian/Pacific Islander, Black, Hispanic&
        ,CityDummy, Income*Dummy, Asian/Pacific Islander*Dummy, Black*Dummy, Hispanic*Dummy'
    
    write(141, 900) 'Violent Crime Rate (cases per 100000 residents), Constant,  Income, Asian/Pacific Islander, Black, Hispanic&
        ,CityDummy, Income*Dummy, Asian/Pacific Islander*Dummy, Black*Dummy, Hispanic*Dummy'
    
    do i = 1,numBS
        
    
    read(530,*) f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11),f(12),grad_la(1),grad_la(2),grad_la(3),grad_la(4),grad_la(5),grad_la(6)
    
    read(531,*) f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11),f(12),grad_sf(1),grad_sf(2),grad_sf(3),grad_sf(4),grad_sf(5),grad_sf(6)
    ! print*, grad_la(1),grad_la(2),grad_la(3),grad_la(4),grad_la(5),grad_la(6),grad_sf(1),grad_sf(2),grad_sf(3),grad_sf(4),grad_sf(5),grad_sf(6)
    
    
    
    ! ###########
    ! Compute estimates
    ! ###########
    beta_data = 0
    alpha2_LB = -20
    alpha2_UB = 20
    Call Timmis_Estimates_Optimization(data_buyers, grad_la, grad_sf, alpha2_LB, alpha2_UB, alpha, alpha2)
    print *, 'final result ', i ,alpha2, alpha(1),alpha(2),alpha(3),alpha(4),alpha(5),alpha(6)

    ! ###########
    ! Write estimates
    ! ###########
    
    write(141, 612, advance = 'no') alpha2,','
    do j = 1,9
        write(141, 612, advance = 'no') alpha(j),','
    end do
    write(141, 613) alpha(10)
    
    
    end do

    
    
900        format(a190)
612        format(f50.8,a1)
613        format(f50.8) 
end subroutine Bootstrapped_MLE

subroutine Timmis_Estimates_Optimization(data_buyers, grad_la, grad_sf, alpha2_LB, alpha2_UB, alpha, alpha2)
    
    USE Global_Data
    
    integer i, j, k, flag

    real(KIND=DOUBLE), DIMENSION(numBUY, numVB), intent(in)   :: data_buyers
    real(KIND=DOUBLE), DIMENSION(6), intent(in)               :: grad_la, grad_sf
    real(KIND=DOUBLE), intent(in)                             :: alpha2_LB, alpha2_UB
    real(KIND=DOUBLE)                                         :: alpha2_left, alpha2_right, gap, loss_left, loss_right, loss_leftleft
    real(KIND=DOUBLE), DIMENSION(10), intent(out)             :: alpha
    real(KIND=DOUBLE), intent(out)                            :: alpha2
    
    flag = 0
    gap = INIT_GAP
    alpha2_left = alpha2_LB
    alpha2_right = alpha2_left + gap
    call Timmis_Estimates_MLE(data_buyers, grad_la, grad_sf, alpha2_left, loss_left, alpha)
    call Timmis_Estimates_MLE(data_buyers, grad_la, grad_sf, alpha2_right, loss_right, alpha)
    do while (gap > 0.0001)

        do while (loss_right < loss_left .AND. alpha2_right < alpha2_UB .OR. flag == 1)
            if (flag == 0) then
                alpha2_leftleft = alpha2_left
                alpha2_left = alpha2_right
            end if 
            alpha2_right = alpha2_left + gap
            call Timmis_Estimates_MLE(data_buyers, grad_la, grad_sf, alpha2_left, loss_left, alpha)
            call Timmis_Estimates_MLE(data_buyers, grad_la, grad_sf, alpha2_right, loss_right, alpha)
            flag = 0
        end do
        
        alpha2 = alpha2_left
        alpha2_left = alpha2_leftleft

        gap = gap * 0.5
        flag = 1
        
    end do

end subroutine Timmis_Estimates_Optimization



subroutine Timmis_Estimates_MLE(data_buyers, grad_la, grad_sf, alpha2, loss, alpha)
    
    USE Global_Data
    
    integer i, j, k
    
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB), intent(in)   :: data_buyers
    real(KIND=DOUBLE), DIMENSION(6), intent(in)               :: grad_la, grad_sf
    real(KIND=DOUBLE), intent(in)                             :: alpha2
    ! The MLE parameter results

    
    
    ! Alpha (Buyer side parameters, expect alpha2)
    real(KIND=DOUBLE), DIMENSION(10), intent(out)             :: alpha
    
    ! Market indicator, buyers characteritics, and their intersection
    real(KIND=DOUBLE), DIMENSION(numBUY, 10)                  :: X
    ! a copy of X, for regression call
    real(KIND=DOUBLE), DIMENSION(numBUY, 10)                  :: A
    ! P'(Z;beta) - alpha * Z
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: Y
    ! a copy of Y, for regression call
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: b
    ! error term
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: nu
    ! OLS std deviation
    real(KIND=DOUBLE)                                         :: sigma
    ! Jacobian
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: joc
    
    real(KIND=DOUBLE)                                         :: loss, likelihood
    real(KIND=DOUBLE), DIMENSION(10)                          :: alpha_max
    
    real(KIND=DOUBLE)                                         :: temp1, pi
    
    ! Parameter for calling regression function
    integer                                                      m, n, nrhs, lda, ldb, rank, lwork, info
    integer, DIMENSION(50)                                    :: jpvt
    real(KIND=DOUBLE), DIMENSION(1000000)                     :: work
    real(KIND=DOUBLE)                                         :: rcond
    
    PI = 4.D0*DATAN(1.D0)
    loss_min =  10000000

    
        X = 0
        Y = 0
        alpha = 0
        nu = 0
        loss = 0
        joc = 0
        A = 0
        b = 0
        ! ###########
        ! Construct X Matrix
        ! ###########
        
        X(:,1) = 1
        X(:,2) = data_buyers(:,_INCOME)
        do i = 1,numBUY
            if (int(data_buyers(i,_RACE)) == 2) then
                X(i,3) = 1
            else if (int(data_buyers(i,_RACE)) == 3) then
                X(i,4) = 1
            else if (int(data_buyers(i,_RACE)) == 4) then
                X(i,5) = 1
            end if
        end do
        !X(:,6) = data_buyers(:,_CITY)
        !X(:,7) = data_buyers(:,_INCOME) * data_buyers(:,_CITY)
        !do i = 1,numBUY
        !    if (int(data_buyers(i,_RACE)) == 2 .AND. int(data_buyers(i,_CITY)) == 1) then
        !        X(i,8) = 1
        !    else if (int(data_buyers(i,_RACE)) == 3 .AND. int(data_buyers(i,_CITY)) == 1) then
        !        X(i,9) = 1
        !    else if (int(data_buyers(i,_RACE)) == 4 .AND. int(data_buyers(i,_CITY)) == 1) then
        !        X(i,10) = 1
        !    end if
        !end do
        
        ! ###########
        ! Construct Y matrix
        ! ###########

        do i = 1,numBUY
            if (data_buyers(i,_city) == 1) then
                Y(i) = grad_la(1) + 2 * grad_la(2) * data_buyers(i,_violent_buyer) + 3 * grad_la(3) * data_buyers(i,_violent_buyer) &
                    + 4 * grad_la(4) * data_buyers(i,_violent_buyer)+ 5 * grad_la(5) * data_buyers(i,_violent_buyer)&
                    + 6 * grad_la(6) * data_buyers(i,_violent_buyer)
            else
                Y(i) = grad_sf(1) + 2 * grad_sf(2) * data_buyers(i,_violent_buyer) + 3 * grad_sf(3) * data_buyers(i,_violent_buyer) &
                    + 4 * grad_sf(4) * data_buyers(i,_violent_buyer)+ 5 * grad_sf(5) * data_buyers(i,_violent_buyer)&
                    + 6 * grad_sf(6) * data_buyers(i,_violent_buyer)
            end if
        end do
        
        do i = 1,numBUY
            Y(i) = Y(i) - alpha2 * data_buyers(i,_violent_buyer)
        end do
        

        
        ! ###########
        ! Do Regression: 
        ! ###########
        
        A = X
        b = Y
        m = numBUY
        n = 5
        nrhs = 1
        lda = numBUY
        ldb = numBUY
        lwork = 1000000
        !do i = 1,5
        !     print *, i, A(436847,i)
        !end do
        call dgelsy(m, n, nrhs, X, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
        

        do i = 1,5
            ! print *, b(i),i
            alpha(i) = b(i)
        end do
        
        ! ###########
        ! Calculate error term
        ! ###########
        
        X(:,1) = 1
        X(:,2) = data_buyers(:,_INCOME)
        do i = 1,numBUY
            if (int(data_buyers(i,_RACE)) == 2) then
                X(i,3) = 1
            else if (int(data_buyers(i,_RACE)) == 3) then
                X(i,4) = 1
            else if (int(data_buyers(i,_RACE)) == 4) then
                X(i,5) = 1
            end if
        end do
        !X(:,6) = data_buyers(:,_CITY)
        !X(:,7) = data_buyers(:,_INCOME) * data_buyers(:,_CITY)
        !do i = 1,numBUY
        !    if (int(data_buyers(i,_RACE)) == 2 .AND. int(data_buyers(i,_CITY)) == 1) then
        !        X(i,8) = 1
        !    else if (int(data_buyers(i,_RACE)) == 3 .AND. int(data_buyers(i,_CITY)) == 1) then
        !        X(i,9) = 1
        !    else if (int(data_buyers(i,_RACE)) == 4 .AND. int(data_buyers(i,_CITY)) == 1) then
        !        X(i,10) = 1
        !    end if
        !end do
        
        do i = 1,numBUY
            if (data_buyers(i,_city) == 1) then
                Y(i) = grad_la(1) + 2 * grad_la(2) * data_buyers(i,_violent_buyer) + 3 * grad_la(3) * data_buyers(i,_violent_buyer) &
                    + 4 * grad_la(4) * data_buyers(i,_violent_buyer)+ 5 * grad_la(5) * data_buyers(i,_violent_buyer)&
                    + 6 * grad_la(6) * data_buyers(i,_violent_buyer)
            else
                Y(i) = grad_sf(1) + 2 * grad_sf(2) * data_buyers(i,_violent_buyer) + 3 * grad_sf(3) * data_buyers(i,_violent_buyer) &
                    + 4 * grad_sf(4) * data_buyers(i,_violent_buyer)+ 5 * grad_sf(5) * data_buyers(i,_violent_buyer)&
                    + 6 * grad_sf(6) * data_buyers(i,_violent_buyer)
            end if
        end do
        
        do i = 1,numBUY
            Y(i) = Y(i) - alpha2 * data_buyers(i,_violent_buyer)
        end do
        
        temp1 = 0
        do i = 1,numBUY
            temp1 = 0
            do j = 1,10
                ! print *,alpha(j), X(i,j)
                temp1 = temp1 + alpha(j)*X(i,j)
            end do
            nu(i) = Y(i) - temp1
            ! print *,nu(i)
        end do
        
        sigma = 0
        do i = 1,numBUY
            ! print *, nu(i)
            sigma = sigma + nu(i)**2
        end do
        sigma = sigma / numBUY
        sigma = sqrt(sigma)
        ! print*, sigma
    
        
        ! #############
        ! MLE Function
        ! #############
        
        ! Jacobian
        do i = 1, numBUY
            if (data_buyers(i,_city) == 1) then
                joc(i) = 2 * grad_la(2) * data_buyers(i,_violent_buyer) + 6 * grad_la(3) * data_buyers(i,_violent_buyer) &
                    + 12 * grad_la(4) * data_buyers(i,_violent_buyer)+ 20 * grad_la(5) * data_buyers(i,_violent_buyer)&
                    + 30 * grad_la(6) * data_buyers(i,_violent_buyer) - alpha2
            else
                joc(i) = 2 * grad_sf(2) * data_buyers(i,_violent_buyer) + 6 * grad_sf(3) * data_buyers(i,_violent_buyer) &
                    + 12 * grad_sf(4) * data_buyers(i,_violent_buyer)+ 20 * grad_sf(5) * data_buyers(i,_violent_buyer)&
                    + 30 * grad_sf(6) * data_buyers(i,_violent_buyer) - alpha2
            end if
        end do
        
        likelihood = 0
        do i = 1, numBUY
            ! print *,  1/( sigma*sqrt(2*pi) ),exp( ( -1/(2*sigma*sigma) ) * nu(i) * nu(i) ),joc(i)
            likelihood = likelihood + ( 1/( sigma*sqrt(2*pi) ) * exp( ( -1/(2*sigma*sigma) ) * nu(i) * nu(i) ) * joc(i) )
        end do
        ! print*, likelihood
        loss = -likelihood
        ! print *, alpha2, loss
    

end subroutine Timmis_Estimates_MLE





subroutine Timmis_Estimates_Results(data_buyers, grad_la, grad_sf, para)

    USE Global_Data
    
    integer i, j, k
    
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB), intent(in)   :: data_buyers
    real(KIND=DOUBLE), DIMENSION(6), intent(in)               :: grad_la, grad_sf
    real(KIND=DOUBLE)                                         :: alpha2
    ! The MLE parameter results
    real(KIND=DOUBLE), DIMENSION(11), intent(out)              :: para
    
    ! Alpha (Buyer side parameters, expect alpha2)
    real(KIND=DOUBLE), DIMENSION(10)                          :: alpha

    
    ! Market indicator, buyers characteritics, and their intersection
    real(KIND=DOUBLE), DIMENSION(numBUY, 10)                  :: X
    ! a copy of X, for regression call
    real(KIND=DOUBLE), DIMENSION(numBUY, 10)                  :: A
    ! P'(Z;beta) - alpha * Z
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: Y
    ! a copy of Y, for regression call
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: b
    ! error term
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: nu
    ! OLS std deviation
    real(KIND=DOUBLE)                                         :: sigma
    ! Jacobian
    real(KIND=DOUBLE), DIMENSION(numBUY)                      :: joc
    
    real(KIND=DOUBLE)                                         :: loss, loss_max, alpha2_max
    real(KIND=DOUBLE), DIMENSION(10)                          :: alpha_max
    
    real(KIND=DOUBLE)                                         :: temp1, pi
    
    ! Parameter for calling regression function
    integer                                                      m, n, nrhs, lda, ldb, rank, lwork, info
    integer, DIMENSION(50)                                    :: jpvt
    real(KIND=DOUBLE), DIMENSION(1000000)                     :: work
    real(KIND=DOUBLE)                                         :: rcond
    
    PI = 4.D0*DATAN(1.D0)
    loss_min =  10000000
    

    
    ! #############
    ! Guess Alpha2
    ! #############
    
    do alpha2 = -100, 100, 1
        X = 0
        Y = 0
        alpha = 0
        nu = 0
        loss = 0
        joc = 0
        A = 0
        b = 0
        ! ###########
        ! Construct X Matrix
        ! ###########
        
        X(:,1) = 1
        X(:,2) = data_buyers(:,_INCOME)
        do i = 1,numBUY
            if (int(data_buyers(i,_RACE)) == 2) then
                X(i,3) = 1
            else if (int(data_buyers(i,_RACE)) == 3) then
                X(i,4) = 1
            else if (int(data_buyers(i,_RACE)) == 4) then
                X(i,5) = 1
            end if
        end do
        X(:,6) = data_buyers(:,_CITY)
        X(:,7) = data_buyers(:,_INCOME) * data_buyers(:,_CITY)
        do i = 1,numBUY
            if (int(data_buyers(i,_RACE)) == 2 .AND. int(data_buyers(i,_CITY)) == 1) then
                X(i,8) = 1
            else if (int(data_buyers(i,_RACE)) == 3 .AND. int(data_buyers(i,_CITY)) == 1) then
                X(i,9) = 1
            else if (int(data_buyers(i,_RACE)) == 4 .AND. int(data_buyers(i,_CITY)) == 1) then
                X(i,10) = 1
            end if
        end do
        
        ! ###########
        ! Construct Y matrix
        ! ###########

        do i = 1,numBUY
            if (data_buyers(i,_city) == 1) then
                Y(i) = grad_la(1) + 2 * grad_la(2) * data_buyers(i,_violent_buyer) + 3 * grad_la(3) * data_buyers(i,_violent_buyer) &
                    + 4 * grad_la(4) * data_buyers(i,_violent_buyer)+ 5 * grad_la(5) * data_buyers(i,_violent_buyer)&
                    + 6 * grad_la(6) * data_buyers(i,_violent_buyer)
            else
                Y(i) = grad_sf(1) + 2 * grad_sf(2) * data_buyers(i,_violent_buyer) + 3 * grad_sf(3) * data_buyers(i,_violent_buyer) &
                    + 4 * grad_sf(4) * data_buyers(i,_violent_buyer)+ 5 * grad_sf(5) * data_buyers(i,_violent_buyer)&
                    + 6 * grad_sf(6) * data_buyers(i,_violent_buyer)
            end if
        end do
        
        do i = 1,numBUY
            Y(i) = Y(i) - alpha2 * data_buyers(i,_violent_buyer)
        end do
        

        
        ! ###########
        ! Do Regression
        ! ###########
        
        A = X
        b = Y
        m = numBUY
        n = 10
        nrhs = 1
        lda = numBUY
        ldb = numBUY
        lwork = 1000000
        !do i = 1,10
        !     print *, i, A(436847,i)
        !end do
        call dgelsy(m, n, nrhs, X, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
        

        do i = 1,10
            ! print *, b(i),i
            alpha(i) = b(i)
        end do
        
        ! ###########
        ! Calculate error term
        ! ###########
        
        X(:,1) = 1
        X(:,2) = data_buyers(:,_INCOME)
        do i = 1,numBUY
            if (int(data_buyers(i,_RACE)) == 2) then
                X(i,3) = 1
            else if (int(data_buyers(i,_RACE)) == 3) then
                X(i,4) = 1
            else if (int(data_buyers(i,_RACE)) == 4) then
                X(i,5) = 1
            end if
        end do
        X(:,6) = data_buyers(:,_CITY)
        X(:,7) = data_buyers(:,_INCOME) * data_buyers(:,_CITY)
        do i = 1,numBUY
            if (int(data_buyers(i,_RACE)) == 2 .AND. int(data_buyers(i,_CITY)) == 1) then
                X(i,8) = 1
            else if (int(data_buyers(i,_RACE)) == 3 .AND. int(data_buyers(i,_CITY)) == 1) then
                X(i,9) = 1
            else if (int(data_buyers(i,_RACE)) == 4 .AND. int(data_buyers(i,_CITY)) == 1) then
                X(i,10) = 1
            end if
        end do
        
        do i = 1,numBUY
            if (data_buyers(i,_city) == 1) then
                Y(i) = grad_la(1) + 2 * grad_la(2) * data_buyers(i,_violent_buyer) + 3 * grad_la(3) * data_buyers(i,_violent_buyer) &
                    + 4 * grad_la(4) * data_buyers(i,_violent_buyer)+ 5 * grad_la(5) * data_buyers(i,_violent_buyer)&
                    + 6 * grad_la(6) * data_buyers(i,_violent_buyer)
            else
                Y(i) = grad_sf(1) + 2 * grad_sf(2) * data_buyers(i,_violent_buyer) + 3 * grad_sf(3) * data_buyers(i,_violent_buyer) &
                    + 4 * grad_sf(4) * data_buyers(i,_violent_buyer)+ 5 * grad_sf(5) * data_buyers(i,_violent_buyer)&
                    + 6 * grad_sf(6) * data_buyers(i,_violent_buyer)
            end if
        end do
        
        do i = 1,numBUY
            Y(i) = Y(i) - alpha2 * data_buyers(i,_violent_buyer)
        end do
        
        temp1 = 0
        do i = 1,numBUY
            temp1 = 0
            do j = 1,10
                ! print *,alpha(j), X(i,j)
                temp1 = temp1 + alpha(j)*X(i,j)
            end do
            nu(i) = Y(i) - temp1
            ! print *,nu(i)
        end do
        
        sigma = 0
        do i = 1,numBUY
            ! print *, nu(i)
            sigma = sigma + nu(i)**2
        end do
        sigma = sigma / numBUY
        sigma = sqrt(sigma)
        ! print*, sigma
    
        
        ! #############
        ! Loss Function
        ! #############
        
        ! Jacobian
        do i = 1, numBUY
            if (data_buyers(i,_city) == 1) then
                joc(i) = 2 * grad_la(2) * data_buyers(i,_violent_buyer) + 6 * grad_la(3) * data_buyers(i,_violent_buyer) &
                    + 12 * grad_la(4) * data_buyers(i,_violent_buyer)+ 20 * grad_la(5) * data_buyers(i,_violent_buyer)&
                    + 30 * grad_la(6) * data_buyers(i,_violent_buyer) - alpha2
            else
                joc(i) = 2 * grad_sf(2) * data_buyers(i,_violent_buyer) + 6 * grad_sf(3) * data_buyers(i,_violent_buyer) &
                    + 12 * grad_sf(4) * data_buyers(i,_violent_buyer)+ 20 * grad_sf(5) * data_buyers(i,_violent_buyer)&
                    + 30 * grad_sf(6) * data_buyers(i,_violent_buyer) - alpha2
            end if
        end do
        
        loss = 0
        do i = 1, numBUY
            ! print *,  1/( sigma*sqrt(2*pi) ),exp( ( -1/(2*sigma*sigma) ) * nu(i) * nu(i) ),joc(i)
            loss = loss + ( 1/( sigma*sqrt(2*pi) ) * exp( ( -1/(2*sigma*sigma) ) * nu(i) * nu(i) ) * joc(i) )
        end do
        ! print*, loss
        
        if (loss > loss_max) then
            loss_max = loss
            alpha2_max = alpha2
            alpha_max = alpha
        end if
        
        print *, alpha2, loss 
    end do
    
    
end subroutine Timmis_Estimates_Results


















subroutine Bootstrapped_Six_Order_Regression_Hedonic_Price_Function(numBS, data_la, data_sf)

    USE Global_Data
    USE Bootstrapped_Hedonic_Price_Function
    
    implicit none
    
    ! Variables
    
    integer i, j, k
    
    ! Data
    real(KIND=DOUBLE), DIMENSION(numLA, numV), intent(in)    :: data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV), intent(in)    :: data_sf
    integer, intent(in)                                      :: numBS
    
    
    real(KIND=DOUBLE), DIMENSION(numLA, numV)                :: boot_data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV)                :: boot_data_sf
    real(KIND=DOUBLE), DIMENSION(37)                         :: beta_la
    real(KIND=DOUBLE), DIMENSION(37)                         :: beta_sf
    
    
    ! File to save regression outputs
    open(unit = 712 , file = trim(dir_mom)//trim(filename_bootstrapped_six_order)//'_LA'//'.csv')
    open(unit = 713 , file = trim(dir_mom)//trim(filename_bootstrapped_six_order)//'_SF'//'.csv')
    
    ! print header

     write(712, 610) 'Constant,# Bathrooms,# Bedrooms,# Stories,Property Crime Rate,&
                    (Property Crime Rate)^2,Year Built,(Year Built)^2,Square Footage,(Square Footage)^2,&
                    # Total Rooms,(# Total Rooms)^2,Violent Crime Rate,(Violent Crime Rate)^2,&
                    (Violent Crime Rate)^3,(Violent Crime Rate)^4,(Violent Crime Rate)^5,(Violent Crime Rate)^6,&
                    1993,1994,1995,1996,1997,1998,2000,2001,2002,2003,2004,2005,2006,2007,2008,Orange,Riverside,San Bernadino,Ventura'
     
     write(713, 611) 'Constant,# Bathrooms,# Bedrooms,# Stories,Property Crime Rate,&
                    (Property Crime Rate)^2,Year Built,(Year Built)^2,Square Footage,(Square Footage)^2,&
                    # Total Rooms,(# Total Rooms)^2,Violent Crime Rate,(Violent Crime Rate)^2,&
                    (Violent Crime Rate)^3,(Violent Crime Rate)^4,(Violent Crime Rate)^5,(Violent Crime Rate)^6,&
                    1993,1994,1995,1996,1997,1998,2000,2001,2002,2003,2004,2005,2006,2007,2008,Contra Costa,San Francisco,San Mateo,Santa Clara'
        
    do i = 1, numBS
        
        boot_data_la = 0
        boot_data_sf = 0
        beta_la = 0
        beta_sf = 0
        ! #######
        ! Acquire the bootstrapped functions
        ! #######
        
        call Bootstrapped(data_la, data_sf, boot_data_la, boot_data_sf)
        
        ! #######
        ! Do linear regression on the bootstrapped dataset
        ! #######
        
        call Regression_Six_Order_Hedonic_Price_Function(boot_data_la, boot_data_sf, beta_la, beta_sf)
        
        ! #######
        ! Write parameters to .csv
        ! #######
        do j = 1,36
            write(712, 612, advance = 'no') beta_la(j),','
        end do
        write(712, 613) beta_la(37)
        
        
        do j = 1,36
            write(713, 612, advance = 'no') beta_sf(j),','
        end do
        write(713, 613) beta_sf(37)
        
        write(*,610, advance = 'no') 'interation: '
        write(*,*) i
    end do
610        format(a509)
611        format(a529)
612        format(f22.8,a1)
613        format(f22.8) 
614        format(a1)           

           
end subroutine Bootstrapped_Six_Order_Regression_Hedonic_Price_Function


    
subroutine Regression_Six_Order_Hedonic_Price_Function(data_la, data_sf, beta_la, beta_sf)

    USE Global_Data
    USE Basic_Algebra
    USE Bootstrapped_Hedonic_Price_Function
    
    implicit none
    
    ! Variables
    
    integer i, j, k
    
    ! Data
    ! In most situations, the inputs are bootstrapped datasets
    real(KIND=DOUBLE), DIMENSION(numLA, numV), intent(in)    :: data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV), intent(in)    :: data_sf
    
    real(KIND=DOUBLE), DIMENSION(37), intent(out)   :: beta_la
    real(KIND=DOUBLE), DIMENSION(37), intent(out)   :: beta_sf
    
    ! Matrix A and b
    real(KIND=DOUBLE), DIMENSION(numLA, 37)                  :: A_la
    real(KIND=DOUBLE), DIMENSION(numSF, 37)                  :: A_sf
    real(KIND=DOUBLE), DIMENSION(numLA)                      :: b_la
    real(KIND=DOUBLE), DIMENSION(numSF)                      :: b_sf
    
    ! Parameter for calling regression function
    integer                                                     m, n, nrhs, lda, ldb, rank, lwork, info
    integer, DIMENSION(50)                                   :: jpvt
    real(KIND=DOUBLE), DIMENSION(1000000)                    :: work
    real(KIND=DOUBLE)                                        :: rcond

    
      
    ! #################
    ! Generate matrix A
    ! #################
    A_la = 0
    A_sf = 0
    b_la = 0
    b_sf = 0
    
    ! ######
    ! LA
    ! ######
    
    A_la(:,1) = 1.
    A_la(:,2) = data_la(:,_BATHROOMS)
    A_la(:,3) = data_la(:,_BEDROOMS)
    A_la(:,4) = data_la(:,_STORIES)
    A_la(:,5) = data_la(:,_PROPERTY)
    A_la(:,6) = data_la(:,_PROPERTY)*data_la(:,_PROPERTY)
    A_la(:,7) = data_la(:,_BUILTYEAR)
    A_la(:,8) = data_la(:,_BUILTYEAR)*data_la(:,_BUILTYEAR)
    A_la(:,9) = data_la(:,_SQUAREFEET)
    A_la(:,10) = data_la(:,_SQUAREFEET)*data_la(:,_SQUAREFEET)
    A_la(:,11) = data_la(:,_ROOMS)
    A_la(:,12) = data_la(:,_ROOMS)*data_la(:,_ROOMS)
    A_la(:,13) = data_la(:,_VIOLENT)
    A_la(:,14) = data_la(:,_VIOLENT)**2
    A_la(:,15) = data_la(:,_VIOLENT)**3
    A_la(:,16) = data_la(:,_VIOLENT)**4
    A_la(:,17) = data_la(:,_VIOLENT)**5
    A_la(:,18) = data_la(:,_VIOLENT)**6
    
    
    ! Vector of Year Dummies (omit 1999)
    ! For example, 1993 -> A(,19)=1, 2006 -> A(,27)=1
    ! print *, data_la(1,_SALEYEAR)
    do i = 1,numLA
        if (data_la(i,_SALEYEAR) < 1998.5) then
            A_la( i, int(data_la(i,_SALEYEAR)-1974) ) = 1
        else if (data_la(i,_SALEYEAR) > 1999.5) then
            A_la( i, int(data_la(i,_SALEYEAR)-1975) ) = 1
        end if
    end do
    
    ! Vector of Dummies for Certain Counties
    ! For Los Angeles, include explicit dummies for counties 59, 65, and 71. For San Francisco, include
    ! dummies for counties 13, 75, 81, and 85.
    ! For example, 59 -> A(,28)=1

    
    do i = 1, numLA
        ! print *, data_la(i,_COUNTYID)
        if ( data_la(i,_COUNTYID) == 59 ) then
            A_la(i, 34) = 1
        else if ( data_la(i,_COUNTYID) == 65 ) then
            A_la(i, 35) = 1
        else if ( data_la(i,_COUNTYID) == 71 ) then
            A_la(i, 36) = 1
        else if ( data_la(i,_COUNTYID) == 111 ) then
            A_la(i, 37) = 1
        end if
    end do
    
    ! ######
    ! SF
    ! ######
    
    A_sf(:,1) = 1.
    A_sf(:,2) = data_sf(:,_BATHROOMS)
    A_sf(:,3) = data_sf(:,_BEDROOMS)
    A_sf(:,4) = data_sf(:,_STORIES)
    A_sf(:,5) = data_sf(:,_PROPERTY)
    A_sf(:,6) = data_sf(:,_PROPERTY)*data_sf(:,_PROPERTY)
    A_sf(:,7) = data_sf(:,_BUILTYEAR)
    A_sf(:,8) = data_sf(:,_BUILTYEAR)*data_sf(:,_BUILTYEAR)
    A_sf(:,9) = data_sf(:,_SQUAREFEET)
    A_sf(:,10) = data_sf(:,_SQUAREFEET)*data_sf(:,_SQUAREFEET)
    A_sf(:,11) = data_sf(:,_ROOMS)
    A_sf(:,12) = data_sf(:,_ROOMS)*data_sf(:,_ROOMS)
    A_sf(:,13) = data_sf(:,_VIOLENT)
    A_sf(:,14) = data_sf(:,_VIOLENT)*data_sf(:,_VIOLENT)
    A_sf(:,15) = data_sf(:,_VIOLENT)**3
    A_sf(:,16) = data_sf(:,_VIOLENT)**4
    A_sf(:,17) = data_sf(:,_VIOLENT)**5
    A_sf(:,18) = data_sf(:,_VIOLENT)**6
    
    ! Vector of Year Dummies (omit 1999)
    ! For example, 1993 -> A(,19)=1, 2008 -> A(,33)=1
    do i = 1,numSF
        if (data_sf(i,_SALEYEAR) < 1998.5) then
            A_sf( i, int(data_sf(i,_SALEYEAR)-1974) ) = 1
        else if (data_sf(i,_SALEYEAR) > 1999.5) then
            A_sf( i, int(data_sf(i,_SALEYEAR)-1975) ) = 1
        end if
    end do
    
    ! Vector of Dummies for Certain Counties
    ! For Los Angels, include explicit dummies for counties 59, 65, and 71. For San Francisco, include
    ! dummies for counties 13, 75, 81, and 85.
    ! For example, 59 -> A_la(,28)=1
    
    do i = 1, numSF
        if ( data_sf(i,_COUNTYID) == 13 ) then
            A_sf(i, 34) = 1
        else if ( data_sf(i,_COUNTYID) == 75 ) then
            A_sf(i, 35) = 1
        else if ( data_sf(i,_COUNTYID) == 81 ) then
            A_sf(i, 36) = 1
        else if ( data_sf(i,_COUNTYID) == 85 ) then
            A_sf(i, 37) = 1
        end if
    end do

    ! #################
    ! Generate vector b
    ! #################
    
    ! LA
    
    b_la = data_la(:,_PRICE)
    
    ! SF
    
    b_sf = data_sf(:,_PRICE)
    
    !Print to check
    !do i = 1,33
    !    print *, A_sf(5000, i)
    !end do
    !print *, b_sf(5000)
    
    ! ##########
    ! Write to bootstrapped_data to .csv file
    ! ##########
    
    ! LA
    
    !open(unit = 101 , file = trim(dir_mom)//'bootstrapped_data_la'//'.csv')
    !do i = 1,numLA
    !    do j = 1,31
    !        write(101,600,advance = 'no') A_la(i, j), ','
    !    end do
    !    write(101,601) b_la(i)
    !end do
    !600 format(f22.8,a1)
    !601 format(f22.8)
    
    ! %%%%%%%%%%%%%%%%%%
    ! Do Regression
    ! %%%%%%%%%%%%%%%%%%
    
    ! LA
    m = numLA
    n = 37
    nrhs = 1
    lda = numLA
    ldb = numLA
    lwork = 1000000
    call dgelsy(m, n, nrhs, A_la, lda, b_la, ldb, jpvt, rcond, rank, work, lwork, info)
    
    do i = 1,37
        beta_la(i) = b_la(i)
    end do
    
    ! SF
    
    m = numSF
    n = 37
    nrhs = 1
    lda = numSF
    ldb = numSF
    lwork = 1000000
    call dgelsy(m, n, nrhs, A_sf, lda, b_sf, ldb, jpvt, rcond, rank, work, lwork, info)
    
    do i = 1,37
        beta_sf(i) = b_sf(i)
    end do

end subroutine Regression_Six_Order_Hedonic_Price_Function

end Module Bishop_Timmins_Estimator