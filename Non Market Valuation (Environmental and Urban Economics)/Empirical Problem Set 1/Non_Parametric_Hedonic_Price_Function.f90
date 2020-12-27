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
    
    
Module Non_Parametric_Hedonic_Price_Function

    Contains
    
subroutine Non_Parametric_Hedonic_Price_Function_Results(data_la)
    
    USE Global_Data
    USE Basic_Algebra
    
    integer chi, i, j, t

    ! Data
    real(KIND=DOUBLE), DIMENSION(numLA, numV), intent(in)         :: data_la
    
    ! K martrix
    real(KIND=DOUBLE), DIMENSION(numLA, 100)                      :: K
    ! the bandwidth
    integer, DIMENSION(4)                                         :: h
    ! Standard Deviation of violence crime
    real(KIND=DOUBLE)                                             :: sigma_VC
    ! PI
    real(KIND=DOUBLE)                                             :: PI
    
    
    
    ! Working matrix for calculation
    real(KIND=DOUBLE), DIMENSION(numLA)                           :: data_VC_la
    real(KIND=DOUBLE), DIMENSION(numLA)                           :: data_price_la
    real(KIND=DOUBLE), DIMENSION(2,2)                             :: A
    real(KIND=DOUBLE), DIMENSION(2,1)                             :: B
    real(KIND=DOUBLE), DIMENSION(2,1)                             :: output_temp
    real(KIND=DOUBLE), DIMENSION(2,2)                             :: A_inv
    
    ! Outputs
    real(KIND=DOUBLE), DIMENSION(4,maxVC)                         :: alpha
    real(KIND=DOUBLE), DIMENSION(4,maxVC)                         :: beta
    
    ! initialize
    
    h = (/1, 3, 10, 1000/)
    PI = 4.D0*DATAN(1.D0)
    print *, PI
    data_VC_la(:) = data_la(:,_VIOLENT)
    data_price_la(:) = data_la(:,_PRICE)
    
    ! Compute the Standard Deviation of violence crime sigma_VC
    
    sigma_VC = sqrt( variance(data_VC_la,numLA,0) )
    print *, sigma_VC
    

    
    do i = 1, 4
        do t = 1, 20
    ! #########
    ! Construct K matrix
    ! #########
        
        
        K = 0
        do chi = 1, 100
            do j = 1, numLA
            
                K(j,chi) = 1/( h(i)*sigma_VC ) * ( 1/sqrt(2*PI) ) * exp( (-0.5) * ( (data_VC_la(j)-((t-1)*100 + chi))/(h(i)*sigma_VC) )**2 )
            
            end do
        end do
        
    ! #########
    ! Compute the maximazation results
    ! #########
        do chi = 1, 100
            A = 0
            B = 0
            do j = 1, numLA
                
                A(1,1) = A(1,1) + data_VC_la(j)**2 * K(j,chi)
                A(1,2) = A(1,2) + data_VC_la(j) * K(j,chi)
                A(2,1) = A(2,1) + data_VC_la(j) * K(j,chi)
                A(2,2) = A(2,2) + K(j,chi)
                B(1,1) = B(1,1) + data_VC_la(j) * K(j,chi) * data_price_la(j)
                B(2,1) = B(2,1) + K(j,chi) * data_price_la(j)
     
            end do 
            output_temp = 0
            A_inv = inverse(A,2)
            ! print *, A(1,1)
            output_temp = matrix_mult( A_inv, B, 2, 2, 1 )

            beta(i,(t-1)*100 + chi) = output_temp(1,1)
            alpha(i,(t-1)*100 + chi) = output_temp(2,1)
        end do
        

        print *, 'step 2: h=', h(i)
        print *, 'step 2: t=', t
        end do
        
    end do
    
    ! ############
    ! write results to csv
    ! ############
    
    ! write alpha
    
    open(unit = 160 , file = trim(dir_mom)//trim(filename_non_paramatric_hedonic_price_gradient)//'.csv')
    
    write(160,950) 'chi, h=1,,h=3,,h=10,,h=1000'
    write(160,951) ',alpha,beta,alpha,beta,alpha,beta,alpha,beta,'
    
    do chi = 1,maxVC
        write(160,952, advance = 'no') chi
        write(160,953) ',',alpha(1,chi),',',beta(1,chi),',',alpha(2,chi),',',beta(2,chi),',',alpha(3,chi),',',beta(3,chi),',',alpha(4,chi),',',beta(4,chi)
    end do
    
950 format(a27)
951 format(a45)
952 format(I5)
953 format(8(a1,f50.8))
    
end subroutine Non_Parametric_Hedonic_Price_Function_Results
    
end Module Non_Parametric_Hedonic_Price_Function