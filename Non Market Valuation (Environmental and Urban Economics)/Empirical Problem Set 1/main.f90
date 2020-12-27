!  main.f90 
!
!  FUNCTIONS:
!  main - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Non_Market_Valuation_Task_1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    
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
program Main

    ! Module defining and setting global variables
    USE Global_Data
    USE Summary_Statistics
    USE Bootstrapped_Hedonic_Price_Function
    USE Basic_Algebra
    USE Bootstrapped_Multimarket_Rosen_Estimates
    USE Non_Parametric_Hedonic_Price_Function
    USE Non_Parametric_MWTP_Estimates
    USE Bishop_Timmins_Estimator

    implicit none
    

    ! Variables
    

    
    integer i, j, numBS

    
    
    
    ! Data
    ! Intergated Data
    real(KIND=DOUBLE), DIMENSION(numLA, numV)   :: data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV)   :: data_sf
    
    real(KIND=DOUBLE), DIMENSION(numBUY, numVB) :: data_buyers
    
    real(KIND=DOUBLE), DIMENSION(numBUYLA, numVBLA) :: data_buyers_la
    
    ! Bootstrapped output
    real(KIND=DOUBLE), DIMENSION(numLA, numV)   :: boot_data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV)   :: boot_data_sf

    
    ! Seperated Data
    ! House ID
    integer           houseID_la(numLA), houseID_sf(numSF)
    ! Price (deflated to year 2000 dollars)
    real(KIND=DOUBLE) price_la(numLA), price_sf(numSF)
    ! County ID # (see below)
    integer           countyID_la(numLA), countyID_sf(numSF)
    ! Year Built
    integer           builtyear_la(numLA), builtyear_sf(numSF)
    ! Square Footage
    real(KIND=DOUBLE) square_la(numLA), square_sf(numSF)
    ! # Bathrooms # Bedrooms # Total Rooms # Stories
    integer           bathrooms_la(numLA), bathrooms_sf(numSF), bedrooms_la(numLA), bedrooms_sf(numSF), &
                      rooms_la(numLA), rooms_sf(numSF), stories_la(numLA), stories_sf(numSF)
    ! Violent Crime Rate (Cases per 100,000)
    real(KIND=DOUBLE) violentcrime_la(numLA), violentcrime_sf(numSF)
    ! Property Crime Rate (Cases Per 100,000)
    real(KIND=DOUBLE) propertycrime_la(numLA), propertycrime_sf(numSF)
    ! Year of Sale (1993 – 2008)
    integer           saleyear_la(numLA), saleyear_sf(numSF)
    
    
    real(KIND=DOUBLE), DIMENSION(numSF, numV)   :: beta_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV)   :: beta_sf
    
    
    ! Body of Non_Market_Valuation_Task_1
    
    
    Call Read_Data_Seperated
    Call Read_Data_Combined
    Call Read_Buyer_Data
    Call Read_Buyer_Data_LA
    Call Compute_Summary_Statistics(data_la, data_sf)
    ! Call Bootstrapped(data_la, data_sf, boot_data_la, boot_data_sf)
    ! Call Regression_Hedonic_Price_Function(data_la, data_sf, beta_la, beta_sf)
    
    numBS = 500
    Call Bootstrapped_Regression_Hedonic_Price_Function(numBS, data_la, data_sf)
    Call Hedonic_Price_Function_Regression_Results(numBS, data_la, data_sf)
    ! print *, data_la(15,_SALEYEAR), saleyear_la(15)
    ! Call Bootstrapped_Rosen_Estimates_Results(data_buyers, numBS)
    Call Non_Parametric_Hedonic_Price_Function_Results(data_la)
    Call Non_Parametric_MWTP_Estimates_Results(data_buyers_la)
    Call Bootstrapped_Six_Order_Regression_Hedonic_Price_Function(numBS, data_la, data_sf)
    Call Bootstrapped_MLE(data_buyers, numBS)
    read(*,*), i

    Contains
    
    
    
Subroutine Read_Data_Seperated

    open(unit=1,file=trim(dir_mom)//'la_data.txt', status='old', action='read')
    do i = 1, numLA
        read(1, *) houseID_la(i), price_la(i), countyID_la(i), &
            builtyear_la(i), square_la(i), bathrooms_la(i), bedrooms_la(i), &
            rooms_la(i), stories_la(i), violentcrime_la(i), propertycrime_la(i), &
            saleyear_la(i)
    end do
    close(1)
    
    open(unit=2,file=trim(dir_mom)//'sf_data.txt', status='old', action='read')
    do i = 1, numSF
        read(2, *) houseID_sf(i), price_sf(i), countyID_sf(i), &
            builtyear_sf(i), square_sf(i), bathrooms_sf(i), bedrooms_sf(i), &
            rooms_sf(i), stories_sf(i), violentcrime_sf(i), propertycrime_sf(i), &
            saleyear_sf(i)
    end do
    close(1)
    
end subroutine Read_Data_Seperated

Subroutine Read_Data_Combined

    open(unit=3,file=trim(dir_mom)//'la_data.txt', status='old', action='read')
    do i = 1, numLA
        read(3, *) data_la(i,:)
    end do
    close(3)
    
    open(unit=4,file=trim(dir_mom)//'sf_data.txt', status='old', action='read')
    do i = 1, numSF
        read(4, *) data_sf(i,:)
    end do
    close(4)


end subroutine Read_Data_Combined

subroutine Read_Buyer_Data
    
    open(unit=5,file=trim(dir_mom)//'buyer_data_sf_la.txt', status='old', action='read')
    do i = 1,numBUY
        read(5, *) data_buyers(i,:)
    end do
    close(5)

end subroutine Read_Buyer_Data

subroutine Read_Buyer_Data_LA
    
    open(unit=6,file=trim(dir_mom)//'buyer_data_la.txt', status='old', action='read')
    do i = 1, numBUYLA
        read(6, *) data_buyers_la(i,:)
    end do

end subroutine Read_Buyer_Data_LA




end program main

    
    
    
    
function Transfer_CountyID(index) result(county)
    INTEGER :: index
    CHARACTER(LEN=120) :: county
    integer i
    
    if (index == 37) then 
        county = 'Los Angeles'
    else if (index == 59) then
        county = 'Orange'
    else if (index == 65) then
        county = 'Riverside'
    else if (index == 71) then
        county = 'San Bernadino'
    else if (index == 111) then
        county = 'Ventura'
    else if (index == 1) then
        county = 'Alameda'
    else if (index == 13) then
        county = 'Contra Costa'
    else if (index == 75) then
        county = 'San Francisco'
    else if (index == 81) then
        county = 'San Mateo'
    else if (index == 85) then
        county = 'Santa Clara'
    end if
end function