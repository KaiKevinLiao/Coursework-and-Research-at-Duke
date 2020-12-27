! Global variables
Module Global_Data

    ! Directory with data
    CHARACTER(LEN=120), parameter :: dir_mom = 'C:\Users\kaike\OneDrive\Desktop\Non_Market_Valuation\Problem Set -- Hedonics\'
    ! Task 1 Summary statistics output filename
    CHARACTER(LEN=120), parameter :: filename_summary = 'Summary_Statistics'
    ! Task 2 Bootstrapped hedonic price function regression outputs filename
    CHARACTER(LEN=120), parameter :: filename_bootstrap = 'Bootstrapped_Regression_Outputs'
    ! Task 2 Hedonic price function regression results, including variance
    CHARACTER(LEN=120), parameter :: filename_hedonic_price_function_regression = 'Hedonic_Price_Function_Regression_Results'
    ! Task 3 Bootstrapped Multimarket Rosen Estimates
    CHARACTER(LEN=120), parameter :: filename_bootstrapped_rosen_estimate = 'Bootstrapped_Rosen_Outputs'
    ! Task 3 Multimarket Rosen Estimates
    CHARACTER(LEN=120), parameter :: filename_rosen_estimate = 'Rosen_Estimate_Regression_Results'
    ! Task 4 Non-paramatric hedonic price gradient
    CHARACTER(LEN=120), parameter :: filename_non_paramatric_hedonic_price_gradient = 'Non_Paramatric_Hedonic_Price_Gradient'
    ! Task 5 Non-paramatric MWTP Estimates
    CHARACTER(LEN=120), parameter :: filename_non_paramatric_MWTP = 'Non_Paramatric_MWTP_results'
    ! Task 6 Bishop-Timmis Bootstrapped Hedoic Price Function
    CHARACTER(LEN=120), parameter :: filename_bootstrapped_six_order = 'Bootstrapped_Six_Order_Hedonic_Price_Function'
    ! Task 6 Bishop-Timmis Bootstrapped Estimates
    CHARACTER(LEN=120), parameter :: filename_bootstrapped_timmis_estimate = 'Bootstrapped_Bishop_Timmis_Outputs'
    ! Task 6 Bishop-Timmis Estimates
    CHARACTER(LEN=120), parameter :: filename_timmis_estimate = 'Bishop_Timmis_Outputs'
    
    integer, parameter :: DOUBLE = SELECTED_REAL_KIND(p=8)
    
    ! Number of observation in LA
    integer, parameter :: numLA = 386063
    ! Number of observation in SF
    integer, parameter :: numSF = 378252
    ! Number of data variables
    integer, parameter :: numV = 12
    ! Number of counties
    integer, parameter :: numC = 5
    ! Number of buyers observations
    integer, parameter :: numBUY = 659548
    ! Number of buyers variables
    integer, parameter :: numVB  = 7
    ! Number of LA buyers observations
    integer, parameter :: numBUYLA = 55498
    ! Number of LA buyers variables
    integer, parameter :: numVBLA = 6
    ! Maximun of VC rate
    integer, parameter :: maxVC = 2000
    integer, dimension(numC), parameter :: County_iden_la = (/37, 59, 65, 71, 111/)
    integer, dimension(numC), parameter :: County_iden_sf = (/1, 13, 75, 81, 85/)
    
    
    real(KIND=DOUBLE), parameter :: alpha2_best = -25, loss_best = 100
    real(KIND=DOUBLE), dimension(10), parameter :: alpha_best = 10000
    integer, parameter :: INIT_GAP = 5

end Module Global_Data