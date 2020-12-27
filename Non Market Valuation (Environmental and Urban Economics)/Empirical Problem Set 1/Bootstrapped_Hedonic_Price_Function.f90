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


Module Bootstrapped_Hedonic_Price_Function

    Contains
    
subroutine Hedonic_Price_Function_Regression_Results(numBS, data_la, data_sf)
    
    USE Global_Data
    USE Basic_Algebra
    
    implicit none
    
    ! Variables
    
    integer i, j, k
    
    ! Data
    real(KIND=DOUBLE), DIMENSION(numLA, numV), intent(in)    :: data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV), intent(in)    :: data_sf
    integer, intent(in)                                      :: numBS
    
    real(KIND=DOUBLE), DIMENSION(33)                         :: beta_la
    real(KIND=DOUBLE), DIMENSION(33)                         :: beta_sf
    
    real(KIND=DOUBLE), DIMENSION(33)                         :: variance_la
    real(KIND=DOUBLE), DIMENSION(33)                         :: variance_sf
    
    real(KIND=DOUBLE), DIMENSION(33)                         :: stddev_la
    real(KIND=DOUBLE), DIMENSION(33)                         :: stddev_sf
    
    real(KIND=DOUBLE), DIMENSION(numBS, 33)                  :: boots_la
    real(KIND=DOUBLE), DIMENSION(numBS, 33)                  :: boots_sf
    
    real(KIND=DOUBLE), DIMENSION(numBS)                      :: A
    

    
    ! File to save regression outputs
    open(unit = 120 , file = trim(dir_mom)//trim(filename_hedonic_price_function_regression)//'_LA'//'.csv')
    open(unit = 121 , file = trim(dir_mom)//trim(filename_hedonic_price_function_regression)//'_SF'//'.csv')
    
    ! print header

     write(120, 610) 'Constant,# Bathrooms,# Bedrooms,# Stories,Property Crime Rate,&
                    (Property Crime Rate)^2,Year Built,(Year Built)^2,Square Footage,(Square Footage)^2,&
                    # Total Rooms,(# Total Rooms)^2,Violent Crime Rate,(Violent Crime Rate)^2,&
                    1993,1994,1995,1996,1997,1998,2000,2001,2002,2003,2004,2005,2006,2007,2008,Orange,Riverside,San Bernadino,Ventura'
     
     write(121, 611) 'Constant,# Bathrooms,# Bedrooms,# Stories,Property Crime Rate,&
                    (Property Crime Rate)^2,Year Built,(Year Built)^2,Square Footage,(Square Footage)^2,&
                    # Total Rooms,(# Total Rooms)^2,Violent Crime Rate,(Violent Crime Rate)^2,&
                    1993,1994,1995,1996,1997,1998,2000,2001,2002,2003,2004,2005,2006,2007,2008,Contra Costa,San Francisco,San Mateo,Santa Clara'
        

        
    ! #######
    ! Compute regression results of original datasets
    ! #######
        
    call Regression_Hedonic_Price_Function(data_la, data_sf, beta_la, beta_sf)
    
    ! #######
    ! Write regression results to .csv
    ! #######
    do j = 1,32
        write(120, 612, advance = 'no') beta_la(j),','
    end do
    write(120, 613) beta_la(33)
        
        
    do j = 1,32
        write(121, 612, advance = 'no') beta_sf(j),','
    end do
    write(121, 613) beta_sf(33)
    
    
    ! #########
    ! Compute parameters variances
    ! #########
    
    ! Read Bootstrapped_Regression_Outputs
    open(unit = 110 , file = trim(dir_mom)//trim(filename_bootstrap)//'_LA'//'.csv')
    open(unit = 111 , file = trim(dir_mom)//trim(filename_bootstrap)//'_SF'//'.csv')
    
    read(110,*)
    do i = 1,numBS
        read(110,*) boots_la(i,:)
    end do
    read(111,*)
    do i = 1,numBS
        read(111,*) boots_sf(i,:)
    end do
    
    ! Calculate the variance
    do i = 1, 33
        A = boots_la(:,i)
        variance_la(i) = variance(A,numBS,0)
    end do
    
    do i = 1, 33
        A = boots_sf(:,i)
        variance_sf(i) = variance(A,numBS,0)
    end do
    
    ! #########
    ! Write parameters variances
    ! #########

    do j = 1,32
        write(120, 612, advance = 'no') variance_la(j),','
    end do
    write(120, 613) beta_la(33)
        
        
    do j = 1,32
        write(121, 612, advance = 'no') variance_sf(j),','
    end do
    write(121, 613) beta_sf(33)
    
    
    ! #########
    ! Compute and write parameters standard deviations
    ! #########
    
    stddev_la(:) = sqrt(variance_la(:))
    stddev_sf(:) = sqrt(variance_sf(:))
    
    do j = 1,32
        write(120, 612, advance = 'no') stddev_la(j),','
    end do
    write(120, 613) beta_la(33)
        
        
    do j = 1,32
        write(121, 612, advance = 'no') stddev_sf(j),','
    end do
    write(121, 613) beta_sf(33)
    

610        format(a333)
611        format(a343)
612        format(f50.8,a1)
613        format(f50.8) 
614        format(a1) 
    
end subroutine Hedonic_Price_Function_Regression_Results
    
subroutine Bootstrapped_Regression_Hedonic_Price_Function(numBS, data_la, data_sf)

    USE Global_Data
    
    implicit none
    
    ! Variables
    
    integer i, j, k
    
    ! Data
    real(KIND=DOUBLE), DIMENSION(numLA, numV), intent(in)    :: data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV), intent(in)    :: data_sf
    integer, intent(in)                                      :: numBS
    
    
    real(KIND=DOUBLE), DIMENSION(numLA, numV)                :: boot_data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV)                :: boot_data_sf
    real(KIND=DOUBLE), DIMENSION(33)                         :: beta_la
    real(KIND=DOUBLE), DIMENSION(33)                         :: beta_sf
    
    
    ! File to save regression outputs
    open(unit = 102 , file = trim(dir_mom)//trim(filename_bootstrap)//'_LA'//'.csv')
    open(unit = 103 , file = trim(dir_mom)//trim(filename_bootstrap)//'_SF'//'.csv')
    
    ! print header

     write(102, 610) 'Constant,# Bathrooms,# Bedrooms,# Stories,Property Crime Rate,&
                    (Property Crime Rate)^2,Year Built,(Year Built)^2,Square Footage,(Square Footage)^2,&
                    # Total Rooms,(# Total Rooms)^2,Violent Crime Rate,(Violent Crime Rate)^2,&
                    1993,1994,1995,1996,1997,1998,2000,2001,2002,2003,2004,2005,2006,2007,2008,Orange,Riverside,San Bernadino'
     
     write(103, 611) 'Constant,# Bathrooms,# Bedrooms,# Stories,Property Crime Rate,&
                    (Property Crime Rate)^2,Year Built,(Year Built)^2,Square Footage,(Square Footage)^2,&
                    # Total Rooms,(# Total Rooms)^2,Violent Crime Rate,(Violent Crime Rate)^2,&
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
        
        call Regression_Hedonic_Price_Function(boot_data_la, boot_data_sf, beta_la, beta_sf)
        
        ! #######
        ! Write parameters to .csv
        ! #######
        do j = 1,32
            write(102, 612, advance = 'no') beta_la(j),','
        end do
        write(102, 613) beta_la(33)
        
        
        do j = 1,32
            write(103, 612, advance = 'no') beta_sf(j),','
        end do
        write(103, 613) beta_sf(33)
        
        write(*,610, advance = 'no') 'interation: '
        write(*,*) i
    end do
610        format(a315)
611        format(a343)
612        format(f22.8,a1)
613        format(f22.8) 
614        format(a1)           

           
end subroutine Bootstrapped_Regression_Hedonic_Price_Function
    




subroutine Regression_Hedonic_Price_Function(data_la, data_sf, beta_la, beta_sf)
    
    USE Global_Data
    
    implicit none
    
    ! Variables
    
    integer i, j, k
    
    ! Data
    ! In most situations, the inputs are bootstrapped datasets
    real(KIND=DOUBLE), DIMENSION(numLA, numV), intent(in)    :: data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV), intent(in)    :: data_sf
    
    real(KIND=DOUBLE), DIMENSION(33), intent(out)   :: beta_la
    real(KIND=DOUBLE), DIMENSION(33), intent(out)   :: beta_sf
    
    ! Matrix A and b
    real(KIND=DOUBLE), DIMENSION(numLA, 33)                  :: A_la
    real(KIND=DOUBLE), DIMENSION(numSF, 33)                  :: A_sf
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
    A_la(:,14) = data_la(:,_VIOLENT)*data_la(:,_VIOLENT)
    
    ! Vector of Year Dummies (omit 1999)
    ! For example, 1993 -> A(,15)=1, 2006 -> A(,27)=1
    ! print *, data_la(1,_SALEYEAR)
    do i = 1,numLA
        if (data_la(i,_SALEYEAR) < 1998.5) then
            A_la( i, int(data_la(i,_SALEYEAR)-1978) ) = 1
        else if (data_la(i,_SALEYEAR) > 1999.5) then
            A_la( i, int(data_la(i,_SALEYEAR)-1979) ) = 1
        end if
    end do
    
    ! Vector of Dummies for Certain Counties
    ! For Los Angeles, include explicit dummies for counties 59, 65, and 71. For San Francisco, include
    ! dummies for counties 13, 75, 81, and 85.
    ! For example, 59 -> A(,28)=1

    
    do i = 1, numLA
        ! print *, data_la(i,_COUNTYID)
        if ( data_la(i,_COUNTYID) == 59 ) then
            A_la(i, 30) = 1
        else if ( data_la(i,_COUNTYID) == 65 ) then
            A_la(i, 31) = 1
        else if ( data_la(i,_COUNTYID) == 71 ) then
            A_la(i, 32) = 1
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
    
    ! Vector of Year Dummies (omit 1999)
    ! For example, 1993 -> A(,15)=1, 2008 -> A(,29)=1
    do i = 1,numSF
        if (data_sf(i,_SALEYEAR) < 1998.5) then
            A_sf( i, int(data_sf(i,_SALEYEAR)-1978) ) = 1
        else if (data_sf(i,_SALEYEAR) > 1999.5) then
            A_sf( i, int(data_sf(i,_SALEYEAR)-1979) ) = 1
        end if
    end do
    
    ! Vector of Dummies for Certain Counties
    ! For Los Angels, include explicit dummies for counties 59, 65, and 71. For San Francisco, include
    ! dummies for counties 13, 75, 81, and 85.
    ! For example, 59 -> A_la(,28)=1
    
    do i = 1, numSF
        if ( data_sf(i,_COUNTYID) == 13 ) then
            A_sf(i, 30) = 1
        else if ( data_sf(i,_COUNTYID) == 75 ) then
            A_sf(i, 31) = 1
        else if ( data_sf(i,_COUNTYID) == 81 ) then
            A_sf(i, 32) = 1
        else if ( data_sf(i,_COUNTYID) == 85 ) then
            A_sf(i, 33) = 1
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
    n = 32
    nrhs = 1
    lda = numLA
    ldb = numLA
    lwork = 1000000
    call dgelsy(m, n, nrhs, A_la, lda, b_la, ldb, jpvt, rcond, rank, work, lwork, info)
    
    do i = 1,32
        beta_la(i) = b_la(i)
    end do
    
    ! SF
    
    m = numSF
    n = 33
    nrhs = 1
    lda = numSF
    ldb = numSF
    lwork = 1000000
    call dgelsy(m, n, nrhs, A_sf, lda, b_sf, ldb, jpvt, rcond, rank, work, lwork, info)
    
    do i = 1,33
        beta_sf(i) = b_sf(i)
    end do

end subroutine  Regression_Hedonic_Price_Function











! #########
! Generate Bootstrapped index (follow the hint)
! #########
    
subroutine Bootstrapped(data_la, data_sf, boot_data_la, boot_data_sf)

    USE Global_Data

    
    implicit none
    
    ! Variables
    
    integer i, j, k
    
    real(KIND=DOUBLE)  u_la(numLA), u_sf(numSF)
    
    ! Data
    ! Intergated Data
    real(KIND=DOUBLE), DIMENSION(numLA, numV), intent(in)   :: data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV), intent(in)   :: data_sf
    ! Bootstrapped output
    real(KIND=DOUBLE), DIMENSION(numLA, numV), intent(out)  :: boot_data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV), intent(out)  :: boot_data_sf
    
    ! Bootstrapped index
    integer   r_la(numLA), r_sf(numSF)
    
    ! Generate random index array
    call init_random_seed()
    
    call RANDOM_NUMBER(u_la)
    call RANDOM_NUMBER(u_sf)
    
    do i = 1,numLA
        ! call random_number(u)
        r_la(i) = 1 + FLOOR(numLA*u_la(i)) 
    end do
    
    do i = 1,numSF
        ! call random_number(u)
        r_sf(i) = 1 + FLOOR(numSF*u_sf(i)) 
    end do
    
    ! Construct bootstrapped dataset
    do i = 1, numLA
        boot_data_la(i,:) = data_la(r_la(i),:)
    end do
    
    do i = 1, numSF
        boot_data_sf(i,:) = data_sf(r_sf(i),:)
    end do
    
end subroutine Bootstrapped

        SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
          END SUBROUTINE
          
    

end Module Bootstrapped_Hedonic_Price_Function