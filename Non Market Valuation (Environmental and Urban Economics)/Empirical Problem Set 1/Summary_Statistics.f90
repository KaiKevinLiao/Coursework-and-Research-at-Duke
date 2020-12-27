
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
    
Module Summary_Statistics

Contains
    
subroutine Compute_Summary_Statistics(data_la, data_sf)
    
    USE Global_Data

    
    implicit none
    
    ! Variables
    
    integer i, j, k
    
    ! Data
    ! Intergated Data
    real(KIND=DOUBLE), DIMENSION(numLA, numV)  :: data_la
    real(KIND=DOUBLE), DIMENSION(numSF, numV)  :: data_sf
    
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
    
    ! Statistics
    ! Means
    real(KIND=DOUBLE) price_mean_la(numC+1), square_mean_la(numC+1), bathrooms_mean_la(numC+1), bedrooms_mean_la(numC+1), &
                      rooms_mean_la(numC+1), stories_mean_la(numC+1), violentcrime_mean_la(numC+1), propertycrime_mean_la(numC+1)
    real(KIND=DOUBLE) price_mean_sf(numC+1), square_mean_sf(numC+1), bathrooms_mean_sf(numC+1), bedrooms_mean_sf(numC+1), &
                      rooms_mean_sf(numC+1), stories_mean_sf(numC+1), violentcrime_mean_sf(numC+1), propertycrime_mean_sf(numC+1)
    
    ! Second Moments
    real(KIND=DOUBLE) price_mom_la(numC+1), square_mom_la(numC+1), bathrooms_mom_la(numC+1), bedrooms_mom_la(numC+1), &
                      rooms_mom_la(numC+1), stories_mom_la(numC+1), violentcrime_mom_la(numC+1), propertycrime_mom_la(numC+1)
    real(KIND=DOUBLE) price_mom_sf(numC+1), square_mom_sf(numC+1), bathrooms_mom_sf(numC+1), bedrooms_mom_sf(numC+1), &
                      rooms_mom_sf(numC+1), stories_mom_sf(numC+1), violentcrime_mom_sf(numC+1), propertycrime_mom_sf(numC+1)
    
    ! Variances
    real(KIND=DOUBLE) price_var_la(numC+1), square_var_la(numC+1), bathrooms_var_la(numC+1), bedrooms_var_la(numC+1), &
                      rooms_var_la(numC+1), stories_var_la(numC+1), violentcrime_var_la(numC+1), propertycrime_var_la(numC+1)
    real(KIND=DOUBLE) price_var_sf(numC+1), square_var_sf(numC+1), bathrooms_var_sf(numC+1), bedrooms_var_sf(numC+1), &
                      rooms_var_sf(numC+1), stories_var_sf(numC+1), violentcrime_var_sf(numC+1), propertycrime_var_sf(numC+1)
    
    ! Function
     CHARACTER(LEN=120) :: Transfer_CountyID



    
    ! %%%%%%%%%%%%%%
    ! Means
    ! %%%%%%%%%%%%%%
    
    ! Mean for all
    ! LA
    do i = 1,numLA
        price_mean_la(numC+1) = price_mean_la(numC+1) + data_la(i, _PRICE)
        square_mean_la(numC+1) = square_mean_la(numC+1) + data_la(i, _SQUAREFEET)
        bathrooms_mean_la(numC+1) = bathrooms_mean_la(numC+1) + data_la(i,_BATHROOMS)
        bedrooms_mean_la(numC+1) = bedrooms_mean_la(numC+1) + data_la(i,_BEDROOMS)
        rooms_mean_la(numC+1) = rooms_mean_la(numC+1) + data_la(i,_ROOMS)
        stories_mean_la(numC+1) = stories_mean_la(numC+1) + data_la(i,_STORIES)
        violentcrime_mean_la(numC+1) = violentcrime_mean_la(numC+1) + data_la(i, _VIOLENT)
        propertycrime_mean_la(numC+1) = propertycrime_mean_la(numC+1) + data_la(i, _PROPERTY)
    end do
    price_mean_la(numC+1) = price_mean_la(numC+1)/numLA
    square_mean_la(numC+1) = square_mean_la(numC+1)/numLA
    bathrooms_mean_la(numC+1) = bathrooms_mean_la(numC+1)/numLA
    bedrooms_mean_la(numC+1) = bedrooms_mean_la(numC+1)/numLA
    rooms_mean_la(numC+1) = rooms_mean_la(numC+1)/numLA
    stories_mean_la(numC+1) = stories_mean_la(numC+1)/numLA
    violentcrime_mean_la(numC+1) = violentcrime_mean_la(numC+1)/numLA
    propertycrime_mean_la(numC+1) = propertycrime_mean_la(numC+1)/numLA
    
    ! SF
    do i = 1,numSF
        price_mean_sf(numC+1) = price_mean_sf(numC+1) + data_sf(i, _PRICE)
        square_mean_sf(numC+1) = square_mean_sf(numC+1) + data_sf(i, _SQUAREFEET)
        bathrooms_mean_sf(numC+1) = bathrooms_mean_sf(numC+1) + data_sf(i,_BATHROOMS)
        bedrooms_mean_sf(numC+1) = bedrooms_mean_sf(numC+1) + data_sf(i,_BEDROOMS)
        rooms_mean_sf(numC+1) = rooms_mean_sf(numC+1) + data_sf(i,_ROOMS)
        stories_mean_sf(numC+1) = stories_mean_sf(numC+1) + data_sf(i,_STORIES)
        violentcrime_mean_sf(numC+1) = violentcrime_mean_sf(numC+1) + data_sf(i, _VIOLENT)
        propertycrime_mean_sf(numC+1) = propertycrime_mean_sf(numC+1) + data_sf(i, _PROPERTY)
    end do
    price_mean_sf(numC+1) = price_mean_sf(numC+1)/numSF
    square_mean_sf(numC+1) = square_mean_sf(numC+1)/numSF
    bathrooms_mean_sf(numC+1) = bathrooms_mean_sf(numC+1)/numSF
    bedrooms_mean_sf(numC+1) = bedrooms_mean_sf(numC+1)/numSF
    rooms_mean_sf(numC+1) = rooms_mean_sf(numC+1)/numSF
    stories_mean_sf(numC+1) = stories_mean_sf(numC+1)/numSF
    violentcrime_mean_sf(numC+1) = violentcrime_mean_sf(numC+1)/numSF
    propertycrime_mean_sf(numC+1) = propertycrime_mean_sf(numC+1)/numSF

    print *, price_mean_la(numC+1)
    
    ! Means for different regions
    ! LA
    do j = 1,numC
        k = 0
        do i = 1, numLA
            if ( County_iden_la(j) .EQ. int(data_la(i,_COUNTYID)) ) then
                price_mean_la(j) = price_mean_la(j) + data_la(i, _PRICE)
                square_mean_la(j) = square_mean_la(j) + data_la(i, _SQUAREFEET)
                bathrooms_mean_la(j) = bathrooms_mean_la(j) + data_la(i,_BATHROOMS)
                bedrooms_mean_la(j) = bedrooms_mean_la(j) + data_la(i,_BEDROOMS)
                rooms_mean_la(j) = rooms_mean_la(j) + data_la(i,_ROOMS)
                stories_mean_la(j) = stories_mean_la(j) + data_la(i,_STORIES)
                violentcrime_mean_la(j) = violentcrime_mean_la(j) + data_la(i, _VIOLENT)
                propertycrime_mean_la(j) = propertycrime_mean_la(j) + data_la(i, _PROPERTY)
                k = k+1
            end if
        end do
        price_mean_la(j) = price_mean_la(j)/k
        square_mean_la(j) = square_mean_la(j)/k
        bathrooms_mean_la(j) = bathrooms_mean_la(j)/k
        bedrooms_mean_la(j) = bedrooms_mean_la(j)/k
        rooms_mean_la(j) = rooms_mean_la(j)/k
        stories_mean_la(j) = stories_mean_la(j)/k
        violentcrime_mean_la(j) = violentcrime_mean_la(j)/k
        propertycrime_mean_la(j) = propertycrime_mean_la(j)/k
    end do
    
    ! SF
    do j = 1,numC
        k = 0
        do i = 1, numSF
            if ( County_iden_sf(j) .EQ. int(data_sf(i,_COUNTYID)) ) then
                price_mean_sf(j) = price_mean_sf(j) + data_sf(i, _PRICE)
                square_mean_sf(j) = square_mean_sf(j) + data_sf(i, _SQUAREFEET)
                bathrooms_mean_sf(j) = bathrooms_mean_sf(j) + data_sf(i,_BATHROOMS)
                bedrooms_mean_sf(j) = bedrooms_mean_sf(j) + data_sf(i,_BEDROOMS)
                rooms_mean_sf(j) = rooms_mean_sf(j) + data_sf(i,_ROOMS)
                stories_mean_sf(j) = stories_mean_sf(j) + data_sf(i,_STORIES)
                violentcrime_mean_sf(j) = violentcrime_mean_sf(j) + data_sf(i, _VIOLENT)
                propertycrime_mean_sf(j) = propertycrime_mean_sf(j) + data_sf(i, _PROPERTY)
                k = k+1
            end if
        end do
        price_mean_sf(j) = price_mean_sf(j)/k
        square_mean_sf(j) = square_mean_sf(j)/k
        bathrooms_mean_sf(j) = bathrooms_mean_sf(j)/k
        bedrooms_mean_sf(j) = bedrooms_mean_sf(j)/k
        rooms_mean_sf(j) = rooms_mean_sf(j)/k
        stories_mean_sf(j) = stories_mean_sf(j)/k
        violentcrime_mean_sf(j) = violentcrime_mean_sf(j)/k
        propertycrime_mean_sf(j) = propertycrime_mean_sf(j)/k
    end do
    
    print *, price_mean_sf(2)
    
    ! %%%%%%%%%%%%%
    ! Second Moments
    ! %%%%%%%%%%%%%
    
    ! Second Moments for all
    ! LA

    do i = 1,numLA
        price_mom_la(numC+1) = price_mom_la(numC+1) + data_la(i, _PRICE)**2
        square_mom_la(numC+1) = square_mom_la(numC+1) + data_la(i, _SQUAREFEET)**2
        bathrooms_mom_la(numC+1) = bathrooms_mom_la(numC+1) + data_la(i,_BATHROOMS)**2
        bedrooms_mom_la(numC+1) = bedrooms_mom_la(numC+1) + data_la(i,_BEDROOMS)**2
        rooms_mom_la(numC+1) = rooms_mom_la(numC+1) + data_la(i,_ROOMS)**2
        stories_mom_la(numC+1) = stories_mom_la(numC+1) + data_la(i,_STORIES)**2
        violentcrime_mom_la(numC+1) = violentcrime_mom_la(numC+1) + data_la(i, _VIOLENT)**2
        propertycrime_mom_la(numC+1) = propertycrime_mom_la(numC+1) + data_la(i, _PROPERTY)**2
    end do
    price_mom_la(numC+1) = price_mom_la(numC+1)/numLA
    square_mom_la(numC+1) = square_mom_la(numC+1)/numLA
    bathrooms_mom_la(numC+1) = bathrooms_mom_la(numC+1)/numLA
    bedrooms_mom_la(numC+1) = bedrooms_mom_la(numC+1)/numLA
    rooms_mom_la(numC+1) = rooms_mom_la(numC+1)/numLA
    stories_mom_la(numC+1) = stories_mom_la(numC+1)/numLA
    violentcrime_mom_la(numC+1) = violentcrime_mom_la(numC+1)/numLA
    propertycrime_mom_la(numC+1) = propertycrime_mom_la(numC+1)/numLA
    
    ! SF
    do i = 1,numSF
        price_mom_sf(numC+1) = price_mom_sf(numC+1) + data_sf(i, _PRICE)**2
        square_mom_sf(numC+1) = square_mom_sf(numC+1) + data_sf(i, _SQUAREFEET)**2
        bathrooms_mom_sf(numC+1) = bathrooms_mom_sf(numC+1) + data_sf(i,_BATHROOMS)**2
        bedrooms_mom_sf(numC+1) = bedrooms_mom_sf(numC+1) + data_sf(i,_BEDROOMS)**2
        rooms_mom_sf(numC+1) = rooms_mom_sf(numC+1) + data_sf(i,_ROOMS)**2
        stories_mom_sf(numC+1) = stories_mom_sf(numC+1) + data_sf(i,_STORIES)**2
        violentcrime_mom_sf(numC+1) = violentcrime_mom_sf(numC+1) + data_sf(i, _VIOLENT)**2
        propertycrime_mom_sf(numC+1) = propertycrime_mom_sf(numC+1) + data_sf(i, _PROPERTY)**2
    end do
    price_mom_sf(numC+1) = price_mom_sf(numC+1)/numSF
    square_mom_sf(numC+1) = square_mom_sf(numC+1)/numSF
    bathrooms_mom_sf(numC+1) = bathrooms_mom_sf(numC+1)/numSF
    bedrooms_mom_sf(numC+1) = bedrooms_mom_sf(numC+1)/numSF
    rooms_mom_sf(numC+1) = rooms_mom_sf(numC+1)/numSF
    stories_mom_sf(numC+1) = stories_mom_sf(numC+1)/numSF
    violentcrime_mom_sf(numC+1) = violentcrime_mom_sf(numC+1)/numSF
    propertycrime_mom_sf(numC+1) = propertycrime_mom_sf(numC+1)/numSF

    ! print *, price_mom_la(numC+1)
    
    ! Second Moments for different regions
    ! LA
    do j = 1,numC
        k = 0
        do i = 1, numLA
            if ( County_iden_la(j) .EQ. int(data_la(i,_COUNTYID)) ) then
                price_mom_la(j) = price_mom_la(j) + data_la(i, _PRICE)**2
                square_mom_la(j) = square_mom_la(j) + data_la(i, _SQUAREFEET)**2
                bathrooms_mom_la(j) = bathrooms_mom_la(j) + data_la(i,_BATHROOMS)**2
                bedrooms_mom_la(j) = bedrooms_mom_la(j) + data_la(i,_BEDROOMS)**2
                rooms_mom_la(j) = rooms_mom_la(j) + data_la(i,_ROOMS)**2
                stories_mom_la(j) = stories_mom_la(j) + data_la(i,_STORIES)**2
                violentcrime_mom_la(j) = violentcrime_mom_la(j) + data_la(i, _VIOLENT)**2
                propertycrime_mom_la(j) = propertycrime_mom_la(j) + data_la(i, _PROPERTY)**2
                k = k+1
            end if
        end do
        price_mom_la(j) = price_mom_la(j)/k
        square_mom_la(j) = square_mom_la(j)/k
        bathrooms_mom_la(j) = bathrooms_mom_la(j)/k
        bedrooms_mom_la(j) = bedrooms_mom_la(j)/k
        rooms_mom_la(j) = rooms_mom_la(j)/k
        stories_mom_la(j) = stories_mom_la(j)/k
        violentcrime_mom_la(j) = violentcrime_mom_la(j)/k
        propertycrime_mom_la(j) = propertycrime_mom_la(j)/k
    end do
    
    ! SF
    do j = 1,numC
        k = 0
        do i = 1, numSF
            if ( County_iden_sf(j) .EQ. int(data_sf(i,_COUNTYID)) ) then
                price_mom_sf(j) = price_mom_sf(j) + data_sf(i, _PRICE)**2
                square_mom_sf(j) = square_mom_sf(j) + data_sf(i, _SQUAREFEET)**2
                bathrooms_mom_sf(j) = bathrooms_mom_sf(j) + data_sf(i,_BATHROOMS)**2
                bedrooms_mom_sf(j) = bedrooms_mom_sf(j) + data_sf(i,_BEDROOMS)**2
                rooms_mom_sf(j) = rooms_mom_sf(j) + data_sf(i,_ROOMS)**2
                stories_mom_sf(j) = stories_mom_sf(j) + data_sf(i,_STORIES)**2
                violentcrime_mom_sf(j) = violentcrime_mom_sf(j) + data_sf(i, _VIOLENT)**2
                propertycrime_mom_sf(j) = propertycrime_mom_sf(j) + data_sf(i, _PROPERTY)**2
                k = k+1
            end if
        end do
        price_mom_sf(j) = price_mom_sf(j)/k
        square_mom_sf(j) = square_mom_sf(j)/k
        bathrooms_mom_sf(j) = bathrooms_mom_sf(j)/k
        bedrooms_mom_sf(j) = bedrooms_mom_sf(j)/k
        rooms_mom_sf(j) = rooms_mom_sf(j)/k
        stories_mom_sf(j) = stories_mom_sf(j)/k
        violentcrime_mom_sf(j) = violentcrime_mom_sf(j)/k
        propertycrime_mom_sf(j) = propertycrime_mom_sf(j)/k
    end do
    
    ! %%%%%%%%%%%%%
    ! Variances
    ! %%%%%%%%%%%%%
    
    ! Second Moments for all
    ! LA
    do j = 1,numC+1
        price_var_la(j) = price_mom_la(j) - price_mean_la(j)**2  
        square_var_la(j) = 	square_mom_la(j) - square_mean_la(j)**2 
        bathrooms_var_la(j) =  bathrooms_mom_la(j) - bathrooms_mean_la(j)**2  
        bedrooms_var_la(j) = bedrooms_mom_la(j) - bedrooms_mean_la(j)**2  
        rooms_var_la(j) = rooms_mom_la(j) - rooms_mean_la(j)**2  
        stories_var_la(j) = stories_mom_la(j) - stories_mean_la(j)**2  
        violentcrime_var_la(j) = violentcrime_mom_la(j) - violentcrime_mean_la(j)**2  
        propertycrime_var_la(j) = propertycrime_mom_la(j) - propertycrime_mean_la(j)**2 
    end do

    ! SF
    do j = 1,numC+1
        price_var_sf(j) = price_mom_sf(j) - price_mean_sf(j)**2  
        square_var_sf(j) = 	square_mom_sf(j) - square_mean_sf(j)**2 
        bathrooms_var_sf(j) =  bathrooms_mom_sf(j) - bathrooms_mean_sf(j)**2  
        bedrooms_var_sf(j) = bedrooms_mom_sf(j) - bedrooms_mean_sf(j)**2  
        rooms_var_sf(j) = rooms_mom_sf(j) - rooms_mean_sf(j)**2  
        stories_var_sf(j) = stories_mom_sf(j) - stories_mean_sf(j)**2  
        violentcrime_var_sf(j) = violentcrime_mom_sf(j) - violentcrime_mean_sf(j)**2  
        propertycrime_var_sf(j) = propertycrime_mom_sf(j) - propertycrime_mean_sf(j)**2 
    end do
    
    
    ! ############
    ! Write summary statistics into a csv file
    ! ############
    
  
    open(unit = 100 , file = trim(dir_mom)//trim(filename_summary)//'.csv')
    
    ! Header
    
    write(100,499) ',,,Price,Square Footage,# Bathrooms,# Bedrooms,# Total Rooms,# Stories,Violent Crime Rate (Cases per 100000),Property Crime Rate (Cases Per 100000)'
    
    ! Means
    
    write(100,500) 'Means,Los Angel,All sample',',',price_mean_la(numC+1), ',', square_mean_la(numC+1),',', bathrooms_mean_la(numC+1),',', bedrooms_mean_la(numC+1),',', &
                      rooms_mean_la(numC+1),',', stories_mean_la(numC+1),',', violentcrime_mean_la(numC+1),',', propertycrime_mean_la(numC+1)
    
    do i = 1,numC
        write(100,501) ',,',Transfer_CountyID(County_iden_la(i)),',',price_mean_la(i), ',', square_mean_la(i),',', bathrooms_mean_la(i),',', bedrooms_mean_la(i),',', &
                      rooms_mean_la(i),',', stories_mean_la(i),',', violentcrime_mean_la(i),',', propertycrime_mean_la(i)
    end do
    
    write(100,500) ',LA,All sample',',',price_mean_sf(numC+1), ',', square_mean_sf(numC+1),',', bathrooms_mean_sf(numC+1),',', bedrooms_mean_sf(numC+1),',', &
                      rooms_mean_sf(numC+1),',', stories_mean_sf(numC+1),',', violentcrime_mean_sf(numC+1),',', propertycrime_mean_sf(numC+1)
    
    do i = 1,numC
        write(100,501) ',,',Transfer_CountyID(County_iden_sf(i)),',',price_mean_sf(i), ',', square_mean_sf(i),',', bathrooms_mean_sf(i),',', bedrooms_mean_sf(i),',', &
                      rooms_mean_sf(i),',', stories_mean_sf(i),',', violentcrime_mean_sf(i),',', propertycrime_mean_sf(i)
    end do
    
    ! Variances
    
    write(100,500) 'Variances,Los Angel,All sample',',',price_var_la(numC+1), ',', square_var_la(numC+1),',', bathrooms_var_la(numC+1),',', bedrooms_var_la(numC+1),',', &
                      rooms_var_la(numC+1),',', stories_var_la(numC+1),',', violentcrime_var_la(numC+1),',', propertycrime_var_la(numC+1)
    
    do i = 1,numC
        write(100,501) ',,',Transfer_CountyID(County_iden_la(i)),',',price_var_la(i), ',', square_var_la(i),',', bathrooms_var_la(i),',', bedrooms_var_la(i),',', &
                      rooms_var_la(i),',', stories_var_la(i),',', violentcrime_var_la(i),',', propertycrime_var_la(i)
    end do
    
    write(100,500) ',SF,All sample',',',price_var_sf(numC+1), ',', square_var_sf(numC+1),',', bathrooms_var_sf(numC+1),',', bedrooms_var_sf(numC+1),',', &
                      rooms_var_sf(numC+1),',', stories_var_sf(numC+1),',', violentcrime_var_sf(numC+1),',', propertycrime_var_sf(numC+1)
    
    do i = 1,numC
        write(100,501) ',,',Transfer_CountyID(County_iden_sf(i)),',',price_var_sf(i), ',', square_var_sf(i),',', bathrooms_var_sf(i),',', bedrooms_var_sf(i),',', &
                      rooms_var_sf(i),',', stories_var_sf(i),',', violentcrime_var_sf(i),',', propertycrime_var_sf(i)
    end do


    
    500 format(a27,8(a1,f22.8))
    499 format(a147)
    501 format(a2,a20,8(a1,f22.8))


    


    
end subroutine Compute_Summary_Statistics

        
        
    
end Module Summary_Statistics