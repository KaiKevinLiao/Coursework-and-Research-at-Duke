program Main

    ! Module defining and setting global variables
    USE Global_Data
    USE LossFunction_MOD
    USE LinearAlgebra_MOD
    USE Nelder_Mead
    USE simulated_anneal
    USE computeNewSteadyState_MOD
    
    implicit none
    
    character(len=1024) folder
    character(len=1024) filename
    
    integer i, k, l, n, o, yr, ind
    integer NP_in, Nobs, mask_country(NP)
    
    real(KIND=DOUBLE) param_vec(NP), param_scaled(NP)
    
    real(KIND=DOUBLE), allocatable, dimension(:) :: param_in, param_out
    
    real(KIND=DOUBLE) F
    
    integer, parameter :: maxiter_Y         = 10000
    real(KIND=DOUBLE), parameter :: toler_Y = 1e-12
    
    integer iter
    real(KIND=DOUBLE) dist_Y, big_Y_new(numK,numC)
    
    ! For Nelder-Mead
    real(KIND=DOUBLE) fvalue, ftol, simp
    integer           MAXFCN, NLOOP, iquad, ifault, iprint
    real(KIND=DOUBLE), allocatable, dimension(:) ::  stepsize, var
    
    ! For Simulated Annealing
    logical maxim
    real(KIND=DOUBLE) rT, eps, Temp, fopt
    real(KIND=DOUBLE), allocatable, dimension(:) :: c(:), &
                                                    VM(:), param_opt(:)
    integer NS, NT, Neps, maxevl, iseed1, iseed2, nacc, nfcnev, nobds, ier
    ! random seed for the simulated annealing procedure
    real(KIND=DOUBLE) random_array(2)
    
    type(SteadyState) SteadyState0, SteadyState1
    type(Shock)       Shock_in
    real(KIND=DOUBLE) big_L0(numK,numC), x_bar0(numK,numC), w_tilde0(numK,numC), big_U0(numK,numC)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Read Fixed Parameters
    open(unit = 1, file ='FixedParams.csv')
        
        do i = 1, NP_fixed
            read(1,*) n, FixedPar(i)
        end do
        
    close(1)
    
    lambda = FixedPar(1) 
    delta  = FixedPar(2)  
    xsi    = FixedPar(3)
    zeta   = FixedPar(4)
    beta   = FixedPar(5)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Read Starting Point for estimation procedure
    open(unit = 1, file ='Starting_Point.csv')
        
        do i = 1, NP
            read(1,*) n, param_vec(i), lb_global(i), ub_global(i), scaling_global(i), mask_global_in(i)
        end do
        
    close(1)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Build an additional variable, mask_country, that will interact with mask_global
    ! To select over what parameters the optimization will be carried out, depending
    ! on the value of COUNTRY_EST (the id of the country we choose to estimate)
    
    mask_country = 0
    
    ! Unemployment Benefit / b
    mask_country(COUNTRY_EST) = 1
    ind = numC
    
    ! Cost Vacancies / kappa
    mask_country((ind+numK*(COUNTRY_EST-1)+1):(ind+numK*(COUNTRY_EST-1)+numK)) = 1
    ind = ind + numK*numC
    
    ! sigma_x
    if (COUNTRY_EST == 1) then
        mask_country(ind+1) = 1
    end if
    ind = ind + 1
    
    ! Exogenous Exit / chi (country component)
    mask_country(ind+COUNTRY_EST) = 1
    ind = ind + numC
    ! Exogenous Exit / chi (sector component)
    if (COUNTRY_EST == 1) then
        mask_country((ind+1):(ind+numK)) = 1
    end if
    ind = ind + numK
    
    ! Compensating Differentials / eta
    mask_country((ind+numK*(COUNTRY_EST-1)+1):(ind+numK*(COUNTRY_EST-1)+numK)) = 1
    ind = ind + numK*numC
    
    ! Mobility Costs / C
    if (COUNTRY_EST == 1) then
        mask_country((ind+1):(ind+numK*numC)) = 1
    end if
    ind = ind+numK*numC
    
    ! Interact mask_global_in with mask_country to focus on COUNTRY_EST specific parameters
    mask_global = mask_global_in*mask_country
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Scale parameters to enter the optimizer
    
    param_scaled = param_vec / scaling_global 

    param_global = param_scaled 

    ! Construct new parameter vector that enters the optimizer -- the remaining parameters,
    ! over which we don't optimize are fixed throughout the procedure
    
    ! number of parameters over which we optimize
    
    NP_in = sum(mask_global) 
    
    allocate(lb_in(NP_in), ub_in(NP_in))
    allocate(param_in(NP_in))
    allocate(param_out(NP_in))

    ! parameter vector that we feed the optimizer with
    
    k = 1 
    do i = 1, NP
        if(mask_global(i) == 1) then
            param_in(k) = param_scaled(i) 
            k = k + 1 
        end if
    end do
    
    k = 1 
    do i = 1, NP
	    if (mask_global(i) == 1) then
    	    lb_in(k) = lb_global(i) / scaling_global(i) 
            ub_in(k) = ub_global(i) / scaling_global(i) 
            k = k + 1 
        end if
    end do
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Read Data Moments
    Call Read_Data
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    !%%%%%%%%%%%%%%%%%%%%%%%%
    ! Step 1: Solve for big_Y
    !%%%%%%%%%%%%%%%%%%%%%%%%
    
    dist_Y = 1
    big_Y  = 1.0
    iter   = 1
    do while (dist_Y > toler_Y .and. iter < maxiter_Y)
        
        big_Y_new = 0.0
        do k = 1, numK
            do o = 1, numC
                do i = 1, numC
                    do l = 1, numK
                        big_Y_new(k,o) = big_Y_new(k,o) + pi_Data(o,i,k)*(mu(k,i)*gamma(l,i) + (1-gamma(l,i))*nu(l,k,i))*big_Y(l,i)
                    end do
                    big_Y_new(k,o) = big_Y_new(k,o) - pi_Data(o,i,k)*mu(k,i)*NX_Data(i,1)
                end do
            end do
        end do
        
        dist_Y = sum( abs((big_Y_new - big_Y)/max(big_Y,0.0000001)) )
        
        big_Y = big_Y_new / sum( big_Y_new )
        
        iter = iter + 1
    
    end do
    
    big_Y = big_Y / sum(big_Y)
    
    ! Print Big_Y to file
    open(unit = 1234 , file = 'big_Y.csv')
        do k = 1, numK
            write(1234,333) (big_Y(k,o), o=1,numC)
        end do
    close(1234)
    
333 format(6(f50.16,','))    
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Read Current Outcomes in "Estimation" Steady State
    
    do i = 1, numC
        folder = ''
        filename = ''
        write (folder, "(A7,I1,A1)") "Country", i, "\"
        write (filename, "(A16,I1)") "Initial_SS_big_L", i
        open(unit = 100 , file = trim(folder)//trim(filename)//".csv")
            do k = 1, numK
                read(100,876) big_L0(k,i)
            end do
        close(100)
    end do
    
    do i = 1, numC
        folder = ''
        filename = ''
        write (folder, "(A7,I1,A1)") "Country", i, "\"
        write (filename, "(A15,I1)") "Initial_SS_xbar", i
        open(unit = 100 , file = trim(folder)//trim(filename)//".csv")
            do k = 1, numK
                read(100,876) x_bar0(k,i)
            end do
        close(100)
    end do
    
    do i = 1, numC
        folder = ''
        filename = ''
        write (folder, "(A7,I1,A1)") "Country", i, "\"
        write (filename, "(A18,I1)") "Initial_SS_w_tilde", i
        open(unit = 100 , file = trim(folder)//trim(filename)//".csv")
            do k = 1, numK
                read(100,876) w_tilde0(k,i)
            end do
        close(100)
    end do
    
    do i = 1, numC
        folder = ''
        filename = ''
        write (folder, "(A7,I1,A1)") "Country", i, "\"
        write (filename, "(A16,I1)") "Initial_SS_big_U", i
        open(unit = 100 , file = trim(folder)//trim(filename)//".csv")
            do k = 1, numK
                read(100,876) big_U0(k,i)
            end do
        close(100)
    end do
    
    SteadyState0 % pi       = pi_Data
    SteadyState0 % big_L    = big_L0
    SteadyState0 % x_bar    = x_bar0
    SteadyState0 % w_tilde  = w_tilde0
    SteadyState0 % big_U    = big_U0
    SteadyState0 % big_Y    = big_Y
    
    Shock_in % d_hat = 1
    Shock_in % A_hat = 1
    Shock_in % NX    = NX_Data
    
    Call NewSteadyState(param_scaled,Shock_in,SteadyState0,SteadyState1)

876 format(f50.16)    
    
    ! %%%%%%%%%%%
    ! Nelder-Mead
    ! %%%%%%%%%%%
    
    if (ALGORITHM == 1) then
    
        FUNCTION_ITERATION = 1
        F_BEST = 1000000000
    
        allocate(stepsize(NP_in))
        allocate(var(NP_in))
    
        ftol     = 0.001
        stepsize = 0.3
        MAXFCN   = 25000
        iprint   = -1
        NLOOP    = 2*NP_in
        iquad    = 1
        simp     = 0.000001                 
       
        write (folder, "(A7,I1,A1)") "Country", COUNTRY_EST, "\"
        write (filename, "(A19,I1)") "FunctionEvaluations", COUNTRY_EST
        
        open(unit = 999 , file = trim(folder)//trim(filename)//'.csv')
    
        Call PrintHeader_FncEval
    
        Call SMM_OBJ_FCN(param_in,F)
    
        Call minim(param_in, stepsize, NP_in, fvalue, MAXFCN, iprint, ftol, NLOOP, iquad,  &
                  simp, var, SMM_OBJ_FCN, ifault) 
    
        close(999)
    
        print*, 'Nelder Mead is done'
        
        deallocate(stepsize)
        deallocate(var)
        
        pause
        
    end if
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! %%%%%%%%%%%%%%%%%%%
    ! Simulated Annealing
    ! %%%%%%%%%%%%%%%%%%%
    
    if (ALGORITHM == 2) then
    
        FUNCTION_ITERATION = 1
        F_BEST = 1000000000
    
        allocate(c(NP_in), VM(NP_in), param_opt(NP_in))
        
        ! maximization?
        maxim = .false.
        ! Temperature reduction factor
        rT = 0.85
        ! tolerance
        eps = 0.0001
        ! Number of cycles. After NS*N function evaluations, each element of VM is adjusted
        NS = 20
        ! random number seeds
        ! Number of iterations before temperature reduction
        NT = max(5*NP_in, 100)
        ! Number of final function values used to decide upon termination
        Neps = 4
        iprint = 0
    
        ! Vector that controls the step length adjustment.
        ! The suggested value for all elements is 2.0.
        c = 2.0
        ! the step length vector
        VM = 0.3
    
        ! Initial temperature 
        !T = ( max_feval - min_feval ) / 10000.0
        Temp = TempSimAnn
        print*
        print*, 'The initial Temperature is: ', Temp
        ! maximum number of evaluations
        maxevl = 100000000
    
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(random_array)
        
        ! resetting the seeds
        ! Random seeds
        iseed1 = int(2200*random_array(1))
        iseed2 = int(2200*random_array(2))
        
        write (folder, "(A7,I1,A1)") "Country", COUNTRY_EST, "\"
        write (filename, "(A19,I1)") "FunctionEvaluations", COUNTRY_EST
        
        open(unit = 999 , file = trim(folder)//trim(filename)//'.csv')
        
        Call PrintHeader_FncEval
    
        Call SMM_OBJ_FCN(param_in,F)
    
        CALL sa(NP_in, param_in, maxim, rT, eps, NS, NT, Neps, maxevl, lb_in, ub_in, c, iprint, iseed1,  &
                iseed2, Temp, VM, param_opt, fopt, nacc, nfcnev, nobds, ier)
    
        close(999)
    
    end if
    
    pause
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    deallocate(lb_in,ub_in)
    deallocate(param_in)
    deallocate(param_out)
    deallocate(c, VM, param_opt)

Contains     

Subroutine Read_Data

    real(KIND=DOUBLE) L_US, Ergodic_Matrix_old(numK+1,numK+1), Ergodic_Matrix_new(numK+1,numK+1), &
                      dist_erg, World_big_Y_Data
    integer iter_erg
    
    !%%%%%%%%%%%
    ! Read gamma
    !%%%%%%%%%%%
    Nobs = numC*numK
    open(unit=1,file=trim(dir_mom)//'gammaMat.csv')
        read(1,*)
        do n = 1, Nobs
            read(1,*) i, k, yr, gamma(k,i)
        end do        
    close(1)
    
    !%%%%%%%%%%%
    ! Read nu
    !%%%%%%%%%%%
    Nobs = numC*numK*numK
    open(unit=2,file=trim(dir_mom)//'nuMat.csv')
        read(2,*)
        do n = 1, Nobs
            read(2,*) i, k, l, yr, nu(k,l,i)
        end do        
    close(2)
    
    do k = 1, numK
        do i = 1, numC
            nu(k,:,i) = nu(k,:,i) / sum(nu(k,:,i))
        end do
    end do
    
    !%%%%%%%%%%%
    ! Read mu
    !%%%%%%%%%%%
    Nobs = numC*numK
    open(unit=3,file=trim(dir_mom)//'muMat.csv')
        read(3,*)
        do n = 1, Nobs
            read(3,*) i, k, yr, mu(k,i)
        end do        
    close(3)
    
    do i = 1, numC
        mu(:,i) = mu(:,i) / sum( mu(:,i) )
    end do
    
    !%%%%%%%%%%%
    ! Read pi
    !%%%%%%%%%%%
    Nobs = numC*numC*numK
    open(unit=4,file=trim(dir_mom)//'piMat.csv')
        read(4,*)
        do n = 1, Nobs
            read(4,*) o, i, k, yr, pi_Data(o,i,k)
        end do        
    close(4)
    
    do k = 1, numK
        do i = 1, numC
            pi_Data(:,i,k) = pi_Data(:,i,k) / sum( pi_Data(:,i,k) )
        end do
    end do
    
    !%%%%%%%%%%%%%
    ! Read NX_Data
    !%%%%%%%%%%%%%
    Nobs = numC
    open(unit=5,file=trim(dir_mom)//'nxMat.csv')
        read(5,*)
        do n = 1, Nobs
            read(5,*) i, yr, NX_Data(i,1)
        end do
    close(5)
    
    !%%%%%%%%%%%%%%%%
    ! Read big_Y_Data
    !%%%%%%%%%%%%%%%%
    Nobs = numK*numC
    open(unit=6,file=trim(dir_mom)//'goMat.csv')
        read(6,*)
        do n = 1, Nobs
            read(6,*) i, k, yr, big_Y_Data(k,i)
        end do        
    close(6)
    
    !%%%%%%%%%%%%%
    ! Read VA_Data
    !%%%%%%%%%%%%%
    !Nobs = numK*numC
    !open(unit=7,file=trim(dir_mom)//'vaMat.csv')
    !    read(7,*)
    !    do n = 1, Nobs
    !        read(7,*) i, k, yr, VA_Data(k,i)
    !    end do        
    !close(7)    
    
    !%%%%%%%%%%%%%%%%
    ! Read Unemp_Data
    !%%%%%%%%%%%%%%%%
    Nobs = numC
    open(unit=8,file=trim(dir_mom)//'unemp_ILO.csv')
        read(8,*)
        do n = 1, Nobs
            read(8,*) i, yr, Unemp_Data(i,1)
            Unemp_Data(i,1) = Unemp_Data(i,1) / 100
        end do        
    close(8) 
    
    !%%%%%%%%%%%%%%%%
    ! Read Emp_Data
    !%%%%%%%%%%%%%%%%
    Nobs = numK*numC
    open(unit=9,file=trim(dir_mom)//'wiot_sea_Employment.csv')
        read(9,*)
        do n = 1, Nobs
            read(9,*) i, k, yr, Emp_Data(k,i)
        end do        
    close(9) 
    
    ! Employment is expressed in thousands of workers
    ! Re-Scale
    Emp_Data = 1000.0*Emp_Data
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct ShareEmp_Data and big_L_bar_Data
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    do k = 1, numK
        do i = 1, numC
            ShareEmp_Data(k,i) = Emp_Data(k,i) / sum( Emp_Data(:,i) )
        end do
    end do
    
    do i = 1, numC
        big_L_bar_Data(i,1) = sum( Emp_Data(:,i) ) / (1 - Unemp_Data(i,1))     
    end do
        
    !%%%%%%%%%%%%%%%%
    ! Read w_bar_Data
    !%%%%%%%%%%%%%%%%
    Nobs = numK*numC
    open(unit=10,file=trim(dir_mom)//'wiot_sea_Wage.csv')
        read(10,*)
        do n = 1, Nobs
            read(10,*) i, k, yr, w_bar_Data(k,i)
        end do        
    close(10) 
    
    !%%%%%%%%%%%%%%%%%%%%%
    ! Read CoeffVar_w_Data
    !%%%%%%%%%%%%%%%%%%%%%
    Nobs = 1
    open(unit=11,file=trim(dir_mom)//'CPS_CoefVar.csv')
        read(11,*)
        do n = 1, Nobs
            read(11,*) l, yr, CoeffVar_w_Data
        end do       
    close(11)
    
    !%%%%%%%%%%%%%%%%%%
    ! Read s_tilde_Data
    !%%%%%%%%%%%%%%%%%%
    Nobs = (numK+1)*(numK+1)
    open(unit=12,file=trim(dir_mom)//'CPS_transitions.csv')
        read(12,*)
        do n = 1, Nobs
            read(12,*) k, l, yr, s_tilde_Data(k+1,l+1)
        end do       
    close(12)
    
    ! Make sure rows of s_tilde_Data sum exactly to 1
    do k = 1, numK+1
        s_tilde_Data(k,:) = s_tilde_Data(k,:) / sum( s_tilde_Data(k,:) )
    end do
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Ergodic Distribution Implied by s_tilde_Data
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ergodic_Matrix_old = 0
    do i = 1, numK+1
        Ergodic_Matrix_old(i,i) = 1
    end do
    iter_erg = 1
    dist_erg = 1
    do while( dist_erg > 0.0000000001 .and. iter_erg <= 2000 )
        Ergodic_Matrix_new = matrix_mult(s_tilde_Data,Ergodic_Matrix_old,numK+1,numK+1,numK+1)
        dist_erg = maxval( (Ergodic_Matrix_new - Ergodic_Matrix_old)/max(abs(Ergodic_Matrix_old),0.00000001) )
        Ergodic_Matrix_old = Ergodic_Matrix_new
        iter_erg = iter_erg + 1
    end do    
    Ergodic_Dist_Data = transpose( Ergodic_Matrix_old(1:1,:) )
    
    ! US unemployment and employment shares are consistent
    ! with the transition matrix
    Unemp_Data(1,1)    = Ergodic_Dist_Data(1,1)
    ShareEmp_Data(:,1) = Ergodic_Dist_Data(2:numK+1,1) / (1-Unemp_Data(1,1))
    
    ! Normalize country sizes -- US has size 1
    L_US           = big_L_bar_Data(1,1)
    big_L_bar(:,1) = big_L_bar_Data(:,1) / L_US
    
    ! World Gross Output is normalized to 1
    ! Normalize other nominal variables accordingly
    World_big_Y_Data = sum( big_Y_Data )
    w_bar_Data       = (w_bar_Data / World_big_Y_Data)*L_US
    NX_Data          = NX_Data / World_big_Y_Data
    big_Y_Data       = big_Y_Data / World_big_Y_Data
    
end subroutine Read_Data 

subroutine PrintHeader_FncEval

    write(999,'(3A)',advance='no') 'F ,'
    write(999,'(9A)',advance='no') 'penalty0,'
    write(999,'(9A)',advance='no') 'penalty1,'
    write(999,'(9A)',advance='no') 'penalty2,'
    write(999,'(9A)',advance='no') 'penalty3,'
    
    do i = 1, numC
        write(999,101,advance='no') 'Unemp. Ben. b(i=', i, '),'
    end do
    
    do i = 1, numC
        do k = 1, numK
            write(999,103,advance='no') 'Cost Vacancies kappa_tilde(k=', k, '-i=', i, '),'
        end do
    end do
    
    write(999,104,advance='no') 'sigma_x,'
    
    do i = 1, numC
        write(999,105,advance='no') 'chi(i=', i, '),'
    end do
    
    do k = 1, numK
        write(999,105,advance='no') 'chi(k=', k, '),'
    end do
    
    do i = 1, numC
        do k = 1, numK
            write(999,106,advance='no') 'eta(k=', k, '-i=', i, '),'
        end do
    end do
    
    do k = 1, numK
        do l = 1, numK
            write(999,107,advance='no') 'C_US(k=', k, '-l=', l, '),'
        end do
    end do
    
    write(999,*)
    
100 format(a15, i1, a3, i1, a2)
101 format(a16, i1, a2)
102 format(a19, i1, a3, i1, a2)
103 format(a29, i1, a3, i1, a2)    
104 format(a8, i1, a2)     
105 format(a6 , i1, a2)
106 format(a6, i1, a3, i1, a2)
107 format(a7, i1, a3, i1, a2)  

end subroutine PrintHeader_FncEval

! FCN2 is defined to be used with
! IMSL optimization routines
SUBROUTINE FCN2(N, X, F)

USE Global_Data
USE LossFunction_MOD

IMPLICIT NONE

INTEGER N
REAL (KIND=DOUBLE) X(N), F

Call SMM_OBJ_FCN(X,F)

RETURN
END SUBROUTINE FCN2
    
end program Main
    