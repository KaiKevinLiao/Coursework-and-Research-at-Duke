Module LossFunction_MOD

Contains
    
Subroutine SMM_OBJ_FCN(param_scaled,F)

USE Global_Data
USE LinearAlgebra_MOD
USE IFPORT

implicit none

integer i, i2, j, k, l
integer ind

real(KIND=DOUBLE), intent(in) :: param_scaled(:)
real(KIND=DOUBLE), intent(out) :: F

type(Parameters) params

real(KIND=DOUBLE) param_vec(NP)

real(KIND=DOUBLE) b(numK,numC), kappa_tilde(numK,numC), sigma_x(numK,numC), &
                  eta(numK,numC), chi(numK,numC), C_US(numK,numK)
                  
real(KIND=DOUBLE) b_cty(numC,1), chi_cty(numC,1), chi_ind(numC,1)
                  
real(KIND=DOUBLE) x_bar(numK), w_tilde(numK), w_bar(numK), Var_w(numK), mean_w, Var_w_tot, &
                  CoeffVar_w, s_tilde(2*numK,2*numK), s_tilde_4(2*numK,2*numK), &
                  s_tilde_annual(numK+1,numK+1), sum_s, ShareEmp(numK), &
                  labor_share(numK), Unemp, u(numK), big_L(numK), &
                  w_bar_norm(numK), w_bar_Data_norm(numK), s_tilde_annual_diag(numK+1,numK+1), &
                  s_tilde_Data_diag(numK+1,numK+1)



integer flag_param

real(KIND=DOUBLE) VacancyCosts(numK), VacancyCosts_Data(numK)

real(KIND=DOUBLE) WeightTr(numK+1,numK+1), WeightU, WeightSh(numK), &
                  WeightW(numK), WeightCoeffVarW, WeightVacCosts(numK)

real(KIND=DOUBLE) F_Tr, F_U, F_Sh, F_W, F_CoeffVarW, F_VacCosts, F_p0, F_p1, F_p2, F_p3
                  
real(KIND=DOUBLE) penalty0(NP), penalty1, penalty2, penalty3, &
                  dist_out, toler_out

integer flag_abort, flag_abort2

real(KIND=DOUBLE) Dev

real(KIND=DOUBLE) q_theta(numK), theta(numK), G_x_bar(numK)

type(MomentsOut) MomentsModel

character(len=1024) folder
character(len=1024) filename

! Elapsed time variables
integer           time_array_0(8), time_array_1(8)
real(KIND=DOUBLE) start_time, finish_time

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Start time
call date_and_time(values=time_array_0)
start_time = time_array_0(5) * 3600 + time_array_0(6) * 60 &
           + time_array_0(7) + 0.001 * time_array_0(8)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! De-scale parameters
! Check if all parameter values are within bounds
! If not, assign bound to parameter value and penalize the size of the violation
! This procedure is necessary when employing unrestricted Nelder-Mead
penalty0 = 0.0
k = 1
do i = 1, NP
    if (mask_global(i) == 1) then
        param_vec(i) = param_scaled(k) * scaling_global(i)
        if (param_vec(i) < lb_global(i)) then
            penalty0(i) = 1000.0*abs(lb_global(i) - param_vec(i))
            param_vec(i) = lb_global(i)
        else if (param_vec(i) > ub_global(i)) then
            penalty0(i) = 1000.0*abs(ub_global(i) - param_vec(i))
            param_vec(i) = ub_global(i)
        end if
        k = k + 1
    else
        param_vec(i) = param_global(i) * scaling_global(i)
    end if
end do

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Construct full numK by numC matrices of parameters
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = 0
do i = 1, numC
    b_cty(i,1) = param_vec(ind+i)
end do
ind = ind + numC
do k = 1, numK
    do i = 1, numC
        b(k,i) = b_cty(i,1)
    end do
end do

do k = 1, numK
    do i = 1, numC
        kappa_tilde(k,i) = param_vec(ind+(i-1)*numK+k)
    end do
end do
ind = ind + numK*numC

sigma_x = param_vec(ind+1)
ind = ind + 1

chi_cty = 0
chi_ind = 0
chi     = 0
do i = 1, numC
    chi_cty(i,1) = param_vec(ind+i)
end do
ind = ind + numC
do k = 1, numK
    chi_ind(k,1) = param_vec(ind+k)
end do
ind = ind + numK
do i = 1, numC
    do k = 1, numK
        chi(k,i) = chi_cty(i,1) + chi_ind(k,1)
    end do
end do

do i = 1, numC
    do k = 1, numK
        eta(k,i) = param_vec(ind+(i-1)*numK+k)
    end do
end do
ind = ind + numK*numC

C_US = 0
do k = 1, numK
    do l = 1, numK
        i2 = (k-1)*numK + l 
        C_US(k,l) = param_vec(ind+i2)
    end do
end do

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Assign parameters to structure params
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params % b           = b 
params % kappa_tilde = kappa_tilde
params % sigma_x     = sigma_x 
params % chi         = chi 
params % eta         = eta 
params % C_US        = C_US

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Check Restrictions on Parameters
! chi's must all be >= 0
! if not, set flag_param = 1, abort and 
! highly penalize that specific parameter
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_param = 0
do k = 1, numK
    do i = 1, numC
        if ( chi(k,i) < 0 ) then
            flag_param = 1
        end if
    end do
end do

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Continue with equilibrium computation
! only if flag_param = 0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (flag_param == 0) then
    
    ! Call routine that solves for the equilibrium
    ! flag_abort = 1 if a given constellation of parameters cannot lead to an equilibrium
    ! Dev is the deviation from equilibrium condition -- positive value of Dev needs to
    ! be penalized in the Loss Function
    Call computeEQ(params,Dev,MomentsModel,flag_abort,flag_abort2)
    
    if (flag_abort .eq. 0 .and. flag_abort2 .eq. 0) then
    
        ! Moments generated by the model
        s_tilde     = MomentsModel % s_tilde
        ShareEmp    = MomentsModel % ShareEmp 
        w_bar       = MomentsModel % w_bar 
        Var_w       = MomentsModel % Var_w  
        Unemp       = MomentsModel % Unemp 
        u           = MomentsModel % u 
        big_L       = MomentsModel % big_L 
        labor_share = MomentsModel % labor_share 
        x_bar       = MomentsModel % x_bar
        w_tilde     = MomentsModel % w_tilde
        q_theta     = MomentsModel % q_theta
        theta       = MomentsModel % theta
        G_x_bar     = MomentsModel % G_x_bar
        dist_out    = MomentsModel % dist_out
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute Economy-Wide Coefficient of Variation
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ! Country-wide Variances
        ! Var(Y) = E(Var(Y|X)) + Var(E(Y|X))
        ! Var(w) = sum(j in 1..numK) Var_w(k)*Prob(Sector=k) + 
        !          sum(j in 1..numK) ((w_bar(k)-mean_w)^2)*Prob(Sector=k)
        mean_w = 0.0
        do k = 1, numK
            mean_w = mean_w + w_bar(k)*ShareEmp(k)
        end do

        Var_w_tot = 0
        do k = 1, numK
            Var_w_tot = Var_w_tot + Var_w(k)*ShareEmp(k) + &
                        ((w_bar(k)-mean_w)**2)*ShareEmp(k)
        end do
        
        CoeffVar_w = sqrt(Var_w_tot)/mean_w
        
        ! Wages are relative to sector 1 wages
        ! We were having difficulty matching the level of wages
        ! Because the model implied value of big_E for the US
        ! was quite different from VA in the data
        do k = 1, numK
            w_bar_norm(k)      = w_bar(k) / w_bar(1)
            w_bar_Data_norm(k) = w_bar_Data(k,COUNTRY_EST) / w_bar_Data(1,COUNTRY_EST)
        end do
    
        !% Annual Transitions

        ! s_tilde to the fourth power
        s_tilde_4(:,:) = 0
        do k = 1, 2*numK
            s_tilde_4(k,k) = 1
        end do
        do j = 1, 4
            s_tilde_4(:,:) = matrix_mult(s_tilde_4(:,:),s_tilde(:,:),2*numK,2*numK,2*numK)
        end do
    
        !% Transitions including unconditional unemployment
        !% unconditional unemployment is "sector 1"
    
        s_tilde_annual = 0
    
        !% Transitions employment-employment
        
        s_tilde_annual(2:numK+1,2:numK+1) = s_tilde_4(numK+1:2*numK,numK+1:2*numK)
    
        !% Transitions to unemployment
    
        do k = 1, numK
            sum_s = 0
            do l = 1, numK
                sum_s = sum_s + big_L(l)*u(l)*s_tilde_4(l,numK+k) 
            end do
            s_tilde_annual(1,1+k) = sum_s / sum(big_L*u) 
        end do
        s_tilde_annual(1,1) = 1 - sum(s_tilde_annual(1,2:numK+1)) 
 
        !% Transitions to unemployment
    
        do k = 1, numK
            sum_s = 0 
            do l = 1, numK
                sum_s = sum_s + s_tilde_4(numK+k,numK+l)
            end do
            s_tilde_annual(1+k,1) = 1 - sum_s 
        end do
    
    end if
    
end if

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Weights on different sets of moments
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (COUNTRY_EST == 1) then
    ! If we estimate the US (COUNTRY_EST=1), then we match US transition rates
    ! and the US coefficient of variation
    WeightTr        = 1
    WeightCoeffVarW = 1
else
    WeightTr        = 0
    WeightCoeffVarW = 0
end if

WeightU         = 1
WeightSh        = 1
WeightW         = 1
WeightVacCosts  = 1

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%
! Loss Function
! %%%%%%%%%%%%%

if (flag_param == 0 .and. flag_abort == 0 .and. flag_abort2 == 0) then
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Penalize Devations from Equilibrium Conditions
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !% Penalize parameter values that do not lead to full convergence, outer loop
    toler_out = 0.0005
    penalty1  = 1000000.0*max(dist_out - toler_out,0.0) 

    !% Penalize deviations from equilibrium, even if convergence is achieved
    penalty2 = 1000000.0*max(Dev - toler_out,0.0)
    
    VacancyCosts_Data = 0.58
    VacancyCosts      = kappa_tilde(:,COUNTRY_EST)*w_tilde / w_bar
    
    s_tilde_Data_diag   = s_tilde_Data
    s_tilde_annual_diag = s_tilde_annual
    do k = 2, numK+1
        s_tilde_Data_diag(k,k)   = 1 - s_tilde_Data(k,k)
        s_tilde_annual_diag(k,k) = 1 - s_tilde_annual(k,k)
    end do
    
    if (COUNTRY_EST == 1) then
        F_Tr = sum( WeightTr*abs( (s_tilde_annual_diag - s_tilde_Data_diag) / abs(s_tilde_Data_diag) ) )
    else
        F_Tr = 0
    end if
    
    if (COUNTRY_EST == 1) then
        F_CoeffVarW = WeightCoeffVarW*abs( (CoeffVar_w - CoeffVar_w_Data) / (abs(CoeffVar_w_Data)) )
    else
        F_CoeffVarW = 0
    end if
    
    F_U  = WeightU*abs( ( Unemp - Unemp_Data(COUNTRY_EST,1) ) / (abs(Unemp_Data(COUNTRY_EST,1))) )

    F_Sh = sum( WeightSh*abs( (ShareEmp - ShareEmp_Data(:,COUNTRY_EST)) / (abs(ShareEmp_Data(:,COUNTRY_EST))) ) ) 
    
    F_W  = sum( WeightW*abs( (w_bar - w_bar_Data(:,COUNTRY_EST)) / (abs(w_bar_Data(:,COUNTRY_EST)) ) ) )
    
    F_VacCosts = sum( WeightVacCosts*abs( VacancyCosts - VacancyCosts_Data ) / &
                                     abs(VacancyCosts_Data) )
    
    F_p0 = sum( penalty0 )
    F_p1 = penalty1
    F_p2 = penalty2
    
    F_p3 = 100000.0*sum( max( (1 - G_x_bar) - 0.99 , 0.0 ) )
    
    F = F_Tr + F_U + F_Sh + F_W + F_CoeffVarW + F_VacCosts + &
        F_p0 + F_p1 + F_p2 + F_p3
    
else
    
    ! we get here if flag_param  = 1 or 
    !                flag_abort  = 0 or 
    !                flag_abort2 = 1
    F = 100000000.0 
    
end if

if ( IsNaN(F) ) then
    F = 100000000.0    
end if

if (flag_abort2 .eq. 1) then
    Call Print_Point_Abort
    Call Abort
end if

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Finish Time
call date_and_time(values=time_array_1)
finish_time = time_array_1(5) * 3600 + time_array_1(6) * 60 &
            + time_array_1(7) + 0.001 * time_array_1(8) 

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( FUNCTION_ITERATION == 1 .or. mod(FUNCTION_ITERATION,1) == 0 ) then
    print*
    print*, '****************************************'
    print*, 'Iteration : ', FUNCTION_ITERATION
    print*, 'Loss Function = ', F
    print*
    print*, 'Iteration Time     :', real(finish_time - start_time,4)
    print*, '****************************************'
    print*
end if

FUNCTION_ITERATION = FUNCTION_ITERATION + 1

if (F < F_BEST) then
    ! Print current point, current F and penalties 
    ! to 'FunctionEvaluations.csv'
    Call PrintCurrentPoint
    ! Save best point so far
    Call Print_BEST_Point
    ! Save comparison between data and model 
    ! moments at the best point so far
    Call PrintMomentComparison
    F_BEST = F
    Call PrintF_BEST
end if

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Contains

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints model vs data moments
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine PrintMomentComparison
    
    folder = ''
    filename = ''
    
    write (folder, "(A7,I1,A1)") "Country", COUNTRY_EST, "\"
    write (filename, "(A16,I1)") "MomentComparison", COUNTRY_EST
        
    open(unit = 100 , file = trim(folder)//trim(filename)//'.csv')
    
    write(100,*) 'Moment, Model, Data, Weight'

    write(100,345) 'Unemployment', ',', Unemp, ',', Unemp_Data(COUNTRY_EST,1), ',', WeightU

    do k = 1, numK
        write(100,346) 'Share Emp Sector ', k, ',', ShareEmp(k), ',', ShareEmp_Data(k,COUNTRY_EST), ',', WeightSh(k)
    end do
    
    do k = 1, numK
        write(100,347) 'Average Wage Sector ', k, ',', w_bar(k), ',', w_bar_Data(k,COUNTRY_EST), ',', WeightW(k)
    end do

    do k = 1, numK
        write(100,348) 'Vacancy Costs Sector ', k, ',', VacancyCosts(k), ',', VacancyCosts_Data(k), ',', WeightVacCosts(k)
    end do
    
    if (COUNTRY_EST == 1) then
        
        write(100,345) 'Coeff Var w', ',', CoeffVar_w, ',', CoeffVar_w_Data, ',', WeightCoeffVarW
        
        do i = 1, (numK+1)
            do j = 1, (numK+1)
                write(100,350) 'Transition Rates US ', i, ' to ', j, ',', s_tilde_annual_diag(i,j), ',' , s_tilde_Data_diag(i,j), ',', WeightTr(i,j)
            end do
        end do
        
    end if
    
    close(100)

345 format(a12,3(a1,f22.8))
346 format(a17,i1,3(a1,f22.8))
347 format(a20,i1,3(a1,f22.8))
348 format(a21,i1,3(a1,f22.8))
349 format(a11,3(a1,f22.8))
350 format(a20,i1,a4,i1,3(a1,f22.8))
    
end subroutine PrintMomentComparison

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prints the current point, current F and penalties to
! 'FunctionEvaluations.csv'
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine PrintCurrentPoint

    write(999,201) F, sum( penalty0 ), penalty1, penalty2, penalty3, (param_vec(i), i = 1, NP)
    
201 format(5(f30.4,','),167(f22.8,','))
    
end subroutine PrintCurrentPoint

subroutine Print_BEST_Point

    folder = ''
    filename = ''
    
    !write (folder, "(A7,I1,A1)") "Country", COUNTRY_EST, "\"
    !write (filename, "(A19,I1)") "Starting_Point_BEST", COUNTRY_EST
    !
    !open(unit = 710 , file =  trim(folder)//trim(filename)//'.csv')

    open(unit = 710 , file =  'Starting_Point_BEST.csv')
    
    do i = 1, NP
        write(710,388) i, param_vec(i), lb_global(i), ub_global(i), scaling_global(i), mask_global_in(i)
    end do
    
    close(710)
    
388 format(i4,',',4(f20.8,','),i4)
    
end subroutine Print_BEST_Point

subroutine Print_Point_Abort

    open(unit = 710 , file = 'Starting_Point_Abort.csv')

    do i = 1, NP
        write(710,388) i, param_vec(i), lb_global(i), ub_global(i), scaling_global(i), mask_global(i)
    end do
    
    close(710)
    
388 format(i4,',',4(f20.8,','),i4)
    
end subroutine Print_Point_Abort

Subroutine PrintF_BEST

write (folder, "(A7,I1,A1)") "Country", COUNTRY_EST, "\"
write (filename, "(A6,I1)") "F_BEST", COUNTRY_EST
open(unit=502,file=trim(folder)//trim(filename)//".csv")

    write(502,*) F
    
close(502)

end Subroutine PrintF_BEST

end subroutine SMM_OBJ_FCN

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This routine computes the equilibrium at parameter
! values given by the structure param
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ComputeEQ(params,Dev,MomentsModel,flag_abort,flag_abort2)

USE Global_Data
USE LinearAlgebra_MOD

implicit none

type(Parameters) , intent(in)  :: params
real(KIND=DOUBLE), intent(out) :: Dev
type(MomentsOut) , intent(out) :: MomentsModel
integer          , intent(out) :: flag_abort, flag_abort2

! maximum # of iterations -- big_U loop
integer, parameter :: maxiter_U   = 50000
! maximum # of iterations -- outer loop
integer, parameter :: maxiter_out = 10000

! convergence tolerance -- big_U loop
real(KIND=DOUBLE), parameter :: toler_U   = 0.00000001
! convergence tolerance -- outer loop
real(KIND=DOUBLE), parameter :: toler_out = 0.000001

! initial steps for updating equilibrium variables
real(KIND=DOUBLE), parameter :: step_L_ini = 0.01
real(KIND=DOUBLE), parameter :: step_x_ini = 0.01

real(KIND=DOUBLE), parameter :: eps = 0.02

real(KIND=DOUBLE) b(numK), chi(numK), kappa_tilde(numK), &
                  sigma_x(numK), eta(numK), C_US(numK,numK), C(numK,numK)

real(KIND=DOUBLE) x_max, x_min, Int_0(numK), varpi(numK), &
                  q_theta_0(numK), test_m(numK), test, x_bar_ub(numK), x_bar0(numK), &
                  big_L0(numK), step_L, step_x, x_bar(numK), x_bar1(numK), big_L(numK), &
                  Int_x_bar(numK), G_x_bar(numK), q_theta(numK), &
                  theta(numK), u(numK), p(numK), Int_2(numK), big_L_tilde(numK), &
                  term1, term2, w_tilde(numK), &
                  big_E_V(numK), big_E_C, big_L_Temp(numK), &
                  lambda_tilde(numK), lambda_tilde_aux

real(KIND=DOUBLE) big_U(numK), max_big_U, big_U_new(numK), C_k(numK), MatAux1(numK)

real(KIND=DOUBLE) s(numK,numK), eye_K(numK,numK), eigenval(numK,1), &
                  eigenvec(numK,numK), AA(numK,numK), y_tilde(numK), phi0

integer           k, i, l, min_loc(1,1), info

real(KIND=DOUBLE) dist_out_vec(maxiter_out,1), &
                  x_bar_best(numK), big_L_best(numK), &
                  dist_out, dist_U, big_L1(numK), &
                  dist_x, dist_L, dist_out_best, mean_dist1, mean_dist2, LHS(numK), RHS(numK)

integer           flag_osc, flag_U, iter_out, iter_U

real(KIND=DOUBLE) Emp_Count(numK), ShareEmp(numK), Unemp, &
                  w_bar(numK), norm_term1(numK), norm_term2(numK), norm_term3(numK), &
                  Var_w(numK), labor_share(numK), s_tilde(2*numK,2*numK)

real(KIND=DOUBLE) spower_old(numK,numK), spower_new(numK,numK), dist_spower
integer           iter_spower

character(len=1024) folder
character(len=1024) filename1, filename2, filename3, filename4, filename5, &
                    filename6, filename7

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Flag to abort the computation of the loss function due inadmissible parameter
flag_abort  = 0 
! Flag to abort the computation of the loss function due to NaN's in big_U
flag_abort2 = 0

! x_max and x_min for the bissection solver finding x_bar_ub
x_max = 100 
x_min = 0 

!%%%%%%%%%%%%%%%%%%%%%%
!% Parameter Assignment
!%%%%%%%%%%%%%%%%%%%%%%

b           = params % b(:,COUNTRY_EST) 
kappa_tilde = params % kappa_tilde(:,COUNTRY_EST)
sigma_x     = params % sigma_x(:,COUNTRY_EST) 
chi         = params % chi(:,COUNTRY_EST) 
eta         = params % eta(:,COUNTRY_EST)
C_US        = params % C_US 

do k = 1, numK
    do l = 1, numK
        C(k,l) = C_US(k,l)
    end do
end do

!%%%%%%%%%%%%%%
!% Step 3
!%%%%%%%%%%%%%%

!% Compute Integral I_ki(x) at x = 0
Int_0     = exp( (sigma_x**2)/2 ) 
varpi     = ( ( 1 - (1-chi)*delta )*kappa_tilde )/( delta*(1-beta) ) 
!% Compute q(theta) evaluated at xmin = 0
q_theta_0 = varpi/Int_0 
!% tests for whether q(theta) < 1-eps for all k,i at xmin = 0
!% If test > 0 then, abort
test_m = -(q_theta_0 >= 1-eps)
test = abs( sum( test_m ) ) 

if (test > 0) then
    flag_abort                 = 1 
    Dev                        = 10000000000.0 
    MomentsModel % s_tilde     = -999 
    MomentsModel % w_bar       = -999 
    MomentsModel % Var_w       = -999 
    MomentsModel % labor_share = -999  
    MomentsModel % Unemp       = -999 
    MomentsModel % ShareEmp    = -999 
    MomentsModel % x_bar       = -999 
end if

if (flag_abort == 0) then
    
!%%%%%%%%%%%%%%
!% Step 4
!%%%%%%%%%%%%%%

!% Find x_bar_ub using bissection method
!% Solve for I_ki(x_bar_ub) = varpi
do k = 1, numK
    x_bar_ub(k) = Find_x_bar_ub(params, varpi, k, x_min, x_max)
end do

!% Initializations for x_bar, bil_L, big_U
if (FUNCTION_ITERATION == 1) then
    x_bar_Guess = x_bar_ub / 10
    big_L_Guess = mu(:,COUNTRY_EST)*big_L_bar(:,1)
end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Step 5: guess big_L x_bar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_bar0 = min( max( x_bar_Guess , 0.000001 ) , x_bar_ub - 0.000001 ) 
big_L0 = big_L_Guess

!% Initialize steps
step_L = step_L_ini 
step_x = step_x_ini 

!% Initialize endogenous variables
x_bar   = x_bar0 
big_L   = big_L0 

dist_out_vec  = 10 
x_bar_best    = x_bar0 
big_L_best    = big_L0 
flag_osc      = 0  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist_out = 1 
iter_out = 1 
do while ( dist_out > toler_out .and. iter_out <= maxiter_out )

    ! If the maximum number of iterations is reached, we
    ! compute all outcomes using the best guess of x_bar / big_L
    if (iter_out == maxiter_out) then
        big_L = big_L_best
        x_bar = x_bar_best
    end if
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 6: Compute Int_x_bar, G_x_bar, theta, u
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !% Compute Integral(x_bar)
    do k = 1, numK
        Int_x_bar(k) = Integral(params,k,x_bar(k)) 
    end do
    
    !% Compute G(x_bar)
    do k = 1, numK
        ! G_x_bar = normcdf(log(x_bar)./sigma_x) 
        call VDCDFNORM(1,log(x_bar(k))/sigma_x(k),G_x_bar(k))
    end do

    !% Compute q(theta)
    q_theta = min(varpi/Int_x_bar,0.99999)

    !% Compute theta
    theta = ( q_theta**(-xsi) - 1 )**(1/xsi) 
    
    !% Compute u
    u = chi / ( theta*q_theta*(1-G_x_bar) + chi )

    !% Compute p
    p = 1 - theta*q_theta*(1-G_x_bar)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 7: Compute big_L_tilde
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    !% Compute big_L_tilde
    do k = 1, numK
        call VDCDFNORM(1,sigma_x(k) - log(x_bar(k))/sigma_x(k),term1)
        call VDCDFNORM(1,-log(x_bar(k))/sigma_x(k),term2)
        Int_2(k) = exp( (sigma_x(k)**2)/2 )*term1 / term2
    end do
    big_L_tilde = big_L*(1-u)*Int_2

    !%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 8: Compute w_tilde
    !%%%%%%%%%%%%%%%%%%%%%%%%%
    
    w_tilde = gamma(:,COUNTRY_EST)*big_Y(:,COUNTRY_EST) / big_L_tilde
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 9: Compute big_E_V
    !%%%%%%%%%%%%%%%%%%%%%%%%%

    big_E_V   = kappa_tilde*w_tilde*theta*u*big_L
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 10: Compute big_E_C
    !%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    big_E_C = sum( gamma(:,COUNTRY_EST)*big_Y(:,COUNTRY_EST) ) - sum( big_E_V ) - NX_Data(COUNTRY_EST,1)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 11: Compute lambda_tilde
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lambda_tilde_aux = big_L_bar(COUNTRY_EST,1) / big_E_C
    lambda_tilde     = lambda_tilde_aux
     
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 12: Solve Bellman Equation
    !% and obtain big_U
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    big_U  = big_U_Guess 
    dist_U = 1 
    iter_U = 1 
    flag_U = 0 
    do while ( dist_U > toler_U .and. iter_U < maxiter_U )

        do i = 1, numC
            max_big_U = maxval(big_U)
        end do
        
        do k = 1, numK
            C_k = C(k,:)
            MatAux1 = exp( ( -C_k + b + theta*kappa_tilde*lambda_tilde*w_tilde*(beta/(1-beta)) + delta*big_U &
                    - delta*max_big_U)/zeta )
            big_U_new(k) = zeta*log( sum(MatAux1) ) + delta*max_big_U
        end do
        
        if ( sum(abs(big_U)) > 0 ) then
            dist_U = maxval( abs( (big_U_new - big_U)/abs(big_U) ) ) 
        else
            dist_U = maxval( abs( big_U_new - big_U ) ) 
        end if

        iter_U = iter_U + 1 

        big_U = big_U_new 

    end do 
    
    ! check if at least one element of big_U is NaN or infinity 
    ! inifinity + 1000.0 is equal to infinity
    if (any(isNaN(big_U)) .or. any(big_U + 1000.0 .eq. big_U)) then
        flag_abort2 = 1
        exit
    else
        big_U_Guess = big_U
    end if
    
    !%%%%%%%%%%%%%%%%%%%%%%
    !% Step 7: Update big_L
    !%%%%%%%%%%%%%%%%%%%%%%

    !% Compute transition rates s
    do k = 1, numK
        C_k = C(k,:)
        MatAux1 = exp( ( -C_k + b + theta*kappa_tilde*lambda_tilde*w_tilde*(beta/(1-beta)) + delta*big_U &
                - delta*max_big_U )/zeta )
        s(k,:) = MatAux1 / sum( MatAux1 )
    end do
    
    eye_K = 0.0
    do i = 1, numK
        eye_K(i,i) = 1
    end do
    
    do i = 1, numC
        AA = eye_K - transpose(s(:,:))
        Call eigen(AA,numK,eigenval,eigenvec,info)
        ! If "info" different from zero, we compute
        ! nullspace of I - s(:,:,i) using the stationary
        ! distribution of s(:,:,i)
        ! info = 1
        if (info == 0) then 
            min_loc(1:1,1) = minloc(abs(eigenval(:,1)))
        else
            min_loc(1:1,1) = 1
        end if
        info = 1
        ! In case there are problems computing eigenvalues and eigenvectors,
        ! that is, if info != 0 or abs(min(eigenval)) is not so small...
        ! We obtain s(:,:,i)^N N-->inf to obtain the stationary distribution
        ! associated with s(:,:,i) 
        if (info == 0 .and. abs(eigenval(min_loc(1,1),1)) < 1e-10) then
            y_tilde(:) = eigenvec(:,min_loc(1,1)) / u(:)
        else
            ! spower_old = identity
            spower_old = 0
            do k = 1, numK
                spower_old(k,k) = 1
            end do
            dist_spower = 1 
            iter_spower = 1 
            do while( dist_spower > 0.000000000001 .and. iter_spower <= 2000 )
                spower_new = matrix_mult(s(:,:),spower_old,numK,numK,numK)
                !dist_spower = norm2(spower_new - spower_old)/max(norm2(spower_old),0.000001)
                dist_spower = maxval( (spower_new - spower_old)/max(abs(spower_old),0.00000001) )
                spower_old = spower_new 
                iter_spower = iter_spower + 1
            end do
            y_tilde(:) = spower_new(1,:) / u(:)
        end if
        phi0   = big_L_bar(COUNTRY_EST,1) / sum( y_tilde )
        big_L1 = phi0*y_tilde
    end do
    
    big_L1 = max(big_L1,0.0)
    
    !%%%%%%%%%%%%%%%%%%%%%%
    !% Step 8: Update x_bar
    !%%%%%%%%%%%%%%%%%%%%%%

    x_bar1 = ( (1-delta)*big_U - eta ) / (lambda_tilde*w_tilde)
    x_bar1 = max( min( x_bar1 , x_bar_ub-0.000001 ) , 0.000001 )
    
    dist_L = maxval( abs(( big_L1 - big_L)/(max(abs(big_L),0.0000001)) ) )
    dist_x = maxval( abs((x_bar1 - x_bar) /(max(abs(x_bar),0.0000001)) ) )
    
    !% error is the max between dist_x and dist_L
    dist_out = max(dist_x,dist_L) 
    
    if (iter_out == 1) then
        dist_out_best = dist_out
    end if
    
    !% Saves full path of dist_out
    dist_out_vec(iter_out,1) = dist_out
    
    !% Saves x_bar, big_L, x_bar1, and big_L1 
    !% leading to the best dist_out
    if (dist_out < dist_out_best) then
        dist_out_best = dist_out 
        x_bar_best    = x_bar 
        big_L_best    = big_L 
    end if
    
    !% Allow 300 iterations before checking for bad behavior
    if (iter_out > 300) then
        !% If dist_out does not improve, or if it improves very little (by less than toler_out/10),
        !% reduce step size and go back to x_bar and L that led to 
        !% the smallest value of dist_out
        if ( flag_osc == 0 .or. iter_out >= flag_osc + 300 ) then
            mean_dist1 = sum(dist_out_vec(iter_out-49:iter_out,1))/50 
            mean_dist2 = sum(dist_out_vec(iter_out-99:iter_out-51,1))/50  
        end if
        if ( ( (mean_dist1 >= mean_dist2) .or. &
             ( (mean_dist1 - mean_dist2 < 0) .and. abs(mean_dist1 - mean_dist2) < toler_out ) ) .and. &
             flag_osc == 0) then
            step_L   = 0.9*step_L 
            step_x   = 0.9*step_x 
            flag_osc = iter_out 
        end if
        !% But we only allow step sizes to be reduced every 300 iterations
        if ( ( (mean_dist1 >= mean_dist2) .or. &
             ( (mean_dist1 - mean_dist2 < 0) .and. abs(mean_dist1 - mean_dist2) < toler_out ) ) .and. &
             iter_out >= flag_osc+300 ) then
            step_L   = 0.9*step_L 
            step_x   = 0.9*step_x 
            flag_osc = iter_out 
        end if
    end if
        
    ! We do NOT allow big changes in big_L
    big_L_Temp = (1-step_L)*big_L + step_L*big_L1 
    big_L      = min(max(big_L_Temp,.9*big_L),1.1*big_L)
    ! We need to renormalize if the bounds above are binding 
    ! Otherwise big_L does not sum to big_L_bar
    do i = 1, numC
       big_L = big_L/sum(big_L)*big_L_bar(COUNTRY_EST,1)
    end do
    
    x_bar = (1-step_x)*x_bar + step_x*x_bar1 
    x_bar = max( min( x_bar , x_bar_ub - 0.000001 ) , 0.000001 )
    
    !% If step size was adjusted after bad behavior, 
    !% restart with x_bar = x_bar_best 
    !% and big_L = big_L_best instead   
    if (flag_osc == iter_out) then
        big_L = big_L_best 
        x_bar = x_bar_best 
    end if
    
    iter_out = iter_out + 1 

end do

if (flag_abort2 == 0) then
    
! Print Initial Steady State

write (folder, "(A7,I1,A1)") "Country", COUNTRY_EST, "\"

write (filename1, "(A15,I1)") "Initial_SS_xbar", COUNTRY_EST
open(unit = 1 , file = trim(folder)//trim(filename1)//'.csv')
    do k = 1, numK
        write(1,101) x_bar(k)
    end do
close(1)

write (filename2, "(A18,I1)") "Initial_SS_w_tilde", COUNTRY_EST
open(unit = 2 , file = trim(folder)//trim(filename2)//'.csv')
    do k = 1, numK
        write(2,101) w_tilde(k)
    end do
close(2)

write (filename3, "(A16,I1)") "Initial_SS_big_L", COUNTRY_EST
open(unit = 3 , file = trim(folder)//trim(filename3)//'.csv')
    do k = 1, numK
        write(3,101) big_L(k)
    end do
close(3)

write (filename4, "(A16,I1)") "Initial_SS_big_U", COUNTRY_EST
open(unit = 4 , file = trim(folder)//trim(filename4)//'.csv')
    do k = 1, numK
        write(4,101) big_U(k)
    end do
close(4)

write (filename5, "(A12,I1)") "Initial_SS_u", COUNTRY_EST
open(unit = 5 , file = trim(folder)//trim(filename5)//'.csv')
    do k = 1, numK
        write(5,101) u(k)
    end do
close(5)

write (filename6, "(A12,I1)") "Initial_SS_s", COUNTRY_EST
open(unit = 6 , file = trim(folder)//trim(filename6)//'.csv')
    do k = 1, numK
        write(6,102) (s(k,l),l=1,numK)
    end do
close(6)

write (filename7, "(A16,I1)") "Initial_SS_theta", COUNTRY_EST
open(unit = 7 , file = trim(folder)//trim(filename7)//'.csv')
    do k = 1, numK
        write(7,101) theta(k)
    end do
close(7)

101 format(6(f50.16,','))
102 format(6(f50.16,','))
    
x_bar_Guess = x_bar
big_L_Guess = big_L

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Check if Equation (94) is met
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LHS = lambda_tilde*w_tilde*x_bar 
RHS = (1-delta)*big_U - eta 

!% Deviations from equilibrium conditions
Dev = maxval( abs(LHS - RHS)/(max(abs(LHS),0.000001)) ) 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Model Simulated Moments
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%
!% Employment Shares
!%%%%%%%%%%%%%%%%%%%

Emp_Count = big_L*(1-u) 
ShareEmp  = Emp_Count / sum(Emp_Count) 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Average unemployment by country
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Unemp = sum( big_L*u ) / sum(big_L)

!%%%%%%%%%%%%%%%
!% Average Wages
!%%%%%%%%%%%%%%%

do k = 1, numK
    Call VDCDFNORM(1,sigma_x(k) - log(x_bar(k))/sigma_x(k),norm_term1(k))
end do

w_bar = (1-beta)*w_tilde*x_bar + beta*w_tilde* &
        exp( (sigma_x**2)/2 )*norm_term1/(1-G_x_bar)

!%%%%%%%%%%%%%%%%%%%
!% Variance of Wages
!%%%%%%%%%%%%%%%%%%%

do k = 1, numK
    Call VDCDFNORM(1,2*sigma_x(k) - log(x_bar(k))/sigma_x(k),norm_term2(k))
    Call VDCDFNORM(1,-log(x_bar(k))/sigma_x(k),norm_term3(k))
end do

Var_w = ( ( beta*w_tilde )**2 )* &
        ( exp(2*sigma_x**2)*norm_term2/norm_term3 - &
          exp(sigma_x**2)*(norm_term1/norm_term3)**2 )


!%%%%%%%%%%%%%
!% Labor Share
!%%%%%%%%%%%%%

labor_share = ( w_bar*Emp_Count )/ big_Y(:,COUNTRY_EST)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% 1-period transition rates
!%%%%%%%%%%%%%%%%%%%%%%%%%%%

!% indices 1 through numK stand for unemployment in 1 through unemployment in K
!% indices K+1 through 2*numK stand for employment in sector 1 to sector K

s_tilde = 0

!% unemployment to unemployment
do l = 1, numK
    do k = 1, numK
        s_tilde(l,k) = s(l,k)*(1 - theta(k)*q_theta(k)*(1-G_x_bar(k)))
    end do
end do

!% unemployment to employment
do l = 1, numK
    do k = 1, numK
        s_tilde(l,numK+k) = s(l,k)*theta(k)*q_theta(k)*(1-G_x_bar(k)) 
    end do 
end do

!% employment to employment
do l = 1, numK
    do k = 1, numK
        if ( l == k ) then
            s_tilde(numK+l,numK+k) = (1-chi(l)) 
        else
            s_tilde(numK+l,numK+k) = 0 
        end if
    end do
end do

!% employment to unemployment
do l = 1, numK
    do k = 1, numK
        if ( l == k ) then
            s_tilde(numK+l,k) = chi(l)
        else
            s_tilde(numK+l,k) = 0
        end if
    end do
end do

MomentsModel % s_tilde      = s_tilde 
MomentsModel % w_bar        = w_bar 
MomentsModel % Var_w        = Var_w 
MomentsModel % labor_share  = labor_share 
MomentsModel % u            = u 
MomentsModel % big_L        = big_L
MomentsModel % Unemp        = Unemp 
MomentsModel % ShareEmp     = ShareEmp 
MomentsModel % lambda_tilde = lambda_tilde 
MomentsModel % x_bar        = x_bar 
MomentsModel % w_tilde      = w_tilde
MomentsModel % q_theta      = q_theta
MomentsModel % theta        = theta
MomentsModel % G_x_bar      = G_x_bar
MomentsModel % dist_out     = dist_out

end if

end if

end subroutine computeEQ

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Defines I_ki(x) function
!%%%%%%%%%%%%%%%%%%%%%%%%%%

function Integral(params,k,x) result(Int)

    USE Global_Data
    
    implicit none
    
    type(Parameters) params
    integer, intent(in) :: k
    real(KIND=DOUBLE), intent(in) :: x
    
    real(KIND=DOUBLE) sigma_x(numK), x_new, term1, term2, Int
    
    sigma_x = params % sigma_x(:,COUNTRY_EST) 
    x_new   = max(x,0.000001) 
    call VDCDFNORM(1,sigma_x(k) - log(x_new)/sigma_x(k),term1)
    call VDCDFNORM(1,-log(x_new)/sigma_x(k),term2)
    Int     = exp( (sigma_x(k)**2) / 2 )*term1 &
            - x_new*term2

end function Integral


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Solves I_ki(x) = varpi using a simple bissection method
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Find_x_bar_ub(params, varpi, k, xmin_in, xmax_in) result(xout)

    USE Global_Data
    USE LinearAlgebra_MOD

    implicit none
    
    type(Parameters), intent(in)   :: params    
    real(KIND=DOUBLE) , intent(in) :: varpi(:)
    integer          , intent(in)  :: k
    real(KIND=DOUBLE), intent(in)  :: xmin_in, xmax_in
    
    real(KIND=DOUBLE) xout

    integer, parameter :: itermax_bisec = 1000 
    integer, parameter :: toler_bisec = 0.000001 
    
    real(KIND=DOUBLE) x, xmin, xmax, f, fmin, fmax
    
    integer iter
    
    xmin = xmin_in 
    xmax = xmax_in 
    
    fmin = Integral(params,k,xmin) - varpi(k) 
    fmax = Integral(params,k,xmax) - varpi(k) 
    
    iter = 1 
    
    !% Adjust xmin and xmax until they lead to opposite signs, if necessary
    do while ( fmin*fmax > 0.0 .and. iter <= itermax_bisec ) 
        
        xmin = xmin / 1.1 
        xmax = xmax * 1.1 
        
        fmin = Integral(params,k,xmin) - varpi(k) 
        fmax = Integral(params,k,xmax) - varpi(k) 
        
        iter = iter + 1 
        
    end do
    
    iter = 1 
    
    if (fmin < 0 .and. fmax > 0) then
        
        do while ( abs(fmin - fmax) > toler_bisec .and. iter <= itermax_bisec )
    
            x = (xmin + xmax) / 2.0 
            f = Integral(params,k,x) - varpi(k) 
            
            if (f > 0) then
                xmax = x 
                fmax = f 
            else
                xmin = x 
                fmin = f 
            end if
            
            iter = iter + 1 
        
        end do
    
    elseif (fmin > 0 .and. fmax < 0) then
    
        do while ( abs(fmin - fmax) > toler_bisec .and. iter <= itermax_bisec )
    
            x = (xmin + xmax) / 2.0 
            f = Integral(params,k,x) - varpi(k) 
        
            if (f > 0) then
                xmin = x 
                fmin = f                
            else
                xmax = x 
                fmax = f 
            end if
            
            iter = iter + 1 
        
        end do
        
    end if
    
    xout = x 

end function Find_x_bar_ub

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end Module LossFunction_MOD