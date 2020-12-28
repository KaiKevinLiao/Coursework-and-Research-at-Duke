Module computeNewSteadyState_MOD

Contains
    
Subroutine NewSteadyState(param_scaled,Shock_in,SteadyState0,SteadyState1)

USE Global_Data
USE LinearAlgebra_MOD

implicit none

real(KIND=DOUBLE), intent(in)  :: param_scaled(:)
type(Shock)      , intent(in)  :: Shock_in
type(SteadyState), intent(in)  :: SteadyState0
type(SteadyState), intent(out) :: SteadyState1

! maximum # of iterations -- big_U loop
integer, parameter :: maxiter_U   = 50000
! maximum # of iterations -- outer loop
integer, parameter :: maxiter_out = 10000
! maximum # of iterations -- w_tilde loop
integer, parameter :: maxiter_w = 10000
! maximum # of iterations -- P_I_hat loop
integer, parameter :: maxiter_P = 10000
! maximum # of iterations -- big_Y1 loop
integer, parameter :: maxiter_Y = 10000

! convergence tolerance -- big_U loop
real(KIND=DOUBLE), parameter :: toler_U   = 1e-12
! convergence tolerance -- outer loop
real(KIND=DOUBLE), parameter :: toler_out = 1e-8
! convergence tolerance -- w_tilde loop
real(KIND=DOUBLE), parameter :: toler_w = 1e-12
! convergence tolerance -- P_I_hat loop
real(KIND=DOUBLE), parameter :: toler_P = 1e-12
! convergence tolerance -- big_Y1 loop
real(KIND=DOUBLE), parameter :: toler_Y = 1e-12

! initial steps for updating equilibrium variables
real(KIND=DOUBLE), parameter :: step_w_ini = 0.1
real(KIND=DOUBLE), parameter :: step_L_ini = 0.025
real(KIND=DOUBLE), parameter :: step_x_ini = 0.025

real(KIND=DOUBLE) param_vec(NP)

real(KIND=DOUBLE) b(numK,numC), chi(numK,numC), kappa_tilde(numK,numC), &
                  sigma_x(numK,numC), eta(numK,numC), C_US(numK,numK), C(numK,numK,numC)

real(KIND=DOUBLE) b_cty(numC,1), chi_cty(numC,1), chi_ind(numC,1)

real(KIND=DOUBLE) w_tilde0(numK,numC), pi0(numC,numC,numK), x_bar0(numK,numC), big_L0(numK,numC), &
                  kappa_tilde0(numK,numC)

real(KIND=DOUBLE) w_tilde1(numK,numC), pi1(numC,numC,numK), x_bar1(numK,numC), big_L1(numK,numC), &
                  big_L1_best(numK,numC), x_bar1_best(numK,numC), w_tilde_hat(numK,numC), &
                  c_hat(numK,numC), P_I_hat(numK,numC), P_I_hat_new(numK,numC), P_F_hat(numC,1), &
                  pi_hat(numC,numC,numK), kappa_tilde1(numK,numC), w_tilde1_new(numK,numC), &
                  big_L1_new(numK,numC), x_bar1_new(numK,numC), log_c_hat(numK,numC), log_P_I_hat(numK,numC)

real(KIND=DOUBLE) Int_x_bar1(numK,numC), G_x_bar1(numK,numC), q_theta1(numK,numC), theta1(numK,numC), &
                  u1(numK,numC), varpi(numK,numC), big_Y1(numK,numC), big_Y1_new(numK,numC), &
                  big_L_tilde1(numK,numC), term1, term2, big_E_V1(numK,numC), big_E_C1(numC,1), &
                  lambda_tilde1(numK,numC), big_U1(numK,numC), Int_2(numK,numC), C_k(numK,numC), &
                  MatAux1(numK,numC), big_U1_new(numK,numC), max_big_U1(1,numC), big_U1_Guess(numK,numC), &
                  lambda_tilde_aux(numC,1), MatAux2(1,numC), gamma_big_Y1(numK,numC), P_I_hat_Guess(numK,numC), &
                  big_Y1_Guess(numK,numC), w_tilde1_Guess(numK,numC)

real(KIND=DOUBLE) s(numK,numK,numC), eye_K(numK,numK), eigenval(numK,1), &
                  eigenvec(numK,numK), AA(numK,numK), y_tilde(numK,numC), &
                  phi0(1,1), y_tilde_tr(1,numK)

real(KIND=DOUBLE) spower_old(numK,numK), spower_new(numK,numK), dist_spower

integer           iter_spower, min_loc(1,1), info

real(KIND=DOUBLE) dist_w_vec(maxiter_w,1), dist_w_best, w_tilde_best(numK,numC), mean_dist_w1, mean_dist_w2, &
                  dist_out_best, dist_out_vec(maxiter_out,1), x_bar_best(numK,numC), big_L_best(numK,numC), &
                  mean_dist1, mean_dist2

integer           flag_osc
        
integer           flag_osc_w

real(KIND=DOUBLE) A_hat(numK,numC), d_hat(numC,numC,numK), NX1(numC,1)

real(KIND=DOUBLE) ones_K(numK,1), ones_K_tr(1,numK)

real(KIND=DOUBLE) dist_L, dist_x, dist_out, dist_w, dist_P, dist_Y, dist_U, &
                  step_L, step_x, step_w

integer iter_out, iter_w, iter_U, iter_P, iter_Y

integer i, j, k, l, o, i2, ind

! Elapsed time variables
integer           time_array_0(8), time_array_1(8)
real(KIND=DOUBLE) start_time, finish_time

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ones_K    = 1
ones_K_tr = transpose(ones_K)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! De-scale parameters
do i = 1, NP
    param_vec(i) = param_global(i) * scaling_global(i)
end do

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

do i = 1, numC
    do k = 1, numK
        do l = 1, numK
            C(k,l,i) = C_US(k,l)
        end do
    end do
end do

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi0          = SteadyState0 % pi
w_tilde0     = SteadyState0 % w_tilde
big_L0       = SteadyState0 % big_L
x_bar0       = SteadyState0 % x_bar

big_U1_Guess  = SteadyState0 % big_U
big_Y1_Guess  = SteadyState0 % big_Y
P_I_hat_Guess = 1.0 

kappa_tilde0 = kappa_tilde

A_hat    = Shock_in % A_hat
d_hat    = SHock_in % d_hat
NX1      = Shock_in % NX

A_hat(:,2) = 2

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Start time
call date_and_time(values=time_array_0)
start_time = time_array_0(5) * 3600 + time_array_0(6) * 60 &
           + time_array_0(7) + 0.001 * time_array_0(8)

! Initializations
big_L1         = big_L0
x_bar1         = x_bar0
w_tilde1_Guess = w_tilde0

dist_out = 1
iter_out = 1 
flag_osc = 0
step_x   = step_x_ini
step_L   = step_L_ini

do while ( dist_out > toler_out .and. iter_out <= maxiter_out )

    ! If the maximum number of iterations is reached, we
    ! compute all outcomes using the best guess of x_bar / big_L
    if (iter_out == maxiter_out) then
        big_L1 = big_L1_best
        x_bar1 = x_bar1_best
    end if
    
    dist_w     = 1
    iter_w     = 1
    flag_osc_w = 0
    step_w     = step_w_ini
    w_tilde1   = w_tilde1_Guess
    do while (dist_w > toler_w .and. iter_w <= maxiter_w)
        
        !open(unit = 999 , file = 'Check_w_tilde1.csv')
        !    do k = 1, numK
        !        write(999,876) (w_tilde1(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 2a: Compute w_tilde_hat and iteratively solve
        ! for P_I_hat and c_hat
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        w_tilde_hat = w_tilde1 / w_tilde0
        
        !open(unit = 999 , file = 'Check_w_tilde_hat.csv')
        !    do k = 1, numK
        !        write(999,876) (w_tilde_hat(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        dist_P      = 1
        iter_P      = 1
        P_I_hat     = P_I_hat_Guess
        log_P_I_hat = log(P_I_hat)
        do while (dist_P > toler_P .and. iter_P <= maxiter_P)
            
            do i = 1, numC
                do k = 1, numK
                    !c_hat(k,i) = w_tilde_hat(k,i)**gamma(k,i)
                    log_c_hat(k,i) = gamma(k,i)*log(w_tilde_hat(k,i))
                    do l = 1, numK
                        !c_hat(k,i) = c_hat(k,i)*( P_I_hat(l,i)**( (1-gamma(k,i))*nu(k,l,i) ) )
                        log_c_hat(k,i) = log_c_hat(k,i) + ( (1-gamma(k,i))*nu(k,l,i) )*log_P_I_hat(l,i)
                    end do
                end do
            end do
            
            c_hat = exp( log_c_hat )
            
            P_I_hat_new = 0
            do i = 1, numC
                do k = 1, numK
                    do o = 1, numC
                        P_I_hat_new(k,i) = P_I_hat_new(k,i) + pi0(o,i,k)*A_hat(k,o)*(c_hat(k,o)*d_hat(o,i,k))**(-lambda)
                    end do
                end do
            end do
            P_I_hat_new = P_I_hat_new**(-1/lambda)
            
            dist_P = maxval( abs( (P_I_hat_new - P_I_hat) / max(abs(P_I_hat),0.000001) ) )
        
            iter_P = iter_P + 1
            
            P_I_hat     = P_I_hat_new
            log_P_I_hat = log( P_I_hat )
            
        end do
        
        P_I_hat_Guess = P_I_hat
        
        !open(unit = 999 , file = 'Check_P_I_hat.csv')
        !    do k = 1, numK
        !        write(999,876) (P_I_hat(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        !open(unit = 999 , file = 'Check_c_hat.csv')
        !    do k = 1, numK
        !        write(999,876) (c_hat(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 2b: Compute P_F_hat
        !%%%%%%%%%%%%%%%%%%%%%%%%%
        
        P_F_hat = 1
        do i = 1, numC
            do k = 1, numK
                P_F_hat(i,1) = P_F_hat(i,1)*( P_I_hat(k,i)**mu(k,i) )
            end do
        end do
        
        !open(unit = 999 , file = 'Check_P_F_hat.csv')
        !    write(999,876) (P_F_hat(i,1), i = 1, numC)
        !close(999)
        
        !%%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 2c: Compute pi_hat
        !%%%%%%%%%%%%%%%%%%%%%%%%
        
        do o = 1, numC
            do i = 1, numC
                do k = 1, numK
                    pi_hat(o,i,k) = A_hat(k,o)*(c_hat(k,o)*d_hat(o,i,k)/P_I_hat(k,i))**(-lambda)
                end do
            end do
        end do
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 2d: Compute pi1 and kappa_tilde1
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        do o = 1, numC
            do i = 1, numC
                do k = 1, numK
                    pi1(o,i,k) = pi0(o,i,k)*pi_hat(o,i,k)
                end do
            end do
        end do
        ! Make sure shares sum to 1
        do i = 1, numC 
            do k = 1, numK
                pi1(:,i,k) = pi1(:,i,k) / sum( pi1(:,i,k) )
            end do
        end do            
        
        do i = 1, numC 
            do k = 1, numK
                kappa_tilde1(k,i) = kappa_tilde0(k,i)*( P_F_hat(i,1)/w_tilde_hat(k,i) )
            end do
        end do
        
        !open(unit = 999 , file = 'Check_pi1.csv')
        !    do k = 1, numK
        !        do o = 1, numC
        !            write(999,876) (pi1(o,i,k), i = 1, numC)
        !        end do
        !    end do
        !close(999)
        
        !open(unit = 999 , file = 'Check_kappa_tilde1.csv')
        !    do k = 1, numK
        !        write(999,876) (kappa_tilde1(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        !%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 3: TO BE COMPLETED
        !%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%%%%%%%%%%%%%%%
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 4: Compute qtheta1, theta1, u1
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        varpi = ( ( 1 - (1-chi)*delta )*kappa_tilde1 )/( delta*(1-beta) ) 
        
        !% Compute Integral(x_bar)
        do i = 1, numC 
            do k = 1, numK
                Int_x_bar1(k,i) = Integral(sigma_x,k,i,x_bar1(k,i)) 
            end do
        end do
    
        !% Compute G(x_bar)
        do i = 1, numC 
            do k = 1, numK
                ! G_x_bar = normcdf(log(x_bar)./sigma_x) 
                call VDCDFNORM(1,log(x_bar1(k,i))/sigma_x(k,i),G_x_bar1(k,i))
            end do
        end do
        
        !% Compute q(theta1)
        q_theta1 = min(varpi/Int_x_bar1,0.99999)

        !% Compute theta
        theta1 = ( q_theta1**(-xsi) - 1 )**(1/xsi) 
    
        !% Compute u
        u1 = chi / ( theta1*q_theta1*(1-G_x_bar1) + chi )
        
        !open(unit = 999 , file = 'Check_u1.csv')
        !    do k = 1, numK
        !        write(999,876) (u1(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        !open(unit = 999 , file = 'Check_theta1.csv')
        !    do k = 1, numK
        !        write(999,876) (theta1(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 5: Solve for big_Y1
        !%%%%%%%%%%%%%%%%%%%%%%%%%
    
        dist_Y = 1
        big_Y1 = big_Y1_Guess
        iter_Y = 1
        do while (dist_Y > toler_Y .and. iter_Y < maxiter_Y)
        
            big_Y1_new = 0.0
            do o = 1, numC
                do k = 1, numK
                    do i = 1, numC
                        do l = 1, numK
                            big_Y1_new(k,o) = big_Y1_new(k,o) + pi1(o,i,k)*(mu(k,i)*gamma(l,i) + (1-gamma(l,i))*nu(l,k,i))*big_Y1(l,i)
                        end do
                        big_Y1_new(k,o) = big_Y1_new(k,o) - pi1(o,i,k)*mu(k,i)*NX1(i,1)
                    end do
                end do
            end do
  
            dist_Y = sum( abs((big_Y1_new - big_Y1)/max(big_Y1,0.0000001)) )
        
            big_Y1 = big_Y1_new / sum( big_Y1_new )
        
            iter_Y = iter_Y + 1
    
        end do
        ! Make sure big_Y1 sums to 1
        big_Y1       = big_Y1 / sum(big_Y1)
        big_Y1_Guess = big_Y1
        
        !open(unit = 999 , file = 'Check_big_Y1.csv')
        !    do k = 1, numK
        !        write(999,876) (big_Y1(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 6: Compute big_L1_tilde
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        !% Compute big_L_tilde
        do i = 1, numC 
            do k = 1, numK
                call VDCDFNORM(1,sigma_x(k,i) - log(x_bar1(k,i))/sigma_x(k,i),term1)
                call VDCDFNORM(1,-log(x_bar1(k,i))/sigma_x(k,i),term2)
                Int_2(k,i) = exp( (sigma_x(k,i)**2)/2 )*term1 / term2
            end do
        end do
        big_L_tilde1 = big_L1*(1-u1)*Int_2
        
        !open(unit = 999 , file = 'Check_big_L_tilde1.csv')
        !    do k = 1, numK
        !        write(999,876) (big_L_tilde1(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        !%%%%%%%%%%%%%%%%%%%%%%%
        ! Step 7: Update w_tilde
        !%%%%%%%%%%%%%%%%%%%%%%%
        
        w_tilde1_new = gamma*big_Y1 / big_L_tilde1
        
        !open(unit = 999 , file = 'Check_w_tilde1_new.csv')
        !    do k = 1, numK
        !        write(999,876) (w_tilde1_new(k,i), i = 1, numC)
        !    end do
        !close(999)
        
        dist_w = maxval( abs( w_tilde1_new - w_tilde1 )/max(abs(w_tilde1),0.0000001) )
        
        dist_w_vec(iter_w,1) = dist_w 
        
        if (iter_w == 1) then
            dist_w_best = dist_w 
        end if
        
        !% Saves w_tilde leading to the best dist_w
        if (dist_w < dist_w_best) then
            dist_w_best  = dist_w 
            w_tilde_best = w_tilde1
        end if
        
        !% Allow 300 iterations before checking for bad behavior
        if (iter_w > 300) then
            !% If dist_out does not improve, or if it improves very little (by less than toler_out/10),
            !% reduce step size and go back to x_bar and L that led to 
            !% the smallest value of dist_out
            if ( flag_osc_w == 0 .or. iter_w >= flag_osc_w + 200 ) then
                mean_dist_w1 = sum(dist_w_vec(iter_w-49:iter_w,1))/50 
                mean_dist_w2 = sum(dist_w_vec(iter_w-99:iter_w-51,1))/50  
            end if
            if ( ( (mean_dist_w1 >= mean_dist_w2) .or. &
                 ( (mean_dist_w1 - mean_dist_w2 < 0) .and. abs(mean_dist_w1 - mean_dist_w2) < toler_w ) ) .and. &
                flag_osc_w == 0) then
                step_w     = 0.9*step_w 
                flag_osc_w = iter_w 
            end if
            !% But we only allow step sizes to be reduced every 200 iterations
            if ( ( (mean_dist_w1 >= mean_dist_w2) .or. &
                ( (mean_dist_w1 - mean_dist_w2 < 0) .and. abs(mean_dist_w1 - mean_dist_w2) < toler_w ) ) .and. &
                iter_w >= flag_osc_w+200 ) then
                step_w     = 0.9*step_w 
                flag_osc_w = iter_w 
            end if
        end if
        
        w_tilde1 = (1-step_w)*w_tilde1 + step_w*w_tilde1_new
        
        !% If step size was adjusted after bad behavior, 
        !% restart with w_tilde = w_tilde_best    
        if (flag_osc_w == iter_w) then
            w_tilde1 = w_tilde_best 
        end if
        
        iter_w = iter_w + 1
        
    end do
    
    w_tilde1_Guess = w_tilde1
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Step 8: Compute big_E_V1 and big_E_C1
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    big_E_V1     = kappa_tilde1*w_tilde1*theta1*u1*big_L1
    gamma_big_Y1 = gamma*big_Y1 
    big_E_C1 = transpose( matrix_mult( ones_K_tr, gamma_big_Y1, 1, numK, numC ) ) - &
               transpose( matrix_mult( ones_K_tr, big_E_V1, 1, numK, numC ) ) - &
               NX1
    
    !open(unit = 999 , file = 'Check_big_E_V1.csv')
    !    do k = 1, numK
    !        write(999,876) (big_E_V1(k,i), i = 1, numC)
    !    end do
    !close(999)
    
    !open(unit = 999 , file = 'Check_big_E_C1.csv')
    !    write(999,876) (big_E_C1(i,1), i = 1, numC)
    !close(999)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Step 9: Compute lambda_tilde1
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lambda_tilde_aux = big_L_bar / big_E_C1
    lambda_tilde1    = repmat_row(lambda_tilde_aux(:,1),numC,numK) 
    
    !open(unit = 999 , file = 'Check_lambda_tilde1.csv')
    !    write(999,876) (lambda_tilde1(1,i), i = 1, numC)
    !close(999)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 10: Solve Bellman Equation
    !% and obtain big_U1
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    big_U1  = big_U1_Guess 
    dist_U = 1 
    iter_U = 1 
    do while ( dist_U > toler_U .and. iter_U < maxiter_U )

        do i = 1, numC
            max_big_U1(1,i) = maxval(big_U1(:,i))
        end do
        
        do k = 1, numK
            C_k = C(k,:,:)
            MatAux1 = exp( ( -C_k + b + theta1*kappa_tilde1*lambda_tilde1*w_tilde1*(beta/(1-beta)) + delta*big_U1 &
                    - delta*repmat_row(max_big_U1(1,:),numC,numK))/zeta )
            big_U1_new(k:k,:) = zeta*log( matrix_mult(ones_K_tr,MatAux1,1,numK,numC) ) + delta*max_big_U1
        end do
        
        if ( sum(abs(big_U1)) > 0 ) then
            dist_U = maxval( abs( (big_U1_new - big_U1)/abs(big_U1) ) ) 
        else
            dist_U = maxval( abs( big_U1_new - big_U1 ) ) 
        end if

        iter_U = iter_U + 1 

        big_U1 = big_U1_new 

    end do 
    ! Update big_U1_Guess to be used in the next iteration
    big_U1_Guess = big_U1
    
    !open(unit = 999 , file = 'Check_big_U1.csv')
    !    do k = 1, numK
    !        write(999,876) (big_U1(k,i), i = 1, numC)
    !    end do
    !close(999)
    
    !%%%%%%%%%%%%%%%%%%%%%%%
    ! Step 11: Update big_L1
    !%%%%%%%%%%%%%%%%%%%%%%%
   
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Step 11a: Compute transition rates s
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    do k = 1, numK
        C_k = C(k,:,:)
        MatAux1 = exp( ( -C_k + b + theta1*kappa_tilde1*lambda_tilde1*w_tilde1*(beta/(1-beta)) + delta*big_U1 &
                - delta*repmat_row(max_big_U1(1,:),numC,numK) )/zeta )
        MatAux2 = matrix_mult(ones_K_tr,MatAux1,1,numK,numC)
        s(k,:,:) = MatAux1 / repmat_row( MatAux2(1,:),numC,numK )
    end do
    
    !open(unit = 999 , file = 'Check_s.csv')
    !    do i = 1, numC
    !        do k = 1, numK
    !            write(999,876) (s(k,l,i), l = 1, numK)
    !        end do
    !    end do
    !close(999)
    
    !%%%%%%%%%%%%%%%%%%
    ! Steps 11b and 11c
    !%%%%%%%%%%%%%%%%%%
    
    eye_K = 0.0
    do i = 1, numK
        eye_K(i,i) = 1
    end do
    
    do i = 1, numC
        AA = eye_K - transpose(s(:,:,i))
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
        if (info == 0 .and. abs(eigenval(min_loc(1,1),1)) < 1e-12) then
            y_tilde(:,i) = eigenvec(:,min_loc(1,1)) / u1(:,i)
        else
            ! spower_old = identity
            spower_old = 0
            do k = 1, numK
                spower_old(k,k) = 1
            end do
            dist_spower = 1 
            iter_spower = 1 
            do while( dist_spower > 0.000000000001 .and. iter_spower <= 2000 )
                spower_new = matrix_mult(s(:,:,i),spower_old,numK,numK,numK)
                !dist_spower = norm2(spower_new - spower_old)/max(norm2(spower_old),0.000001)
                dist_spower = maxval( (spower_new - spower_old)/max(abs(spower_old),0.00000001) )
                spower_old = spower_new 
                iter_spower = iter_spower + 1
            end do
            y_tilde(:,i:i) = transpose(spower_new(1:1,:)) / u1(:,i:i)
        end if
        y_tilde_tr  = transpose(y_tilde(:,i:i))
        phi0        = big_L_bar(i,1) / ( matrix_mult(y_tilde_tr,ones_K,1,numK,1 ) ) 
        big_L1_new(:,i) = phi0(1,1)*y_tilde(:,i) 
    end do
    
    big_L1_new = max(big_L1_new,0.0)
    
    !open(unit = 999 , file = 'Check_big_L1_new.csv')
    !    do k = 1, numK
    !        write(999,876) (big_L1_new(k,i), i = 1, numC)
    !    end do
    !close(999)
    
    !%%%%%%%%%%%%%%%%%%%%%%%
    !% Step 12: Update x_bar
    !%%%%%%%%%%%%%%%%%%%%%%%

    x_bar1_new = ( (1-delta)*big_U1 - eta ) / (lambda_tilde1*w_tilde1)
    x_bar1_new = max( x_bar1_new , 0.000001 )
    
    !open(unit = 999 , file = 'x_bar1_new.csv')
    !    do k = 1, numK
    !        write(999,876) (x_bar1_new(k,i), i = 1, numC)
    !    end do
    !close(999)
    
    dist_L = maxval( abs(( big_L1_new - big_L1)/(max(abs(big_L1),0.0000001)) ) )
    dist_x = maxval( abs((x_bar1_new - x_bar1) /(max(abs(x_bar1),0.0000001)) ) )
    
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
        x_bar_best    = x_bar1 
        big_L_best    = big_L1
    end if
    
    !% Allow 300 iterations before checking for bad behavior
    if (iter_out > 300) then
        !% If dist_out does not improve, or if it improves very little (by less than toler_out/10),
        !% reduce step size and go back to x_bar and L that led to 
        !% the smallest value of dist_out
        if ( flag_osc == 0 .or. iter_out >= flag_osc + 200 ) then
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
             iter_out >= flag_osc + 200 ) then
            step_L   = 0.9*step_L 
            step_x   = 0.9*step_x 
            flag_osc = iter_out 
        end if
    end if
    
    big_L1 = (1-step_L)*big_L1 + step_L*big_L1_new 
    
    x_bar1 = (1-step_x)*x_bar1 + step_x*x_bar1_new 
    x_bar1 = max( x_bar1 , 0.000001 )
    
    !% If step size was adjusted after bad behavior, 
    !% restart with x_bar = x_bar_best 
    !% and big_L = big_L_best instead   
    if (flag_osc == iter_out) then
        big_L1 = big_L_best 
        x_bar1 = x_bar_best 
    end if
    
    iter_out = iter_out + 1
    
    write(*,*) 
    write(*,*) 'iter_out =', iter_out
    write(*,*) 'dist_out =', dist_out
    write(*,*)
    
end do

! Finish Time
call date_and_time(values=time_array_1)
finish_time = time_array_1(5) * 3600 + time_array_1(6) * 60 &
            + time_array_1(7) + 0.001 * time_array_1(8) 

print*
print*, 'Total Time     :', real(finish_time - start_time,4)
print*

pause

876 format(6(f50.16,','))

end subroutine NewSteadyState

!%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Defines I_ki(x) function
!%%%%%%%%%%%%%%%%%%%%%%%%%%

function Integral(sigma_x,k,i,x) result(Int)

    USE Global_Data
    
    implicit none
    
    real(KIND=DOUBLE), intent(in) :: sigma_x(:,:)
    integer, intent(in) :: k, i
    real(KIND=DOUBLE), intent(in) :: x
    
    real(KIND=DOUBLE) x_new, term1, term2, Int
    
    x_new   = max(x,0.000001) 
    call VDCDFNORM(1,sigma_x(k,i) - log(x_new)/sigma_x(k,i),term1)
    call VDCDFNORM(1,-log(x_new)/sigma_x(k,i),term2)
    Int     = exp( (sigma_x(k,i)**2) / 2 )*term1 &
            - x_new*term2

end function Integral

end Module computeNewSteadyState_MOD