% Problem Set 6, Question 2
%
% linearized solution of a Ramsey model with an unaticipated and 
% permanent shock to depreciation
% CRRA utility function
% u(c) = (c^(1 - \sigma) - 1)/(1 - \sigma),
% Cobb-Douglas production function
% ----------------------------------------------
clear all;   % some housekeeping

beta = 0.95; % discount factor
theta = 1.0/beta - 1.0; % discount rate
sigma = 2.0; % utility curvature
alpha = 0.3; % exponent on capital
Omega = 1; % productivity

delta_L = 0.1; % old depreciation rate
delta_H = 0.2; % new depreciation rate

T = 50; % number of time periods to simulate


% COMPUTING THE INITIAL STEADY STATE
% ===================================
Kss_R1 = ( (1-beta+delta_L*beta) / (alpha*beta*Omega) )^(1/(alpha-1));
Css_R1 = Omega*Kss_R1^alpha - delta_L*Kss_R1;

rSS_R1 = 1/beta - 1;
RSS_R1 = rSS_R1 + delta_L;
wSS_R1 = Kss_R1^alpha - alpha*Kss_R1^(alpha-1)*Kss_R1;
% COMPUTING THE NEW STEADY STATE AND ASSOCIATED DYNAMICS
% ======================================================

% computing the new steady state
% ==============================
Kss_R2 = ( (1-beta+delta_H*beta) / (alpha*beta*Omega) )^(1/(alpha-1));
Css_R2 = Omega*Kss_R2^alpha - delta_H*Kss_R2;

rSS_R2 = 1/beta - 1;
RSS_R2 = rSS_R2 + delta_H;
wSS_R2 = Kss_R2^alpha - alpha*Kss_R2^(alpha-1)*Kss_R2;

% computing the deviation on which the new policy rule will be applied
% ====================================================================
K_dev_R2 = Kss_R1 - Kss_R2;


% getting entries of the A matrix for regime 2,
% using steady state values and corresponding derivatives
% =======================================================
coeff1_R2 = -(1.0/sigma)*Css_R2*alpha*(alpha-1.0)*Omega*Kss_R2^(alpha - 2.0);

A_R2 = NaN*zeros(2,2);

A_R2(1,1) = 1.0 + beta*coeff1_R2;

A_R2(1,2) = -coeff1_R2;

A_R2(2,1) = -1.0;

A_R2(2,2) = 1.0/beta;

[V_R2,D_R2] = eig(A_R2); % returns eigenvectors collected in the columns of V,
                % and the corresponding eigenvalues on the diagonal of D

% finding out at which position the stable eigenvalue appears
for i_D = 1:2
    if abs(D_R2(i_D,i_D)) < 1.0
       stab_col_ind_R2 = i_D;
    end % of if      
end % of for loop running over the diagonal of the matrix D


% Simulating dynamics from the old to the new steady state
x_t(2,1) = K_dev_R2;
x_t(1,1) = (V_R2(1,stab_col_ind_R2)/V_R2(2,stab_col_ind_R2))*x_t(2,1);

for t = 2:T
    x_t(2,t) = A_R2(2,:)*x_t(:,t-1);
    x_t(1,t) = (V_R2(1,stab_col_ind_R2)/V_R2(2,stab_col_ind_R2))*x_t(2,t);
end

% Computing levels
c_sim = x_t(1,:) + Css_R2;
k_sim = x_t(2,:) + Kss_R2;

% Compute prices.

% Levels.
R_sim_level = alpha.*k_sim.^(alpha-1);                                      % Rental price of capital.
r_sim_level = R_sim_level - delta_H;                                        % Real interest rate.
w_sim_level = k_sim.^alpha - alpha.*k_sim.^(alpha-1).*k_sim;                % Wage rate.

% Deviations (from new steady state).
R_sim_dev = R_sim_level  - RSS_R2;                  
r_sim_dev = r_sim_level - rSS_R2;
w_sim_dev = w_sim_level - wSS_R2;

% Compute labor income share.
lab_inc_share_sim = w_sim_level./(k_sim.^alpha);

figure(1);
plot(1:T,x_t,'LineWidth',3);
xlabel('time','FontSize',20);
legend('consumption dev.','capital dev.');

figure(2);
plot(1:T,[c_sim;k_sim],'LineWidth',3);
xlabel('time','FontSize',20);
legend('consumption level','capital level');

figure(3);
subplot(1,2,1)
plot(k_sim/Kss_R1,'LineWidth',3)
xlabel('Time period t','FontSize',20);
ylabel('Capital k','FontSize',20);
title('Transition path capital k','FontSize',20);

subplot(1,2,2)
plot(c_sim/Css_R1,'LineWidth',3)
xlabel('Time period t','FontSize',20);
ylabel('Consumtion c','FontSize',20);
title('Transition path consumption c','FontSize',20);

figure(4)

subplot(1,2,1)
plot(1:T,[r_sim_level;R_sim_level],'LineWidth',3);
xlabel('time','FontSize',20);
ylabel('levels','FontSize',20);
legend('Real interest rate, level','Rental price of capital, level');

subplot(1,2,2)
plot(1:T,w_sim_level,'LineWidth',3);
xlabel('time','FontSize',20);
ylabel('levels','FontSize',20);
legend('Wage rate, level');

figure(5)

subplot(1,2,1)
plot(1:T,[r_sim_dev;R_sim_dev],'LineWidth',3);
xlabel('time','FontSize',20);
ylabel('Deviations from new s.s.','FontSize',20);
legend('Real interest rate','Rental price of capital');

subplot(1,2,2)
plot(1:T,w_sim_dev,'LineWidth',3);
xlabel('time','FontSize',20);
ylabel('Deviations from new s.s.','FontSize',20);
legend('Wage rate');

figure(6)

plot(1:T,lab_inc_share_sim,'LineWidth',3);
xlabel('time','FontSize',20);
ylabel('levels','FontSize',20);
ylim([0 1])
legend('Labor Income Share');

figure(7)
plot(1:T,k_sim.^alpha,'LineWidth',3);
xlabel('time','FontSize',20);
ylabel('levels','FontSize',20);
legend('Output');