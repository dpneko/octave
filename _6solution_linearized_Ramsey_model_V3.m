% linearized solution of a Ramsey model
% CRRA utility function
% u(c) = (c^(1 - \sigma) - 1)/(1 - \sigma),
% Cobb-Douglas production function
% Thomas Hintermaier, hinterma@uni-bonn.de
% November 3, 2019
% ----------------------------------------------
clear all;   % some housekeeping

beta = 0.95; % discount factor

theta = 1.0/beta - 1.0; % discount rate

sigma = 2.0; % utility curvature

alpha = 0.3; % exponent on capital

delta = 0.1; % depreciation rate

T = 50; % number of time periods to simulate

% computing the steady state
% ==========================
Kss = ( (1-beta+delta*beta) / (alpha*beta) )^(1/(alpha-1));

Css = Kss^alpha - delta*Kss;

% getting entries of the A matrix,
% using steady state values and corresponding derivatives
% =======================================================
coeff1 = -(1.0/sigma)*Css*alpha*(alpha-1.0)*Kss^(alpha - 2.0);

A = NaN*zeros(2,2);

A(1,1) = 1.0 + beta*coeff1;

A(1,2) = -coeff1;

A(2,1) = -1.0;

A(2,2) = 1.0/beta;

[V,D] = eig(A); % returns eigenvectors collected in the columns of V,
                % and the corresponding eigenvalues on the diagonal of D

% finding out at which position the stable eigenvalue appears
for i_D = 1:2;
    if abs(D(i_D,i_D)) < 1.0;
       stab_col_ind = i_D;
    end; % of if      
end; % of for loop running over the diagonal of the matrix D

% Defining the grid range and number of grid points
K_range = -Kss*0.6:Kss/30:Kss*0.6;

% The linearly approximated consumption policy function
C_dev_lin = (V(1,stab_col_ind)/V(2,stab_col_ind))*K_range;

figure(11);
plot(K_range,C_dev_lin);
set(gca,"fontsize", 24)
xlabel('Capital today, deviation from st. st.','FontSize',24);
ylabel('Consumption today, deviation from st. st.','FontSize',24);
title('Linearly approximated consumption policy function','FontSize',24);

x_t = NaN*zeros(2,T); % initializing the array for the simulated variables
% Simulating dynamics for some variables
x_t(2,1) = -1;
x_t(1,1) = (V(1,stab_col_ind)/V(2,stab_col_ind))*x_t(2,1);

for t = 2:T;
    
    x_t(2,t) = A(2,:)*x_t(:,t-1);
    x_t(1,t) = (V(1,stab_col_ind)/V(2,stab_col_ind))*x_t(2,t);
    
end;

figure(21);
plot(1:T,x_t);
set(gca,"fontsize", 24)
xlabel('time','FontSize',24);
legend('consumption dev.','capital dev.');

figure(22);
plot(x_t(2,:),x_t(1,:),'-x');
set(gca,"fontsize", 24)
ylabel('Consumption deviation','FontSize',24);
xlabel('Capital deviation','FontSize',24);












