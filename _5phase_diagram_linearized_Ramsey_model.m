% linearized dynamics of the Ramsey model
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
% % % % % % T = 500; % ! note the difference if instead using solution_linearized_Ramsey_model_.m

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

x_t = NaN*zeros(2,T); % initializing the array for the simulated variables
% starting values
% ---------------
% capital STATE is in the second row of x_t
x_t(2,1) = -1; 
% consumption BEHAVIOR/CHOICE is in the first row of x_t,
% putting it on the stable eigenvector ensures convergence
x_t(1,1) = (V(1,stab_col_ind)/V(2,stab_col_ind))*x_t(2,1);

% iterating on the dynamics from some starting value for T periods
% Note: If you don't choose the two coordinates of the starting vector
% to be (proportional to) the eigenvector related to an eigenvalue < 1,
% then the trajecotry of those applications of x_t+1 = A*x_t will NOT lead to the steady state.

for t = 2:T;
    
    x_t(:,t) = A*x_t(:,t-1);
    
end;

figure(1);
plot(1:T,x_t);
set(gca,"fontsize", 24)
xlabel('time','FontSize',24);
legend('consumption dev.','capital dev.');

figure(2);
plot(x_t(2,:),x_t(1,:),'-x');
set(gca,"fontsize", 24)
ylabel('Consumption deviation','FontSize',24);
xlabel('Capital deviation','FontSize',24);

% additional objects to be shown in the plot are created below
% ------------------------------------------------------------
% setting the range for plotting
K_range = -Kss*0.6:Kss/30:Kss*0.6;
C_range = -Css*0.4:Css/30:Css*0.4;

[Km, Cm] = meshgrid(K_range,C_range);

% calculating one-period changes of variables (in each direction)
% as implied by the (linearized) dynamics,
% if starting from some position in the (c_t,k_t)-deviations space today

K_diff = NaN*zeros(size(Km));
C_diff = NaN*zeros(size(Cm));

I_2 = eye(2);
A_I = A - I_2;

for i = 1:size(Km,1);
    for j = 1:size(Km,2);
        C_diff(i,j) = A_I(1,:)*[Cm(i,j);Km(i,j)];  
        K_diff(i,j) = A_I(2,:)*[Cm(i,j);Km(i,j)];              
    end;
end;

K_diff(K_diff.^2 > 0.2*mean(mean(K_diff.^2))) = NaN;
C_diff(C_diff.^2 > 0.2*mean(mean(C_diff.^2))) = NaN;

% Note: The selection above of values to be plotted later
% is just there to make the plots look nicer, it is thus an example of a detail
% which  - unlike other parts of this code - is not central to the concepts I have emphasized in class.

figure(27);
quiver(Km,Cm,K_diff,C_diff)
hold on;

% finding out at which position the unstable eigenvalue appears
for i_D = 1:2;
    if abs(D(i_D,i_D)) > 1;
       unstab_col_ind = i_D;
    end; % of if      
end; % of for loop running over the diagonal of the matrix D

plot(K_range,K_range*(V(1,stab_col_ind)/V(2,stab_col_ind)),'b','LineWidth',2);
plot(K_range,K_range*(V(1,unstab_col_ind)/V(2,unstab_col_ind)),'m','LineWidth',2);

plot(K_range,-K_range*(A_I(2,2)/A_I(2,1)),'g','LineWidth',2);
plot(K_range,-K_range*(A_I(1,2)/A_I(1,1)),'r','LineWidth',2);

% xlim([min(K_range),max(K_range)]);
% ylim([min(C_range),max(C_range)]);

xlim([min(min(Km(real(~(isnan(K_diff))).*real(~(isnan(C_diff))) > 0.5))),max(max(Km(real(~(isnan(K_diff))).*real(~(isnan(C_diff))) > 0.5)))]);
ylim([min(min(Cm(real(~(isnan(K_diff))).*real(~(isnan(C_diff))) > 0.5))),max(max(Cm(real(~(isnan(K_diff))).*real(~(isnan(C_diff))) > 0.5)))]);
set(gca,"fontsize", 24)
xlabel('capital deviation from steady state','FontSize',24);
ylabel('consumption deviation from steady state','FontSize',24);
title('Dynamics of LINEARIZED system of deviations from steady state','FontSize',24);
% axis equal;
hold off;









