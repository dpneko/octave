% Value function iteration
% to perform optimization for a Ramsey model,
% CRRA utility function, Cobb-Douglas production
% Thomas Hintermaier, hinterma@uni-bonn.de
% January 15, 2018
% ----------------------------------------------

% Setting paramters
beta = 0.95;
sigma = 2.0;
alpha = 0.3;
delta = 0.1;

% We get steady states just for making 
% reasonable choices for the grid range 
% of the state variable (capital)
Kss = ( (1-beta+delta*beta) / (alpha*beta) )^(1/(alpha-1));
Css = Kss^alpha - delta*Kss;

% Defining the grid range and number of grid points
K_grid = linspace(0.8*Kss,1.2*Kss,400);

% We maintain the linearly approximated solution,
% obtained from the analysis of the phase diagram in the program linearized_quiver_CRRA.m
% just in case we later want to compare it to the non-linear policy functions obtained from dynamic programming.
% For the linearly approximated solution, 
% we need to get entries of the A matrix, using steady state values and corresponding derivatives

coeff1 = -(1.0/sigma)*Css*alpha*(alpha-1.0)*Kss^(alpha - 2.0);
A = NaN*zeros(2,2);
A(1,1) = 1.0 + beta*coeff1;
A(1,2) = -coeff1;
A(2,1) = -1.0;
A(2,2) = 1.0/beta;

[Q,Lambda] = eig(A);

[trash_one,   stab_col_ind] = max(abs(diag(Lambda) < 1.0));
[trash_one, unstab_col_ind] = max(abs(diag(Lambda) > 1.0));

% The linearly approximated consumption policy function
C_lin = (Q(1,stab_col_ind)/Q(2,stab_col_ind))*(K_grid - Kss) + Css;

figure(11);
plot(K_grid,C_lin);
set(gca,"fontsize",20);
xlabel('Capital today','FontSize',24);
ylabel('Consumption today','FontSize',24);
title('Linearly approximated consumption policy function','FontSize',24);

% Performing iteration on the value function
% ------------------------------------------

% Get current-period utility payoff
% implied by moving from some K_now (today) to some K_next (next period)
% Note that it is convenient (to save on computations) to compute this upfront, 
% before we enter the iteration on the value function,
% so we can rely on this (constant-over-iterations) object later in the iterations.

U_now_K_next_K_now = -Inf*zeros(length(K_grid),length(K_grid));

for j_now = 1:length(K_grid);
    for i_next = 1:length(K_grid);
    
    K_now = K_grid(j_now);
    K_next = K_grid(i_next);
    C_now = K_now^alpha + (1 - delta)*K_now - K_next;
    if C_now > 0;
    U_now_K_next_K_now(i_next,j_now) = (C_now^(1 - sigma) - 1)/(1 - sigma);
    end; % of if
    % Note that the period-utility calculation
    % relies on a CRRA utility function, and uses a Cobb-Douglas production function
    
    end; % of for over i_next
 end; % of for over j_now    

% Initialize the value function
V_old = zeros(1,length(K_grid));
V = NaN*zeros(1,length(K_grid));
index_K_policy = NaN*zeros(1,length(K_grid));
K_policy = NaN*zeros(1,length(K_grid));

crit = 0.000001; % used in the criterion for checking convergence
max_iter = 500;  % maximum number of iterations, relevant only if convergence is not achieved before (which should be!)

% Just for illustration, not needed for algorithm,
% storing the value functions obtained in all iterations
V_store_iter = NaN*zeros(max_iter,length(K_grid)); 
 
% This is the key iteration (Bellman equation)
% --------------------------------------------
for i_iter = 1:max_iter;

   % Bellman equation
   for j_now = 1:length(K_grid);
   [V(j_now), index_K_policy(j_now)] = max(U_now_K_next_K_now(:,j_now) + beta*V_old');
   end;

   max_diff = max(abs(V - V_old)); 
   if max_diff < crit;
      break; % quit iterations, when convergence is achieved
   end;

   % displaying on screen the current stage of iterations, to have some feedback during computations
   if mod(i_iter,10)== 0; % just at specific iterations, to avoid lengthy output on screen
   disp(['Iteration number ',num2str(i_iter)]);
   end;
   
   V_store_iter(i_iter,:) = V;

   V_old = V;
end; % of for-loop over iterations on value function

if max_diff >= crit;
   disp('Did not achieve convergence!');
end;

% Getting capital levels corresponding to indexes of the maximizers
K_policy = K_grid(index_K_policy);

% Getting consumption policy from capital policy
C_policy = K_grid.^alpha + (1 - delta)*K_grid - K_policy; 

figure(20);
plot(K_grid,V);
set(gca,"fontsize",20);
xlabel('Capital today','FontSize',24);
ylabel('Value function','FontSize',24);
title('Value function','FontSize',24);

figure(22);
plot(K_grid,K_policy);
set(gca,"fontsize",20);
hold on;
plot(K_grid,K_grid,'k');
hold off;
xlabel('Capital today','FontSize',24);
ylabel('Capital tomorrow','FontSize',24);
title('Capital policy function and 45 degree line','FontSize',24);

figure(21);
plot(K_grid,C_policy);
set(gca,"fontsize",20);

% the following three lines are to also include the linearized policy
% comment out, if just showing the results of dynamic programming
hold on;
plot(K_grid,C_lin);
hold off;

xlabel('Capital today','FontSize',24);
ylabel('Consumption today','FontSize',24);
title('Consumption policy function','FontSize',24);

figure(30);
plot(K_grid,V_store_iter(1:i_iter,:));
set(gca,"fontsize",20);
xlabel('Capital today','FontSize',24);
ylabel('Value functions','FontSize',24);
title('Value functions during iterations towards convergence','FontSize',24);













