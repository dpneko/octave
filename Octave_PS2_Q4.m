% Problem Set 2, Question 4

clear all;
clc;

% Setting up the model, parameters
% ================================
alpha  = 0.5;       % capital share
delta  = 0.1;       % depreciation rate
sigma  = 2;         % utility curvature
beta   = 0.95;      % discount factor

time_vec = 1:30;

K_vec = zeros(size(time_vec));
C_vec = zeros(size(time_vec));

% b) Initial values for part b)
%  K_vec(1) = 1.0;
%  C_vec(1) = 1.0;

% d) Initial values for part d)
K_vec(1) = 1.0;
C_vec(1) = 0.627582223701021;

for t = 1:max(size(time_vec)) 
    % Dynamic relations for the two variables as derived in part a) of the
    % question.
    K_vec(t+1) = K_vec(t)^alpha + (1 - delta)*K_vec(t) - C_vec(t);
%%%%%     C_vec(t+1) = C_vec(t)*(beta*( alpha*(K_vec(t)^alpha + (1 - delta)*K_vec(t) - C_vec(t))^(alpha-1) + 1 - delta ) )^(1/sigma);
    C_vec(t+1) = C_vec(t)*(beta*( alpha*K_vec(t+1)^(alpha-1) + 1 - delta ) )^(1/sigma);
    if K_vec(t+1) < 0
        K_vec(t+1) = NaN;
        C_vec(t+1) = NaN;
        fprintf('Infeasible choice - Negative capital in period %d!\n', t+1)
    end
end

figure;
plot(K_vec,C_vec,'rd');
xlabel('Capital','FontSize',25)
ylabel('Consumption','FontSize',25)
title('Sequence for capital and consumption','FontSize',25)
