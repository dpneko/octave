% This implements the solution of a
% multivariate rational expectations model,
% using the example of a basic RBC model,
% i.e., the model economy used in the lectures
% to explain the general solution procedure.
% The computations performed in the following code therefore
% correspond to pages 41 - 51 of the book
% by Roger E. A. Farmer,
% "Macroeconomics of Self-fulfilling Prophecies", 2nd edition, MIT Press, 1999.
% The comments refer to page numbers and equation numbers from the book.
% Thomas Hintermaier, hinterma@uni-bonn.de, December 16, 2019

% Specify paramter values used in the model
% (Remark: Notation is in line with symbols in the book.)
beta  = 0.95; % discount factor
alpha = 0.3;  % capital share
delta = 0.1;  % depreciation rate
rho   = 0.9;  % autocorrelation of productivity

% Getting non-stochastic steady-state values
% book: p. 46; lecture slides: Part 9, Slide 13
k_bar = ((1/beta - 1 + delta)/alpha)^(1/(alpha - 1));
c_bar = k_bar^alpha  - delta*k_bar;
s_bar = 1;

% Describing relationships underlying linear coefficients
% book: p. 46; lecture slides: Part 9, Slide 14.
a1 = beta*alpha*(alpha - 1)*k_bar^(alpha - 1);
a2 = beta*alpha*k_bar^(alpha - 1);
b1 = 1 - delta + alpha*k_bar^(alpha - 1);
b2 = k_bar^(alpha - 1);
b3 = -c_bar/k_bar;

% Writing dynamic equilibrium conditions into matrix form

% m371 stands for the 1st matrix showing up in equation (3.7),
% i.e. the matrix multiplying period-t variables on the LHS of (3.7).
% book: p. 47; lecture slides: Part 9, Slide 15.

m371 = zeros(3,3);
m371(1,1) = -1;
m371(2,1) = b3;
m371(2,2) = b1;
m371(2,3) = b2;
m371(3,3) = rho;

% m372 stands for the 2nd matrix showing up in equation (3.7),
% i.e. the matrix multiplying period-t variables on the LHS of (3.7).
% book: p. 47; lecture slides: Part 9, Slide 15.
m372 = zeros(3,3);
m372(1,1) = -1;
m372(1,2) = a1;
m372(1,3) = a2;
m372(2,2) = 1;
m372(3,3) = 1;

% m38 stands for the inverse of m371
% book: equation (3.8), p. 47; lecture slides: Part 9, Slide 16 (top).
m38 = inv(m371);

% Matrix A
% book: equation (3.9), p. 47; lecture slides: Part 9, Slide 16.
% Note that this encodes all the relevant dynamics for the solution.
A = m38*m372;

% Calculate Eigenvectors and Eigenvalues of A
[Q_temp,Lam_temp] = eig(A); % These are labeled "temp", 
                            % for being temporary in the sense that we still
                            % need to make sure that the ORDERING of eigenvalues
                            % puts the forward-stable one in the first position

% The next lines are there to identify the appropriate eigenvalue,
% which is forward stable, and therefore can be used to derive
% a restriction to express a free variable as an equilibrium function of predetermined variables
diagonal_Lam_temp = diag(Lam_temp);
select_i_stable = abs(diagonal_Lam_temp) < 1;
mark_unstable = not(select_i_stable);

% Reordering to make sure the stable eigenvalue and the corresponding eigenvector
% are in the first positions in the matrix of eigenvalues and the matrix of eigenvectors.
Lam = diag([diagonal_Lam_temp(select_i_stable);diagonal_Lam_temp(mark_unstable)]);
Q = zeros(3,3);
Q(:,1) =     Q_temp(:,select_i_stable);
Q(:,2:end) = Q_temp(:,mark_unstable);
% The eigenvalue less than one in absolute value is in the first position now
% book: p. 50, equation (3.12); lecture slides: Part 10, Slide 7.

% The inverse of Q, as used in the decomposition A = Q*Lam*inv(Q)
% book, p. 50 (top), lecture slides: Part 10, Slide 4.
Q_inv = inv(Q);

% The key restriction: getting a free variable (consumption)
% as an equilibrium function of predetermined variables (capital, productivity-state).
% book: p. 51, lecture slides: Part 10, Slide 9.
disp('coefficients of rational expectations equilibrium function');
c_coeff_k = -Q_inv(1,2)/Q_inv(1,1)
c_coeff_s = -Q_inv(1,3)/Q_inv(1,1)

% % disp('Q');
% % Q
% % 
% % disp('Lam');
% % Lam
% % 
% % disp('inverse of Q');
% % Q_inv

var_v = 0.007^2;
% T = 10;
T = 100;
% T = 1000;
% T = 100000;

randn('seed',123); % Seeding the random number generator, to be able to control and fix randomness for various runs (experiments)
% Getting draws of the shock realizations
v_t = zeros(T,1);
% ---------------------------------------------------
% Specification to obtain impulse-response functions:
% v_t(2) = sqrt(var_v)*1; % make sure to comment out random shocks for later periods
% ---------------------------------------------------
% Specification for stochastic realizations:
v_t(2:end)=sqrt(var_v)*randn(length(v_t) - 1,1);
% ---------------------------------------------------

% Generating the sequence of autocorrelated productivity states
s_t = NaN*zeros(T,1);
s_t(1) = 0;
for t = 2:T
 
    s_t(t) = rho*s_t(t-1) + v_t(t);
    
end;

% % plot(s_t);
% % s_t(1:10)

% Generating the time series of variables, which are equilibria for the given shock realizations
c_t = NaN*zeros(T,1);
k_t = NaN*zeros(T+1,1);
k_t(1) = 0;

for t = 1:T;
    % Imposing the key equilibrium restriction from rational expectations:
    % The free variable (consumption) is expressed as a function of
    % the predetermined variables (capital and productivity state)
    % book: p. 51, lecture slides: Part 10, Slide 9.
    % This is the recursive structure described in Part 10, Slide 12.
    c_t(t) = -Q_inv(1,2)/Q_inv(1,1)*k_t(t) -Q_inv(1,3)/Q_inv(1,1)*s_t(t);
    k_t(t+1) = b1*k_t(t) + b2*s_t(t) + b3*c_t(t);   
    
end;

% % figure(1);
% % plot(s_t);
% % title('productivity state');
% % 
% % figure(2);
% % plot(c_t);
% % title('consumption');
% % 
% % figure(3);
% % plot(k_t);
% % title('capital stock');

figure(10);
plot(s_t,'r-x');
hold on;
plot(c_t,'b-x');
plot(k_t,'g-x');
set(gca,"fontsize", 18)
title('productivity state, consumption, capital','FontSize',24);
legend('productivity state','consumption','capital');
hold off;

% % disp('s_t        c_t         k_t+1');
% % [s_t,c_t,k_t(2:end)]













