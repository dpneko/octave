% Problem Set 1, Question 5
% Golden Rule dynamics for CES function: F(K,L) = ( alpha*K^s + (1-alpha)*L^s )^(1/s)
% L = Lbar = constant

% Parameters
alpha = 1/3;
delta = 0.5;
Lbar  = 1;
s     = 0.1;

T = 20; % number of periods simulated

% Golden rule steady state values
ksharp = ((delta^(s/(1-s))-alpha^(1/(1-s)))/((1-alpha)*alpha^(s/(1-s))*Lbar^s))^(-1/s);
csharp = (alpha*ksharp^s+(1-alpha)*Lbar^s)^(1/s) - delta*ksharp;

kgrid = linspace(0,2*ksharp,1000);

% Law of motion for capital
kplus = (alpha*kgrid.^s+(1-alpha)*Lbar^s).^(1/s) + (1 - delta)*kgrid - csharp;

figure(1);
plot(kgrid,kplus,'LineWidth',2)
hold on;
plot(kgrid,kgrid,'r','LineWidth',2)
plot(ksharp,ksharp,'s','MarkerSize',10)
xlim([min(kgrid) max(kgrid)]);
ylim([min(kgrid) max(kgrid)]);
xlabel('Capital today','FontSize',25);
ylabel('Capital tomorrow','FontSize',25);
title('Dynamics under the Golden Rule','FontSize',25);


initial_multiple = 2; % Situation (a)

% initial_multiple = 1/2; % Situation (b)

k_t = NaN(1,T);

k_t(1) = ksharp*initial_multiple; % initial value of capital
for i = 1:T
    k_t(i+1) = (alpha*k_t(i)^s+(1-alpha)*Lbar^s)^(1/s)+ (1 - delta)*k_t(i) - csharp;
    if k_t(i+1) < 0
        k_t(i+1) = NaN;
        break
    end
end


plot(k_t(1:end-1),k_t(2:end),'d','MarkerSize',8,'MarkerFaceColor', 'b');
xlabel('Capital today','FontSize',25);
ylabel('Capital tomorrow','FontSize',25);
hold off;

figure(2);
plot(1:T,k_t(1:T),'d','MarkerSize',8,'MarkerFaceColor', 'b');
xlabel('Time','FontSize',25);
ylabel('Capital','FontSize',25);

