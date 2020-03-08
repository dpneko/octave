
% Golden Rule dynamics for the production function: F(K) = K^alpha


% parameters
alpha = 0.33;
delta = 0.5;


T = 20; % number of periods simulated

% golden rule levels of capital stock and consumption

ksharp = (delta/alpha)^(1/(alpha-1));
csharp = ksharp^alpha - delta*ksharp;

kgrid = linspace(0,2*ksharp,1000); % (grid) range considered

% law of motion for capital, under specific assumption
kplus = kgrid.^alpha + (1 - delta)*kgrid - csharp;

figure(1);
plot(kgrid,kplus,'LineWidth',2)
hold on;
plot(kgrid,kgrid,'r','LineWidth',2)
plot(ksharp,ksharp,'s','MarkerSize',10)
xlim([min(kgrid) max(kgrid)]);
ylim([min(kgrid) max(kgrid)]);
axis("square") % This is just to make the 45 degree line really 45 degree, avoiding distortions from the proportions of the plot.



initial_multiple = 2; % Situation (a)

% initial_multiple = 0.5; % Situation (b)

k_t = NaN(1,T);
k_t(1) = ksharp*initial_multiple; % initial value of capital
for i = 1:T
    k_t(i+1) = k_t(i).^alpha + (1 - delta)*k_t(i) - csharp;
    if k_t(i+1) < 0
        k_t(i+1) = NaN;
        break
    end
end


plot(k_t(1:end-1),k_t(2:end),'d','MarkerSize',8,'MarkerFaceColor', 'b');
xlabel('Capital today','FontSize',25);
ylabel('Capital tomorrow','FontSize',25);
title('Dynamics under the Golden Rule','FontSize',25);
hold off;

figure(2);
plot(1:T,k_t(1:T),'d','MarkerSize',8,'MarkerFaceColor', 'b');
xlabel('Time','FontSize',25);
ylabel('Capital','FontSize',25);

