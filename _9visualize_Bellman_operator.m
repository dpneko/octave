% visualize_Bellman_operator.m
% Thomas Hintermaier, hinterma@uni-bonn.de
% January 15, 2018

% This script visualizes what happens in the operation
% specified by the Bellman equation.
% This follows the explanation of objects involved,
% as developed on the blackboard in class.

% Run this script AFTER you have run the script val_iter.m,
% making sure that all variables from val_iter.m are still in memory (Workspace).

% Pick some grid points of today's capital state,
% for which the components in the Bellman operator are visualized.

pick_some_index1 = round(0.3*numel(K_grid));
% The 0.3 on the previous line means that a grid point at about 30% into the grid is picked.
% Change 0.3 to any other relative position (between 0 and 1) to look at other states of capital.

pick_some_index2 = round(0.4*numel(K_grid));
% Another level of the state of capital today,
% so you get to see how the components vary in the level of the capital stock.

% Plot slices of the current-period utility consequence,
% for starting from the state of capital specified above,
% and considering the variation from choosing various levels of next-period capital.

figure(401);
plot(K_grid,U_now_K_next_K_now(:,pick_some_index1), 'r', 'LineWidth',3,...
     K_grid,U_now_K_next_K_now(:,pick_some_index2), 'g', 'LineWidth',3);
set(gca,"fontsize", 24);
xlabel('Choice of next-period capital, k_{t+1}','FontSize',24);
ylabel('U(c_t)','FontSize',24);
title('Current-period utility from consumption, if choosing some level of next-period capital','FontSize',24);
legend(['capital state, k_t = ',num2str(K_grid(pick_some_index1),4)],
       ['capital state, k_t = ',num2str(K_grid(pick_some_index2),4)],"location","southwest");

% Plot the discounted continuation value,
% considering the variation from choosing various levels of next-period capital.

figure(402);
plot(K_grid,beta*V_store_iter(i_iter-1,:)','y','LineWidth',3);
set(gca,"fontsize", 24);
xlabel('Choice of next-period capital, k_{t+1}','FontSize',24);
ylabel('beta*V_{t+1}(k_{t+1})','FontSize',24);
title('discounted continuation value, if choosing some level of next-period capital','FontSize',24);

% Plot the total of the objective on the right-hand-side of Bellman equation,
% for starting from the state of capital specified above,
% and considering the variation from choosing various levels of next-period capital.
figure(403);
plot(K_grid,U_now_K_next_K_now(:,pick_some_index1) + beta*V_store_iter(i_iter-1,:)', 'r', 'LineWidth',3,...
     K_grid,U_now_K_next_K_now(:,pick_some_index2) + beta*V_store_iter(i_iter-1,:)', 'g', 'LineWidth',3);
set(gca,"fontsize", 24);
xlabel('Choice of next-period capital, k_{t+1}','FontSize',24);
ylabel('U(c_t) + beta*V_{t+1}(k_{t+1})','FontSize',24);
title('total of objective on right-hand-side of Bellman equation, if choosing some level of next-period capital','FontSize',24);
legend(['capital state, k_t = ',num2str(K_grid(pick_some_index1),4)],
       ['capital state, k_t = ',num2str(K_grid(pick_some_index2),4)],"location","southwest");

 







