%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% This code determines the parameter constraints (bounds) for fitting   %
% the complex model. It does so manually, by considering any biological %
% constraints and trying to get the trajectories to cover a range at    %
% least as wide as the region allowable for a VP.                       %
%                                                                       %
% Updated: 6/4/2025                                                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read and data and visualize
clear all; clc; close all;
format shortG
load data_AD.mat
tf = days(end); % stopping time
rng(1); 
NVpops = 500;

%% Model parameters/ICs to fix
p_fixed = [1 0.35 tumor_mean(1) 0 10 0]; % dI, dT, U0, I0, V0, T0

%% Model parameters to fit
% r constrained by fit to individual control data
% delta_V constrained by experimental measurements
% alpha constrained by experimental measurements
%     r,    beta*1000,  delta_V,  alpha/1000, gamma_I,  epsilon,  lambda
lb = [0.28, 0.4,        0.1,      1,          10^(-5),  0,         0]; 
ub = [0.37, 4,          5,        5,          10^(-2),  1,         1];        
Nparams = length(lb);

%% ADD TO PLOT: THE NORMAL DISTRIBUTION PLOTS MADE IN evaluate_distributions
mean_params = zeros(1,Nparams);
std_params = zeros(1,Nparams);
Ngen = NVpops*100;
plausible_params_normal = zeros(Ngen,Nparams);
for i = 1:Nparams
    mean_params(i) = (lb(i)+ub(i))/2; % mean set to center of region
    std_params(i) = (ub(i)-mean_params(i))/3; % set bounds to be 3 std from mean
    plausible_params_normal(:,i) = std_params(i).*randn(Ngen,1) + mean_params(i); 
end

%% Define region solution trajectory must stay within for a virtual patient
tumor_STD = std(tumorsizes,0,2);
traj_std = 3; % Can change if want a tighter or wider distribution of trajectories
curve1_Nstd = tumor_mean' + traj_std*tumor_STD'; % must be size 1xnumDays
curve2_Nstd = tumor_mean' - traj_std*tumor_STD';
x_Nstd = [days', fliplr(days')]; % time pt vector must be of size 1xnumDays
inBetween_Nstd = [curve1_Nstd, fliplr(curve2_Nstd)];
grayColor = [.7 .7 .7];

params_lower = zeros(Nparams,Nparams);
for i = 1:Nparams
    params_lower(i,:) = mean_params;
    params_lower(i,i) = lb(i);
end

params_upper = zeros(Nparams,Nparams);
for i = 1:Nparams
    params_upper(i,:) = mean_params;
    params_upper(i,i) = ub(i);
end

Nplot1 = 4;
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.7, 0.8]);
sgtitle('Complex Model (1 of 2)','FontSize',16,'FontWeight','bold')
for i = 1:Nplot1
    % Solve model at lower bounds
    [sol1_lower, sol2_lower, sol3_lower] = solve_model(p_fixed, params_lower(i,:) ,tf);
    sol_lower.x = [sol1_lower.x    sol2_lower.x    sol3_lower.x];
    sol_lower.y = [sol1_lower.y    sol2_lower.y    sol3_lower.y];

    % Solve model at upper bounds
    [sol1_upper, sol2_upper, sol3_upper] = solve_model(p_fixed, params_upper(i,:) ,tf);
    sol_upper.x = [sol1_upper.x    sol2_upper.x    sol3_upper.x];
    sol_upper.y = [sol1_upper.y    sol2_upper.y    sol3_upper.y];

    %% Plot time-course of all VPs in comparison to mean and num_STD*STD bounds
    subplot(2,4,i)
    h1 = fill(x_Nstd, inBetween_Nstd, grayColor,'DisplayName','VP Range'); 
    hold on;
    h2 = plot(days,tumor_mean,'b','LineWidth',3,'DisplayName','Mean of Data'); 
    h3 = plot(sol_lower.x, sol_lower.y(1,:)+sol_lower.y(2,:),'r',...
        'LineWidth',2,'DisplayName','Min Param Traj'); 
    h4 = plot(sol_upper.x, sol_upper.y(1,:)+sol_upper.y(2,:),'k',...
        'LineWidth',2,'DisplayName','Max Param Traj'); 
    hold off;
    ax = gca;
    ax.FontSize = 14;
    ylim([0,inf]);
    %ylim([0,6366]);
    xlim([0,tf]);
    xlabel('Time (days)','FontSize',16)
    ylabel('Tumor Volume (mm^3)','FontSize',16);
    legend([h1 h2 h3 h4],'FontSize',16,'Location','NorthWest'); 
    if i == 1
        title('r','FontSize',14)
        fprintf('r: [min best max]: ');
        %ylim([0,6200])
    elseif i == 2
        title('\beta','FontSize',14)
        fprintf('beta: [min best max]: ');
    elseif i == 3
        title('\delta_V','FontSize',14)
        fprintf('delta_V: [min best max]: ');
    elseif i == 4
        title('\alpha','FontSize',14)
        fprintf('alpha: [min best max]: ');
    end
    subtitle(['min = ' num2str(lb(i)) ', median = ' num2str(mean_params(i)) ...
        ', max = ' num2str(ub(i))]);
    disp([lb(i) mean_params(i) ub(i)]);

    subplot(2,4,i+Nplot1)
    hold on;
    if i == 1
        histogram(plausible_params_normal(:,i),'Normalization','pdf'); 
        xlabel('r','FontSize',16)
    elseif i == 2
        histogram(plausible_params_normal(:,i)/1000,'Normalization','pdf'); 
        xlabel('\beta','FontSize',16)
    elseif i == 3
        histogram(plausible_params_normal(:,i),'Normalization','pdf'); 
        xlabel('\delta_V','FontSize',16)
    elseif i == 4
        histogram(plausible_params_normal(:,i)*1000,'Normalization','pdf'); 
        xlabel('\alpha','FontSize',16)
    end
    if i == 2
        xline([lb(i) ub(i)]/1000,'--r','LineWidth',2);
        xline(mean_params(i)/1000,'-r','LineWidth',2);
    elseif i == 4
        xline([lb(i) ub(i)]*1000,'--r','LineWidth',2);
        xline(mean_params(i)*1000,'-r','LineWidth',2);
    else
        xline([lb(i) ub(i)],'--r','LineWidth',2);
        xline(mean_params(i),'-r','LineWidth',2);
    end
    hold off;
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,'param_bounds_complex1.fig');
exportgraphics(gcf,'param_bounds_complex1.png');

Nplot2 = Nparams - Nplot1;
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.7, 0.8]);
sgtitle('Complex Model (2 of 2)','FontSize',16,'FontWeight','bold')
for i = Nplot1+1:Nparams
    % Solve model at lower bounds
    [sol1_lower, sol2_lower, sol3_lower] = solve_model(p_fixed,...
        params_lower(i,:) ,tf);
    sol_lower.x = [sol1_lower.x    sol2_lower.x    sol3_lower.x];
    sol_lower.y = [sol1_lower.y    sol2_lower.y    sol3_lower.y];

    % Solve model at upper bounds
    [sol1_upper, sol2_upper, sol3_upper] = solve_model(p_fixed, ...
        params_upper(i,:) ,tf);
    sol_upper.x = [sol1_upper.x    sol2_upper.x    sol3_upper.x];
    sol_upper.y = [sol1_upper.y    sol2_upper.y    sol3_upper.y];

    %% Plot time-course of all VPs in comparison to mean and num_STD*STD bounds
    subplot(2,3,i-Nplot1)
    h1 = fill(x_Nstd, inBetween_Nstd, grayColor,'DisplayName','VP Range'); 
    hold on;
    h2 = plot(days,tumor_mean,'b','LineWidth',3,'DisplayName','Mean of Data'); 
    h3 = plot(sol_lower.x, sol_lower.y(1,:)+sol_lower.y(2,:),'r',...
        'LineWidth',2,'DisplayName','Min Param Traj'); 
    h4 = plot(sol_upper.x, sol_upper.y(1,:)+sol_upper.y(2,:),'k',...
        'LineWidth',2,'DisplayName','Max Param Traj'); 
    hold off;
    ax = gca;
    ax.FontSize = 14;
    ylim([0,inf]);
    xlim([0,tf]);
    xlabel('Time (days)','FontSize',16)
    ylabel('Tumor Volume (mm^3)','FontSize',16);
    legend([h1 h2 h3 h4],'FontSize',16,'Location','NorthWest'); 
    if i == 5
        title('\gamma_I','FontSize',14)
        fprintf('gamma_I: [min best max]: ');
    elseif i == 6
        title('\epsilon (\gamma_U = \epsilon\times \gamma_I)','FontSize',14)
        fprintf('epsilon: [min best max]: ');
    elseif i == 7
        title('\lambda','FontSize',14)
        fprintf('lambda: [min best max]: ');
    end
    subtitle(['min = ' num2str(lb(i)) ', best = ' num2str(mean_params(i)) ...
        ', max = ' num2str(ub(i))]);
    disp([lb(i) mean_params(i) ub(i)]);

    subplot(2,3,i-Nplot1+Nplot2)
    hold on;
    if i == 5
        histogram(plausible_params_normal(:,i),'Normalization','pdf'); 
        xlabel('\gamma_I','FontSize',16)
    elseif i == 6
        histogram(plausible_params_normal(:,i),'Normalization','pdf'); 
        xlabel('\epsilon (\gamma_U = \epsilon\times \gamma_I)','FontSize',16)
    elseif i == 7
        histogram(plausible_params_normal(:,i),'Normalization','pdf'); 
        xlabel('\lambda','FontSize',16)
    end
    xline([lb(i) ub(i)],'--r','LineWidth',2);
    xline(mean_params(i),'-r','LineWidth',2);

    hold off;
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,'param_bounds_complex2.fig');
exportgraphics(gcf,'param_bounds_complex2.png');

save param_bounds_complex.mat lb ub mean_params std_params p_fixed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions needed to execute main code                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OV model with dosing protocol hard-wired in
function [sol1, sol2, sol3] = solve_model(p_fixed,p_fit,tf)
    ode_opts = odeset('RelTol',5e-3,'AbsTol',5e-5); 
    % Fixed
    U0 = p_fixed(3); 
    I0 = p_fixed(4);
    V0 = p_fixed(5);
    T0 = p_fixed(6);

    % First virus dose
    ICs = [U0    I0    V0   T0];
    sol1 = ode15s(@(t,x) OV_model(t,x,p_fixed,p_fit), [0 2], ICs,ode_opts);

    % Second virus dose
    ICs = sol1.y(:,end)' + [0    0    V0    0];
    sol2 = ode15s(@(t,x) OV_model(t,x,p_fixed,p_fit), [2 4], ICs,ode_opts);

    % Third virus dose
    ICs = sol2.y(:,end)' + [0    0    V0    0];
    sol3 = ode15s(@(t,x) OV_model(t,x,p_fixed,p_fit), [4 tf], ICs,ode_opts);
end

%% OV model
function dydt = OV_model(t,y,p_fixed,p_fit)
    % Fixed
    d_I = p_fixed(1);
    d_T = p_fixed(2);

    % To fit
    r = p_fit(1);
    beta = p_fit(2);
    d_V = p_fit(3);
    alpha = p_fit(4);
    g_I = p_fit(5);
    epsilon = p_fit(6);
    g_U = epsilon*g_I;
    lambda = p_fit(7);
    
    % Variables: enforce non-negativity
    U  = max(y(1),0); 
    I  = max(y(2),0);
    V  = max(y(3),0);
    T  = max(y(4),0);

    % Frequency-dependent killing: avoid division by 0 on freq-dep terms
    N = U + I;
    dU  = r*U - max(0,beta*U*V/N) - max(0,g_U*U*T/N);
    dI  = max(0,beta*U*V/N) - d_I*I - max(0,g_I*I*T/N);
    dV  = alpha*d_I*I - d_V*V; 
    dT  = lambda*I - d_T*T;
    dydt = [dU; dI; dV; dT]; 
end


%% Define cost to minimize
function cost = objective(p_fit,p_fixed,tumor_mean,tumor_STD,tf) 
    [sol1,sol2,sol3] = solve_model(p_fixed,p_fit,tf);

    % Evaluate numerical solution to ODE at experimental time points
    modeldata = [deval(sol1, 0:1, 1)   deval(sol2, 2:3, 1)    deval(sol3, 4:tf, 1)];
    
    % Goodness of fit: divide by 2 to match likelihood formula 
    %SSE = sum((tumor_mean(:)-modeldata(:)).^2);
    cost = 0.5*sum(((tumor_mean(:)-modeldata(:)).^2)./(tumor_STD.^2));
end
