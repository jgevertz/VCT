%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% This code determines the parameter constraints (bounds) for the       %
% simple model. It does so manually, by considering any biological      %
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
p_fixed = [tumor_mean(1) 10]; % T0, V0

%% Model parameters to fit
%     r,        beta,    delta_V,
lb = [0.28, 0.001,   0.1]; % r constrained by fit to individual control data
ub = [0.37, 0.1,      5]; % delta_V constrained by experimental measurements
Nparams = length(lb);

%% ADD TO PLOT: THE NORMAL DISTRIBUTION PLOTS MADE IN evaluate_distributions
mean_params = zeros(1,Nparams);
std_params = zeros(1,Nparams);
Ngen = NVpops*100;
plausible_params_normal = zeros(Ngen,Nparams);
for i = 1:Nparams
    mean_params(i) = (lb(i)+ub(i))/2; % mean set to center of region
    %param_std = std_scale*mean_params(i); % STD is a multiple of mean
    std_params(i) = (ub(i)-mean_params(i))/3; % set bounds to be 3 std from mean
    plausible_params_normal(:,i) = ...
        std_params(i).*randn(Ngen,1) + mean_params(i); 
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

figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.7, 0.8]);
sgtitle('Simple Model','FontSize',16,'FontWeight','bold')
for i = 1:Nparams
    % Solve model at lower bounds
    [sol1_lower, sol2_lower, sol3_lower] = solve_model(p_fixed, params_lower(i,:) ,tf);
    sol_lower.x = [sol1_lower.x    sol2_lower.x    sol3_lower.x];
    sol_lower.y = [sol1_lower.y    sol2_lower.y    sol3_lower.y];

    % Solve model at upper bounds
    [sol1_upper, sol2_upper, sol3_upper] = solve_model(p_fixed, params_upper(i,:) ,tf);
    sol_upper.x = [sol1_upper.x    sol2_upper.x    sol3_upper.x];
    sol_upper.y = [sol1_upper.y    sol2_upper.y    sol3_upper.y];

    %% Plot time-course of all VPs in comparison to mean and num_STD*STD bounds
    subplot(2,3,i)
    h1 = fill(x_Nstd, inBetween_Nstd, grayColor,'DisplayName','VP Range'); 
    hold on;
    h2 = plot(days,tumor_mean,'b','LineWidth',3,'DisplayName','Mean of Data'); 
    h3 = plot(sol_lower.x, sol_lower.y(1,:),'r','LineWidth',2,...
        'DisplayName','Min Param Traj'); 
    h4 = plot(sol_upper.x, sol_upper.y(1,:),'k','LineWidth',2,...
        'DisplayName','Max Param Traj'); 
    hold off;
    ax = gca;
    ax.FontSize = 14;
    ylim([0,6366]);
    xlim([0,tf]);
    xlabel('Time (days)','FontSize',16)
    ylabel('Tumor Volume (mm^3)','FontSize',16);
    legend([h1 h2 h3 h4],'FontSize',16,'Location','NorthWest'); 
    if i == 1
        title('r','FontSize',14)
        fprintf('r: [min best max]: ');
        ylim([0,6200])
    elseif i == 2
        title('\beta','FontSize',14)
        fprintf('beta: [min best max]: ');
    elseif i == 3
        title('\delta_V','FontSize',14)
        fprintf('delta_V: [min best max]: ');
    end
    subtitle(['min = ' num2str(lb(i)) ', median = ' num2str(mean_params(i)) ...
        ', max = ' num2str(ub(i))]);
    disp([lb(i) mean_params(i) ub(i)]);

    subplot(2,3,i+Nparams)
    hold on;
    if i == 1
        histogram(plausible_params_normal(:,i),'Normalization','pdf'); 
        xlabel('r','FontSize',16)
    elseif i == 2
        histogram(plausible_params_normal(:,i),'Normalization','pdf'); 
        xlabel('\beta','FontSize',16)
    elseif i == 3
        histogram(plausible_params_normal(:,i),'Normalization','pdf'); 
        xlabel('\delta_V','FontSize',16)
    end
    xline([lb(i) ub(i)],'--r','LineWidth',2);
    xline(mean_params(i),'-r','LineWidth',2);
    hold off;
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,'param_bounds_simple.fig');
exportgraphics(gcf,'param_bounds_simple.png');

save param_bounds_simple.mat lb ub mean_params std_params p_fixed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions needed to execute main code                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OV model with dosing protocol hard-wired in
function [sol1, sol2, sol3] = solve_model(p_fixed,p_fit,tf)
    % Fixed
    T0 = p_fixed(1); 
    V0 = p_fixed(2);

    % First virus dose
    ICs = [T0    V0];
    sol1 = ode23s(@(t,x) OV_model(t,x,p_fixed,p_fit), [0 2], ICs);
    % Second virus dose
    ICs = sol1.y(:,end)' + [0   V0];
    sol2 = ode23s(@(t,x) OV_model(t,x,p_fixed,p_fit), [2 4], ICs);
    % Third virus dose
    ICs = sol2.y(:,end)' + [0    V0];
    sol3 = ode23s(@(t,x) OV_model(t,x,p_fixed,p_fit), [4 tf], ICs);
end

%% OV model
function dydt = OV_model(t,y,p_fixed,p_fit)
    % To fit
    r = p_fit(1);
    beta = p_fit(2); 
    d_V = p_fit(3); 
    
    % Variables: enforce non-negativity
    T  = max(y(1),0); 
    V  = max(y(2),0);

    dT  = r*T - beta*T*V;  % changed simple model to be well-defined
                           % (cannot have density-dependence)
    dV  = - d_V*V; 
    dydt = [dT; dV]; 
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
