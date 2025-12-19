%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                       %
% This codes finds the best fit parameters in the SIMPLE OV model:      %
%   T' = r*T - beta*T*V,   dV' = -d_V*V                                 %
% to the average of the experimental data.                              %
% - Requires the actual OV experimental data in data_Ad.mat             %
% - Also requires the fixed parameters and initial conditions specified %
%   in p_fixed and the lower bounds (lb) and upper bounds (ub) for each %
%   parameter to be fit. These are read in from the param_bounds.mat    %
%   file generated from running the OV_param_bnds_v2.m code.            %
% - Implements multi-start fitting from numSobols starting points.      %
%   numSobols is computed based on dimensionality of parameter space.   %
% - fmincon gets called from each of the numSobols starting parameter   %
%   guesses                                                             %
% - Best of the best is selected as the optimal model parameterization, %
%   and the near-optimal parameters have a goodness-of-fit within 5% of %
%   optimal, though that threshold can be changed.                      %
%                                                                       %   
% Output using numSobols = 125: Best fit has cost of = 0.386148 with    %
% 	r = 0.280000                                                        %
% 	beta = 0.007629                                                     %
% 	delta_V = 0.594669                                                  %
% Elapsed time is 63.750703 seconds.                                    %
%                                                                       %
% Updated: 6/4/2025                                                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load in the time and data
clear all; close all; tic; clc; 
rng(1); % fixed seed for debugging
options = optimoptions(@fmincon,'Display','notify'); %fmincon options
load data_AD.mat
tumor_STD = std(tumorsizes,0,2);
tf = days(end); % stopping time

%% Model ICs, fixed parameters, bounds for fit parameters
load param_bounds_simple.mat 
% Contains lb, ub, p_fixed: % T0, V0
% Parameters to fit in order: r, beta, delta_V,

%% Multi-start fitting algorithm setup
numParams = length(lb);
numSobols = 5^numParams; % how many random points to sample
n_skip = 1000; n_leap = 0; % Parameters needed by sobolset 
uniform_sobol = sobolset(numParams,'Skip',n_skip,'Leap',n_leap);
uniform_sobol = net(uniform_sobol,numSobols);
% Rescale parameters to be in range [lb,ub]
p_fit = zeros(numSobols,numParams);
for i = 1:numParams
    p_fit(:,i) = (ub(i) - lb(i))*uniform_sobol(:,i)+lb(i);
end
fun = @(z)objective(z,p_fixed,tumor_mean,tumor_STD,tf); % goodness-of-fit function

%% Multi-start minimization across Mpool parallel pools
param = zeros(numSobols,numParams);
fits = zeros(1,numSobols);
exitflag = zeros(1,numSobols);
for i = 1:numSobols
    %call fmincon from each starting parameter guess
    fprintf('Up to multi-start #%d of %d\n',i,numSobols)
    [param(i,:), fits(i), exitflag(i)] = fmincon(fun,...
        p_fit(i,:),[],[],[],[],lb,ub,[],options);
end
params_fits = [param fits'];

%% Sort best fit per multistart from best to worst 
[fit_sorted, index] = sort(fits); 
best_fit_objective = fit_sorted(1);
best_fit_params  = param(index(1),:);
best_fit_objective_sorted = [];
best_fit_params_sorted = [];
best_fit_objective_sorted(1) = fit_sorted(1);
best_fit_params_sorted(1,:)  = param(index(1),:);
fprintf('Best fit has cost of = %f with parameters\n',best_fit_objective);
for j = 1:numParams
    if j == 1
        fprintf('\tr = ')
	    fprintf('%f\n',best_fit_params(j));
    elseif j == 2
        fprintf('\tbeta = ')
	    fprintf('%f\n',best_fit_params(j));
    elseif j == 3
        fprintf('\tdelta_V = ')
	    fprintf('%f\n',best_fit_params(j));
    end
end

%% View solution at top 5 parameters
for i = 1:5
    [sol1, sol2, sol3] = solve_model(p_fixed, param(index(i),:) ,tf);
    sol.x = [sol1.x    sol2.x    sol3.x];
    sol.y = [sol1.y    sol2.y    sol3.y];

    figure;
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.75, 0.85]);
    sgtitle(['r = ' num2str(param(index(i),1)) ', \beta = ' ...
        num2str(param(index(i),2)),', \delta_V = ' ...
        num2str(param(index(i),3))],'FontSize',16,'FontWeight','bold'); 
    subplot(1,2,1)
    plot(sol.x, sol.y(1,:),'b','LineWidth',2)
    hold on;
    errorbar(days, tumor_mean, tumor_STD,'or')
    xlim([0,max(sol.x)+1])
    xlabel('Time (days)','FontSize',14)
    ylabel('Tumor Volume (mm^3)','FontSize',14)
    legend('Model','Data','Location','NorthWest','FontSize',14)
    hold off;
        
    subplot(1,2,2)
    plot(sol.x, sol.y(2,:),'b','LineWidth',2)
    xlim([0,max(sol.x)+1])
    xlabel('Time (days)','FontSize',14)
    ylabel('Virus Volume (mm^3)','FontSize',14)
        
    exportgraphics(gcf,['best_fit' num2str(i) '.png']);
    saveas(gcf,['best_fit' num2str(i) '.fig']);
end

%% Finding all fits within (perc_opt*100)% of the best fit
perc_opt = 0.05;
perc_diff = zeros(numSobols,1);
best_counter=2;
for i=2:length(fit_sorted)
    %percent difference from best objective function
    perc_diff(i)=(fit_sorted(i)-fit_sorted(1))/fit_sorted(1);

    if perc_diff(i) < perc_opt % within (100*perc_opt)% from optimal
        best_fit_objective_sorted(best_counter) = fit_sorted(i);
        best_fit_params_sorted(best_counter,:)  = param(index(i),:);
        best_counter = best_counter+1;
    end
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.75, 0.85]);
sgtitle(['Distributions within ' num2str(100*perc_opt) '% of Optimal'],...
    'FontSize',16,'FontWeight','bold')          
for j =1:numParams
    subplot(2,2,j)
    histogram(best_fit_params_sorted(:,j));
    if j == 1
        xlabel('r','FontSize',14)
    elseif j == 2
        xlabel('\beta','FontSize',14)
    elseif j == 3
        xlabel('\delta_V','FontSize',14)
    end
    title(['Range: [' num2str( min(best_fit_params_sorted(:,j)) ) ','...
        num2str( max(best_fit_params_sorted(:,j)) ) ']']); 
end
sgtitle('Parameter Distributions','FontSize',16,'FontWeight','bold')      
exportgraphics(gcf,'best_param_distributions.png');
saveas(gcf,'best_param_distributions.fig');
toc;

save best_fit_simple.mat best_fit_params best_fit_objective p_fixed ...
    fits param fit_sorted index params_fits best_fit_objective_sorted ...
    best_fit_params_sorted lb ub

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



