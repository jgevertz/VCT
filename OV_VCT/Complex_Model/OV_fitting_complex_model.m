%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                       %
% This codes finds the best fit parameters in the complex OV model:     %
%   U' = r*U - beta*U*V/N g_U*U*T/N, N = U + I,                         %
%   I' = beta*U*V/N - d_I*I - epsilon*g_U*I*T/N,                        %
%   V' = alpha*d_I*I - d_V*V, T' = lambda*I - d_T*T                     %
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
% - Model is solved using ode15s with all variables forced to be >=0.   %
%   ode15s was required as some parameterizations resulted in a system  %
%   that was too stiff for ode23s to handle. Further, the relative and  %
%   absolute tolerance had to be adjusted to handle other integration   %
%   issues:                                                             %
%    RelTol = 5e-3 (instead of the default of 1e-3)                     %
%    AbsTol = 5e-5 (instead of the default of 1e-5)                     %
%                                                                       %   
% Output using numSobols = 78125: Best fit has cost of = 0.270134 with  %
%   r = 0.280000                                                        %
%   beta = 0.002273                                                     %
%   delta_V = 4.999816                                                  %
%   alpha = 1848.526415                                                 %
%   gamma_I = 0.000454                                                  %
%   epsilon = 0.502501 so gamma_U = 0.000228                            %
%   lambda = 0.289367                                                   %
% Elapsed time is 13223.720585 seconds.                                 %
%                                                                       %
% Updated: 6/5/2025                                                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load in the time and data
clear all; close all; clc; tic;
rng(1); % fixed seed for debugging
options = optimoptions(@fmincon,'Display','notify'); %fmincon options
load data_AD.mat
tumor_STD = std(tumorsizes,0,2);
tf = days(end); % stopping time

%% Model ICs, fixed parameters, bounds for fit parameters
load param_bounds_complex.mat 
% Contains lb, ub, p_fixed: % dI, U0, I0, V0, T0
% Parameters to fit in order: r, beta, delta_V, alpha, gamma_I, epsilon,  lambda

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

%% Multi-start minimization across npools parallel pools
param = zeros(numSobols,numParams);
fits = zeros(1,numSobols);
exitflag = zeros(1,numSobols);
npools = 4; % number of parallel pools
parpool('local', npools); % Open distributed processing 
poolobj = gcp('nocreate'); % If no pool, don’t create 
if isempty(poolobj)
    poolsize = 0;
else
    poolsize = poolobj.NumWorkers;
end
parfor i = 1:numSobols
%for i = 1:numSobols
    %call fmincon from each starting parameter guess
    fprintf('Up to multi-start #%d of %d with parameters:\n',i,numSobols)
    disp(p_fit(i,:));
    [sol1, sol2, sol3] = solve_model(p_fixed,p_fit(i,:),tf);
    [param(i,:), fits(i), exitflag(i)] = fmincon(fun,...
        p_fit(i,:),[],[],[],[],lb,ub,[],options);
end
delete(gcp('nocreate'));
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
        fprintf('\tr = %f\n',best_fit_params(j));
    elseif j == 2
        fprintf('\tbeta = %f\n',best_fit_params(j)/1000);
    elseif j == 3
        fprintf('\tdelta_V = %f\n',best_fit_params(j));
    elseif j == 4
        fprintf('\talpha = %f\n',1000*best_fit_params(j));
    elseif j == 5
        fprintf('\tgamma_I = %f\n',best_fit_params(j));
    elseif j == 6
        fprintf('\tepsilon = %f so gamma_U = %f\n',best_fit_params(j),...
            best_fit_params(j)*best_fit_params(j-1));
    elseif j == 7
        fprintf('\tlambda = %f\n',best_fit_params(j));
    end
end

% View solution at top N parameters
for i = 1:5
    [sol1, sol2, sol3] = solve_model(p_fixed, param(index(i),:) ,tf);
    sol.x = [sol1.x    sol2.x    sol3.x];
    sol.y = [sol1.y    sol2.y    sol3.y];

    figure;
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.95, 0.95]);
    sgtitle(['r = ' num2str(param(index(i),1)) ', \beta = ' ...
        num2str(param(index(i),2)/1000) ', \delta_V = ' ...
        num2str(param(index(i),3)) ', \alpha = ' ...
        num2str(param(index(i),4)*1000) ', \gamma_I = ' ...
        num2str(param(index(i),5)) ', \epsilon = ' ...
        num2str(param(index(i),6))  ', (so \gamma_U = \epsilon\times \gamma_I = '...
        num2str(param(index(i),5)*param(index(i),6)) '), \lambda = ' ...
        num2str(param(index(i),7))],'FontSize',16,'FontWeight','bold'); 
    subplot(3,2,1)
    plot(sol.x, sol.y(1,:)+sol.y(2,:),'b','LineWidth',2)
    hold on;
    errorbar(days, tumor_mean, tumor_STD,'or')
    xlim([0,max(sol.x)+1])
    xlabel('Time (days)','FontSize',14)
    ylabel('Tumor Volume (mm^3)','FontSize',14)
    legend('Model','Data','Location','NorthWest','FontSize',14)
    hold off;
    
    subplot(3,2,2)
    plot(sol.x, sol.y(1,:),'b','LineWidth',2)
    xlim([0,max(sol.x)+1])
    xlabel('Time (days)','FontSize',14)
    ylabel('Uninfected Volume (mm^3)','FontSize',14)
    
    subplot(3,2,3)
    plot(sol.x, sol.y(2,:),'b','LineWidth',2)
    xlim([0,max(sol.x)+1])
    xlabel('Time (days)','FontSize',14)
    ylabel('Infected Volume (mm^3)','FontSize',14)
    
    subplot(3,2,4)
    plot(sol.x, sol.y(3,:),'b','LineWidth',2)
    xlim([0,max(sol.x)+1])
    xlabel('Time (days)','FontSize',14)
    ylabel('Virus Volume (mm^3)','FontSize',14)

    subplot(3,2,5)
    plot(sol.x, sol.y(4,:),'b','LineWidth',2)
    xlim([0,max(sol.x)+1])
    xlabel('Time (days)','FontSize',14)
    ylabel('T Cell Volume (mm^3)','FontSize',14)
    
    exportgraphics(gcf,['best_fit' num2str(i) '.png']);
    saveas(gcf,['best_fit' num2str(i) '.fig']);
end

%% Finding all fits within 5% of the best fit
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

%% Near optimal parameter histograms
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.95, 0.95]);
sgtitle(['Distributions within ' num2str(100*perc_opt) '% of Optimal'],...
    'FontSize',16,'FontWeight','bold')      
for j =1:numParams
    subplot(2,4,j)
    if j == 2
        histogram(best_fit_params_sorted(:,j)/1000);
        title(['Range: [' num2str( min(best_fit_params_sorted(:,j))/1000 ) ...
            ',' num2str( max(best_fit_params_sorted(:,j))/1000 ) ']']); 
    elseif j == 4
        histogram(best_fit_params_sorted(:,j)*1000);
        title(['Range: [' num2str( min(best_fit_params_sorted(:,j))*1000 )...
            ',' num2str( max(best_fit_params_sorted(:,j))*1000 ) ']']); 
    else
        histogram(best_fit_params_sorted(:,j));
        title(['Range: [' num2str( min(best_fit_params_sorted(:,j)) ) ','...
        num2str( max(best_fit_params_sorted(:,j)) ) ']']); 
    end

    if j == 1
        xlabel('r','FontSize',14)
    elseif j == 2
        xlabel('\beta','FontSize',14)
    elseif j == 3
        xlabel('\delta_V','FontSize',14)
    elseif j == 4
        xlabel('\alpha','FontSize',14)
    elseif j == 5
        xlabel('\gamma_I','FontSize',14)
    elseif j == 6
        xlabel('\epsilon (\gamma_U=\epsilon\times \gamma_I)','FontSize',14)
    elseif j == 7
        xlabel('\lambda','FontSize',14)
    end
end
exportgraphics(gcf,'best_param_distributions.png');
saveas(gcf,'best_param_distributions.fig');
toc;

save best_fit_complex.mat best_fit_params best_fit_objective p_fixed ...
    fits param fit_sorted index params_fits best_fit_objective_sorted ...
    best_fit_params_sorted lb ub


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
    modeldata_U = [deval(sol1, 0:1, 1)   deval(sol2, 2:3, 1)  ...
        deval(sol3, 4:tf, 1)];
    modeldata_I = [deval(sol1, 0:1, 2)   deval(sol2, 2:3, 2)  ...
        deval(sol3, 4:tf, 2)];    
    modeldata = (modeldata_U + modeldata_I)'; % Total tumor volume
    
    % Goodness of fit: divide by 2 to match likelihood formula 
    %SSE = sum((tumor_mean(:)-modeldata(:)).^2);
    cost = 0.5*sum(((tumor_mean(:)-modeldata(:)).^2)./(tumor_STD.^2));
end



