%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                       %
% This codes generates a virtual population with NVPops virtual         %
% patients for the SIMPLE OV model:                                     %
%   T' = r*T - beta*T*V,   dV' = -d_V*V                                 %
% using either the accept-or-reject (enter 1 at prompt) or the          %
% accept-or-perturb (enter 2 at prompt) method.                         %
% - Requires the actual OV experimental data in data_Ad.mat             %
% - Also requires the fixed parameters and initial conditions specified %
%   in p_fixed and the lower bounds (lb) and upper bounds (ub) for each %
%   parameter to be fit. These are read in from the param_bounds.mat    %
%   file generated from running the OV_param_bnds_v2.m code.            %
%                                                                       %
% If you run accept-or-reject, you'll further be prompted to say how    %
% you want to impose the constraints on your parameters specified in lb %
% and ub.                                                               %
% - The first option is to just sample from a normal prior, and after   %
%   sampling, reject any parameterization where the parameter           %
%   constraints do not hold.                                            %
% - The second option is to truncate the normal prior so no parameters  %
%   are sampled outside the specified lower and upper bounds.           %
%                                                                       %
% If you run accept-or-perturb, the lower and upper bounds for each     %
% uniform parameter distribution are set to match the bounds set in lb  %
% and ub, respectively.                                                 %
%                                                                       %
% A solution trajectory meets the constrains of a VP provided it stays  %
% within traj_std*sigma from the mean of the experimental data at each  %
% time point, where both the mean and standard deviation (sigma) are    %
% computed from the experimental data.                                  %
%  - traj_std is set at 3, but can be made smaller for a narrow         %
%    distribution or larger for a wider distribution.                   %
%                                                                       %
% For Accept-or-Reject: Each parameter is assumed to be normally        %
% distributed about its best-fit value.                                 %
%  - NVPops = 500 required generating 945 plausible patients if we      %
%    truncate the parameters after generating the prior.                %
%    Run time = 12.415544 seconds.                                      %
%  - NVPops = 500 required generating 948 plausible patients if we      %
%    truncate the parameters after generating the prior.                %
%    Run time = 12.866473 seconds.                                      %
%  - Note: now that the parameter bounds well match the normal          %
%    distribution extent, it really doesn't matter if we truncate the   %
%    prior or simply reject paramters that don't satisfy the bounds.    %
%                                                                       %
% For Accept-or-Perturb with NVPops = 500 and bounds as specified in    %
% corresponding fit file. Run time = 1463.043060 seconds.               %
%                                                                       %
% Last update: 6/4/2025                                                 %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close all; tic;
rng(1); % fixed seed for replicability

%% Read and data and visualize
load data_AD.mat
tf = days(end); % stopping time
load param_bounds_simple.mat 
% Contains lb, ub, p_fixed = [T0, V0]
% VP params in order = [r, beta, delta_V, alpha]

%% Define region solution trajectory must stay within for a virtual patient
tumor_STD = std(tumorsizes,0,2);
traj_std = 3; % Can change if want a tighter or wider distribution of trajectories
curve1_Nstd = tumor_mean' + traj_std*tumor_STD'; % must be size 1xnumDays
curve2_Nstd = tumor_mean' - traj_std*tumor_STD';
x_Nstd = [days', fliplr(days')]; % time pt vector must be of size 1xnumDays
inBetween_Nstd = [curve1_Nstd, fliplr(curve2_Nstd)];
grayColor = [.7 .7 .7];

%% Query for VP approach
fprintf('Do you want to generate VPs using accept-or-reject (enter 1) '); 
prompt = "or accept-or-perturb (enter 2)? ";
VP_method = input(prompt); % Stores value answered to prompt
NVpops = 50; %0;
if VP_method == 1 % accept-or-reject 
    fprintf('Do you truncate the parameters after generating the prior (enter 1)\n'); 
    prompt2 = "or before generating the prior (enter 2)? ";
    distribution_method = input(prompt2); 

    if distribution_method == 1 % truncate after generation
        % Generate plausible patients from a normal distribution
        plausible_params = normal_distributions(NVpops,mean_params,...
            std_params); 
        
        % Now sample the distributions to find NVPops virtual patients that 
        % satisfy both parameter and trajectory constraints
        [Vpop,Traj] = accept_or_reject(NVpops,curve1_Nstd,...
            curve2_Nstd,plausible_params,p_fixed,tf,lb,ub);
        fname_fig_params = 'VPop_distributions_acc_rej';
        fname_fig_trajects = 'VPop_trajectories_acc_rej';
        title_text = 'Accept-or-Reject';
        fsave = 'Vpops_acc_rej.mat';

    elseif distribution_method == 2 % truncate before generation
        % Generate plausible patients from a truncated normal distribution
        % that satisfies parameter constraints
        plausible_params = normal_distributions_trunc(NVpops,...
            mean_params,lb,ub,std_params); 

        % Now sample the distributions to find NVPops virtual patients that 
        % satisfy only the trajectory constraints
        [Vpop,Traj] = accept_or_reject_trunc(NVpops,curve1_Nstd,...
            curve2_Nstd,plausible_params,p_fixed,tf);
        fname_fig_params = 'VPop_distributions_acc_rej_trunc';
        fname_fig_trajects = 'VPop_trajectories_acc_rej_trunc';
        title_text = 'Accept-or-Reject (Truncate Priors)';
        fsave = 'Vpops_acc_rej_trunc.mat';
    end
    
elseif VP_method == 2 % accept-or-perturb
    % Generate plausible patients from a uniform distribution
    plausible_params = uniform_distributions(NVpops,mean_params,lb,ub);

    % Now sample the distributions to find NVPops virtual patients
    [Vpop,Traj] = accept_or_perturb(NVpops,curve1_Nstd,...
        curve2_Nstd,plausible_params,p_fixed,tf,lb,ub);

    fname_fig_params = 'VPop_distributions_acc_perturb';
    fname_fig_trajects = 'VPop_trajectories_acc_perturb';
    title_text = 'Accept-or-Perturb';
    fsave = 'Vpops_acc_pert.mat';

else
    fprinf('Must enter 1 or 2 to generate your virtual population\n');
end

%% Plot initial distributions and distributions of VP parameters
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.6]);
sgtitle([title_text ' Parameter Distributions'],'fontSize',16,'fontweight','bold');
for i = 1:length(lb)
    subplot(1,3,i)
    histogram(plausible_params(:,i),'Normalization','pdf'); hold on; 
    histogram(Vpop(:,i),'Normalization','pdf'); hold off;
    if i == 1
        xlabel('r','FontSize',16)
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
        legend('Sampling','VPop','FontSize',16,'location','northwest');
    end
    ax = gca;
    ax.FontSize = 14;
    title(['VPop range: [' num2str(min(Vpop(:,i))) ',' num2str(max(Vpop(:,i))) ...
        ']'],'FontSize',14);
end
% saveas(gcf,[fname_fig_params,'.fig']);
% exportgraphics(gcf,[fname_fig_params,'.png']);

%% Plot time-course of all VPs in comparison to mean and num_STD*STD bounds
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.3, 0.6]);
h1 = fill(x_Nstd, inBetween_Nstd, grayColor,'DisplayName','VP Range'); 
hold on;
h2 = plot(days,Traj','o-'); 
h3 = plot(days,tumor_mean,'k','LineWidth',3,'DisplayName','Mean of Data'); 
hold off;
ax = gca;
ax.FontSize = 14;
ylim([0,inf]);
xlabel('Time (days)','FontSize',16)
ylabel('Tumor Volume (mm^3)','FontSize',16)
title([title_text ' VPop Trajectories within 3\sigma of Mean'],'FontSize',16)
legend([h1 h3],'FontSize',16,'Location','NorthWest'); 
% saveas(gcf,[fname_fig_trajects,'.fig']);
% exportgraphics(gcf,[fname_fig_trajects,'.png']);

toc;
% save(fsave,'curve1_Nstd','curve2_Nstd','plausible_params','Vpop','Traj',...
%     'traj_std');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Functions%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plausible_params = normal_distributions(NVpops,mean_params,std_params)
    Nparams = length(mean_params);
    Ngen = NVpops*100; % generate many more that needed due to rejections
    plausible_params = zeros(Ngen,Nparams);
    for i = 1:Nparams
        param_mean = mean_params(i); 
        param_std = std_params(i); 
        plausible_params(:,i) = param_std.*randn(Ngen,1) + param_mean; 
    end
end

% Truncated normal
function plausible_params = normal_distributions_trunc(NVpops,mean_params,...
                                                lb,ub,std_params)
    Nparams = length(mean_params);
    Ngen = NVpops*10;  % generate more that needed due to rejections
    plausible_params = zeros(Ngen,Nparams);
    for i = 1:Nparams
        param_mean = mean_params(i); 
        param_std = std_params(i); 
        pd = makedist('Normal','mu',param_mean,'sigma',param_std);
        t = truncate(pd,lb(i),ub(i))
        plausible_params(:,i) = random(t,Ngen,1);
    end
end

function plausible_params = uniform_distributions(NVpops,mean_params,lb,ub)
    Nparams = length(mean_params);
    Ngen = NVpops*2; % need much less VPs in accept-or-perturb
    plausible_params = zeros(Ngen,Nparams);
    for i = 1:Nparams
        param_min = lb(i);
        param_max = ub(i);
        plausible_params(:,i) = param_min + (param_max-param_min)*rand(Ngen,1);
    end
end

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

function [Vpop,Traj] = accept_or_reject(NVpops,curve1,curve2,VP_params,...
                                                          p_fixed,tf,lb,ub)   
    Nparams = size(VP_params,2);
    Vpop = zeros(NVpops,Nparams); % 3 parameters to represent a VP
    Traj = zeros(NVpops,tf+1);                                        
    i = 1;
    count_VPs = 0;
    VP_attempts = 0; 
    while count_VPs < NVpops
        plausible_params = VP_params(i,:);
        fprintf('Plausible patient number %d has r = %f, beta = %f, delta_V = %f\n',...
            i,VP_params(i,1), VP_params(i,2), VP_params(i,3)); 
        VP_attempts = VP_attempts + 1;
        
        lb_holds = plausible_params - lb; % Need >= so params >=lb
        ub_holds = ub - plausible_params; % Need >= so ub >= params

        if all(plausible_params>=0) && all(lb_holds>=0) && all(ub_holds>=0) 
            [sol1, sol2, sol3] = solve_model(p_fixed,plausible_params,tf);
            modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) ...
                deval(sol3, 4:tf, 1)];
            modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) ...
                deval(sol3, 4:tf, 2)];    
            traj_sample = (modeldata_U + modeldata_I)'; % Total tumor volume

            accept = 0; 
            for j = 1:length(traj_sample)
                if (traj_sample(j)<curve2(j))||(traj_sample(j)>curve1(j)) 
                    accept = -1;
                end
            end
            
            if accept == 0 % add VP
                count_VPs = count_VPs+1; 
                fprintf('\t%d. Accepting virtual patient so now have %d VPs\n',...
                    i,count_VPs);
                Vpop(count_VPs,:) = plausible_params;
                Traj(count_VPs,:) = traj_sample;
            else
                fprintf('\t%d. Outside trajectory bounds: rejecting virtual patient, ',i);
                fprintf('so still have %d VPs\n',count_VPs)
            end
        else
            fprintf('\t%d. Outside parameter bounds: rejecting virtual patient, ',i);
            fprintf('so still have %d VPs\n',count_VPs)
        end
        i = i+1;
    end                                                
end

function [Vpop,Traj] = accept_or_reject_trunc(NVpops,curve1,curve2,...
                                                     VP_params,p_fixed,tf)   
    Nparams = size(VP_params,2);
    Vpop = zeros(NVpops,Nparams); % 3 parameters to represent a VP
    Traj = zeros(NVpops,tf+1);                                        
    i = 1;
    count_VPs = 0;
    VP_attempts = 0; 
    while count_VPs < NVpops
        plausible_params = VP_params(i,:);
        fprintf('Plausible patient number %d has r = %f, beta = %f, delta_V = %f\n',...
            i,VP_params(i,1), VP_params(i,2), VP_params(i,3)); 
        VP_attempts = VP_attempts + 1;
        
        [sol1, sol2, sol3] = solve_model(p_fixed,plausible_params,tf);
        modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) ...
            deval(sol3, 4:tf, 1)];
        modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) ...
            deval(sol3, 4:tf, 2)];    
        traj_sample = (modeldata_U + modeldata_I)'; % Total tumor volume

        accept = 0; 
        for j = 1:length(traj_sample)
            if (traj_sample(j)<curve2(j))||(traj_sample(j)>curve1(j)) 
                accept = -1;
            end
        end
        
        if accept == 0 % add VP
            count_VPs = count_VPs+1; 
            fprintf('\t%d. Accepting virtual patient so now have %d VPs\n',...
                i,count_VPs);
            Vpop(count_VPs,:) = plausible_params;
            Traj(count_VPs,:) = traj_sample;
        else
            fprintf('\t%d. Outside trajectory bounds: rejecting virtual patient, ',i);
            fprintf('so still have %d VPs\n',count_VPs)
        end

        i = i+1;
    end                                                
end

function [Vpop,Traj] = accept_or_perturb(NVpops,curve1,curve2,VP_params,...
                                                          p_fixed,tf,lb,ub)
    Nparams = size(VP_params,2);
    Vpop = zeros(NVpops,Nparams); % 3 parameters to represent a VP
    Traj = zeros(NVpops,tf+1); 
    i = 1;
    count_VPs = 0;
    VP_attempts = 0; 
    num_perturbed = 0; 
    while count_VPs < NVpops
        plausible_params = VP_params(i,:);
        fprintf('Plausible patient number %d has r = %f, beta = %f, delta_V = %f\n',...
            i,VP_params(i,1), VP_params(i,2), VP_params(i,3)); 
        VP_attempts = VP_attempts + 1;

        [sol1, sol2, sol3] = solve_model(p_fixed,plausible_params,tf);
        modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) ...
            deval(sol3, 4:tf, 1)];
        modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) ...
            deval(sol3, 4:tf, 2)];    
        traj_sample = (modeldata_U + modeldata_I)'; % Total tumor volume
      
        accept = 0; 
        for j = 1:length(traj_sample)
            if (traj_sample(j)<curve2(j))||(traj_sample(j)>curve1(j)) 
                accept = -1;
            end
        end

        if accept == 0 % add VP
            count_VPs = count_VPs+1; 
            fprintf('\t%d. Accepting virtual patient so now have %d VPs\n',...
                    i,count_VPs);
            Vpop(count_VPs,:) = plausible_params;
            Traj(count_VPs,:) = traj_sample;
        else
            %% Instead of rejecting, any ones out of bounds get moved into
            %% bounds by finding r,d minimizing g(r,d) = as defined in chapter
            fprintf('\t%d. Perturbing parameters for virtual patient\n',i);
            num_perturbed = num_perturbed + 1;

            fun = @(z)objective(z,curve2,curve1,p_fixed,tf);
            VP_sim_ann = simulannealbnd(fun,plausible_params,lb,ub);

            [sol1, sol2, sol3] = solve_model(p_fixed,VP_sim_ann,tf);
            modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) ...
                deval(sol3, 4:tf, 1)];
            modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) ...
                deval(sol3, 4:tf, 2)];  
            traj_sim_ann  = (modeldata_U + modeldata_I)'; % Total tumor volume

            accept = 0; 
            for j = 1:length(traj_sim_ann)
                if (traj_sim_ann(j)<curve2(j))||(traj_sim_ann(j)>curve1(j)) 
                    accept = -1;
                end
            end
            if accept == 0 % add VP
                count_VPs = count_VPs+1; 
                fprintf('\t\t%d. Accepting perturbed virtual patient so now have %d VPs\n',...
                    i,count_VPs);
                Vpop(count_VPs,:) = VP_sim_ann;
                Traj(count_VPs,:) = traj_sim_ann;
            else
                fprintf('\t\t%d. Rejecting perturbed virtual patient, ',i);
                fprintf('so still have %d VPs\n',count_VPs)
            end
        end
        i = i+1;
    end
end

function value = objective(VP,lb,ub,p_fixed,tf)
    [sol1, sol2, sol3] = solve_model(p_fixed,VP,tf);
    modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) ...
        deval(sol3, 4:tf, 1)];
    modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) ...
        deval(sol3, 4:tf, 2)];    
    VP_sample  = (modeldata_U + modeldata_I)';

    sum = 0;
    for i = 1:length(VP_sample)
        term1 = (lb(i)+ub(i))/2; 
        term2 = (VP_sample(i)-term1)^2; 
        term3 = ( (ub(i)/2) - (lb(i)/2) )^2; 
        sum_term = (term2 - term3);
        sum = sum + max(sum_term,0); 
    end
    value = max(sum, 0);
end