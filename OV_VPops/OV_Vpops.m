%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                       %
% This codes generates virtual populations using the accept-or-reject   %
% method for the OV model in which r, beta and U0 vary across Vpops.    %
% Each parameter is assumed to be normally distributed about its best-  %
% fit value, and the standard deviation of the distribution is set to a %
% quarter of the mean. A Vpop is accepted if its trajectory falls       %
% within 3 standard deviations of the mean of the experimental data.    %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clearvars; close all; clc; 
rng(1);

%% Read and display data, and fix treatment protocol 
load data_AD.mat
tf = days(end); % stopping time
I0 = 0; % no infected cells initially
V0  = 10; % initial number of viruses, scaled to volume

%% Best fits in scenario where r, beta, U0 fixed, all others fixed
% EMD_Randstad/Identifiability_for_Sampling/OV_Data/Fit_beta_U0_r/OV_fitting.m
p.r = 0.267632968652492; 
p.beta = 0.659286153982597;
p.d_I = 1; 
p.alpha = 3; 
p.d_V = 2.3; 
p.U0 =57.249205506759170; 
in_cond = [p.U0 I0 V0]; 
[sol1, sol2, sol3] = OV_model(p,in_cond,tf);
% Combine data
sol.x = [sol1.x    sol2.x    sol3.x];
sol.y = [sol1.y    sol2.y    sol3.y];
% Plot
figure;
errorbar(days, tumor_mean, tumor_SE,'o','LineWidth',2);
hold on;
plot(sol.x, sol.y(1,:)+sol.y(2,:),'LineWidth',2);
%plot(days, tumorsizes(:,1),'o');
hold off;
xlim([0,max(sol.x)+1])
xlabel('Time (days)','FontSize',16);
ylabel('Volume (mm^3)','FontSize',16);
legend('Data','Model best fit','FontSize',16','location','northwest')

%% Define region solution trajectory must stay within for a virtual patient
tumor_STD = std(tumorsizes,0,2);
Nstd_dev = 3;
curve1_Nstd = tumor_mean' + Nstd_dev*tumor_STD'; % must be size 1xnumDays
curve2_Nstd = tumor_mean' - Nstd_dev*tumor_STD';
x_Nstd = [days', fliplr(days')]; % time pt vector must be of size 1xnumDays
inBetween_Nstd = [curve1_Nstd, fliplr(curve2_Nstd)];
grayColor = [.7 .7 .7];

%% Generate Vpop using accept-or-reject method
N_Vpops = 100;
r_mean = p.r; 
r_std = 0.25*r_mean; 
r_norm = r_std.*randn(5*N_Vpops,1) + r_mean; % generate 5x more to ensure I have enough data after rejection
beta_mean = p.beta; 
beta_std = 0.25*beta_mean; 
beta_norm = beta_std.*randn(5*N_Vpops,1) + beta_mean;  % generate 5x more to ensure I have enough data after rejection
U0_mean = p.U0; 
U0_std = 0.25*U0_mean; 
U0_norm = U0_std.*randn(5*N_Vpops,1) + U0_mean; 
[Vpop_AccRej,Traj_AccRej] = accept_or_reject(N_Vpops,curve1_Nstd,...
    curve2_Nstd,p, r_norm,beta_norm,U0_norm,in_cond,tf);


%% Plot time-course of all VPs in comparison to mean and num_STD*STD bounds
figure; 
fill(x_Nstd, inBetween_Nstd, grayColor); hold on;
plot(days,Traj_AccRej','o-'); 
plot(days,tumor_mean,'k','LineWidth',3); hold off;
ylim([0,inf]);
xlabel('Time (days)','FontSize',16)
ylabel('Tumor Volume (mm^3)','FontSize',16)
title('Accept-or-Reject Vpops within 3 STD of Mean','FontSize',16)

save Vpops.mat curve1_Nstd curve2_Nstd r_norm beta_norm U0_norm Vpop_AccRej Traj_AccRej
    

%%% Functions %%%
function [sol1, sol2, sol3] = OV_model(p,in_cond,tf)
    %% Experimental protocol is to give 3 doses of OV
    V0 = in_cond(3);
    % First virus dose
    sol1 = ode23s(@odemodel,[0 2], in_cond);
    % Second virus dose
    in_cond = sol1.y(:,end)' + [0    0    V0];
    sol2 = ode23s(@odemodel,[2 4], in_cond);
    % Third virus dose
    in_cond = sol2.y(:,end)' + [0    0    V0];
    sol3 = ode23s(@odemodel,[4 tf], in_cond);
    
    %% OV model
    function dydt = odemodel(t,y)
        U  = y(1); 
        I  = y(2);
        V  = y(3);
        % Frequency-dependent killing
        N = U + I;
        dU  = p.r*U - p.beta*U*V/N;
        dI  = p.beta*U*V/N - p.d_I*I;
        dV  = p.alpha*p.d_I*I - p.d_V*V; 
        dydt = [dU; dI; dV]; 
    end
end

function [Vpop_AccRej,Traj_AccRej] = accept_or_reject(N_Vpops,curve1,curve2,p,...
                                        r_norm,beta_norm,U0_norm,in_cond,tf)                                          
    i = 1;
    count_VPs = 0;
    r_VP = [];
    beta_VP = []; 
    U0_VP = [];
    Nvirt = length(r_norm);
    VP_sample = zeros(Nvirt,length(0:1:tf)); 
    VP_trajectory = [];
    VP_attempts = zeros(1,2); 
    while count_VPs < N_Vpops
        fprintf('Virtual patient number %d has r = %f, beta = %f, U0 = %f\n',...
            i,r_norm(i), beta_norm(i),U0_norm(i)); 
        VP_attempts(1) = VP_attempts(1) + 1;
        
        if((r_norm(i)>=0)&&(beta_norm(i)>=0)&&(U0_norm(i)>=0))
            in_cond(1) = U0_norm(i);
            p.r = r_norm(i);
            p.beta = beta_norm(i);
            p.U0 = U0_norm(i);
            [sol1, sol2, sol3] = OV_model(p,in_cond,tf);
            % Combine data
            sol.x = [sol1.x    sol2.x    sol3.x];
            sol.y = [sol1.y    sol2.y    sol3.y]; 
            modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) deval(sol3, 4:tf, 1)];
            modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) deval(sol3, 4:tf, 2)];    
            VP_sample = (modeldata_U + modeldata_I)'; % Total tumor volume

            accept = 0; 
            for j = 1:length(VP_sample)
                if (VP_sample(j)<curve2(j))||(VP_sample(j)>curve1(j)) 
                    accept = -1;
                end
            end
            params = [p.r p.beta p.U0];
            if accept == 0 % add VP
                count_VPs = count_VPs+1; 
                fprintf('\t%d. Accepting virtual patient\n',i);
                Vpop_AccRej(count_VPs,:) = params;
                Traj_AccRej(count_VPs,:) = VP_sample;
            else
                fprintf('\t%d. Outside bounds: rejecting virtual patient with:\n',i);
                fprintf('\t\tr = %f, beta = %f, U0 = %f\n',params(1),params(2),params(3));
            end
        end
        i = i+1;
    end                                                
end

function [Vpop_AccPert,Traj_AccPert] = accept_or_perturb(curve1,curve2,lb,ub,...
                                                    p,Vpop_params,in_cond,tf)
    count_VPs = 0;
    Vpop_AccPert = []; Traj_AccPert =[];
    for i = 1:size(Vpop_params,2)
        p.r = Vpop_params(1,i);
        p.beta = Vpop_params(2,i);
        p.d_I = Vpop_params(3,i);
        %p.alpha is set globally
        p.d_V = Vpop_params(4,i);
        p.U0 = Vpop_params(5,i);
        in_cond(1) = p.U0;
        [sol1, sol2, sol3] = OV_model(p,in_cond,tf);
        % Combine data
        sol.x = [sol1.x    sol2.x    sol3.x];
        sol.y = [sol1.y    sol2.y    sol3.y]; 
        modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) deval(sol3, 4:tf, 1)];
        modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) deval(sol3, 4:tf, 2)];    
        VP_sample = (modeldata_U + modeldata_I)'; % Total tumor volume
                
        accept = 0; 
        for j = 1:length(VP_sample)
            if (VP_sample(j)<curve2(j))||(VP_sample(j)>curve1(j)) 
                accept = -1;
            end
        end
        params = [p.r p.beta p.d_I p.alpha p.d_V p.U0];
        if accept == 0 % add VP
            count_VPs = count_VPs+1; 
            fprintf('\t%d. Accepting virtual patient\n',i);
            Vpop_AccPert(count_VPs,:) = params;
            Traj_AccPert(count_VPs,:) = VP_sample;
        else
            %% Instead of rejecting, any ones out of bounds get moved into
            %% bounds by finding r,d minimizing g(r,d) = as defined in chapter
            fprintf('\t%d. Perturbing parameters for virtual patient\n',i);
            p0 = [p.r p.beta p.d_I p.d_V p.U0]; % parameters without alpha
            fun = @(z)objective(z,curve2,curve1,p,in_cond,tf);
            p_simAnn_bnds = simulannealbnd(fun,p0,lb,ub) %,options); % remove options if don't want to see progression
            ptemp = zeros(1,length(p_simAnn_bnds)+1); 
            for j = 1:length(p_simAnn_bnds)
                if j<=3
                    ptemp(j) = p_simAnn_bnds(j); 
                else
                    ptemp(j+1) = p_simAnn_bnds(j); 
                end
            end
            ptemp(4) = p.alpha; 
            p.r = ptemp(1);
            p.beta = ptemp(2);
            p.d_I = ptemp(3);
            p.d_V = ptemp(4);
            p.U0 = ptemp(6);
            in_cond(1) = p.U0;
            [sol1, sol2, sol3] = OV_model(p,in_cond,tf);
            % Combine data
            sol.x = [sol1.x    sol2.x    sol3.x];
            sol.y = [sol1.y    sol2.y    sol3.y]; 
            modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) deval(sol3, 4:tf, 1)];
            modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) deval(sol3, 4:tf, 2)];  
            VP_sample  = (modeldata_U + modeldata_I)'; % Total tumor volume
            
            accept = 0; 
            for j = 1:length(VP_sample)
                if (VP_sample(j)<curve2(j))||(VP_sample(j)>curve1(j)) 
                    accept = -1;
                end
            end
            if accept == 0 % add VP
                count_VPs = count_VPs+1; 
                fprintf('\t\t%d. Accepting perturbed virtual patient\n',i);
                Vpop_AccPert(count_VPs,:) = ptemp;
                Traj_AccPert(count_VPs,:) = VP_sample;
            else
                fprintf('\t\t%d. Rejecting perturbed virtual patient\n',i);
            end
        end
    end
end

function value = objective(p0,l,u,p,in_cond,tf)
    p.r = p0(1);
    p.beta = p0(2);
    p.d_I = p0(3);
    p.d_V = p0(4);
    p.U0 =p0(5);
    in_cond(1) = p.U0;
    [sol1, sol2, sol3] = OV_model(p,in_cond,tf);
    sol.x = [sol1.x    sol2.x    sol3.x];
    sol.y = [sol1.y    sol2.y    sol3.y]; 
    modeldata_U = [deval(sol1, 0:1, 1) deval(sol2, 2:3, 1) deval(sol3, 4:tf, 1)];
    modeldata_I = [deval(sol1, 0:1, 2) deval(sol2, 2:3, 2) deval(sol3, 4:tf, 2)];    
    VP_sample  = (modeldata_U + modeldata_I)';
    
    sum = 0;
    for i = 1:length(VP_sample)
        term1 = (l(i)+u(i))/2; 
        term2 = (VP_sample(i)-term1)^2; 
        term3 = ( (u(i)/2) - (l(i)/2) )^2; 
        sum_term = (term2 - term3);
        sum = sum + max(sum_term,0); 
    end
    value = max(sum, 0);
end
