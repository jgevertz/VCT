%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Implementation of two commonly-used virtual population methods with   %
% with uniform and normal prior on toy model:                           %
% x' = rx(1-x)-dx, x(0) = x0.                                           %
% - Code reads in simulated patient data generated from cellular        %
%   automaton model developed by K. Storie et al (Reference:            %
%   https://link.springer.com/chapter/10.1007/978-3-030-57129-0_8)      %
% - First it fits toy model to average of the data to establish a best- %
%   fit value for the growth rate r, kill rate d, and IC x0             %
% - VP Method 1 ('Accept-or-Reject'): Generate value of r and d from    %
%   uniform distribution with min = param_std below mean (or            %
%   0 if this is negative) and max = param_std above mean Keep IC fixed %
%   at best-fit value of x0. Then with a normal prior with mean = mean  %
%   of r in fitting patient data and standard deviation = param_std*mean%
%   VP only accepted if tumor trajectory falls within Nstd of mean      %
%   trajectory in data.                                                 %
%                                                                       %
% - VP Method 2 ('Accept-or-Perturb'): Generate value of r and d from   %
%   uniform and normal prior distributions.                             %
%   A VP is automatically accepted if its trajectory falls within Nstd  %
%   of mean trajectory.                                                 %
%   If it does not, simulated annealing is used to evolve               %
%   the parameters (according to a specific objective function) so that %
%   the trajectory does fall within Nstd of mean trajectory in data.    %
% - Figures 4, 5 and 6 for                                              %
%    Assessing the Role of Patient Generation Techniques in             %
%    Virtual Clinical Trial Outcomes                                    %
% - Figure 7 data generated using:                                      %
%    (a) beta = param_std = [0.25, 0.75, 1.25, 1.75, 2.25], N = 100,    % 
%        alpha = Nstd = 3, 10 replicates                                %
%    (b) alpha = Nstd = [1,2,3,4], N = 100, beta = 0.75, 10 replicates  %
%    (d) alpha 3, beta = 0.75, NVirt = [100, 500, 1000, 5000, 10000]    %
%                                                                       %
%   Authors: Jana L Gevertz and Joanna R Wares                          %
%   Last Updated: 5/29/2024                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc;
set(0, 'DefaultAxesFontSize', 14);
rng(1); % fix seed
dfile ='Toy_VCT_flip.txt';
if exist(dfile, 'file')
    delete(dfile); 
end
diary(dfile)
diary on
tic

%% Features of Vpop generation to control
Nvirt = 1000; % target number of virtual patients
Nstd = 3; % number of standard deviations from mean for accepted trajectory (alpha in paper)
param_std = 0.75; % standard deviation of parameter is computed as param_std*param_mean (beta in paper)


%% Read in simulated data from cellular automaton model
load CA_tumors.mat
month = 1:31;
time_month = time(month)-time(1); % adjust time so t = 0 is start of treatment
[max_volume, max_idx] = max(mean(volume,2));
vol_sample = volume(month+max_idx,:);  % further adjust so t = 0 aligns with beginning of treatment working
tf = time_month(end); 
vol_mean = mean(vol_sample,2); 
vol_std = std(vol_sample,0,2);

%% Visualize data that has been read in 
figure; hold on;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.4, 0.65]);
count = 1;
for i = idx'
    plot(time_month,vol_sample(:,i),'o-','LineWidth',2,....
        'Color',cmap(count,:),'DisplayName',...
        sprintf('Dose = %.2f',dose_contin_vary(i)) );
    count = count+1;
end
errorbar(time_month,vol_mean,vol_std,'o-','LineWidth',2,'Color',[0 0 0 0.5]); 
hold off; 
h = colorbar;
set(h, 'XTick', 0.01:0.01:0.1)
set(gca, 'CLim', [0.01,0.1])
ylim([0,inf])
xlabel('Time (days)','FontSize',16)
ylabel('Relative volume','FontSize',16)
title([num2str(totRep) ' simulated tumors with dose uniformly distributed in [' ...
    num2str(dose_min) ', ' num2str(dose_max) ']'],'FontSize',16);
fname_fig = 'simulated_tumors_truncated';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);


%% Fit to average of simulated data from CA model
p0 = [0.1 0.1 0.1]; % Guess for starting values: r, d, x0
lb = [0 0 0];
ub = [1 1 1]; 
options2 = optimoptions(@simulannealbnd,'MaxFunctionEvaluations', 3*10000); %,...
   % 'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping});
fun = @(z)objective2(z,time_month,vol_mean,vol_std);
[best_fit, best_fit_objective] = simulannealbnd(fun,p0,lb,ub,options2);
fprintf('Best fit of %f with r = %f, d = %f, x0 = %f\n',best_fit_objective,...
    best_fit(1),best_fit(2),best_fit(3)); 

%% Count number of synthetic datasets for which tumor grows/shrinks
data_traj_start = vol_sample(1,:);
data_traj_end = vol_sample(end,:);
data_shrink = sum(data_traj_end<data_traj_start);
data_grow = sum(data_traj_end>data_traj_start);
fprintf('Data: Of %d synthetic tumors, %d grow and %d shrink\n',size(vol_sample,2),...
    data_grow,data_shrink); 

%% Display best fit
sol_avg = ode23s(@(t,x) logistic(t,x,best_fit(1),best_fit(2)),time_month,best_fit(3));
figure; errorbar(time_month,vol_mean,vol_std,'o','LineWidth',2); hold on;
plot(sol_avg.x,sol_avg.y,'LineWidth',3); hold off
xlabel('Time (days)','FontSize',16);
ylabel('Relative tumor size','FontSize',16);
legend('Average data','Model best-fit','Location','NorthEast','FontSize',16); 
fname_fig = 'best_fit_to_avg';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Visualize all trajectories within 1 standard deviation of mean
curve1 = vol_mean' + vol_std'; % must be size 1xnumDays
curve2 = vol_mean' - vol_std';
x2 = [time_month', fliplr(time_month')];  % time pt vector must be of size 1xnumDays
inBetween = [curve1, fliplr(curve2)];
grayColor = [.7 .7 .7];
figure;
fill(x2, inBetween, grayColor); hold on;
plot(time_month,vol_sample,'o'); 
plot(time_month,vol_mean,'k','LineWidth',3); hold off;
ylim([0,1]);
xlabel('Time (days)','FontSize',16)
ylabel('Relative tumor size','FontSize',16)
title('Patient Data (Colors), Average (Black), 1 STD (Grey)','FontSize',16)
fname_fig = 'patient_data_1STD_fromMean';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Visualize all trajectories within 3 standard deviations of mean
curve1_Nstd = vol_mean' + Nstd*vol_std'; % must be size 1xnumDays
curve2_Nstd = vol_mean' - Nstd*vol_std';
x_Nstd = [time_month', fliplr(time_month')]; % time pt vector must be of size 1xnumDays
inBetween_Nstd = [curve1_Nstd, fliplr(curve2_Nstd)];
figure; 
fill(x_Nstd, inBetween_Nstd, grayColor); hold on;
plot(time_month,vol_sample,'o'); 
plot(time_month,vol_mean,'k','LineWidth',3); hold off;
ylim([0,1]);
xlabel('Time (days)','FontSize',16)
ylabel('Relative tumor size','FontSize',16)
title('Average (Black), 3 STD (Grey)','FontSize',16)
fname_fig = 'patient_data_3STD_fromMean';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png']);

%% Generate virtual patient parameters
r_mean = best_fit(1); %mean(r);
d_mean = best_fit(2); %mean(d); 
r_std = param_std*r_mean; %std(r);
d_std = param_std*d_mean; %std(d);

r_min = max(0,r_mean-3*r_std); % within 3 standard deviations to "match" normal
r_max = r_mean+3*r_std;        % within 3 standard deviations to "match" normal
d_min = max(0,d_mean-3*d_std); % within 3 standard deviations to "match" normal
d_max = d_mean+3*d_std;        % within 3 standard deviations to "match" normal
lb = [r_min,d_min];
ub = [r_max,d_max]; 

%% Generate normal parameters
r_norm1 = r_std.*randn(5*Nvirt,1) + r_mean; % generate 5x more to ensure I have enough data after rejection
d_norm1 = d_std.*randn(5*Nvirt,1) + d_mean;  % generate 5x more to ensure I have enough data after rejection
r_norm2 = r_std.*randn(5*Nvirt,1) + r_mean; % generate 5x more to ensure I have enough data after rejection
d_norm2 = d_std.*randn(5*Nvirt,1) + d_mean;  % generate 5x more to ensure I have enough data after rejection

%% Generate uniformly distributed parameters
r_uniform1 = r_min + (r_max-r_min)*rand(5*Nvirt,1); 
d_uniform1 = d_min + (d_max-d_min)*rand(5*Nvirt,1); 
r_uniform2 = r_min + (r_max-r_min)*rand(5*Nvirt,1); 
d_uniform2 = d_min + (d_max-d_min)*rand(5*Nvirt,1); 
x0 = best_fit(3); % fix initial condition

%% Generate virtual patients using method 1 with FLIPPED distribution :
%% 1) assume normal distribution with std = 0.75*mean (instead of 0.2*mean as per chapter)
%% 2) provided no parameters negative, check if trajectory lies within 3 std

%% method 1 = ar
%% method 2 = ap
VP_attempts = zeros(1,4); 

%% create data for each method for each prior
for exp_ind = 1:4
    if exp_ind == 1
        %% accept-or-reject with uniform
        r_params1 = r_uniform1
        d_params1 = d_uniform1

    elseif exp_ind == 2
        %% accept-or-reject with normal       
        r_params2 = r_norm1
        d_params2 = d_norm1

    elseif exp_ind == 3
        %% accept-or-perturb with uniform
        r_params3 = r_uniform2
        d_params3 = d_uniform2
        
    elseif exp_ind == 4
        %% accept-or-perturb with normal
        r_params4 = r_norm2
        d_params4 = d_norm2
    end

    %% Accept-or-Reject uniform
    if exp_ind == 1 
        i = 1;
        count_VPs = 0;
        r_VP1 = zeros(Nvirt,1);
        d_VP1 = zeros(Nvirt,1);
        VP_sample = zeros(Nvirt,length(time_month)); 

        %% traj 1 is AR with uniform
        VP_trajectory1 =zeros(Nvirt,length(time_month)); 
        while count_VPs < Nvirt
                fprintf('Virtual patient number %d has r = %f, d = %f\n',i,r_params1(i),d_params1(i)); 
                VP_attempts(1) = VP_attempts(1) + 1;
                if((r_params1(i)>=0)&&(d_params1(i)>=0))
                    sol_VP = ode23s(@(t,x) logistic(t,x,r_params1(i),d_params1(i)),time_month,x0);
                    VP_sample(i,:) = deval(sol_VP,time_month,1); 
                    reject = 0; 
                    for j = 1:length(VP_sample(i,:))
                        %if (VP_sample(i,j)<curve2(j))||(VP_sample(i,j)>curve1(j)) % 1 STD
                        if (VP_sample(i,j)<curve2_Nstd(j))||(VP_sample(i,j)>curve1_Nstd(j)) % 3 STD
                            reject = -1;
                        end
                    end
                    if reject == 0 % add VP
                        count_VPs = count_VPs+1; 
                        fprintf('\t%d. Accepting virtual patient\n',count_VPs); 
                        r_VP1(count_VPs) = r_params1(i);
                        d_VP1(count_VPs) = d_params1(i); 
                        VP_trajectory1(count_VPs,:) = VP_sample(i,:);
                    else
                        fprintf('\tOutside bounds: rejecting virtual patient with r = %f, d = %f\n',...
                            r_params1(i),d_params1(i)); 
                    end
                else
                    fprintf('\tNegative parameter: rejecting virtual patient with r = %f, d = %f\n',...
                        r_params1(i),d_params1(i)); 
                end     
                i = i+1;
        end
        
        %% Plot time-course of all VPs in comparison to mean and num_STD*STD bounds
        figure;
        fill(x_Nstd, inBetween_Nstd, grayColor); hold on; % 3 std
        plot(time_month,VP_trajectory1,'o');
        plot(time_month,vol_mean,'k','LineWidth',3); hold off;
        ylim([0,1]);
        xlabel('Time (days)','FontSize',16)
        ylabel('Relative tumor size','FontSize',16)
        title('Method 1: VPs (Colors), Average (Black), 3 STD (Grey)','FontSize',16)
        fname_fig = 'VP_trajectories_AccRej_Uniform';
        saveas(gcf,[fname_fig,'.fig'])
        saveas(gcf,[fname_fig,'.png']);

        %% Count number of VPs for which tumor grow or shrank
        %VP_traj_end = VP_trajectory1(:,end); 
    end 

    %% Accept-or-Reject normal
    if exp_ind == 2 
        i = 1;
        count_VPs = 0;
        r_VP2 = zeros(Nvirt,1);
        d_VP2 = zeros(Nvirt,1); 
        VP_sample2 = zeros(Nvirt,length(time_month)); 

        %% traj 2 is AR with normal
        VP_trajectory2 = zeros(Nvirt,length(time_month)); 
        count_reject = 0;
        while count_VPs < Nvirt
            fprintf('Virtual patient number %d has r = %f, d = %f\n',i,r_params2(i),d_params2(i)); 
            VP_attempts(2) = VP_attempts(2) + 1;
            
            if((r_params2(i)>=0)&&(d_params2(i)>=0))

                %% solve with the toy model to get the trajectory of tumor growth

                sol_VP = ode23s(@(t,x) logistic(t,x,r_params2(i),d_params2(i)),time_month,x0);
                VP_sample2(i,:) = deval(sol_VP,time_month,1); 
                reject = 0; 
                for j = 1:length(VP_sample2(i,:))
                    %if (VP_sample2(i,j)<curve2(j))||(VP_sample2(i,j)>curve1(j)) % 1 STD
                    if (VP_sample2(i,j)<curve2_Nstd(j))||(VP_sample2(i,j)>curve1_Nstd(j)) % 3 STD
                        reject = -1;                        
                    end
                end
                if reject == 0 % add VP
                    count_VPs = count_VPs+1; 
                    fprintf('\t%d. Accepting virtual patient\n',count_VPs); 

                    %% data with a 2 is AR with normal
                    r_VP2(count_VPs) = r_params2(i);
                    d_VP2(count_VPs) = d_params2(i); 
                    VP_trajectory2(count_VPs,:) = VP_sample2(i,:);
                else
                    fprintf('\tOutside bounds: rejecting virtual patient with r = %f, d = %f\n',...
                        r_params2(i),d_params2(i)); 
                    count_reject = count_reject + 1;
                end
            else
                fprintf('\tNegative parameter: rejecting virtual patient with r = %f, d = %f\n',...
                    r_params2(i),d_params2(i)); 
                count_reject = count_reject + 1;
            end     
            i = i+1;
        end
        
        %% Plot time-course of all VPs in comparison to mean and num_STD*STD bounds
        figure;
        fill(x_Nstd, inBetween_Nstd, grayColor); hold on; % 3 std
        plot(time_month,VP_trajectory2,'o');
        plot(time_month,vol_mean,'k','LineWidth',3); hold off;
        ylim([0,1]);
        xlabel('Time (days)','FontSize',16)
        ylabel('Relative tumor size','FontSize',16)
        title('Method 1 (Norm): VPs (Colors), Average (Black), 3 STD (Grey)','FontSize',16)
        fname_fig = 'VP_trajectories_AccRej_Uniform';
        saveas(gcf,[fname_fig,'.fig'])
        saveas(gcf,[fname_fig,'.png']);

        %% Count number of VPs for which tumor grow or shrank
        VP_traj_end2 = VP_trajectory2(:,end);
    end 
    
    %% Generate virtual patients using method 2 with FLIPPED distribution
    %% 1) assume normal distribution with std 0.75*mean not = 0.2*mean as per chapter
    %% 2) provided no parameters negative, check if trajectory lies within 3 std
    
    %% Accept-or-perturb with a uniform prior is labelled with a 3
    if exp_ind == 3 
        i = 1;
        count_VPs = 0;
        r_VP3 = zeros(Nvirt,1);
        d_VP3 = zeros(Nvirt,1);
        VP_sample = zeros(Nvirt,length(time_month)); 
        VP_trajectory3 = zeros(Nvirt,length(time_month)); 
        num_perturbed1 = 0; 
        while count_VPs < Nvirt
            fprintf('Virtual patient number %d has r = %f, d = %f\n',i,r_params3(i),d_params3(i)); 
            VP_attempts(3) = VP_attempts(3)+1; 
            if((r_params3(i)>=0)&&(d_params3(i)>=0))
                sol_VP = ode23s(@(t,x) logistic(t,x,r_params3(i),d_params3(i)),time_month,x0);
                VP_sample(i,:) = deval(sol_VP,time_month,1); 
                reject = 0; 
                for j = 1:length(VP_sample(i,:))
                    %if (VP_sample(i,j)<curve(j))||(VP_sample(i,j)>curve1(j)) % 1 STD
                    if (VP_sample(i,j)<curve2_Nstd(j))||(VP_sample(i,j)>curve1_Nstd(j)) % 3 STD
                        reject = -1;
                    end
                end
                if reject == 0 % add VP
                    count_VPs = count_VPs+1; 
                    fprintf('\t%d. Accepting virtual patient\n',count_VPs); 
                    r_VP3(count_VPs) = r_params3(i);
                    d_VP3(count_VPs) = d_params3(i); 
                    VP_trajectory3(count_VPs,:) = VP_sample(i,:); 
                else    
                    %% Instead of rejecting, any ones out of bounds get moved into 
                    %% bounds by finding r,d minimizing g(r,d) = as defined in chapter
                    %% Figure shows how the trajectory was perturbed - can remove plotting
                    %% if this data isn't needed (it is a lot of plots)
                    num_perturbed1 = num_perturbed1 + 1;

                    %% If you want to plot out of bounds trajectory
                   % if Nvirt <= 100
                    %    figure;
                     %   fill(x_Nstd, inBetween_Nstd, grayColor); hold on; % 3 std
                      %  plot(time_month,VP_sample(i,:),'s');
                    %end
                    p0 = [r_params3(i),d_params3(i)];
                    fun = @(z)objective(z,curve2_Nstd,curve1_Nstd,time_month,x0);
                    x_simAnn_bnds = simulannealbnd(fun,p0,[0,0],...
                        [r_mean+3*r_std,d_mean+3*d_std]); %,options); % remove options if don't want to see progression
                    fprintf('\tsimulannealbnd with lower and upper bound: r = %f and d = %f\n',...
                        x_simAnn_bnds(1),x_simAnn_bnds(2));
                    count_VPs = count_VPs+1; 
                    r_VP3(count_VPs) = x_simAnn_bnds(1); 
                    d_VP3(count_VPs) = x_simAnn_bnds(2);
                    sol_VP = ode23s(@(t,x) logistic(t,x,r_VP3(count_VPs),d_VP3(count_VPs)),time_month,x0);
                    VP_sample(i,:) = deval(sol_VP,time_month,1); 
                    VP_trajectory3(count_VPs,:) = VP_sample(i,:);
        
                    %% If you want to plot the out of range parameter trajectories
                    %if Nvirt <= 100
                     %   plot(time_month,VP_trajectory(count_VPs,:),'o');
                      %  plot(time_month,vol_mean,'k','LineWidth',3); hold off;
                       % ylim([0,1]);
                        %xlabel('Time','FontSize',16)
                        %ylabel('Relative tumor size','FontSize',16)
                        %title('Adjusting out-of-range VP','FontSize',16) 
                        %legend('3 STD','Out-of-range VP','Adjusted VP','Location','Northwest'); 
                    %end
                end
            else
               fprintf('\tNegative parameter: rejecting virtual patient with r = %f, d = %f\n',...
                   r_params3(i),d_params3(i));  
            end
            i = i+1;
        end
    end %% if method == 3

    %% Accept-or-perturb with a normal prior labelled with a 4
    if exp_ind == 4
        i = 1;
        count_VPs = 0;
        r_VP4 = zeros(Nvirt,1);
        d_VP4 = zeros(Nvirt,1);
        VP_sample = zeros(Nvirt,length(time_month)); 
        VP_trajectory4 = zeros(Nvirt,length(time_month)); 
        num_perturbed2 = 0; 
        while count_VPs < Nvirt
            fprintf('Virtual patient number %d has r = %f, d = %f\n',i,r_params4(i),d_params4(i)); 
            VP_attempts(4) = VP_attempts(4)+1; 
            if ((r_params4(i)>=0)&&(d_params4(i)>=0))
                sol_VP = ode23s(@(t,x) logistic(t,x,r_params4(i),d_params4(i)),time_month,x0);
                VP_sample(i,:) = deval(sol_VP,time_month,1); 
                reject = 0; 
                for j = 1:length(VP_sample(i,:))
                    %if (VP_sample(i,j)<curve2(j))||(VP_sample(i,j)>curve1(j)) % 1 STD
                    if (VP_sample(i,j)<curve2_Nstd(j))||(VP_sample(i,j)>curve1_Nstd(j)) % 3 STD
                        reject = -1;
                    end
                end
                if reject == 0 % add VP
                    count_VPs = count_VPs+1; 
                    fprintf('\t%d. Accepting virtual patient\n',count_VPs); 
                    r_VP4(count_VPs) = r_params4(i);
                    d_VP4(count_VPs) = d_params4(i); 
                    VP_trajectory4(count_VPs,:) = VP_sample(i,:);     
                else    
                    %% Instead of rejecting, any ones out of bounds get moved into 
                    %% bounds by finding r,d minimizing g(r,d) = as defined in chapter
                    %% Figure shows how the trajectory was perturbed - can remove plotting
                    %% if this data isn't needed (it is a lot of plots)
                    num_perturbed2 = num_perturbed2 + 1;
                    
                    %%Plot out of bounds trajectories
                    %if Nvirt <= 100
                     %   figure;
                      %  fill(x_Nstd, inBetween_Nstd, grayColor); hold on; % 3 std
                       % plot(time_month,VP_sample2(i,:),'s');
                    %end
                    p0 = [r_params4(i),d_params4(i)];
                    fun = @(z)objective(z,curve2_Nstd,curve1_Nstd,time_month,x0);
                    x_simAnn_bnds = simulannealbnd(fun,p0,[0,0],...
                        [r_mean+3*r_std,d_mean+3*d_std]); %,options); % remove options if don't want to see progression
                    fprintf('\tsimulannealbnd with lower and upper bound: r = %f and d = %f\n',...
                        x_simAnn_bnds(1),x_simAnn_bnds(2));
                    count_VPs = count_VPs+1; 
                    r_VP4(count_VPs) = x_simAnn_bnds(1); 
                    d_VP4(count_VPs) = x_simAnn_bnds(2);
                    sol_VP = ode23s(@(t,x) logistic(t,x,r_VP4(count_VPs),d_VP4(count_VPs)),time_month,x0);
                    VP_sample(i,:) = deval(sol_VP,time_month,1); 
                    VP_trajectory4(count_VPs,:) = VP_sample(i,:);
        
                    %% If you want to plot the out of range parameter trajectories
                    %if Nvirt <= 100
                     %   plot(time_month,VP_trajectory2(count_VPs,:),'o');
                      %  plot(time_month,vol_mean,'k','LineWidth',3); hold off;
                       % ylim([0,1]);
                        %xlabel('Time','FontSize',16)
                        %ylabel('Relative tumor size','FontSize',16)
                        %title('Adjusting out-of-range VP','FontSize',16) 
                        %legend('3 STD','Out-of-range VP','Adjusted VP','Location','Northwest'); 
                    %end
                end
            else
               fprintf('\tNegative parameter: rejecting virtual patient with r = %f, d = %f\n',...
                   r_params4(i),d_params4(i));  
            end
            i = i+1;
        end
    end 

    switch exp_ind
        case 1
            %% Histogram of distribution sampled per parameter
            x1 = -0.1:0.001:ub(1)+0.1;
            y1 = makedist('Uniform','lower',lb(1),'upper',ub(1));
            pdf1 = pdf(y1,x1);
            x2 = -0.1:0.001:ub(2)+0.1;
            y2 = makedist('Uniform','lower',lb(2),'upper',ub(2));
            pdf2 = pdf(y2,x2);
    
            %% Histogram of distribution sampled per parameter
            range = max([x1,x2]);
            x3 = r_mean-range:.001:r_mean+range;
            y3 = normpdf(x3,r_mean,r_std);
            x4 = d_mean-range:.001:d_mean+range;
            y4 = normpdf(x4,d_mean,d_std);
        
        case 2
            %% Histogram of distribution sampled per parameter
            x1 = -0.1:0.001:ub(1)+0.1;
            y1 = makedist('Uniform','lower',lb(1),'upper',ub(1));
            pdf1 = pdf(y1,x1);
            x2 = -0.1:0.001:ub(2)+0.1;
            y2 = makedist('Uniform','lower',lb(2),'upper',ub(2));
            pdf2 = pdf(y2,x2);
    
            %% Histogram of distribution sampled per parameter
            range = max([x1,x2]);
            x3 = r_mean-range:.001:r_mean+range;
            y3 = normpdf(x3,r_mean,r_std);
            x4 = d_mean-range:.001:d_mean+range;
            y4 = normpdf(x4,d_mean,d_std);
        
        case 3
            %% Histogram of distribution sampled per parameter
            x1 = -0.1:0.001:ub(1)+0.1;
            y1 = makedist('Uniform','lower',lb(1),'upper',ub(1));
            pdf1 = pdf(y1,x1);
            x3 = x1;
            y3 = pdf1;
    
            x2 = -0.1:0.001:ub(2)+0.1;
            y2 = makedist('Uniform','lower',lb(2),'upper',ub(2));
            pdf2 = pdf(y2,x2);
            x4 = x2;
            y4 = pdf2;
    
        case 4
            %% Histogram of distribution sampled per parameter
            range = max([x1,x2]);
            x1 = r_mean-range:.001:r_mean+range;
            y1 = normpdf(x1,r_mean,r_std);
            pdf1 = y1;
            x2 = d_mean-range:.001:d_mean+range;
            y2 = normpdf(x2,d_mean,d_std);
            pdf2 = y2;
    
            x3 = x1;
            y3 = y1;
            x4 = x2;
            y4 = y2;
    end
end %% end for loop for the exp_ind, choosing which experiment to run
  
VP_shrink = zeros(1,4);
VP_grow = zeros(1,4); 
%% Count number of VPs for which tumor grow or shrank
%% AR with uniform
VP_traj1_end = VP_trajectory1(:,end);
VP_shrink(1) = sum(VP_traj1_end<x0);
VP_grow(1) = sum(VP_traj1_end>x0);
fprintf('Method 1: Of %d VPs, %d grow and %d shrink\n',length(VP_traj1_end),...
    VP_grow(1),VP_shrink(1)); 

%% Count number of VPs for which tumor grow or shrank
%% AR with normal
VP_traj2_end = VP_trajectory2(:,end);
VP_shrink(2) = sum(VP_traj2_end<x0);
VP_grow(2) = sum(VP_traj2_end>x0);
fprintf('Method 2: Of %d VPs, %d grow and %d shrink\n',length(VP_traj2_end),...
    VP_grow(2),VP_shrink(2)); 

%% Count number of VPs for which tumor grow or shrank
%% AP with uniform
VP_traj3_end = VP_trajectory3(:,end);
VP_shrink(3) = sum(VP_traj3_end<x0);
VP_grow(3) = sum(VP_traj3_end>x0);
fprintf('Method 3: Of %d VPs, %d grow and %d shrink\n',length(VP_traj3_end),...
    VP_grow(3),VP_shrink(3)); 

%% Count number of VPs for which tumor grow or shrank
%% AP with normal
VP_traj4_end = VP_trajectory4(:,end);
VP_shrink(4) = sum(VP_traj4_end<x0);
VP_grow(4) = sum(VP_traj4_end>x0);
fprintf('Method 4: Of %d VPs, %d grow and %d shrink\n',length(VP_traj4_end),...
    VP_grow(4),VP_shrink(4)); 


%% Figure 4
%% Prior distribution plots
figure;
subplot(3,2,1);
%% r uniform 
x1 = -0.1:0.001:ub(1)+0.1;
y1 = makedist('Uniform','lower',lb(1),'upper',ub(1));
pdf1 = pdf(y1,x1);
%% d uniform
x2 = -0.1:0.001:ub(2)+0.1;
y2 = makedist('Uniform','lower',lb(2),'upper',ub(2));
pdf2 = pdf(y2,x2);
range = max([x1,x2]);
%% r normal
x3 = r_mean-range:.001:r_mean+range;
y3 = normpdf(x3,r_mean,r_std);
%% d normal
x4 = d_mean-range:.001:d_mean+range;
y4 = normpdf(x4,d_mean,d_std);
plot(x1,pdf1,'LineWidth',2); hold on;
plot(x3,y3,'LineWidth',2); hold off; 
xlim([0,1.5]);
xlabel('r','FontSize',16)

subplot(3,2,2);
plot(x2,pdf2,'LineWidth',2); hold on;
plot(x4,y4,'LineWidth',2); hold off; 
xlim([0,1.5]);
xlabel('d','FontSize',16)

%% Posterior Distributions
%% Accept-or-reject with uniform and normal priors for r
subplot(3,2,3)
histogram(r_VP1,0:0.05:1.5,'Normalization','pdf');  hold on; 
histogram(r_VP2,0:0.05:1.5,'Normalization','pdf'); hold off; 
xlabel('r','FontSize',16); 
legend('A-R(u)','A-R(n)','FontSize',16,'location','NorthEast');

%% Accept-or-reject with uniform and normal priors for d
subplot(3,2,4)
histogram(d_VP1,0:0.05:1.5,'Normalization','pdf');  hold on; 
histogram(d_VP2,0:0.05:1.5,'Normalization','pdf'); hold off; 
xlabel('r','FontSize',16); 
legend('A-R(u)','A-R(n)','FontSize',16,'location','NorthEast');

%% Accept-or-perturb with uniform and normal priors for r
subplot(3,2,5)
histogram(r_VP3,0:0.05:1.5,'Normalization','pdf');  hold on; 
histogram(r_VP4,0:0.05:1.5,'Normalization','pdf'); hold off; 
xlabel('r','FontSize',16); 
legend('A-P(u)','A-P(n)','FontSize',16,'location','NorthEast');

%% Accept-or-perturb with uniform and normal priors for d
subplot(3,2,6)
histogram(d_VP3,0:0.05:1.5,'Normalization','pdf');  hold on; 
histogram(d_VP4,0:0.05:1.5,'Normalization','pdf'); hold off; 
xlabel('d','FontSize',16); 
legend('A-P(u)','A-P(n)','FontSize',16,'location','NorthEast');

fname_fig = '6_dists_hists_1000';
saveas(gcf,[fname_fig,'.fig'])


%% Figure 5
figure;
%% Histogram AR and AP with uniform priors
subplot(3,2,1)
histogram(r_VP1,0:0.05:1.5,'Normalization','pdf');  hold on; 
histogram(r_VP3,0:0.05:1.5,'Normalization','pdf'); hold off; 
xlabel('r','FontSize',16); 
legend('A-R(u)','A-P(u)','FontSize',16,'location','NorthEast');

%% Scatter plot of AP and AR with uniform priors
subplot(3,2,2)
scatter(r_VP1,d_VP1,15,'filled'); hold on; 
scatter(r_VP3,d_VP3,15,'s','filled'); hold off; 
legend('A-R(u)','A-P(u)','FontSize',16,'location','SouthEast');
xlabel('r','FontSize',16); 
ylabel('d','FontSize',16); 

%% Historgram AR and AP with normal priors
subplot(3,2,3)
histogram(r_VP2,0:0.05:1.5,'Normalization','pdf');  hold on; 
histogram(r_VP4,0:0.05:1.5,'Normalization','pdf'); hold off; 
xlabel('r','FontSize',16); 
legend('A-R(n)','A-P(n)','FontSize',16,'location','NorthEast');

%% Scatter plot of AP and AR with normal priors
subplot(3,2,4)
scatter(r_VP2,d_VP2,15,'filled'); hold on; 
scatter(r_VP4,d_VP4,15,'s','filled'); hold off; 
legend('A-R(n)','A-P(n)','FontSize',16,'location','SouthEast');
xlabel('r','FontSize',16); 
ylabel('d','FontSize',16);

%% Histogram of AR with a normal prior and AP with a uniform prior
subplot(3,2,5)
histogram(r_VP2,0:0.05:1.5,'Normalization','pdf');  hold on; 
histogram(r_VP3,0:0.05:1.5,'Normalization','pdf'); hold off; 
xlabel('r','FontSize',16); 
legend('A-R(n)','A-P(u)','FontSize',16,'location','NorthEast');

%% Scatter plot of AR with a normal prior and AP with a uniform prior
subplot(3,2,6)
scatter(r_VP2,d_VP2,15,'filled'); hold on; 
scatter(r_VP3,d_VP3,15,'s','filled'); hold off; 
legend('A-R(n)','A-P(u)','FontSize',16,'location','SouthEast');
xlabel('r','FontSize',16); 
ylabel('d','FontSize',16);
fname_fig = '6_hist_scatter_1000';
saveas(gcf,[fname_fig,'.fig'])


%% Figure 6 
%% bar plot of responders
figure;
VP_predictions = VP_shrink/Nvirt;
b = bar(diag(VP_predictions),'stacked');
xtips1 = b(1).XEndPoints;
ytips1 = b(4).YEndPoints;
labels1 = string(ytips1);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xlabel('','FontSize',16)
ylabel('Probability of Tumor Shrinking','FontSize',16)
ylim([0,1])
xticklabels({'A-R (u)','A-R (n)','A-P (u)','A-P (n)'})
fname_fig = 'bar_1000';
saveas(gcf,[fname_fig,'.fig'])

toc
diary off
%% FIX to save num_perturbed if you want
save output.mat time_month vol_sample vol_mean vol_std lb ub best_fit ...
    best_fit_objective r_VP1 d_VP1 VP_trajectory1 VP_traj1_end VP_grow ...
    VP_shrink r_VP2 d_VP2 VP_trajectory2 VP_traj2_end ...
     r_VP3 d_VP3 VP_trajectory3 VP_traj3_end ...
     r_VP4 d_VP4 VP_trajectory4 VP_traj4_end ...
    VP_attempts %num_perturbed


%% Functions %%
%% Toy model
function xp = logistic(t,x,r,d)
 xp = r*x*(1-x)-d*x;
end

%% For VP Method #2
function value = objective(p0,l,u,time,x0)    
    sol_VP = ode23s(@(t,x) logistic(t,x,p0(1),p0(2)),time,x0);
    VP_sample = deval(sol_VP,time,1); 
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

%% Fit function (to be minimized)
function SSE_divideby_Var = objective2(p,time,vol_mean,vol_std)
    sol = ode23s(@(t,x) logistic(t,x,p(1),p(2)),time,p(3));
    
    % Evaluate numerical solution to ODE at experimental time points
    modeldata = deval(sol,time,1);

    % Goodness of fit: divide by 2 to match likelihood formula 
    SSE_divideby_Var = 0.5*sum(((vol_mean-modeldata').^2)./(vol_std.^2));
end


