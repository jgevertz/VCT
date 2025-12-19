%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                       %
% This codes reads in the virtual populations with NVPops virtual       %
% patients for the simple OV model:                                     %
%   T' = r*T - beta*T*V,   dV' = -d_V*V                                 %
% generated using both the accept-or-reject and the accept-or-perturb   %
% methods.                                                              %
%                                                                       %
% Conducts a virtual clinical trial that looks at the impact of varying %
% dose of drug. VPs are then classified using the RECIST criteria:      %
% - Complete responder (CR): tumor undetectable                         %
% - Partial responder (PR): 30% decrease compared to baseline           %
% - Stable disease (SD): do not meet criteria of PR or PD               %
% - Progressive disease (PD): >=20% increase compared to baseline       %
%                                                                       %
% Output:                                                               %
% - Figure 1 displays the results in 2 subplots: one for accept-or-     %
%   reject, the other accept-or-perturb.                                %
% - Figure 2 displays the same results in 4 subplots, where each        %
%   subplot shows one response for the two VP methods.                  %
% - The label assigned to each VP is also saved in classify_acc_rej     %
%   (for accept-or-reject) and classify_acc_pert (for accept-or-        %
%   perturb), which is needed if we want to try to understand what      %
%   features explain treatment response.                                %
%                                                                       %
% Elapsed time is 127.631905 seconds (for dose_scaling = 1:1:40)        %
% Updated: 6/6/2025                                                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc; tic;
load data_AD.mat
tf = days(end); % stopping time
load best_fit_simple.mat % Contains lb, ub, p_fixed: % T0, V0
load Vpops_acc_rej_flip.mat  % Contains VP params in order
                             % r, beta, delta_V,
undetectable_threshold = 1; % small
T0 = p_fixed(1);
V0 = p_fixed(2); % experimental dose     
NVPops = size(Vpop,1);
dose_scaling = 1:1:40; 

%% Accept-or-reject VPs first
VPop_acc_rej = Vpop;
count_CR_acc_rej = zeros(size(dose_scaling));
count_PR_acc_rej = zeros(size(dose_scaling));
count_SD_acc_rej = zeros(size(dose_scaling));
count_PD_acc_rej = zeros(size(dose_scaling));
classify_acc_rej = zeros(NVPops,length(dose_scaling)); % 1 = CR, 2 = PR
                                                       % 3 = SD, 4 = PD
for j = 1:length(dose_scaling)
    p_fixed(2) = V0*dose_scaling(j); % scale dose
    %fprintf('For dose of %d\n',p_fixed(2));
    %figure; hold on;
    for i=1:NVPops
        [sol1, sol2, sol3] = solve_model(p_fixed, VPop_acc_rej(i,:) ,tf);
        sol.x = [sol1.x    sol2.x    sol3.x];
        sol.y = [sol1.y    sol2.y    sol3.y];
        %plot(sol.x, sol.y(1,:)+sol.y(2,:),'LineWidth',2);
        
        Tf = sol.y(1,end);
        %fprintf('\tVP %d has Tf = %f compared to T0 = %f => ',i,Tf,T0);
        % final volume  >=20% increase compared to baseline
        if Tf >= 1.2*T0
            count_PD_acc_rej(j) = count_PD_acc_rej(j) + 1;
            classify_acc_rej(i,j) = 4; % PD = 4
            %fprintf('PD\n');
        elseif Tf <= undetectable_threshold % tumor undetectable
            count_CR_acc_rej(j) = count_CR_acc_rej(j) + 1;
            classify_acc_rej(i,j) = 1; % CR = 1
            %fprintf('CR\n');
        elseif Tf <= 0.7*T0 % at least 30% decrease compared to baseline
            count_PR_acc_rej(j) = count_PR_acc_rej(j) + 1;
            classify_acc_rej(i,j) = 2; % PR = 2
            %fprintf('PR\n');
        else
            count_SD_acc_rej(j) = count_SD_acc_rej(j) + 1;
            classify_acc_rej(i,j) = 3; % SD = 3
            %fprintf('SD\n');
        end 
    end
    %hold off;
    fprintf('Accept-or-reject: For dose of %d:\n',p_fixed(2));
    fprintf('\tCR fraction = %f\n',count_CR_acc_rej(j)/NVPops);
    fprintf('\tPR fraction = %f\n',count_PR_acc_rej(j)/NVPops);
    fprintf('\tSD fraction = %f\n',count_SD_acc_rej(j)/NVPops);
    fprintf('\tPD fraction = %f\n',count_PD_acc_rej(j)/NVPops);
end


%% Accept-or-perturb VPs second
clear VPop;
load Vpops_acc_pert_flip.mat  % Contains VP params in order
                              % r, beta, delta_V
VPop_acc_pert = Vpop;
count_CR_acc_pert = zeros(size(dose_scaling));
count_PR_acc_pert = zeros(size(dose_scaling));
count_SD_acc_pert = zeros(size(dose_scaling));
count_PD_acc_pert = zeros(size(dose_scaling));
classify_acc_pert = zeros(NVPops,length(dose_scaling)); % 1 = CR, 2 = PR
                                                        % 3 = SD, 4 = PD
for j = 1:length(dose_scaling)
    p_fixed(2) = V0*dose_scaling(j); % scale dose
    %fprintf('For dose of %d\n',p_fixed(2));
    %figure; hold on;
    for i=1:NVPops
        [sol1, sol2, sol3] = solve_model(p_fixed, VPop_acc_pert(i,:) ,tf);
        sol.x = [sol1.x    sol2.x    sol3.x];
        sol.y = [sol1.y    sol2.y    sol3.y];
        %plot(sol.x, sol.y(1,:)+sol.y(2,:),'LineWidth',2);
        
        Tf = sol.y(1,end);
        %fprintf('\tVP %d has Tf = %f compared to T0 = %f => ',i,Tf,T0);
        % final volume  >=20% increase compared to baseline
        if Tf >= 1.2*T0
            count_PD_acc_pert(j) = count_PD_acc_pert(j) + 1;
            classify_acc_pert(i,j) = 4; % PD = 4
            %fprintf('PD\n');
        elseif Tf <= undetectable_threshold % tumor undetectable
            count_CR_acc_pert(j) = count_CR_acc_pert(j) + 1;
            classify_acc_pert(i,j) = 1; % CR = 1
            %fprintf('CR\n');
        elseif Tf <= 0.7*T0 % at least 30% decrease compared to baseline
            count_PR_acc_pert(j) = count_PR_acc_pert(j) + 1;
            classify_acc_pert(i,j) = 2; % PR = 2
            %fprintf('PR\n');
        else
            count_SD_acc_pert(j) = count_SD_acc_pert(j) + 1;
            classify_acc_pert(i,j) = 3; % SD = 3
            %fprintf('SD\n');
        end 
    end
    %hold off;
    fprintf('Accept-or-perturb: For dose of %d:\n',p_fixed(2));
    fprintf('\tCR fraction = %f\n',count_CR_acc_pert(j)/NVPops);
    fprintf('\tPR fraction = %f\n',count_PR_acc_pert(j)/NVPops);
    fprintf('\tSD fraction = %f\n',count_SD_acc_pert(j)/NVPops);
    fprintf('\tPD fraction = %f\n',count_PD_acc_pert(j)/NVPops);
end
toc;

figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.6, 0.6]);
subplot(1,2,1) % accept-or-reject
hold on;
plot(V0*dose_scaling,count_CR_acc_rej/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PR_acc_rej/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_SD_acc_rej/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PD_acc_rej/NVPops,'o-','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Accept-or-Reject','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(1,2,2) % accept-or-perturb
hold on;
plot(V0*dose_scaling,count_CR_acc_pert/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PR_acc_pert/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_SD_acc_pert/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PD_acc_pert/NVPops,'o-','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Accept-or-Perturb','FontSize',16)
legend('CR','PR','SD','PD','FontSize',16','location','best');
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;
saveas(gcf,'VCT_by_method_flip.fig');
exportgraphics(gcf,'VCT_by_method_flip.png');

figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.9]);
subplot(2,2,1) % CR
hold on;
plot(V0*dose_scaling,count_CR_acc_rej/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_CR_acc_pert/NVPops,'s--','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Complete Responders (CR)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,2) % PR
hold on;
plot(V0*dose_scaling,count_PR_acc_rej/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PR_acc_pert/NVPops,'s--','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Partial Responders (PR)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,3) % SD
hold on;
plot(V0*dose_scaling,count_SD_acc_rej/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_SD_acc_pert/NVPops,'s--','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Stable Disease (SD)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,4) % PD
hold on;
plot(V0*dose_scaling,count_PD_acc_rej/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PD_acc_pert/NVPops,'s--','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Progressive Disease (PD)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;
legend('Accept-or-Reject','Accept-or-Perturb','FontSize',16',...
    'location','NorthEast');
saveas(gcf,'VCT_by_classification_flip.fig');
exportgraphics(gcf,'VCT_by_classification_flip.png');

save VCT_flip.mat dose_scaling count_CR_acc_rej count_PR_acc_rej ...
    count_SD_acc_rej count_PD_acc_rej count_CR_acc_pert count_PR_acc_pert ...
    count_SD_acc_pert count_PD_acc_pert classify_acc_rej ...
    classify_acc_pert p_fixed VPop_acc_rej VPop_acc_pert

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
