%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% This code reads in the output from all models and creates plots that  %
% allow us to compare the output across the models.                     %
% - Figure 1: Just the VCT output from accept-or-reject VPs across      %
%   models.                                                             %
% - Figure 2: Just the VCT output from accept-or-perturb VPs across     %
%   models.                                                             %
% - Figure 3: VCT output for all VPs (generated from both methods)      %
%   across models.                                                      %
% - Figure 4: VP parameters for simple model across generation methods  %
%   (accept-or-reject versus accept-or-perturb).                        %
% - Figure 5: VP parameters for intermediate model across generation    %
%   methods (accept-or-reject versus accept-or-perturb).                %
% - Figure 6: VP parameters for complex model across generation methods %
%   (accept-or-reject versus accept-or-perturb).                        %
%                                                                       %
% Updated: 6/6/2025                                                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

load Simple_Model/data_AD.mat
tf = days(end); % stopping time
tumor_STD = std(tumorsizes,0,2);

load Simple_Model/best_fit_simple.mat
best_fit_params_simple = best_fit_params;
params_fixed_simple = p_fixed;
best_fit_simple = best_fit_objective; 
[sol_simple1, sol_simple2, sol_simple3] = ...
    solve_model_simple(params_fixed_simple, ...
    best_fit_params_simple, tf);
sol_simple.x = [sol_simple1.x  sol_simple2.x  sol_simple3.x];
sol_simple.y = [sol_simple1.y  sol_simple2.y  sol_simple3.y];

load Intermediate_Model/best_fit_intermediate.mat
best_fit_params_int = best_fit_params;
params_fixed_int = p_fixed;
best_fit_int = best_fit_objective; 
[sol_int1, sol_int2, sol_int3] = ...
    solve_model_int(params_fixed_int, ...
    best_fit_params_int, tf);
sol_int.x = [sol_int1.x  sol_int2.x  sol_int3.x];
sol_int.y = [sol_int1.y  sol_int2.y  sol_int3.y];

load Complex_Model/best_fit_complex.mat
best_fit_params_complex = best_fit_params;
params_fixed_complex = p_fixed;
best_fit_complex = best_fit_objective; 
[sol_complex1, sol_complex2, sol_complex3] = ...
    solve_model_complex(params_fixed_complex, ...
    best_fit_params_complex, tf);
sol_complex.x = [sol_complex1.x  sol_complex2.x  sol_complex3.x];
sol_complex.y = [sol_complex1.y  sol_complex2.y  sol_complex3.y];

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.85, 0.95]);
subplot(2,2,1)
plot(days,tumorsizes,'-o','LineWidth',2)
xlim([0,tf])
yl = ylim; 
title('Individual Mouse Trajectories','FontSize',16)
xlabel('Time (days)','FontSize',14)
ylabel('Tumor Volume (mm^3)','FontSize',14)

subplot(2,2,2)
errorbar(days, tumor_mean, tumor_STD,'or','LineWidth',2)
hold on;
plot(sol_simple.x, sol_simple.y(1,:),'b','LineWidth',3);
hold off;
xlim([0,tf]);
ylim([0,yl(end)]);
title(['Simple Model: Optimal Cost = ' num2str(best_fit_simple,'%.2f')],...
    'FontSize',16)
xlabel('Time (days)','FontSize',14)
ylabel('Virus Volume (mm^3)','FontSize',14)
legend('Model','Data','Location','NorthWest','FontSize',14)

subplot(2,2,3)
errorbar(days, tumor_mean, tumor_STD,'or','LineWidth',2)
hold on;
plot(sol_int.x, sol_int.y(1,:),'b','LineWidth',3);
hold off;
xlim([0,tf]);
ylim([0,yl(end)]);
title(['Intermediate Model: Optimal Cost = ' num2str(best_fit_int,'%.2f')],...
    'FontSize',16)
xlabel('Time (days)','FontSize',14)
ylabel('Virus Volume (mm^3)','FontSize',14)
legend('Model','Data','Location','NorthWest','FontSize',14)    

subplot(2,2,4)
errorbar(days, tumor_mean, tumor_STD,'or','LineWidth',2)
hold on;
plot(sol_complex.x, sol_complex.y(1,:),'b','LineWidth',3);
hold off;
xlim([0,tf]);
ylim([0,yl(end)]);
title(['Complex Model: Optimal Cost = ' num2str(best_fit_complex,'%.2f')],...
    'FontSize',16)
xlabel('Time (days)','FontSize',14)
ylabel('Virus Volume (mm^3)','FontSize',14)
legend('Model','Data','Location','NorthWest','FontSize',14)    

exportgraphics(gcf,'best_fit_all.png');
saveas(gcf,'best_fit_all.fig');

V0 = 10;
NVPops = 500;
load Simple_Model/VCT.mat
count_CR_acc_rej_simple = count_CR_acc_rej;
count_PR_acc_rej_simple = count_PR_acc_rej;
count_SD_acc_rej_simple = count_SD_acc_rej;
count_PD_acc_rej_simple = count_PD_acc_rej;
count_CR_acc_pert_simple = count_CR_acc_pert;
count_PR_acc_pert_simple = count_PR_acc_pert;
count_SD_acc_pert_simple = count_SD_acc_pert;
count_PD_acc_pert_simple = count_PD_acc_pert;

load Intermediate_Model/VCT.mat
count_CR_acc_rej_inter = count_CR_acc_rej;
count_PR_acc_rej_inter = count_PR_acc_rej;
count_SD_acc_rej_inter = count_SD_acc_rej;
count_PD_acc_rej_inter = count_PD_acc_rej;
count_CR_acc_pert_inter = count_CR_acc_pert;
count_PR_acc_pert_inter = count_PR_acc_pert;
count_SD_acc_pert_inter = count_SD_acc_pert;
count_PD_acc_pert_inter = count_PD_acc_pert;

load Complex_Model/VCT.mat
count_CR_acc_rej_complex = count_CR_acc_rej;
count_PR_acc_rej_complex = count_PR_acc_rej;
count_SD_acc_rej_complex = count_SD_acc_rej;
count_PD_acc_rej_complex = count_PD_acc_rej;
count_CR_acc_pert_complex = count_CR_acc_pert;
count_PR_acc_pert_complex = count_PR_acc_pert;
count_SD_acc_pert_complex = count_SD_acc_pert;
count_PD_acc_pert_complex = count_PD_acc_pert;


%% Accept-or-reject across models
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.9]);
sgtitle('Compare VCT Outcome by Model Complexity (Accept-or-Reject)',...
    'FontSize',16,'FontWeight','bold');
subplot(2,2,1) % CR
hold on;
plot(V0*dose_scaling,count_CR_acc_rej_simple/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_CR_acc_rej_inter/NVPops,'s--','LineWidth',2); 
plot(V0*dose_scaling,count_CR_acc_rej_complex/NVPops,'^-.','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Complete Responders (CR)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,2) % PR
hold on;
plot(V0*dose_scaling,count_PR_acc_rej_simple/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PR_acc_rej_inter/NVPops,'s--','LineWidth',2); 
plot(V0*dose_scaling,count_PR_acc_rej_complex/NVPops,'^-.','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Partial Responders (PR)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,3) % SD
hold on;
plot(V0*dose_scaling,count_SD_acc_rej_simple/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_SD_acc_rej_inter/NVPops,'s--','LineWidth',2); 
plot(V0*dose_scaling,count_SD_acc_rej_complex/NVPops,'^-.','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Stable Disease (SD)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,4) % PD
hold on;
plot(V0*dose_scaling,count_PD_acc_rej_simple/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PD_acc_rej_inter/NVPops,'s--','LineWidth',2); 
plot(V0*dose_scaling,count_PD_acc_rej_complex/NVPops,'^-.','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Progressive Disease (PD)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;
legend('Simple','Intermediate','Complex','FontSize',16',...
    'location','NorthEast');
saveas(gcf,'VCT_by_model_accept_reject.fig');
exportgraphics(gcf,'VCT_by_model_accept_reject.png');

%% Accept-or-perturb across models
figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.9]);
sgtitle('Compare VCT Outcome by Model Complexity (Accept-or-Perturb)',...
    'FontSize',16,'FontWeight','bold');
subplot(2,2,1) % CR
hold on;
plot(V0*dose_scaling,count_CR_acc_pert_simple/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_CR_acc_pert_inter/NVPops,'s--','LineWidth',2); 
plot(V0*dose_scaling,count_CR_acc_pert_complex/NVPops,'^-.','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Complete Responders (CR)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,2) % PR
hold on;
plot(V0*dose_scaling,count_PR_acc_pert_simple/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PR_acc_pert_inter/NVPops,'s--','LineWidth',2); 
plot(V0*dose_scaling,count_PR_acc_pert_complex/NVPops,'^-.','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Partial Responders (PR)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,3) % SD
hold on;
plot(V0*dose_scaling,count_SD_acc_pert_simple/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_SD_acc_pert_inter/NVPops,'s--','LineWidth',2); 
plot(V0*dose_scaling,count_SD_acc_pert_complex/NVPops,'^-.','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Stable Disease (SD)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,4) % PD
hold on;
plot(V0*dose_scaling,count_PD_acc_pert_simple/NVPops,'o-','LineWidth',2); 
plot(V0*dose_scaling,count_PD_acc_pert_inter/NVPops,'s--','LineWidth',2); 
plot(V0*dose_scaling,count_PD_acc_pert_complex/NVPops,'^-.','LineWidth',2); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Progressive Disease (PD)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;
legend('Simple','Intermediate','Complex','FontSize',16',...
    'location','NorthEast');
saveas(gcf,'VCT_by_model_accept_perturb.fig');
exportgraphics(gcf,'VCT_by_model_accept_perturb.png');

%% All together now
cmap = lines(6); 
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.9]);
sgtitle('Compare VCT Outcome by VPop Method and Model Complexity',...
    'FontSize',16,'FontWeight','bold');
subplot(2,2,1) % CR
hold on;
plot(V0*dose_scaling,count_CR_acc_rej_simple/NVPops,'o--','LineWidth',2,...
    'Color',cmap(1,:)); 
plot(V0*dose_scaling,count_CR_acc_rej_inter/NVPops,'s-.','LineWidth',2,...
    'Color',cmap(2,:)); 
plot(V0*dose_scaling,count_CR_acc_rej_complex/NVPops,'^:','LineWidth',2,...
    'Color',cmap(3,:)); 
plot(V0*dose_scaling,count_CR_acc_pert_simple/NVPops,'*--','LineWidth',2,...
    'Color',cmap(4,:)); 
plot(V0*dose_scaling,count_CR_acc_pert_inter/NVPops,'+-.','LineWidth',2,...
    'Color',cmap(5,:)); 
plot(V0*dose_scaling,count_CR_acc_pert_complex/NVPops,'x:','LineWidth',2,...
    'Color',cmap(6,:)); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Complete Responders (CR)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,2) % PR
hold on;
plot(V0*dose_scaling,count_PR_acc_rej_simple/NVPops,'o--','LineWidth',2,...
    'Color',cmap(1,:)); 
plot(V0*dose_scaling,count_PR_acc_rej_inter/NVPops,'s-.','LineWidth',2,...
    'Color',cmap(2,:)); 
plot(V0*dose_scaling,count_PR_acc_rej_complex/NVPops,'^:','LineWidth',2,...
    'Color',cmap(3,:)); 
plot(V0*dose_scaling,count_PR_acc_pert_simple/NVPops,'*--','LineWidth',2,...
    'Color',cmap(4,:)); 
plot(V0*dose_scaling,count_PR_acc_pert_inter/NVPops,'+-.','LineWidth',2,...
    'Color',cmap(5,:)); 
plot(V0*dose_scaling,count_PR_acc_pert_complex/NVPops,'x:','LineWidth',2,...
    'Color',cmap(6,:)); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Partial Responders (PR)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,3) % SD
hold on;
plot(V0*dose_scaling,count_SD_acc_rej_simple/NVPops,'o--','LineWidth',2,...
    'Color',cmap(1,:)); 
plot(V0*dose_scaling,count_SD_acc_rej_inter/NVPops,'s-.','LineWidth',2,...
    'Color',cmap(2,:)); 
plot(V0*dose_scaling,count_SD_acc_rej_complex/NVPops,'^:','LineWidth',2,...
    'Color',cmap(3,:)); 
plot(V0*dose_scaling,count_SD_acc_pert_simple/NVPops,'*--','LineWidth',2,...
    'Color',cmap(4,:)); 
plot(V0*dose_scaling,count_SD_acc_pert_inter/NVPops,'+-.','LineWidth',2,...
    'Color',cmap(5,:)); 
plot(V0*dose_scaling,count_SD_acc_pert_complex/NVPops,'x:','LineWidth',2,...
    'Color',cmap(6,:)); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Stable Disease (SD)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;

subplot(2,2,4) % PD
hold on;
plot(V0*dose_scaling,count_PD_acc_rej_simple/NVPops,'o--','LineWidth',2,...
    'Color',cmap(1,:)); 
plot(V0*dose_scaling,count_PD_acc_rej_inter/NVPops,'s-.','LineWidth',2,...
    'Color',cmap(2,:)); 
plot(V0*dose_scaling,count_PD_acc_rej_complex/NVPops,'^:','LineWidth',2,...
    'Color',cmap(3,:)); 
plot(V0*dose_scaling,count_PD_acc_pert_simple/NVPops,'*--','LineWidth',2,...
    'Color',cmap(4,:)); 
plot(V0*dose_scaling,count_PD_acc_pert_inter/NVPops,'+-.','LineWidth',2,...
    'Color',cmap(5,:)); 
plot(V0*dose_scaling,count_PD_acc_pert_complex/NVPops,'x:','LineWidth',2,...
    'Color',cmap(6,:)); 
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Fraction of Virtual Population','FontSize',16);
title('Progressive Disease (PD)','FontSize',16)
ylim([0,1.02]);
ax = gca;
ax.FontSize = 14;
legend('Accept-or-Reject: Simple','Accept-or-Reject: Intermediate',...
    'Accept-or-Reject: Complex','Accept-or-Perturb: Simple',...
    'Accept-or-Perturb: Intermediate','Accept-or-Perturb: Complex',...
    'FontSize',16, 'location','NorthEast');
saveas(gcf,'VCT_by_model_all.fig');
exportgraphics(gcf,'VCT_by_model_all.png');

%%%%
%% Now plot VP parameters for the same model, across methods
%%%%
load Simple_Model/Vpops_acc_rej.mat
Vpop_simple_acc_rej = Vpop; 
clear Vpop
load Simple_Model/Vpops_acc_pert.mat
Vpop_simple_acc_pert = Vpop; 
clear Vpop

%% Simple model
Nparams = size(Vpop_simple_acc_rej,2);
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.6]);
sgtitle('Simple Model: Parameter Distributions','fontSize',16,...
    'fontweight','bold');
for i = 1:Nparams
    subplot(1,3,i)
    histogram(Vpop_simple_acc_rej(:,i),'Normalization','pdf'); hold on; 
    histogram(Vpop_simple_acc_pert(:,i),'Normalization','pdf'); hold off;
    if i == 1
        xlabel('r','FontSize',16)
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
        legend('Accept-or-Reject','Accept-or-Perturb','FontSize',16,...
            'location','northwest');
    end
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,'compare_VPops_simple.fig');
exportgraphics(gcf,'compare_VPops_simple.png');

%% Intermediate model
load Intermediate_Model/Vpops_acc_rej.mat
Vpop_inter_acc_rej = Vpop; 
clear Vpop
load Intermediate_Model/Vpops_acc_pert.mat
Vpop_inter_acc_pert = Vpop; 
clear Vpop

Nparams = size(Vpop_inter_acc_rej,2);
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.55, 0.8]);
sgtitle('Intermediate Model: Parameter Distributions','fontSize',16,...
    'fontweight','bold');
for i = 1:Nparams
    subplot(2,2,i)
    hold on; 
    if i == 2
        histogram(Vpop_inter_acc_rej(:,i)/1000,'Normalization','pdf'); 
        histogram(Vpop_inter_acc_pert(:,i)/1000,'Normalization','pdf');
    elseif i == 4
        histogram(Vpop_inter_acc_rej(:,i)*1000,'Normalization','pdf'); 
        histogram(Vpop_inter_acc_pert(:,i)*1000,'Normalization','pdf');
    else
        histogram(Vpop_inter_acc_rej(:,i),'Normalization','pdf'); 
        histogram(Vpop_inter_acc_pert(:,i),'Normalization','pdf');
    end

    hold off;
    if i == 1
        xlabel('r','FontSize',16)
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    elseif i == 4
        xlabel('\alpha','FontSize',16)
        legend('Accept-or-Reject','Accept-or-Perturb','FontSize',16,...
            'location','northeast');
    end
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,'compare_VPops_intermediate.fig');
exportgraphics(gcf,'compare_VPops_intermediate.png');

%% Complex model
load Complex_Model/Vpops_acc_rej.mat
Vpop_complex_acc_rej = Vpop; 
clear Vpop
load Complex_Model/Vpops_acc_pert.mat
Vpop_complex_acc_pert = Vpop; 
clear Vpop

Nparams = size(Vpop_complex_acc_rej,2);
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle('Complex Model: Parameter Distributions','fontSize',16,...
    'fontweight','bold');
for i = 1:Nparams
    subplot(2,4,i)
    hold on; 
    if i == 2
        histogram(Vpop_complex_acc_rej(:,i)/1000,'Normalization','pdf'); 
        histogram(Vpop_complex_acc_pert(:,i)/1000,'Normalization','pdf');
    elseif i == 4
        histogram(Vpop_complex_acc_rej(:,i)*1000,'Normalization','pdf'); 
        histogram(Vpop_complex_acc_pert(:,i)*1000,'Normalization','pdf');
    else
        histogram(Vpop_complex_acc_rej(:,i),'Normalization','pdf'); 
        histogram(Vpop_complex_acc_pert(:,i),'Normalization','pdf');
    end

    hold off;
    if i == 1
        xlabel('r','FontSize',16)
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    elseif i == 4
        xlabel('\alpha','FontSize',16)
    elseif i == 5
        xlabel('\gamma_I','FontSize',16)
    elseif i == 6
        xlabel('\epsilon (\gamma_U = \epsilon\times \gamma_I)','FontSize',16)
    elseif i == 7
        xlabel('\lambda','FontSize',16)
        legend('Accept-or-Reject','Accept-or-Perturb','FontSize',16,...
            'location','northwest');
    end
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,'compare_VPops_complex.fig');
exportgraphics(gcf,'compare_VPops_complex.png');


function [sol1, sol2, sol3] = solve_model_simple(p_fixed,p_fit,tf)
    % Fixed
    T0 = p_fixed(1); 
    V0 = p_fixed(2);

    % First virus dose
    ICs = [T0    V0];
    sol1 = ode23s(@(t,x) OV_model_simple(t,x,p_fixed,p_fit), [0 2], ICs);
    % Second virus dose
    ICs = sol1.y(:,end)' + [0   V0];
    sol2 = ode23s(@(t,x) OV_model_simple(t,x,p_fixed,p_fit), [2 4], ICs);
    % Third virus dose
    ICs = sol2.y(:,end)' + [0    V0];
    sol3 = ode23s(@(t,x) OV_model_simple(t,x,p_fixed,p_fit), [4 tf], ICs);
end

function dydt = OV_model_simple(t,y,p_fixed,p_fit)
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

function [sol1, sol2, sol3] = solve_model_int(p_fixed,p_fit,tf)
    % Fixed
    U0 = p_fixed(2); 
    I0 = p_fixed(3);
    V0 = p_fixed(4);

    % First virus dose
    ICs = [U0    I0    V0];
    sol1 = ode23s(@(t,x) OV_model_int(t,x,p_fixed,p_fit), [0 2], ICs);

    % Second virus dose
    ICs = sol1.y(:,end)' + [0    0    V0];
    sol2 = ode23s(@(t,x) OV_model_int(t,x,p_fixed,p_fit), [2 4], ICs);

    % Third virus dose
    ICs = sol2.y(:,end)' + [0    0    V0];
    sol3 = ode23s(@(t,x) OV_model_int(t,x,p_fixed,p_fit), [4 tf], ICs);
end

%% OV model
function dydt = OV_model_int(t,y,p_fixed,p_fit)
    % Fixed
    d_I = p_fixed(1);

    % To fit
    r = p_fit(1);
    beta = p_fit(2);
    d_V = p_fit(3);
    alpha = p_fit(4);
    
    % Variables: enforce non-negativity
    U  = max(y(1),0); 
    I  = max(y(2),0);
    V  = max(y(3),0);

    % Frequency-dependent killing: avoid division by 0 on freq-dep terms
    N = U + I;
    dU  = r*U - max(0,beta*U*V/N); 
    dI  = max(0,beta*U*V/N) - d_I*I;
    dV  = alpha*d_I*I - d_V*V; 
    dydt = [dU; dI; dV]; 
end

function [sol1, sol2, sol3] = solve_model_complex(p_fixed,p_fit,tf)
    ode_opts = odeset('RelTol',5e-3,'AbsTol',5e-5); 
    % Fixed
    U0 = p_fixed(3); 
    I0 = p_fixed(4);
    V0 = p_fixed(5);
    T0 = p_fixed(6);

    % First virus dose
    ICs = [U0    I0    V0   T0];
    sol1 = ode15s(@(t,x) OV_model_complex(t,x,p_fixed,p_fit), [0 2], ICs,ode_opts);

    % Second virus dose
    ICs = sol1.y(:,end)' + [0    0    V0    0];
    sol2 = ode15s(@(t,x) OV_model_complex(t,x,p_fixed,p_fit), [2 4], ICs,ode_opts);

    % Third virus dose
    ICs = sol2.y(:,end)' + [0    0    V0    0];
    sol3 = ode15s(@(t,x) OV_model_complex(t,x,p_fixed,p_fit), [4 tf], ICs,ode_opts);
end

%% OV model
function dydt = OV_model_complex(t,y,p_fixed,p_fit)
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