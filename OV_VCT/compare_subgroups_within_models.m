%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% This code reads in the output from all models and looks at VP         %
% parameter distributions based on how the VP was classified. This is   %
% only done for one dose selected by the dose_idx variable. Can change  %
% this dose if want to analyze another dose. I tried to choose a dose   %
% for each where there was substantial variability, in patient          %
% response.                                                             %
% - Figure 1: Simple model with VP classified into one of four groups:  %
%   complete responder (CR), partial responder (PR), stable disease     %
%   (SD), or progressive disease (PD).                                  %
% - Figure 2: Simple model with VP classified as either non-progressive %
%   (defined as CR, PR or SD) or progressive (PD).                      %
% - Figure 3: Intermediate model with VP classified into one of four    %
%   groups: CR, PR, SD, or PD.                                          %
% - Figure 4: Intermediate model with VP classified as either           %
%   non-progressive or progressive (PD).                                %
% - Figure 5+6: Complex model with VP classified as either              %
%   non-progressive or progressive (PD).                                %
% - Figure 7+8: Complex model with VP classified as either              %
%   non-progressive or progressive (PD).                                %
%                                                                       %
% Updated: 6/6/2025                                                     %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
V0 = 10;

%% Simple model
load Simple_Model/VCT.mat
% Pick dose of 100 for comparison
dose_idx = 5; 
fprintf('Will compare simple model distributions at dose = %d\n',...
    V0*dose_scaling(dose_idx))
Vpop_simple_acc_rej = VPop_acc_rej; 
Vpop_simple_acc_pert = VPop_acc_pert; 
classify_simple_acc_rej = classify_acc_rej(:,dose_idx);
classify_simple_acc_pert = classify_acc_pert(:,dose_idx);
clear Vpop_acc_rej Vpop_acc_pert classify_acc_rej classify_acc_pert
NVPops = size(Vpop_simple_acc_rej,1);
Nparams = size(Vpop_simple_acc_rej,2);

% Compare all four subgroups
Nsubgroups = 4; % Key: 1 = CR, 2 = PR, 3 = SD, 4 = PD
VP_class_4sub_simple_acc_rej = cell(1,Nsubgroups); 
VP_params_4sub_simple_acc_rej = cell(1,Nsubgroups); 
VP_class_4sub_simple_acc_pert = cell(1,Nsubgroups); 
VP_params_4sub_simple_acc_pert = cell(1,Nsubgroups); 
for i = 1:Nsubgroups
    % Accept-or-reject
    idx = find(classify_simple_acc_rej==i);
    VP_class_4sub_simple_acc_rej{i} = idx;
    VP_params_4sub_simple_acc_rej{i} = Vpop_simple_acc_rej(idx,:);

    % Accept-or-perturb
    idx2 = find(classify_simple_acc_pert==i);
    VP_class_4sub_simple_acc_pert{i} = idx2;
    VP_params_4sub_simple_acc_pert{i} = Vpop_simple_acc_pert(idx2,:);
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle(['Simple Model: 4 Subgroup Parameter Distributions '...
    'at Dose of ' num2str(V0*dose_scaling(dose_idx)) ...
    'x10^9 Virions'],'fontSize',16,'fontweight','bold');
for i = 1:Nparams
    subplot(2,3,i)
    hold on; 
    for j = 1:Nsubgroups
        histogram(VP_params_4sub_simple_acc_rej{j}(:,i),...
            'Normalization','pdf'); 
    end
    hold off;
    if i == 1
        xlabel('r','FontSize',16)
        legend('CR','PR','SD','PD','FontSize',16,'location','northeast');
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    end
    title('Accept-or-Reject','FontSize',16);
    ax = gca;
    ax.FontSize = 14;

    subplot(2,3,i+Nparams)
    hold on; 
    for j = 1:Nsubgroups
        histogram(VP_params_4sub_simple_acc_pert{j}(:,i),...
            'Normalization','pdf'); 
    end
    hold off;
    if i == 1
        xlabel('r','FontSize',16)
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    end
    title('Accept-or-Perturb','FontSize',16);
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,['VPops_4sub_simple' num2str(dose_idx) '.fig']);
exportgraphics(gcf,['VPops_4sub_simple' num2str(dose_idx) '.png']);

% If we only want to look at progressive disease or not
Nsubgroups = 2; % Key: 1 = CR, 2 = PR, 3 = SD, 4 = PD
VP_class_2sub_simple_acc_rej = cell(1,Nsubgroups); 
VP_params_2sub_simple_acc_rej = cell(1,Nsubgroups); 
VP_class_2sub_simple_acc_pert = cell(1,Nsubgroups); 
VP_params_2sub_simple_acc_pert = cell(1,Nsubgroups); 

% Accept-or-reject
idxYes = find(classify_simple_acc_rej<4); % CR, PR, SD
idxNo = find(classify_simple_acc_rej==4); % PD
VP_class_2sub_simple_acc_rej{1} = idxYes;
VP_params_2sub_simple_acc_rej{1} = Vpop_simple_acc_rej(idxYes,:);
VP_class_2sub_simple_acc_rej{2} = idxNo;
VP_params_2sub_simple_acc_rej{2} = Vpop_simple_acc_rej(idxNo,:);

% Accept-or-perturb
idx2Yes = find(classify_simple_acc_pert<4); % CR, PR, SD
idx2No = find(classify_simple_acc_pert==4); % PD
VP_class_2sub_simple_acc_pert{1} = idx2Yes;
VP_params_2sub_simple_acc_pert{1} = Vpop_simple_acc_pert(idx2Yes,:); 
VP_class_2sub_simple_acc_pert{2} = idx2No;
VP_params_2sub_simple_acc_pert{2} = Vpop_simple_acc_pert(idx2No,:);

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle(['Simple Model: Binary Subgroup Parameter Distributions '...
    'at Dose of ' num2str(V0*dose_scaling(dose_idx)) ...
    'x10^9 Virions'],'fontSize',16,'fontweight','bold');
for i = 1:Nparams
    subplot(2,3,i)
    hold on; 
    for j = 1:Nsubgroups
        histogram(VP_params_2sub_simple_acc_rej{j}(:,i), ...
            'Normalization','pdf'); 
    end
    hold off;
    if i == 1
        xlabel('r','FontSize',16)
        legend('Non-PD','PD','FontSize',16,'location','northeast');
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    end
    title('Accept-or-Reject','FontSize',16);
    ax = gca;
    ax.FontSize = 14;

    subplot(2,3,i+Nparams)
    hold on; 
    for j = 1:Nsubgroups
        histogram(VP_params_2sub_simple_acc_pert{j}(:,i),...
            'Normalization','pdf'); 
    end
    hold off;
    if i == 1
        xlabel('r','FontSize',16)
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    end
    title('Accept-or-Perturb','FontSize',16);
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,['VPops_2sub_simple' num2str(dose_idx) '.fig']);
exportgraphics(gcf,['VPops_2sub_simple' num2str(dose_idx) '.png']);

%% Intermediate model
load Intermediate_Model/VCT.mat
% Pick dose of 50 for comparison
dose_idx = 5; 
fprintf('Will compare simple model distributions at dose = %d\n',...
    V0*dose_scaling(dose_idx))
Vpop_inter_acc_rej = VPop_acc_rej; 
Vpop_inter_acc_pert = VPop_acc_pert; 
classify_inter_acc_rej = classify_acc_rej(:,dose_idx);
classify_inter_acc_pert = classify_acc_pert(:,dose_idx);
clear Vpop_acc_rej Vpop_acc_pert classify_acc_rej classify_acc_pert
NVPops = size(Vpop_inter_acc_rej,1);
Nparams = size(Vpop_inter_acc_rej,2);

% Compare all four subgroups
Nsubgroups = 4; % Key: 1 = CR, 2 = PR, 3 = SD, 4 = PD
VP_class_4sub_inter_acc_rej = cell(1,Nsubgroups); 
VP_params_4sub_inter_acc_rej = cell(1,Nsubgroups); 
VP_class_4sub_inter_acc_pert = cell(1,Nsubgroups); 
VP_params_4sub_inter_acc_pert = cell(1,Nsubgroups); 
for i = 1:Nsubgroups
    % Accept-or-reject
    idx = find(classify_inter_acc_rej==i);
    VP_class_4sub_inter_acc_rej{i} = idx;
    VP_params_4sub_inter_acc_rej{i} = Vpop_inter_acc_rej(idx,:);

    % Accept-or-perturb
    idx2 = find(classify_inter_acc_pert==i);
    VP_class_4sub_inter_acc_pert{i} = idx2;
    VP_params_4sub_inter_acc_pert{i} = Vpop_inter_acc_pert(idx2,:);
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle(['Intermediate Model: 4 Subgroup Parameter Distributions '...
    'at Dose of ' num2str(V0*dose_scaling(dose_idx)) ...
    'x10^9 Virions'],'fontSize',16,'fontweight','bold');
for i = 1:Nparams
    subplot(2,4,i)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_4sub_inter_acc_rej{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_4sub_inter_acc_rej{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_4sub_inter_acc_rej{j}(:,i),...
                'Normalization','pdf'); 
        end
    end
    hold off;
    if i == 1
        xlabel('r','FontSize',16)
        legend('CR','PR','SD','PD','FontSize',16,'location','northeast');
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    elseif i == 4
        xlabel('\alpha','FontSize',16)
    end
    title('Accept-or-Reject','FontSize',16);
    ax = gca;
    ax.FontSize = 14;

    subplot(2,4,i+Nparams)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_4sub_inter_acc_pert{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_4sub_inter_acc_pert{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_4sub_inter_acc_pert{j}(:,i),...
                'Normalization','pdf'); 
        end
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
    end
    title('Accept-or-Perturb','FontSize',16);
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,['VPops_4sub_intermediate' num2str(dose_idx) '.fig']);
exportgraphics(gcf,['VPops_4sub_intermediate' num2str(dose_idx) '.png']);

% If we only want to look at progressive disease or not
Nsubgroups = 2; % Key: 1 = CR, 2 = PR, 3 = SD, 4 = PD
VP_class_2sub_inter_acc_rej = cell(1,Nsubgroups); 
VP_params_2sub_inter_acc_rej = cell(1,Nsubgroups); 
VP_class_2sub_inter_acc_pert = cell(1,Nsubgroups); 
VP_params_2sub_inter_acc_pert = cell(1,Nsubgroups); 

% Accept-or-reject
idxYes = find(classify_inter_acc_rej<4); % CR, PR, SD
idxNo = find(classify_inter_acc_rej==4); % PD
VP_class_2sub_inter_acc_rej{1} = idxYes;
VP_params_2sub_inter_acc_rej{1} = Vpop_inter_acc_rej(idxYes,:);
VP_class_2sub_inter_acc_rej{2} = idxNo;
VP_params_2sub_inter_acc_rej{2} = Vpop_inter_acc_rej(idxNo,:);

% Accept-or-perturb
idx2Yes = find(classify_inter_acc_pert<4); % CR, PR, SD
idx2No = find(classify_inter_acc_pert==4); % PD
VP_class_2sub_inter_acc_pert{1} = idx2Yes;
VP_params_2sub_inter_acc_pert{1} = Vpop_inter_acc_pert(idx2Yes,:); 
VP_class_2sub_inter_acc_pert{2} = idx2No;
VP_params_2sub_inter_acc_pert{2} = Vpop_inter_acc_pert(idx2No,:);

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle(['Intermediate Model: Binary Subgroup Parameter Distributions '...
    'at Dose of ' num2str(V0*dose_scaling(dose_idx)) ...
    'x10^9 Virions'],'fontSize',16,'fontweight','bold');
for i = 1:Nparams
    subplot(2,4,i)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_2sub_inter_acc_rej{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_2sub_inter_acc_rej{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_2sub_inter_acc_rej{j}(:,i),...
                'Normalization','pdf'); 
        end
    end
    hold off;
    if i == 1
        xlabel('r','FontSize',16)
        legend('Non-PD','PD','FontSize',16,'location','northeast');
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    elseif i == 4
        xlabel('\alpha','FontSize',16)
    end
    title('Accept-or-Reject','FontSize',16);
    ax = gca;
    ax.FontSize = 14;

    subplot(2,4,i+Nparams)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_2sub_inter_acc_pert{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_2sub_inter_acc_pert{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_2sub_inter_acc_pert{j}(:,i),...
                'Normalization','pdf'); 
        end
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
    end
    title('Accept-or-Perturb','FontSize',16);
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,['VPops_2sub_intermediate' num2str(dose_idx) '.fig']);
exportgraphics(gcf,['VPops_2sub_intermediate' num2str(dose_idx) '.png']);


%% Complex model
load Complex_Model/VCT.mat
% Pick dose of 50 for comparison
dose_idx = 5; 
fprintf('Will compare simple model distributions at dose = %d\n',...
    V0*dose_scaling(dose_idx))
Vpop_complex_acc_rej = VPop_acc_rej; 
Vpop_complex_acc_pert = VPop_acc_pert; 
classify_complex_acc_rej = classify_acc_rej(:,dose_idx);
classify_complex_acc_pert = classify_acc_pert(:,dose_idx);
clear Vpop_acc_rej Vpop_acc_pert classify_acc_rej classify_acc_pert
NVPops = size(Vpop_complex_acc_rej,1);
Nparams = size(Vpop_complex_acc_rej,2);

% Compare all four subgroups
Nsubgroups = 4; % Key: 1 = CR, 2 = PR, 3 = SD, 4 = PD
VP_class_4sub_complex_acc_rej = cell(1,Nsubgroups); 
VP_params_4sub_complex_acc_rej = cell(1,Nsubgroups); 
VP_class_4sub_complex_acc_pert = cell(1,Nsubgroups); 
VP_params_4sub_complex_acc_pert = cell(1,Nsubgroups); 
for i = 1:Nsubgroups
    % Accept-or-reject
    idx = find(classify_complex_acc_rej==i);
    VP_class_4sub_complex_acc_rej{i} = idx;
    VP_params_4sub_complex_acc_rej{i} = Vpop_complex_acc_rej(idx,:);

    % Accept-or-perturb
    idx2 = find(classify_complex_acc_pert==i);
    VP_class_4sub_complex_acc_pert{i} = idx2;
    VP_params_4sub_complex_acc_pert{i} = Vpop_complex_acc_pert(idx2,:);
end

Nplot1 = 4;
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle(['Complex Model (1 of 2): 4 Subgroup Parameter Distributions '...
    'at Dose of ' num2str(V0*dose_scaling(dose_idx)) ...
    'x10^9 Virions'],'fontSize',16,'fontweight','bold');
for i = 1:Nplot1
    subplot(2,4,i)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_4sub_complex_acc_rej{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_4sub_complex_acc_rej{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_4sub_complex_acc_rej{j}(:,i),...
                'Normalization','pdf'); 
        end
    end
    hold off;
    if i == 1
        xlabel('r','FontSize',16)
        legend('CR','PR','SD','PD','FontSize',16,'location','northeast');
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    elseif i == 4
        xlabel('\alpha','FontSize',16)
    end
    title('Accept-or-Reject','FontSize',16);
    ax = gca;
    ax.FontSize = 14;

    subplot(2,4,i+Nplot1)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_4sub_complex_acc_pert{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_4sub_complex_acc_pert{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_4sub_complex_acc_pert{j}(:,i),...
                'Normalization','pdf'); 
        end
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
    end
    title('Accept-or-Perturb','FontSize',16);
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,['VPops1_4sub_complex' num2str(dose_idx) '.fig']);
exportgraphics(gcf,['VPops1_4sub_complex' num2str(dose_idx) '.png']);


Nplot2 = Nparams - Nplot1;
figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle(['Complex Model (2 of 2): 4 Subgroup Parameter Distributions '...
    'at Dose of ' num2str(V0*dose_scaling(dose_idx)) ...
    'x10^9 Virions'],'fontSize',16,'fontweight','bold');
for i = Nplot1+1:Nparams
    subplot(2,3,i-Nplot1)
    hold on; 
    for j = 1:Nsubgroups
        histogram(VP_params_4sub_complex_acc_rej{j}(:,i),...
            'Normalization','pdf'); 
    end
    hold off;
    if i == 5
        xlabel('\gamma_I','FontSize',16)
        legend('CR','PR','SD','PD','FontSize',16,'location','northeast');
    elseif i == 6
        xlabel('\epsilon (\gamma_U = \epsilon\times \gamma_I)','FontSize',16)
    elseif i == 7
        xlabel('\lambda','FontSize',16)
    end
    title('Accept-or-Reject','FontSize',16);
    ax = gca;
    ax.FontSize = 14;

    subplot(2,3,i-Nplot1+Nplot2)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_4sub_complex_acc_pert{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_4sub_complex_acc_pert{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_4sub_complex_acc_pert{j}(:,i),...
                'Normalization','pdf'); 
        end
    end
    hold off;
    if i == 5
        xlabel('\gamma_I','FontSize',16)
    elseif i == 6
        xlabel('\epsilon (\gamma_U = \epsilon\times \gamma_I)','FontSize',16)
    elseif i == 7
        xlabel('\lambda','FontSize',16)
    end
    title('Accept-or-Perturb','FontSize',16);
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,['VPops2_4sub_complex' num2str(dose_idx) '.fig']);
exportgraphics(gcf,['VPops2_4sub_complex' num2str(dose_idx) '.png']);


% If we only want to look at progressive disease or not
Nsubgroups = 2; % Key: 1 = CR, 2 = PR, 3 = SD, 4 = PD
VP_class_2sub_complex_acc_rej = cell(1,Nsubgroups); 
VP_params_2sub_complex_acc_rej = cell(1,Nsubgroups); 
VP_class_2sub_complex_acc_pert = cell(1,Nsubgroups); 
VP_params_2sub_complex_acc_pert = cell(1,Nsubgroups); 

% Accept-or-reject
idxYes = find(classify_complex_acc_rej<4); % CR, PR, SD
idxNo = find(classify_complex_acc_rej==4); % PD
VP_class_2sub_complex_acc_rej{1} = idxYes;
VP_params_2sub_complex_acc_rej{1} = Vpop_complex_acc_rej(idxYes,:);
VP_class_2sub_complex_acc_rej{2} = idxNo;
VP_params_2sub_complex_acc_rej{2} = Vpop_complex_acc_rej(idxNo,:);

% Accept-or-perturb
idx2Yes = find(classify_complex_acc_pert<4); % CR, PR, SD
idx2No = find(classify_complex_acc_pert==4); % PD
VP_class_2sub_complex_acc_pert{1} = idx2Yes;
VP_params_2sub_complex_acc_pert{1} = Vpop_complex_acc_pert(idx2Yes,:); 
VP_class_2sub_complex_acc_pert{2} = idx2No;
VP_params_2sub_complex_acc_pert{2} = Vpop_complex_acc_pert(idx2No,:);

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle(['Complex Model (1 of 2): Binary Subgroup Parameter Distributions '...
    'at Dose of ' num2str(V0*dose_scaling(dose_idx)) ...
    'x10^9 Virions'],'fontSize',16,'fontweight','bold');
for i = 1:Nplot1
    subplot(2,4,i)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_2sub_complex_acc_rej{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_2sub_complex_acc_rej{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_2sub_complex_acc_rej{j}(:,i),...
                'Normalization','pdf'); 
        end
    end
    hold off;
    if i == 1
        xlabel('r','FontSize',16)
        legend('Non-PD','PD','FontSize',16,'location','northeast');
    elseif i == 2
        xlabel('\beta','FontSize',16)
    elseif i == 3
        xlabel('\delta_V','FontSize',16)
    elseif i == 4
        xlabel('\alpha','FontSize',16)
    end
    title('Accept-or-Reject','FontSize',16);
    ax = gca;
    ax.FontSize = 14;

    subplot(2,4,i+Nplot1)
    hold on; 
    for j = 1:Nsubgroups
        if i == 2
            histogram(VP_params_2sub_complex_acc_pert{j}(:,i)/1000,...
                'Normalization','pdf'); 
        elseif i == 4
            histogram(VP_params_2sub_complex_acc_pert{j}(:,i)*1000,...
                'Normalization','pdf'); 
        else
            histogram(VP_params_2sub_complex_acc_pert{j}(:,i),...
                'Normalization','pdf'); 
        end
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
    end
    title('Accept-or-Perturb','FontSize',16);
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,['VPops1_2sub_complex' num2str(dose_idx) '.fig']);
exportgraphics(gcf,['VPops1_2sub_complex' num2str(dose_idx) '.png']);

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.8]);
sgtitle(['Complex Model (2 of 2): Binary Subgroup Parameter Distributions '...
    'at Dose of ' num2str(V0*dose_scaling(dose_idx)) ...
    'x10^9 Virions'],'fontSize',16,'fontweight','bold');
for i = Nplot1+1:Nparams
    subplot(2,3,i-Nplot1)
    hold on; 
    for j = 1:Nsubgroups
        histogram(VP_params_2sub_complex_acc_rej{j}(:,i),...
            'Normalization','pdf'); 
    end
    hold off;
    if i == 5
        xlabel('\gamma_I','FontSize',16)
        legend('Non-PD','PD','FontSize',16,'location','northeast');
    elseif i == 6
        xlabel('\epsilon (\gamma_U = \epsilon\times \gamma_I)','FontSize',16)
    elseif i == 7
        xlabel('\lambda','FontSize',16)
    end
    title('Accept-or-Reject','FontSize',16);
    ax = gca;
    ax.FontSize = 14;

    subplot(2,3,i-Nplot1+Nplot2)
    hold on; 
    for j = 1:Nsubgroups
        histogram(VP_params_2sub_complex_acc_pert{j}(:,i),...
            'Normalization','pdf'); 
    end
    hold off;
    if i == 5
        xlabel('\gamma_I','FontSize',16)
        legend('Non-PD','PD','FontSize',16,'location','northeast');
    elseif i == 6
        xlabel('\epsilon (\gamma_U = \epsilon\times \gamma_I)','FontSize',16)
    elseif i == 7
        xlabel('\lambda','FontSize',16)
    end
    title('Accept-or-Perturb','FontSize',16);
    ax = gca;
    ax.FontSize = 14;
end
saveas(gcf,['VPops2_2sub_complex' num2str(dose_idx) '.fig']);
exportgraphics(gcf,['VPops2_2sub_complex' num2str(dose_idx) '.png']);


