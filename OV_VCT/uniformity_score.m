%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Reads in VPops and computes the uniformity score (normalized Earth    %
% Movers Distance) for all model types (simple, intermediate, complex), %
% all priors (normal, uniform) and all inclusion/exclusion criteria     %
% (accept-or-reject and accept-or-perturb).                             %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test metric on simple data
% x1 = [0.1 0.25 0.5 0.75 0.9]; 
% D1 = uniformity_score(x1,0,1)
% 
% x2 = 2*[0.1 0.25 0.5 0.75 0.9]; 
% D2 = uniformity_score(x2,0,2)
% 
% x3 = 0:0.0025:1;
% D3 = uniformity_score(x3,0,1)
% 
% x4 = [0.01 0.001 0.002 0.03 0.004];
% D4 = uniformity_score(x4,0,1)

clear all; close all; clc;
load Simple_Model/Vpops_acc_rej.mat
lb = curve2_Nstd(end);
ub = curve1_Nstd(end);
simple_acc_rej = Traj(:,end);
D_simple_acc_rej = emd(simple_acc_rej,lb,ub)

load Simple_Model/Vpops_acc_pert.mat
simple_acc_pert = Traj(:,end);
D_simple_acc_pert = emd(simple_acc_pert,lb,ub)

load Intermediate_Model/Vpops_acc_rej.mat
inter_acc_rej = Traj(:,end);
D_inter_acc_rej = emd(inter_acc_rej,lb,ub)

load Intermediate_Model/Vpops_acc_pert.mat
inter_acc_pert = Traj(:,end);
D_inter_acc_pert = emd(inter_acc_pert,lb,ub)


load Complex_Model/Vpops_acc_rej.mat
complex_acc_rej = Traj(:,end);
D_complex_acc_rej = emd(complex_acc_rej,lb,ub)

load Complex_Model/Vpops_acc_pert.mat
complex_acc_pert = Traj(:,end);
D_complex_acc_pert = emd(complex_acc_pert,lb,ub)

D = [D_simple_acc_rej D_simple_acc_pert; ...
    D_inter_acc_rej D_inter_acc_pert; ...
    D_complex_acc_rej D_complex_acc_pert];
figure; 
b = bar(D);
str = {'Simple','Intermediate','Complex'};
set(gca,'XTickLabel',str,'XTick',1:numel(str),'FontSize',14)
xlabel('Model','FontSize',16)
ylabel('Uniformity Score','FontSize',16)
legend('Accept-or-Reject','Accept-or-Perturb','Location','NorthWest',...
    'FontSize',16)

%%%%%
% FLIPPED priors
%%%%%
load Simple_Model/Vpops_acc_rej_flip.mat
simple_acc_rej_flip = Traj(:,end);
D_simple_acc_rej_flip = emd(simple_acc_rej_flip,lb,ub)

load Simple_Model/Vpops_acc_pert_flip.mat
simple_acc_pert_flip = Traj(:,end);
D_simple_acc_pert_flip = emd(simple_acc_pert_flip,lb,ub)

load Intermediate_Model/Vpops_acc_rej_flip.mat
inter_acc_rej_flip = Traj(:,end);
D_inter_acc_rej_flip = emd(inter_acc_rej_flip,lb,ub)

load Intermediate_Model/Vpops_acc_pert_flip.mat
inter_acc_pert_flip = Traj(:,end);
D_inter_acc_pert_flip = emd(inter_acc_pert_flip,lb,ub)


load Complex_Model/Vpops_acc_rej_flip.mat
complex_acc_rej_flip = Traj(:,end);
D_complex_acc_rej_flip = emd(complex_acc_rej_flip,lb,ub)

load Complex_Model/Vpops_acc_pert_flip.mat
complex_acc_pert_flip = Traj(:,end);
D_complex_acc_pert_flip = emd(complex_acc_pert_flip,lb,ub)

D_flip = [D_simple_acc_rej_flip D_simple_acc_pert_flip; ...
    D_inter_acc_rej_flip D_inter_acc_pert_flip; ...
    D_complex_acc_rej_flip D_complex_acc_pert_flip];

figure; 
subplot(1,2,1)
bar(D);
str = {'Simple','Intermediate','Complex'};
set(gca,'XTickLabel',str,'XTick',1:numel(str),'FontSize',14)
xlabel('Model','FontSize',16)
ylabel('Uniformity Score','FontSize',16)
legend('Accept-or-Reject','Accept-or-Perturb','Location','NorthWest',...
    'FontSize',16)
subtitle('Standard Priors','FontSize',16)

subplot(1,2,2)
bar(D_flip);
str = {'Simple','Intermediate','Complex'};
set(gca,'XTickLabel',str,'XTick',1:numel(str),'FontSize',14)
xlabel('Model','FontSize',16)
ylabel('Uniformity Score','FontSize',16)
legend('Accept-or-Reject','Accept-or-Perturb','Location','NorthWest',...
    'FontSize',16)
subtitle('Flipped Priors','FontSize',16)

D_all = [D_simple_acc_rej D_simple_acc_pert D_simple_acc_rej_flip D_simple_acc_pert_flip; ...
    D_inter_acc_rej D_inter_acc_pert D_inter_acc_rej_flip D_inter_acc_pert_flip; ...
    D_complex_acc_rej D_complex_acc_pert D_complex_acc_rej_flip D_complex_acc_pert_flip];
figure;
b = bar(D_all);
str = {'Simple','Intermediate','Complex'};
set(gca,'XTickLabel',str,'XTick',1:numel(str),'FontSize',14)

% Will include numerical height of graph for dataset 1
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = compose('%.2f', b(1).YData); 
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',14); 
% Will include numerical height of graph for dataset 1
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = compose('%.2f', b(2).YData); 
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',14);
% Will include numerical height of graph for dataset 3
xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = compose('%.2f', b(3).YData); 
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',14);
% Will include numerical height of graph for dataset 4
xtips4 = b(4).XEndPoints;
ytips4 = b(4).YEndPoints;
labels4 = compose('%.2f', b(4).YData);  
text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',14);

xlabel('Model','FontSize',16)
ylabel('Uniformity Score,\it{D}','FontSize',16,'Interpreter','tex')
legend('Accept-or-Reject (Normal)','Accept-or-Perturb (Uniform)',...
    'Accept-or-Reject (Uniform)','Accept-or-Perturb (Normal)',...
    'Location','NorthWest','FontSize',16)


function D = emd(x,a,b)
% Inputs:
%   x : vector of point locations
%   a : interval start
%   b : interval end
%
% Output:
%   D : uniformity score in [0,1]
%       D = 1   --> perfectly evenly spaced
%       D ~ 0   --> points highly clustered

    % Rescale to [0,1]
    y = sort((x - a) / (b - a));
    n = numel(y);

    % Empirical CDF jumps at i/n for y(i)
    % Integral of |F_n(t)-t| over [0,1]:
    % closed-form: sum over intervals between sorted y's
    W1 = 0;
    y_ext = [0; y(:); 1];  % add boundaries
    for i = 1:n+1
        % On [y_ext(i), y_ext(i+1)], F_n(t) = (i-1)/n
        left = y_ext(i);
        right = y_ext(i+1);
        Fn_val = (i-1)/n;

        % Integrate |Fn_val - t| dt over [left,right]
        % Split case where t=Fn_val crosses interval
        c = Fn_val;
        if c <= left || c >= right
            % sign is constant
            W1 = W1 + abs(Fn_val - (left+right)/2) * (right-left);
        else
            % split into two sub-intervals
            W1 = W1 + ...
                (integrate_abs(Fn_val,left,c)) + ...
                (integrate_abs(Fn_val,c,right));
        end
    end

    % Uniformity score
    D = 1 - 2*W1;
end

function val = integrate_abs(c,left,right)
    % integrate |c - t| dt from left to right, with c in [left,right]
    val = (c-left)^2/2 + (right-c)^2/2;
end