%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% This code reads in the VCT output from all models and computes the    %
% Shannon entropy associated with the outcomes of the VCT.              %
% a uniform prior and accept-or-perturb using a normal prior.           %
% - Figure 1: Shannon entropy across OV doses using standard pairing of %
%   prior and inclusion/exclusion criteria. That is, normal with        %
%   accept-or-reject and uniform with accept-or-perturb.                %
% - Figure 2: Same result for "flipped" priors.                         %
%                                                                       %
% Updated: 6/14/2025                                                    %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
V0 = 10;
NVPops = 500;

%% Simple Model
load Simple_Model/VCT.mat
count_CR_acc_rej_simple = count_CR_acc_rej;
count_PR_acc_rej_simple = count_PR_acc_rej;
count_SD_acc_rej_simple = count_SD_acc_rej;
count_PD_acc_rej_simple = count_PD_acc_rej;
count_CR_acc_pert_simple = count_CR_acc_pert;
count_PR_acc_pert_simple = count_PR_acc_pert;
count_SD_acc_pert_simple = count_SD_acc_pert;
count_PD_acc_pert_simple = count_PD_acc_pert;

Ncat = 4; 
Ndoses = size(count_CR_acc_rej_simple,2);

% Accept-or-reject
p_acc_rej_simple = zeros(Ncat, Ndoses);
p_acc_rej_simple(1,:) = count_CR_acc_rej_simple/NVPops;
p_acc_rej_simple(2,:) = count_PR_acc_rej_simple/NVPops;
p_acc_rej_simple(3,:) = count_SD_acc_rej_simple/NVPops;
p_acc_rej_simple(4,:) = count_PD_acc_rej_simple/NVPops;

H_acc_rej_simple = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_rej_simple(i,j);
        if p ~= 0
            H_acc_rej_simple(1,j) = H_acc_rej_simple(1,j) - p*log10(p);
        end
    end
end
H_acc_rej_simple = H_acc_rej_simple/log10(Ncat); 

% Accept-or-perturb
p_acc_pert_simple = zeros(Ncat, Ndoses);
p_acc_pert_simple(1,:) = count_CR_acc_pert_simple/NVPops;
p_acc_pert_simple(2,:) = count_PR_acc_pert_simple/NVPops;
p_acc_pert_simple(3,:) = count_SD_acc_pert_simple/NVPops;
p_acc_pert_simple(4,:) = count_PD_acc_pert_simple/NVPops;

H_acc_pert_simple = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_pert_simple(i,j);
        if p ~= 0
            H_acc_pert_simple(1,j) = H_acc_pert_simple(1,j) - p*log10(p);
        end
    end
end
H_acc_pert_simple = H_acc_pert_simple/log10(Ncat); 

%% Intermediate Model
load Intermediate_Model/VCT.mat
count_CR_acc_rej_inter = count_CR_acc_rej;
count_PR_acc_rej_inter = count_PR_acc_rej;
count_SD_acc_rej_inter = count_SD_acc_rej;
count_PD_acc_rej_inter = count_PD_acc_rej;
count_CR_acc_pert_inter = count_CR_acc_pert;
count_PR_acc_pert_inter = count_PR_acc_pert;
count_SD_acc_pert_inter = count_SD_acc_pert;
count_PD_acc_pert_inter = count_PD_acc_pert;

% Accept-or-reject
p_acc_rej_inter = zeros(Ncat, Ndoses);
p_acc_rej_inter(1,:) = count_CR_acc_rej_inter/NVPops;
p_acc_rej_inter(2,:) = count_PR_acc_rej_inter/NVPops;
p_acc_rej_inter(3,:) = count_SD_acc_rej_inter/NVPops;
p_acc_rej_inter(4,:) = count_PD_acc_rej_inter/NVPops;

H_acc_rej_inter = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_rej_inter(i,j);
        if p ~= 0
            H_acc_rej_inter(1,j) = H_acc_rej_inter(1,j) - p*log10(p);
        end
    end
end
H_acc_rej_inter = H_acc_rej_inter/log10(Ncat); 

% Accept-or-perturb
p_acc_pert_inter = zeros(Ncat, Ndoses);
p_acc_pert_inter(1,:) = count_CR_acc_pert_inter/NVPops;
p_acc_pert_inter(2,:) = count_PR_acc_pert_inter/NVPops;
p_acc_pert_inter(3,:) = count_SD_acc_pert_inter/NVPops;
p_acc_pert_inter(4,:) = count_PD_acc_pert_inter/NVPops;

H_acc_pert_inter = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_pert_inter(i,j);
        if p ~= 0
            H_acc_pert_inter(1,j) = H_acc_pert_inter(1,j) - p*log10(p);
        end
    end
end
H_acc_pert_inter = H_acc_pert_inter/log10(Ncat); 

%% Complex Model
load Complex_Model/VCT.mat
count_CR_acc_rej_complex = count_CR_acc_rej;
count_PR_acc_rej_complex = count_PR_acc_rej;
count_SD_acc_rej_complex = count_SD_acc_rej;
count_PD_acc_rej_complex = count_PD_acc_rej;
count_CR_acc_pert_complex = count_CR_acc_pert;
count_PR_acc_pert_complex = count_PR_acc_pert;
count_SD_acc_pert_complex = count_SD_acc_pert;
count_PD_acc_pert_complex = count_PD_acc_pert;

% Accept-or-reject
p_acc_rej_complex = zeros(Ncat, Ndoses);
p_acc_rej_complex(1,:) = count_CR_acc_rej_complex/NVPops;
p_acc_rej_complex(2,:) = count_PR_acc_rej_complex/NVPops;
p_acc_rej_complex(3,:) = count_SD_acc_rej_complex/NVPops;
p_acc_rej_complex(4,:) = count_PD_acc_rej_complex/NVPops;

H_acc_rej_complex = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_rej_complex(i,j);
        if p ~= 0
            H_acc_rej_complex(1,j) = H_acc_rej_complex(1,j) - p*log10(p);
        end
    end
end
H_acc_rej_complex = H_acc_rej_complex/log10(Ncat); 

% Accept-or-perturb
p_acc_pert_complex = zeros(Ncat, Ndoses);
p_acc_pert_complex(1,:) = count_CR_acc_pert_complex/NVPops;
p_acc_pert_complex(2,:) = count_PR_acc_pert_complex/NVPops;
p_acc_pert_complex(3,:) = count_SD_acc_pert_complex/NVPops;
p_acc_pert_complex(4,:) = count_PD_acc_pert_complex/NVPops;

H_acc_pert_complex = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_pert_complex(i,j);
        if p ~= 0
            H_acc_pert_complex(1,j) = H_acc_pert_complex(1,j) - p*log10(p);
        end
    end
end
H_acc_pert_complex = H_acc_pert_complex/log10(Ncat); 

figure; 
hold on;
h1 = plot(V0*dose_scaling,H_acc_rej_simple,'o-','LineWidth',2,...
    'DisplayName','Simple Model: A-R');  
h2 = plot(V0*dose_scaling,H_acc_pert_simple,'o-','LineWidth',2,...
    'DisplayName','Simple Model: A-P');  
h3 = plot(V0*dose_scaling,H_acc_rej_inter,'s-','LineWidth',2,...
    'DisplayName','Intermediate Model: A-R');  
h4 = plot(V0*dose_scaling,H_acc_pert_inter,'s-','LineWidth',2,...
    'DisplayName','Intermediate Model: A-P');  
h5 = plot(V0*dose_scaling,H_acc_rej_complex,'x-','LineWidth',2,...
    'DisplayName','Complex Model: A-R');  
h6 = plot(V0*dose_scaling,H_acc_pert_complex,'x-','LineWidth',2,...
    'DisplayName','Complex Model: A-P');  
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Normalized Shannon Entropy (H)','FontSize',16);
ylim([0,1]);
ax = gca;
ax.FontSize = 14;
legend([h1 h2 h3 h4 h5 h6],'Location','NorthEast','FontSize',16); 

%%%%%
% FLIPPED Priors
%%%%%

%% Simple Model
load Simple_Model/VCT_flip.mat
count_CR_acc_rej_simple_flip = count_CR_acc_rej;
count_PR_acc_rej_simple_flip = count_PR_acc_rej;
count_SD_acc_rej_simple_flip = count_SD_acc_rej;
count_PD_acc_rej_simple_flip = count_PD_acc_rej;
count_CR_acc_pert_simple_flip = count_CR_acc_pert;
count_PR_acc_pert_simple_flip = count_PR_acc_pert;
count_SD_acc_pert_simple_flip = count_SD_acc_pert;
count_PD_acc_pert_simple_flip = count_PD_acc_pert;

% Accept-or-reject
p_acc_rej_simple_flip = zeros(Ncat, Ndoses);
p_acc_rej_simple_flip(1,:) = count_CR_acc_rej_simple_flip/NVPops;
p_acc_rej_simple_flip(2,:) = count_PR_acc_rej_simple_flip/NVPops;
p_acc_rej_simple_flip(3,:) = count_SD_acc_rej_simple_flip/NVPops;
p_acc_rej_simple_flip(4,:) = count_PD_acc_rej_simple_flip/NVPops;

H_acc_rej_simple_flip = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_rej_simple_flip(i,j);
        if p ~= 0
            H_acc_rej_simple_flip(1,j) = H_acc_rej_simple_flip(1,j) - p*log10(p);
        end
    end
end
H_acc_rej_simple_flip = H_acc_rej_simple_flip/log10(Ncat); 

% Accept-or-perturb
p_acc_pert_simple_flip = zeros(Ncat, Ndoses);
p_acc_pert_simple_flip(1,:) = count_CR_acc_pert_simple_flip/NVPops;
p_acc_pert_simple_flip(2,:) = count_PR_acc_pert_simple_flip/NVPops;
p_acc_pert_simple_flip(3,:) = count_SD_acc_pert_simple_flip/NVPops;
p_acc_pert_simple_flip(4,:) = count_PD_acc_pert_simple_flip/NVPops;

H_acc_pert_simple_flip = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_pert_simple_flip(i,j);
        if p ~= 0
            H_acc_pert_simple_flip(1,j) = H_acc_pert_simple_flip(1,j) - p*log10(p);
        end
    end
end
H_acc_pert_simple_flip = H_acc_pert_simple_flip/log10(Ncat); 

%% Intermediate Model
load Intermediate_Model/VCT_flip.mat
count_CR_acc_rej_inter_flip = count_CR_acc_rej;
count_PR_acc_rej_inter_flip = count_PR_acc_rej;
count_SD_acc_rej_inter_flip = count_SD_acc_rej;
count_PD_acc_rej_inter_flip = count_PD_acc_rej;
count_CR_acc_pert_inter_flip = count_CR_acc_pert;
count_PR_acc_pert_inter_flip = count_PR_acc_pert;
count_SD_acc_pert_inter_flip = count_SD_acc_pert;
count_PD_acc_pert_inter_flip = count_PD_acc_pert;

% Accept-or-reject
p_acc_rej_inter_flip = zeros(Ncat, Ndoses);
p_acc_rej_inter_flip(1,:) = count_CR_acc_rej_inter_flip/NVPops;
p_acc_rej_inter_flip(2,:) = count_PR_acc_rej_inter_flip/NVPops;
p_acc_rej_inter_flip(3,:) = count_SD_acc_rej_inter_flip/NVPops;
p_acc_rej_inter_flip(4,:) = count_PD_acc_rej_inter_flip/NVPops;

H_acc_rej_inter_flip = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_rej_inter_flip(i,j);
        if p ~= 0
            H_acc_rej_inter_flip(1,j) = H_acc_rej_inter_flip(1,j) - p*log10(p);
        end
    end
end
H_acc_rej_inter_flip = H_acc_rej_inter_flip/log10(Ncat); 

% Accept-or-perturb
p_acc_pert_inter_flip = zeros(Ncat, Ndoses);
p_acc_pert_inter_flip(1,:) = count_CR_acc_pert_inter_flip/NVPops;
p_acc_pert_inter_flip(2,:) = count_PR_acc_pert_inter_flip/NVPops;
p_acc_pert_inter_flip(3,:) = count_SD_acc_pert_inter_flip/NVPops;
p_acc_pert_inter_flip(4,:) = count_PD_acc_pert_inter_flip/NVPops;

H_acc_pert_inter_flip = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_pert_inter_flip(i,j);
        if p ~= 0
            H_acc_pert_inter_flip(1,j) = H_acc_pert_inter_flip(1,j) - p*log10(p);
        end
    end
end
H_acc_pert_inter_flip = H_acc_pert_inter_flip/log10(Ncat); 

%% Complex Model
load Complex_Model/VCT_flip.mat
count_CR_acc_rej_complex_flip = count_CR_acc_rej;
count_PR_acc_rej_complex_flip = count_PR_acc_rej;
count_SD_acc_rej_complex_flip = count_SD_acc_rej;
count_PD_acc_rej_complex_flip = count_PD_acc_rej;
count_CR_acc_pert_complex_flip = count_CR_acc_pert;
count_PR_acc_pert_complex_flip = count_PR_acc_pert;
count_SD_acc_pert_complex_flip = count_SD_acc_pert;
count_PD_acc_pert_complex_flip = count_PD_acc_pert;

% Accept-or-reject
p_acc_rej_complex_flip = zeros(Ncat, Ndoses);
p_acc_rej_complex_flip(1,:) = count_CR_acc_rej_complex_flip/NVPops;
p_acc_rej_complex_flip(2,:) = count_PR_acc_rej_complex_flip/NVPops;
p_acc_rej_complex_flip(3,:) = count_SD_acc_rej_complex_flip/NVPops;
p_acc_rej_complex_flip(4,:) = count_PD_acc_rej_complex_flip/NVPops;

H_acc_rej_complex_flip = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_rej_complex_flip(i,j);
        if p ~= 0
            H_acc_rej_complex_flip(1,j) = H_acc_rej_complex_flip(1,j) - p*log10(p);
        end
    end
end
H_acc_rej_complex_flip = H_acc_rej_complex_flip/log10(Ncat); 

% Accept-or-perturb
p_acc_pert_complex_flip = zeros(Ncat, Ndoses);
p_acc_pert_complex_flip(1,:) = count_CR_acc_pert_complex_flip/NVPops;
p_acc_pert_complex_flip(2,:) = count_PR_acc_pert_complex_flip/NVPops;
p_acc_pert_complex_flip(3,:) = count_SD_acc_pert_complex_flip/NVPops;
p_acc_pert_complex_flip(4,:) = count_PD_acc_pert_complex_flip/NVPops;

H_acc_pert_complex_flip = zeros(1,Ndoses);
for j = 1:Ndoses
    for i = 1:Ncat
        p = p_acc_pert_complex_flip(i,j);
        if p ~= 0
            H_acc_pert_complex_flip(1,j) = H_acc_pert_complex_flip(1,j) - p*log10(p);
        end
    end
end
H_acc_pert_complex_flip = H_acc_pert_complex_flip/log10(Ncat); 

figure; 
hold on;
h1 = plot(V0*dose_scaling,H_acc_rej_simple_flip,'o-','LineWidth',2,...
    'DisplayName','Simple Model: A-R');  
h2 = plot(V0*dose_scaling,H_acc_pert_simple_flip,'o-','LineWidth',2,...
    'DisplayName','Simple Model: A-P');  
h3 = plot(V0*dose_scaling,H_acc_rej_inter_flip,'s-','LineWidth',2,...
    'DisplayName','Intermediate Model: A-R');  
h4 = plot(V0*dose_scaling,H_acc_pert_inter_flip,'s-','LineWidth',2,...
    'DisplayName','Intermediate Model: A-P');  
h5 = plot(V0*dose_scaling,H_acc_rej_complex_flip,'x-','LineWidth',2,...
    'DisplayName','Complex Model: A-R');  
h6 = plot(V0*dose_scaling,H_acc_pert_complex_flip,'x-','LineWidth',2,...
    'DisplayName','Complex Model: A-P');  
hold off;
xlabel('Virions (\times 10^9)','FontSize',16)
ylabel('Normalized Shannon Entropy (H)','FontSize',16);
ylim([0,1]);
ax = gca;
ax.FontSize = 14;
legend([h1 h2 h3 h4 h5 h6],'Location','NorthEast','FontSize',16); 


