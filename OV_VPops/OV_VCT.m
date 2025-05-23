%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                                                       %
% Reads in parametrization of each Vpop generated in OV_Vpops and       %
% conducts a virtual clinical trial that looks at the impact of varying %
% dose of drug.                                                         %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc; 
load data_AD.mat
tf = days(end); % stopping time
I0 = 0; % no infected cells initially
V0  = 10; % initial number of viruses, scaled to volume
p.d_I = 1; 
p.alpha = 3; 
p.d_V = 2.3; 

load Vpops.mat
dose_scaling = 1:1:90; 
count_shrink = zeros(size(dose_scaling));
for j = 1:length(dose_scaling)
    %figure; hold on;
    for i=1:size(Vpop_AccRej,1)
        p.r = Vpop_AccRej(i,1); 
        p.beta = Vpop_AccRej(i,2); 
        p.U0 = Vpop_AccRej(i,3);
        in_cond = [p.U0 I0 V0]; 
        
        [sol1, sol2, sol3] = OV_model_VCT(dose_scaling(j),p,in_cond,tf);
        % Combine data
        sol.x = [sol1.x  sol2.x  sol3.x];
        sol.y = [sol1.y  sol2.y  sol3.y];
        
        %plot(sol.x, sol.y(1,:)+sol.y(2,:),'LineWidth',2);
        Vf = sol.y(1,end)+sol.y(2,end);
        if Vf<p.U0
            count_shrink(j) = count_shrink(j) + 1;
        end 
    end
    %hold off;
    fprintf('Dose scaling = %d: %d of %d Vpops have tumor shrink\n',...
        dose_scaling(j),count_shrink(j),size(Vpop_AccRej,1));
end

figure; 
plot(V0*dose_scaling,count_shrink,'o-','LineWidth',2); 
xlabel('OV dose','FontSize',16)
ylabel('Percent of Vpops in which tumor shrinks','FontSize',16);


%%% Functions %%%
function [sol1, sol2, sol3] = OV_model_VCT(dose_scaling,p,in_cond,tf)
    V0 = dose_scaling*in_cond(3);
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
% 
% function [sol1, sol2, sol3, sol4, sol5, sol6] = OV_model_VCT(p,in_cond,tf)
%     %% VCT asks if virtual patient will be effectively treated with one additional dose
%     V0 = 3*in_cond(3);
%     % First virus dose
%     sol1 = ode23s(@odemodel,[0 2], in_cond);
%     % Second virus dose
%     in_cond = sol1.y(:,end)' + [0    0    V0];
%     sol2 = ode23s(@odemodel,[2 4], in_cond);
%     % Third virus dose
%     in_cond = sol2.y(:,end)' + [0    0    V0];
%     sol3 = ode23s(@odemodel,[4 6], in_cond);
%     % Fourth virus dose
%     in_cond = sol3.y(:,end)' + [0    0    V0];
%     sol4 = ode23s(@odemodel,[6 8], in_cond);
%     % Fifth virus dose
%     in_cond = sol4.y(:,end)' + [0    0    V0];
%     sol5 = ode23s(@odemodel,[8 10], in_cond);
%     % Sixth virus dose
%     in_cond = sol5.y(:,end)' + [0    0    V0];
%     sol6 = ode23s(@odemodel,[10 tf], in_cond);
%     
%     %% OV model
%     function dydt = odemodel(t,y)
%         U  = y(1); 
%         I  = y(2);
%         V  = y(3);
%         % Frequency-dependent killing
%         N = U + I;
%         dU  = p.r*U - p.beta*U*V/N;
%         dI  = p.beta*U*V/N - p.d_I*I;
%         dV  = p.alpha*p.d_I*I - p.d_V*V; 
%         dydt = [dU; dI; dV]; 
%     end
% end

