clear all; close all; clc;
options = optimoptions(@fmincon,'Display','notify'); %fmincon options
T = readtable("Control_Data.xlsx",'VariableNamingRule','preserve');
A = table2array(T); 
time = A(:,1); % time in first column
N = length(time); % number of time points we have data for
data = A(:,2:16); % data in second column

t0 = 0;            % start time
tf = time(end);    % end time determined by getting end time in data

% Fit linear model for each mouse
param_guess = 0.3;  % starting guess for r
lb = 0; % lower bound on each parameter is 0 
ub = 10; % upper bound unknown for each
A = []; b = [];  % Linear inequality constraints of form Ax <= b
Aeq = []; beq = []; % Linear equality constrains of form Aeq*x = beq
best_params = zeros(size(data,2),1);
fits = zeros(size(data,2),1);
t_best = cell(size(data,2),1);
tum_best = cell(size(data,2),1);
x0 = mean(data(1,:));

cmap = parula(size(data,2));
h = [];
figure;
subplot(1,2,1)
hold on;
for i = 1:size(data,2)
    %x0 = data(1,i);
    time_data = [time data(:,i)]; % save in a matrix
    fun = @(z)goodness_of_fit(z,time_data,t0,tf,x0); % function to minimize
    [best_params(i), fits(i)] = fmincon(fun,param_guess,A,b,Aeq,beq,lb,...
        ub,[],options);
    fprintf('Best fit for Mouse %d: r = %f with cost = %f\n',i,...
        best_params(i),fits(i));

    % To visualize fit
    [t_best{i}, tum_best{i}] = ode45(@(t,x) DE_model(t,x,best_params(i)),...
        [t0,tf],x0);
    %figure;
    plot(t_best{i}, tum_best{i},'Color',cmap(i,:),'LineWidth',2); 
    %hold on;
    h{i} = plot(time,data(:,i),'o','Color',cmap(i,:),'LineWidth', 2,...
        'DisplayName', sprintf('Mouse #%d',i));
    %legend('Model','Data','Location','NorthWest','FontSize',14)
    xlabel('Time (days)','FontSize',14)
    ylabel('Tumor Volume (mm^3)','FontSize',14)
end
hold off;
ax = gca;
ax.FontSize = 14;
legend([h{1} h{2} h{3} h{4} h{5} h{6} h{7} h{8} h{9} h{10} ...
    h{11} h{12} h{13} h{14} h{15}],'Location','NorthWest','FontSize',14,...
    'NumColumns',3);

subplot(1,2,2)
hold on;
for i = 1:size(data,2)
    plot(i,best_params(i),'o','MarkerFaceColor',cmap(i,:),'MarkerSize',9,...
        'MarkerEdgeColor','none'); 
end
hold off;
xlabel('Mouse #','FontSize',14)
ylabel('Best fit \it{r}','FontSize',14,'Interpreter','tex')
ax = gca;
ax.FontSize = 14;
% exportgraphics(gcf,'r_histogram.png');
% saveas(gcf,'r_histogram.fig');

% 
% for i = 1:size(data,2)
%     figure;
%     plot(t_best{i}, tum_best{i},'Color',cmap(i,:),'LineWidth',2); 
%     hold on;
%     plot(time,data(:,i),'o','Color',cmap(i,:),'DisplayName',...
%         sprintf('Mouse #%d',i));
%     hold off;
%     xlabel('Time (days)','FontSize',14)
%     ylabel('Tumor Volume (mm^3)','FontSize',14)
%     title(['Mouse #' num2str(i)],'FontSize',16)
% end

save best_fit_control_avg_x0.mat best_params fits

function SSE = goodness_of_fit(params,data,t0,tf,x0)
    t_data = data(:,1);
    x_data = data(:,2);
    sol = ode45(@(t,x) DE_model(t,x,params),[t0,tf],x0);
    x_model = deval(sol,t_data)';
    SSE = sum((x_model-x_data).^2); 
end

function dxdt = DE_model(t,x,params)
    r = params(1);
    dxdt = r*x;
end
