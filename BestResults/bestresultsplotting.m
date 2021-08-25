%% bestresultsplotting.m
% Levi Manring, Duke University
% 2021
%
% This m-file plots the results contained in the "BestResults" folder as
% well as the desired reward function and the friction coefficient

%%

% Initialization
clear all
close all
clc

% Add the relevant paths to the directory
addpath(genpath('/VehicleDynamicsModelRWD'));

% Load the model for comparison
load('CurveFitModel');

%% RWD NOSLIP results
% Either load PILCO or DDPG results
% load('PILCO_RWD_NOSLIP');
load('DDPG_RWD_NOSLIP');
tbetween = params.dt/2;

% get complete initial conditions
[~, yAx, ~, ~, ~] = y_eval_fcn(states(1,1),params);
phiA = sqrt(yAx^2+1)/params.dim.R;
theta_A = 0;
thetadot_A = phiA*states(1,2);
init_conds = [states(1,1), states(1,2), theta_A, thetadot_A];

% Simulate the RWD SLIP model for comparison
slip_cond_start = curves_startslipcheck_fun(init_conds,action(1,1),params,curves);
stored_states = zeros(length(times),5);
stored_cond = zeros(length(times),1);
stored_cond(1,1) = slip_cond_start;
stored_states(1,1:4) = init_conds;
refine = 20;
% Step through integration
for k = 1:(length(times) - 1)
    % Get time and action
    t0 = times(k,1);
    current_Action = action(k,1);
    
    % Integrate over the timestep
    [slipstates, slip_cond_end, event_num, failed] = curves_stickslip_nr(init_conds,slip_cond_start,current_Action,[t0, t0 + params.dt],params,curves,refine);
    
    % Update slip condition and initial conditions
    slip_cond_start = slip_cond_end;
    
    stored_cond(k,1) = slip_cond_start;
    init_conds = slipstates(1,1:4);
    stored_states(k+1,:) = slipstates(1,:);
    
    % Report progress
    fprintf('RWD SLIP Comparison: k # %3.0f, Failed? %3.0f \n', k, failed);
end

% Find when the RWD SLIP model slips
[rows, ~] = find(stored_cond > 0);

% find bounds on the slipping areas, and locate the times
bounds = bounds_fun(rows);
tslip = times(bounds);
if size(tslip,2) < 2
    tslip = tslip';
end

% Create a comparison plot
basevalue1 = -10;
basevalue2 = -800;
figure;
h1 = subplot(2,1,1); hold on
h2 = subplot(2,1,2); hold on
% Create shaded areas denoting slippage
for k = 1:length(tslip)
    if tslip(k,1) == tslip(k,2)
        tslip(k,2) = tslip(k,2) + tbetween;
    end
    a1 = area(h1,[tslip(k,1), tslip(k,2)], [-basevalue1; -basevalue1],basevalue1,'FaceColor',[220/256 220/256, 220/256],'LineStyle','none');
    a2 = area(h2,[tslip(k,1), tslip(k,2)], [-basevalue2; -basevalue2],basevalue2,'FaceColor',[220/256 220/256, 220/256],'LineStyle','none');
end
p1 = plot(h1,times,stored_states(:,1),'r-','linewidth',2);
p2 = plot(h1,times,states(:,1),'b-','linewidth',2);
p3 = plot(h1,[times(1) times(end)],[params.goal params.goal],'k--','linewidth',3);
p4 = plot(h2,times,action,'k-','linewidth',2);
set(h1,'TickLength',[0, 0],'TickLabelInterpreter','latex','fontsize',16);
set(h2,'TickLength',[0, 0],'TickLabelInterpreter','latex','fontsize',16);
xlabel(h1,'$t$ (s)','Interpreter','latex','fontsize',24);
ylabel(h1,'$x$ (m)','Interpreter','latex','fontsize',24);
xlabel(h2,'$t$ (s)','Interpreter','latex','fontsize',24);
ylabel(h2,'$\tau_A$ (N-m)','Interpreter','latex','fontsize',24);
legend(h1,[p2 p1 a1 p3],'Case 1 Trajectory','Case 1/Case 3 Trajectory','Wheel $A$ Slip Region','Target','location','northwest','interpreter','latex','fontsize',18);
legend(h2,[p4 a2],'Control','Wheel $A$ Slip Region','location','northeast','interpreter','latex','fontsize',18);
ylim(h2,[basevalue2, -basevalue2]);
xlim(h1,[0, 22]);
xlim(h2,[0, 22]);
set(gcf,'color','none');

%% RWD SLIP RESULTS
clear all;
% Choose one of two results files to load
load('DDPG_RWD_SLIP_1');
% load('DDPG_RWD_SLIP_2');

stored_cond = states(:,4);
tbetween = params.dt/2;
[rows, ~] = find(stored_cond > 0);

% find bounds on the slipping areas
bounds = bounds_fun(rows);
tslip = times(bounds);
if size(tslip,2) < 2
    tslip = tslip';
end

% Create a comparison plot
basevalue1 = -10;
basevalue2 = -800;
figure;
h1 = subplot(2,1,1); hold on
h2 = subplot(2,1,2); hold on
% Create shaded areas denoting slippage
for k = 1:length(tslip)
    if tslip(k,1) == tslip(k,2)
        tslip(k,2) = tslip(k,2) + tbetween;
    end
    a1 = area(h1,[tslip(k,1), tslip(k,2)], [-basevalue1; -basevalue1],basevalue1,'FaceColor',[220/256 220/256, 220/256],'LineStyle','none');
    a2 = area(h2,[tslip(k,1), tslip(k,2)], [-basevalue2; -basevalue2],basevalue2,'FaceColor',[220/256 220/256, 220/256],'LineStyle','none');
end
p2 = plot(h1,times,states(:,1),'b-','linewidth',2);
p3 = plot(h1,[times(1) times(end)],[params.goal params.goal],'k--','linewidth',3);
p4 = plot(h2,times,action,'k-','linewidth',2);
set(h1,'TickLength',[0, 0],'TickLabelInterpreter','latex','fontsize',16);
set(h2,'TickLength',[0, 0],'TickLabelInterpreter','latex','fontsize',16);
xlabel(h1,'$t$ (s)','Interpreter','latex','fontsize',24);
ylabel(h1,'$x$ (m)','Interpreter','latex','fontsize',24);
xlabel(h2,'$t$ (s)','Interpreter','latex','fontsize',24);
ylabel(h2,'$\tau_A$ (N-m)','Interpreter','latex','fontsize',24);
legend(h1,[p2 a1 p3],'Case 1/Case 3 Trajectory','Wheel $A$ Slip Region','Target','location','east','interpreter','latex','fontsize',18);
legend(h2,[p4 a2],'Control','Wheel $A$ Slip Region','location','east','interpreter','latex','fontsize',18);
ylim(h2,[basevalue2, -basevalue2]);
xlim(h1,[0, 22]);
xlim(h2,[0, 22]);
set(gcf,'color','none');

%% AWD Results
% run through the AWD results
for k = 1:5
    % load the AWD results
    load(strcat('DDPG_AWD_',num2str(k)));
    tbetween = params.dt/2;
    stored_cond = states(:,5);
    % find the regions for each slipping condition
    [rows2, ~] = find(stored_cond == 2);
    [rows3, ~] = find(stored_cond == 3);
    [rows4, ~] = find(stored_cond == 4);
    
    basevalue1 = -10;
    basevalue2 = -800;
    case2color = [0.6 1 0.6];
    case3color = [1 1 0.6];
    case4color = [1 0.6 1];
    
    figure;
    h1 = subplot(2,1,1); hold on
    h2 = subplot(2,1,2); hold on
    if ~isempty(rows2)
        bounds2 = bounds_fun(rows2);
        tslip2 = times(bounds2);
        if size(tslip2,2) < 2
            tslip2 = tslip2';
        end
        for k = 1:size(tslip2,1)
            if tslip2(k,1) == tslip2(k,2)
                tslip2(k,2) = tslip2(k,2) + tbetween;
            end
            a1_2 = area(h1,[tslip2(k,1), tslip2(k,2)], [-basevalue1; -basevalue1],basevalue1,'FaceColor',case2color,'LineStyle','none');
            a2_2 = area(h2,[tslip2(k,1), tslip2(k,2)], [-basevalue2; -basevalue2],basevalue2,'FaceColor',case2color,'LineStyle','none');
        end
    end
    if ~isempty(rows3)
        bounds3 = bounds_fun(rows3);
        tslip3 = times(bounds3);
        if size(tslip3,2) < 2
            tslip3 = tslip3';
        end
        for k = 1:size(tslip3,1)
            if tslip3(k,1) == tslip3(k,2)
                tslip3(k,2) = tslip3(k,2) + tbetween;
            end
            a1_3 = area(h1,[tslip3(k,1), tslip3(k,2)], [-basevalue1; -basevalue1],basevalue1,'FaceColor',case3color,'LineStyle','none');
            a2_3 = area(h2,[tslip3(k,1), tslip3(k,2)], [-basevalue2; -basevalue2],basevalue2,'FaceColor',case3color,'LineStyle','none');
        end
    end
    
    if ~isempty(rows4)
        bounds4 = bounds_fun(rows4);
        tvrel4 = times(bounds4);
        if size(tvrel4,2) < 2
            tvrel4 = tvrel4';
        end
        for k = 1:size(tvrel4,1)
            if tvrel4(k,1) == tvrel4(k,2)
                tvrel4(k,2) = tvrel4(k,2) + tbetween;
            end
            a1_4 = area(h1,[tvrel4(k,1), tvrel4(k,2)], [-basevalue1; -basevalue1],basevalue1,'FaceColor',case4color,'LineStyle','none');
            a2_4 = area(h2,[tvrel4(k,1), tvrel4(k,2)], [-basevalue2; -basevalue2],basevalue2,'FaceColor',case4color,'LineStyle','none');
        end
    end
    p2 = plot(h1,times,states(:,1),'b-','linewidth',2);
    p3 = plot(h1,[times(1) times(end)],[params.goal params.goal],'k--','linewidth',3);
    p4 = plot(h2,times,action(:,1),'k-','linewidth',2);
    p5 = plot(h2,times,action(:,2),'ro--','linewidth',2);
    set(h1,'TickLength',[0, 0],'TickLabelInterpreter','latex','fontsize',16);
    set(h2,'TickLength',[0, 0],'TickLabelInterpreter','latex','fontsize',16);
    xlabel(h1,'$t$ (s)','Interpreter','latex','fontsize',24);
    ylabel(h1,'$x$ (m)','Interpreter','latex','fontsize',24);
    xlabel(h2,'$t$ (s)','Interpreter','latex','fontsize',24);
    ylabel(h2,'$\tau_K$ (N-m)','Interpreter','latex','fontsize',24);
    if isempty(rows4)
        legend(h1,[p2 a1_3 a1_2 p3],'Case 1-4 Trajectory','Wheel $A$ Slip Region',...
            'Wheel $A$ \& $B$ Slip Region','Target','location','east','interpreter','latex','fontsize',18);
        legend(h2,[p4 p5 a2_3 a2_2],'Control $\tau_A$','Control $\tau_B$','Wheel $A$ Slip Region',...
            'Wheel $A$ \& $B$ Slip Region','location','east','interpreter','latex','fontsize',18);
    else
        legend(h1,[p2 a1_3 a1_4 a1_2 p3],'Case 1-4 Trajectory','Wheel $A$ Slip Region',...
            'Wheel $B$ Slip Region','Wheel $A$ \& $B$ Slip Region','Target','location','east','interpreter','latex','fontsize',18);
        legend(h2,[p4 p5 a2_3 a2_4 a2_2],'Control $\tau_A$','Control $\tau_B$','Wheel $A$ Slip Region',...
            'Wheel $B$ Slip Region','Wheel $A$ \& $B$ Slip Region','location','east','interpreter','latex','fontsize',18);
    end
    ylim(h2,[basevalue2, -basevalue2]);
    xlim(h1,[0, 14]);
    xlim(h2,[0, 14]);
    set(gcf,'color','none');
end

%% Reward Function Plotting
% addpath to get to the desired reward function
addpath('/RL_DDPG_RWD');

params.goal = 9.5;
params.goaltolerance = 0.1;
params.statestolerance = 0.1;
e_xvec = (-params.goal:0.1:5)';
e_xdotvec = (-8:0.1:8)';
[e_xmesh, e_xdotmesh] = meshgrid(e_xvec, e_xdotvec);
Rewardk = zeros(size(e_xmesh,1),size(e_xmesh,2));
colorvalue = Rewardk;
for k = 1:size(e_xmesh,1)
    for p = 1:size(e_xmesh,2)
        x = e_xmesh(k,p) + params.goal;
        xdot = e_xdotmesh(k,p);
        Rewardk(k,p) = RWDrewardfun_augmented(x,xdot,0,0,params);
        colorvalue(k,p) = RWDrewardfun_augmented(x,xdot,0,0,params);
    end
end

figure;
surf(e_xmesh,e_xdotmesh,Rewardk,colorvalue); zlim([0, 25]);
xlim([e_xvec(1) e_xvec(end)]); ylim([e_xdotvec(1) e_xvec(end)]);
ylabel('$e_{\dot{x}}$','Interpreter','latex','fontsize',24);
xlabel('$e_x$','Interpreter','latex','fontsize',24);
zlabel('$r_s$','Interpreter','latex','fontsize',24);
set(gca,'TickLabelInterpreter','latex','fontsize',20);
set(gcf,'color','none');

%% Friction coefficient plotting
addpath('/VehicleDynamicsModelRWD/OtherFilesNeeded');
v_rel = (-3:0.0001:3)';
mu = zeros(size(v_rel));
for k = 1:length(v_rel)
    mu(k,1) = mu_fcn(v_rel(k,1),params.dim.mu_s,2);
end
figure;
p1 = subplot(1,1,1);
plot(v_rel,mu,'k-','linewidth',2);
set(p1,'TickLabelInterpreter','latex','fontsize',16);
ylabel(p1,'$\mu_K$','Interpreter','latex','fontsize',24);
xlabel(p1,'$v_{r,K}$','Interpreter','latex','fontsize',24);
set(gcf,'color','none');
