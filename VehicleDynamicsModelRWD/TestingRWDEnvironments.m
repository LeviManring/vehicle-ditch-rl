%% TestingRWDEvironments.m
% Levi Manring, Duke University
% 2021
%
% This m-file tests 5 different model/solver combinations for a vehicle 
% moving on a user-defined surface profile:
%   1. A curve-fit model incorporating no-slip/slip dynamics and event
%   detection using Matlab's ODE event detection capability
%   2. A curve-fit model incorporating no-slip/slip dynamics and event
%   detection using a user-coded Newton-Raphson routine
%   3. An exact solution model incorporating stick/slip dynamics and event
%   detection using Matlab's ODE event detection capability
%   4. A curve-fit model incorporating JUST no-slip dynamics
%   5. An exact solution model incorporating JUST no-slip dynamics
%
% Run this m-file to compare the accuracy/results of using a curve-fit
% model versus an exact solution model. Also, note the importance of
% incorporating the slip dynamics for different dynamic behavior.
%
% In judging these estimates: Scenarios 1, 2, and 3 should be close, with 1
% and 2 being nearly identical and 3 having some slight integration
% difference compared to 1 or 2. Scenarios 4 and 5 should be nearly identical
% as well. The reason to use the curvefit model is for speed in
% evaluating the dynamic model(s) as the exact method is computationally more
% demanding.
%
% There are two different control/torque profiles we have included, one is
% generated using PILCO, the other is generated using DDPG in Matlab. These
% two scenarios were trained on the no-slip dynamics model. We have saved
% the trained agents that work with the discontinuous model, but since
% agents work as feedback controllers, controlling the complete
% discontinuous model in an open-loop simulation like this does not work as
% well as. If you want to test the trained DDPG agent for the no-slip/slip 
% dynamics, see the README in the RL_Car_RWD folder.
%
% To see a video plot of the results, simply uncomment the video plotting
% line in each section.

%% Initialize
clear all; close all; clc;

% Generate the paths for the RWD vehicle dynamics model (modify this line
% for your file directory)
addpath(genpath('/Users/levimanring/Documents/Duke University/Research/RL Car/VehicleDynamicsModelRWD'));

%% Load or calculate VehicleParameters if these are not loaded previously

% % if calculating uncomment these lines
% syms x
% f = (exp(x^2*(-16/225))*3-3)*(x^2/1225); % user can define f as any continuous function of x
% params = f_params(f,x);
% GetCurves;

% % otherwise load these lines
load('VehicleParameters');
load('CurveFitModel');

%% Load the desired scenario
% load('PILCOFinalProfile'); % torque profile 1, PILCO trained no-slip
load('RL_NS_FinalProfile'); % torque profile 2, DDPG trained no-slip

%% Set parameters for each scenario
params.dt = tvec(2) - tvec(1);
refine = 20;
framerate = 10; % framerate for video recording
record = 0; % choosing to record or not
guess_Dx = params.dim.l;
params.dim.goal = 9.5;
[phi, ~] = phi_fun(x(1,1),'A',guess_Dx,params);
global_init_conds = [x(1,1), xdot(1,1), 0, xdot(1,1)*phi];

%% 1. Curve-fit model w/ no-slip/slip dynamics and Matlab's ODE event detection
% initialization
init_conds = global_init_conds;
slip_cond_start = curves_startslipcheck_fun(init_conds,torque(1,1),params,curves);
stored_states_1 = zeros(length(tvec),5);

% Step through integration
for k = 1:(length(tvec) - 1)
    % Get time and torque
    t0 = tvec(k,1);
    current_Action = torque(k,1);
    
    % Integrate over the timestep
    [states, slip_cond_end, ~] = curves_stickslip_ode(init_conds,slip_cond_start,current_Action,[t0, t0 + params.dt],params,refine,curves);
    
    % Update slip condition and initial conditions
    slip_cond_start = slip_cond_end;
    init_conds = states(1,1:4);
    stored_states_1(k+1,:) = states(1,:);
    
    % Report progress
    fprintf('Scenario 1: k # %3.0f \n', k);
end

% % Uncomment to visualize a video plot
close all;
titledescription = 'Scenario 1';
videoname = 'Scenario1';
thetaM_1 = curves.tM(stored_states_1(:,1));
video_fun(stored_states_1(:,1),thetaM_1,torque,stored_states_1(:,5),titledescription,params,record,framerate,videoname);

%% 2. Curve-fit model w/ no-slip/slip dynamics and event detection w/ Newton-Raphson routine
% Initialization

init_conds = global_init_conds;
slip_cond_start = curves_startslipcheck_fun(init_conds,torque(1,1),params,curves);
stored_states_2 = zeros(length(tvec),5);

% Step through integration
for k = 1:(length(tvec) - 1)
    % Get time and torque
    t0 = tvec(k,1);
    current_Action = torque(k,1);
    
    % Integrate over the timestep
    [states, slip_cond_end, event_num, failed] = curves_stickslip_nr(init_conds,slip_cond_start,current_Action,[t0, t0 + params.dt],params,curves,refine);
    
    % Update slip condition and initial conditions
    slip_cond_start = slip_cond_end;
    init_conds = states(1,1:4);
    stored_states_2(k+1,:) = states(1,:);
    
    % Report progress
    fprintf('Scenario 2: k # %3.0f, Failed? %3.0f \n', k, failed);
end

% % Uncomment to visualize a video plot
% close all;
% titledescription = 'Scenario 2';
% videoname = 'Scenario2';
% thetaM_2 = curves.tM(stored_states_2(:,1));
% video_fun(stored_states_2(:,1),thetaM_2,torque,stored_states_2(:,5),titledescription,params,record,framerate,videoname);

%% 3. Exact solution model w/ stick/slip dynamics and Matlab's ODE event detection
% Initialization
init_conds = global_init_conds;
slip_cond_start = startslipcheck_fun(init_conds,torque(1,1),guess_Dx,params);
stored_states_3 = zeros(length(tvec),5);

% Step through integration
for k = 1:(length(tvec) - 1)
    % Get time and torque
    t0 = tvec(k,1);
    current_Action = torque(k,1);
    
    % Integrate over the timestep
    [states, slip_cond_end] = stickslip_sim_ode(init_conds,slip_cond_start,current_Action,[t0, t0+params.dt],guess_Dx,params,refine);
    
    % Update slip condition and initial conditions
    slip_cond_start = slip_cond_end;
    init_conds = states(1,1:4);
    stored_states_3(k+1,:) = states(1,:);
    
    % Report progress
    fprintf('Scenario 3: k # %3.0f \n', k);
end

% % Uncomment to visualize a video plot
% close all;
% titledescription = 'Scenario 3';
% videoname = 'Scenario3';
% thetaM_3 = curves.tM(stored_states_3(:,1)); % using curves here to get vehicle angle
% video_fun(stored_states_3(:,1),thetaM_3,torque,stored_states_3(:,5),titledescription,params,record,framerate,videoname);

%% 4. Curve-fit model w/ JUST no-slip dynamics
% Initialization
init_conds = global_init_conds;
stored_states_4 = zeros(length(tvec),1);

% Step through integration
for k = 1:(length(tvec) - 1)
    % Get time and torque
    t0 = tvec(k,1);
    current_Action = torque(k,1);
    
    % Integrate over the timestep
    [tns_c, yns_c] = ode45(@(t,y) curves_no_slip_ode(t,y,current_Action,params,curves), [t0, t0+params.dt], init_conds);
    
    % Update initial conditions
    init_conds = yns_c(end,:);
    stored_states_4(k+1,:) = yns_c(end,1);
    
    % Report progress
    fprintf('Scenario 4: k # %3.0f \n', k);
end

% % Uncomment to visualize a video plot
% close all;
% titledescription = 'Scenario 4';
% videoname = 'Scenario4';
% thetaM_4 = curves.tM(stored_states_4(:,1));
% video_fun(stored_states_4(:,1),thetaM_4,torque,zeros(length(thetaM_4)),titledescription,params,record,framerate,videoname);

%% 5. Exact solution model w/ JUST no-slip dynamics
% Initialization
init_conds = global_init_conds;
stored_states_5 = zeros(length(tvec),1);

% Step through integration
for k = 1:(length(tvec) - 1)
    % Get time and torque
    t0 = tvec(k,1);
    current_Action = torque(k,1);
    
    % Integrate over the timestep
    [tns, yns] = ode45(@(t,y) no_slip_ode(t,y,current_Action,guess_Dx,params), [t0, t0+params.dt], init_conds);
    
    % Update initial conditions
    init_conds = yns(end,:);
    stored_states_5(k+1,:) = yns(end,1);
    
    % Report progress
    fprintf('Scenario 5: k # %3.0f \n', k);
end

% % Uncomment to visualize a video plot
% close all;
% titledescription = 'Scenario 5';
% videoname = 'Scenario5';
% thetaM_5 = curves.tM(stored_states_5(:,1)); % using curves here to get vehicle angle
% video_fun(stored_states_5(:,1),thetaM_5,torque,zeros(length(thetaM_5)),titledescription,params,record,framerate,videoname);

%% Plot and compare the results from the different scenarios

% Comparing Scenarios 1, 2, and 3
figure;
subplot(3,1,1);
plot(tvec,stored_states_1(:,1),'bo','MarkerFaceColor','b'); hold on;
plot(tvec,stored_states_2(:,1),'r-','linewidth',2);
plot(tvec,stored_states_3(:,1),'--','linewidth',2);
yline(params.dim.goal,'k--','LineWidth',2); xlim([0 max(tvec)]);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$x$ (m)','Interpreter','Latex');
title('Position','Interpreter','Latex');
legend('Scenario 1','Scenario 2','Scenario 3','Goal','Interpreter','Latex');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');
subplot(3,1,2);
plot(tvec,stored_states_1(:,5),'bo','MarkerFaceColor','b');hold on;
plot(tvec,stored_states_2(:,5),'r-','linewidth',2);
plot(tvec,stored_states_3(:,5),'--','linewidth',2); xlim([0 max(tvec)]);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$v_{rel}$ (m/s)','Interpreter','Latex');
title('Relative Velocity','Interpreter','Latex');
legend('Scenario 1','Scenario 2','Scenario 3','Interpreter','Latex');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');
subplot(3,1,3);
plot(tvec,torque,'k-','linewidth',2); xlim([0 max(tvec)]);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$\tau_{A}$ (N-m)','Interpreter','Latex');
title('Torque','Interpreter','Latex');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');


% Comparing Scenarios 4 and 5
figure;
subplot(2,1,1);
plot(tvec,stored_states_4(:,1),'b-','linewidth',2); hold on;
plot(tvec,stored_states_5(:,1),'r--','linewidth',2); 
yline(params.dim.goal,'k--','LineWidth',2); ylim([-15 15]); xlim([0 max(tvec)]);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$x$ (m)','Interpreter','Latex');
title('Position','Interpreter','Latex');
legend('Scenario 4','Scenario 5','Goal','Interpreter','Latex');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');
subplot(2,1,2);
plot(tvec,torque,'k-','linewidth',2); xlim([0 max(tvec)]);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$\tau_{A}$ (N-m)','Interpreter','Latex');
title('Torque','Interpreter','Latex');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');

