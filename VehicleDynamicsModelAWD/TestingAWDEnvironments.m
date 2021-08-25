%% TestingAWDEvironments.m
% Levi Manring, Duke University
% 2021
%
% This m-file tests 2 different solver options for a vehicle 
% moving on a user-defined surface profile. In particular, this tests a planar
% model with torque control applied to both the rear wheel A and front
% wheel B. It is assumed that both front and rear wheel can slip, and thus
% we have four possible dynamic scenarios for this vehicle:
%   1 -> noslip, 2 -> allslip, 3 -> rearslipfrontstick, 4 -> rearstickfrontslip
% A Newton-Raphson solver is used to detect an event and transition between
% these four possible dynamic scenarios.
%
% The two solvers take two different approaches to integrating the dynamics
% model. The ExactSolutionModelAWD seeks to numerically integrate the exact
% dynamics at each time step. The CurvesModelAWD uses curvefit stored
% interpolated functions to interpolate dynamic parameters for the
% solution. The curves strategy can take a long time to initially calculate
% the curvefit interpolated functions, thus it is more convenient to
% calculate these initially, and then save/store them to load in the
% future.
%
% The ExactSolutionMOdelAWD appears to outperform the CurvesModel slightly
% in computational performance, and it is slightly more accurate. A
% comparison of the two models will show nearly identical performance.
%
% To see a video plot of the results, simply uncomment the video plotting
% line in each section.
%
% NOTICE: To visually compare these models WITHOUT running this
% computationally time-intensive simulation of both models, simply go to
% the plotting section and uncomment the lines allowing the user to load
% the presimulated data. The data recorded in the m-file
% 'SimulatedAWDModels' shows nearly identical performance for these two
% models. The visual results can be seen in folder 'AWDExample'

%% Initialize
clear all; close all; clc;

% Generate the paths for the RWD vehicle dynamics model (modify this line
% for your file directory)
addpath(genpath('/Users/levimanring/Documents/Duke University/Research/RL Car/VehicleDynamicsModelAWD'));

%% Load or calculate VehicleParameters if these are not loaded previously

% % if calculating uncomment these lines
% syms x
% f = (exp(x^2*(-16/225))*3-3)*(x^2/1225); % user can define f as any continuous function of x
% params = f_params(f,x);
% GetCurvesAWD;

% % otherwise load these lines
load('VehicleParameters');
% load('CurvesAllModels');

% % % NOTE: if loading or generating the curves model takes too long, you
% can always start by just simulating the exact solution model

%% Initialize the torque profiles
dt = 0.1;
params.dt = dt;
tvec = (0:dt:20)';
global_init_conds = [0, 0, 0, 0, 0, 0];
cond_start = 1; % assume not slipping initially, this will be confirmed later
refine = 20;
event_numtotal = 0;

% create sinusoidal profiles for the torques at wheel A and wheel B
omegaA = 6/(2*pi);
omegaB = 7/(2*pi);
phaseA = rand;
phaseB = rand;
tauAvec = 1*params.dim.maxT*sin(omegaA*tvec + phaseA*2*pi);
tauBvec = 1*params.dim.maxT*sin(omegaB*tvec + phaseB*2*pi);

%% Simulate the exact solution model
totalstates_exactsoln = zeros(length(tvec),8);
init_conds = global_init_conds;
tic % check how long this takes just for reference
guess_Dx = params.dim.l;
for k = 1:(length(tvec)-1)
    
    tauA = tauAvec(k,1);
    tauB = tauBvec(k,1);
    t_int = [tvec(k,1), tvec(k+1,1)];
    
    % Simulate the timestep
    [states, cond_end, event_num, failed] = exactsoln_dynamicsAWD_nr(init_conds,cond_start,tauA,tauB,t_int,guess_Dx,params,refine);
    % Give a progress report
    fprintf('Exact Soln Progress:- k: %3.0f, Cond: %3.0f, Eventnum: %3.0f, Failed? %3.0f \n',k,cond_end,event_num, failed);
    
    % Store the new states and update the initial conditions and slipping
    % condition
    totalstates_exactsoln(k+1,:) = states;
    init_conds = states(1:6);
    cond_start = cond_end;
    event_numtotal = event_numtotal + event_num;
    
end
toc

%% Simulate the curvefit model
totalstates_curves = zeros(length(tvec),8);
init_conds = global_init_conds;
cond_start = 1;
tic % check how long this takes just for reference
guess_Dx = params.dim.l;
for k = 1:(length(tvec)-1)
    
    tauA = tauAvec(k,1);
    tauB = tauBvec(k,1);
    t_int = [tvec(k,1), tvec(k+1,1)];
    
    % Simulate the timestep
    [states, cond_end, event_num, failed] = curves_dynamicsAWD_nr(init_conds,cond_start,tauA,tauB,t_int,params,curves,refine);
    % Give a progress report
    fprintf('Curves Progress:- k: %3.0f, Cond: %3.0f, Eventnum: %3.0f, Failed? %3.0f \n',k,cond_end,event_num, failed);
    
    % Store the new states and update the initial conditions and slipping
    % condition
    totalstates_curves(k+1,:) = states;
    init_conds = states(1:6);
    cond_start = cond_end;
    event_numtotal = event_numtotal + event_num;
    
end
toc

%% Plot the results
% To avoid simulation (see note in the title) uncomment this line
% load('SimulatedAWDModels');

% Get the position and relative velocity results from exactsoln
x_es = totalstates_exactsoln(:,1);
vrelA_es = totalstates_exactsoln(:,7);
vrelB_es = totalstates_exactsoln(:,8);

% Get the position and relative velocity results from curves
x_c = totalstates_curves(:,1);
vrelA_c = totalstates_curves(:,7);
vrelB_c = totalstates_curves(:,8);

% Create Comparison Plot
figure; 
subplot(2,2,1);
plot(tvec,x_es,'b-','linewidth',2); hold on;
plot(tvec,x_c,'r--','linewidth',2);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$x$ (m)','Interpreter','Latex');
title('Vehicle Position','Interpreter','Latex');
legend('$x$ Exact Soln','$x$ Curves','Interpreter','Latex','Location','Northwest');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');
subplot(2,2,2);
plot(tvec,tauAvec,'c-','linewidth',2); hold on;
plot(tvec,tauBvec,'m-','linewidth',2);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$\tau$ (N-m)','Interpreter','Latex');
title('Applied Torques','Interpreter','Latex');
legend('$\tau_A$','$\tau_B$','Interpreter','Latex');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');
subplot(2,2,3);
plot(tvec,vrelA_es,'b-','linewidth',2); hold on;
plot(tvec,vrelA_c,'r--','linewidth',2); ylim([-max(abs(vrelA_es)), max(abs(vrelA_es))]);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$v_{rel}$ (m/s)','Interpreter','Latex');
title('Relative Velocity at Wheel A, $v_{rel,A}$','Interpreter','Latex');
legend('$v_{rel,A}$ Exact Soln','$v_{rel,A}$ Curves','Interpreter','Latex','Location','Northwest');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');
subplot(2,2,4);
plot(tvec,vrelB_es,'b-','linewidth',2); hold on;
plot(tvec,vrelB_c,'r--','linewidth',2); ylim([-max(abs(vrelB_es)), max(abs(vrelB_es))]);
xlabel('$t$ (s)','Interpreter','Latex');
ylabel('$v_{rel}$ (m/s)','Interpreter','Latex');
title('Relative Velocity at Wheel B, $v_{rel,B}$','Interpreter','Latex');
legend('$v_{rel,B}$ Exact Soln','$v_{rel,B}$ Curves','Interpreter','Latex','Location','Northwest');
set(gca,'fontsize',18,'TickLabelInterpreter','Latex');

% Create Video
record = 0; % choose to record or not, 1 means record, 0 means don't record
framerate = 10;

videoname = 'Exact Solution';
titledescription = 'Exact Solution';
video_fun_AWD(x_es,curves.tM(x_es),tauAvec,tauBvec,titledescription,params,record,framerate,videoname);

videoname = 'Curves Solution';
titledescription = 'Curves Solution';
video_fun_AWD(x_c,curves.tM(x_c),tauAvec,tauBvec,titledescription,params,record,framerate,videoname);
