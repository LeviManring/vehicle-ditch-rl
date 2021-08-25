function [InitialObservation,LoggedSignals] = myRWDResetFunction(params)
%% myAWDResetFunction
% Levi Manring, Duke University
% 2021
%
% This function is used to determine the initial states at the start of a
% RL training episode
%
% Inputs:
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%   InitialObservation: 4x1 array containing the initial states for the
%           training episode
% 	LoggedSignals: a structure of the logged signals

%% Setup the initial conditions
rand_perturb_upper = params.init_state_mean + params.init_state_bound;
rand_perturb_lower = params.init_state_mean - params.init_state_bound;
rand_vec = (rand_perturb_upper - rand_perturb_lower).*rand(2,1) + rand_perturb_lower;
x = rand_vec(1);
xdot = rand_vec(2);

[~, yAx, ~, ~, ~] = y_eval_fcn(x,params);
phiA = sqrt(yAx^2+1)/params.dim.R;
theta_A = 0;
thetadot_A = phiA*xdot;
complex_init_conds = [x, xdot, theta_A, thetadot_A];

v_rel = v_rel_fcn(complex_init_conds,params);

LoggedSignals.Time = 0;
LoggedSignals.State = [complex_init_conds, v_rel];
LoggedSignals.TotalReward = 0;
LoggedSignals.Action = 0;
LoggedSignals.SlipCondEnd = 0;
LoggedSignals.TotalEventNum = 0;
LoggedSignals.StepCount = 0;

InitialObservation = [x, xdot, v_rel, 0]';
