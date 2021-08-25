function [NextObs,Reward,IsDone,LoggedSignals] = myNoSlipStepFunction(Action,LoggedSignals,params,curves)
%% myNoSlipStepFunction
% Levi Manring, Duke University
% 2021
%
% This function updates the environment states given the applied action
%
% Inputs:
%   Action: action applied to the rear wheel A
% 	LoggedSignals: a structure of the logged signals
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs:
%   NextObs: 2x1 array containing the ending states for the time-step
%   Reward: 1x1 double indicating the reward achieved this time-step
%   IsDone: 1x1 boulean indicating whether or not the episode is terminated
% 	LoggedSignals: a structure of the logged signals


%% Define Initial States
States = LoggedSignals.State;
init_conds = [States(1), States(2)];

current_Action = double(Action);

start_time = LoggedSignals.Time;

[~, ysim] = ode45(@(t,y) curves_no_slipminimum_ode(t,y,current_Action,params,curves), [start_time, start_time+params.dt],init_conds);

xStates = ysim(end,:);

% Update the logged signals
LoggedSignals.Time = LoggedSignals.Time + params.dt;
LoggedSignals.State = xStates;
LoggedSignals.Action = Action;
LoggedSignals.StepCount = LoggedSignals.StepCount + 1;
fprintf('Step %3.0f \n',LoggedSignals.StepCount);

% Update the observation
NextObs = xStates';

%% Determine Reward
x = NextObs(1);
xdot = NextObs(2);

Reward = RWDrewardfun(x,xdot,params);
LoggedSignals.TotalReward = Reward + LoggedSignals.TotalReward;

%% Determine Terminal Condition
if abs(x) > 17
    IsDone = 1;
else
    IsDone = 0;
end



