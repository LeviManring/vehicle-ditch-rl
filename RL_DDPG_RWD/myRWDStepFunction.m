function [NextObs,Reward,IsDone,LoggedSignals] = myRWDStepFunction(Action,LoggedSignals,params,curves)
%% myRWDStepFunction
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
%   NextObs: 4x1 array containing the ending states for the time-step
%   Reward: 1x1 double indicating the reward achieved this time-step
%   IsDone: 1x1 boulean indicating whether or not the episode is terminated
% 	LoggedSignals: a structure of the logged signals


%% Define Initial States
States = LoggedSignals.State;
init_conds = [States(1), States(2), States(3), States(4)];

% extract the action
current_Action = double(Action);

start_time = LoggedSignals.Time;
slip_cond_start = LoggedSignals.SlipCondEnd;
refine = 20;

% simulate the RWD dynamic model
try
    [xStates, slip_cond_end, event_num, failure] = curves_stickslip_nr(init_conds,slip_cond_start,current_Action,[start_time, start_time + params.dt],params,curves,refine);
    if any(isnan(xStates(:)))
        failure = 1;
    end
catch
    xStates = States;
    slip_cond_end = slip_cond_start;
    event_num = 0;
    failure = 1;
end

% Update the logged signals
LoggedSignals.TotalEventNum = LoggedSignals.TotalEventNum + event_num;
LoggedSignals.SlipCondEnd = slip_cond_end;
LoggedSignals.Time = LoggedSignals.Time + params.dt;
LoggedSignals.State = xStates;
LoggedSignals.Action = Action;
LoggedSignals.StepCount = LoggedSignals.StepCount + 1;
fprintf('Step: %3.0f, x: %3.2f, TauA: %3.4f, VrelA: %3.4f, SlipCond: %3.0f, Event #: %3.0f, Failure?: %3.0f \n',...
    LoggedSignals.StepCount,xStates(1),current_Action,xStates(5),slip_cond_end,event_num,failure);

% Update the observation
NextObs = [xStates(1), xStates(2), xStates(5), slip_cond_end]';

%% Determine Reward
x = NextObs(1);
xdot = NextObs(2);
if slip_cond_end == 0
    vrelA = 0;
else
    vrelA = NextObs(3);
end

if failure == 0
    try
        Reward = RWDrewardfun_augmented(x,xdot,vrelA,slip_cond_end,params);
    catch
        Reward = 0;
        failure = 1;
    end
else
    Reward = 0;
end

LoggedSignals.TotalReward = Reward + LoggedSignals.TotalReward;

%% Determine Terminal Condition
if abs(x) > 17 || failure == 1
    IsDone = 1;
    if failure ~= 1
        if x < 0
            fprintf('Car exited out the back \n');
        else
            fprintf('Car sped out the ditch but didnt stop \n');
        end
    end
else
    IsDone = 0;
end



