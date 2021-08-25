function [NextObs,Reward,IsDone,LoggedSignals] = myAWDStepFunction(Action,LoggedSignals,params)
%% myAWDStepFunction
% Levi Manring, Duke University
% 2021
%
% This function updates the environment states given the applied action
%
% Inputs:
%   Action: an action vector containing the torques applied to front and
%           rear wheels
% 	LoggedSignals: a structure of the logged signals
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%   NextObs: 5x1 array containing the ending states for the time-step
%   Reward: 1x1 double indicating the reward achieved this time-step
%   IsDone: 1x1 boulean indicating whether or not the episode is terminated
% 	LoggedSignals: a structure of the logged signals


%% Define Initial States 
States = LoggedSignals.State;
init_conds = [States(1), States(2), States(3), States(4), States(5), States(6)];

% extract the action
tauA = double(Action(1));
tauB = double(Action(2));

start_time = LoggedSignals.Time;
cond_start = LoggedSignals.CondEnd;

t_int = [start_time, start_time + params.dt];
refine = 20;
guess_Dx = params.dim.l;

% simulate the AWD dynamic model
try
    [xStates, cond_end, event_num, failure] = exactsoln_dynamicsAWD_nr(init_conds,cond_start,tauA,tauB,t_int,guess_Dx,params,refine);
catch
    xStates = States;
    cond_end = cond_start;
    event_num = 0;
    failure = 1;
end

% Update the logged signals
LoggedSignals.TotalEventNum = LoggedSignals.TotalEventNum + event_num;
LoggedSignals.CondEnd = cond_end;
LoggedSignals.Time = t_int(2);
LoggedSignals.State = [xStates, cond_end];
LoggedSignals.Action = Action;
LoggedSignals.StepCount = LoggedSignals.StepCount + 1;
fprintf('Step: %3.0f, x: %3.2f, TauA: %3.4f, TauB: %3.4f, VrelA: %3.4f, VrelB: %3.4f, Cond: %3.0f, Event #: %3.0f, Failure?: %3.0f \n',...
    LoggedSignals.StepCount,xStates(1),tauA,tauB,xStates(7),xStates(8),cond_end,event_num,failure);

% Update the observation
NextObs = [xStates(1), xStates(2), xStates(7), xStates(8), cond_end]';

%% Determine Reward
x = NextObs(1);
xdot = NextObs(2);
if cond_end == 1 || cond_end == 4
    vrelA = 0;
else
    vrelA = NextObs(3);
end
if cond_end == 1 || cond_end == 3
    vrelB = 0;
else
    vrelB = NextObs(4);
end

if failure == 0
    try
        Reward = AWDrewardfun_augmented(x,xdot,vrelA,vrelB,cond_end,params);
    catch
        Reward = 0;
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



