%% RWDlearn.m
% Levi Manring, Duke University
% 2021
%
% This m-file is used to train a RWD vehicle dynamics model to achieve
% exit from a ditch. The environment, agent, and training parameters are
% all defined here as well as how to simulate the trained agent.

%% Initialization
clear all
close all
clc

% add the relevant paths to the directory
addpath(genpath('/VehicleDynamicsModelRWD'));
load('VehicleParameters');
load('CurveFitModel');

%% Setting up the environment
% Define the observations
ObservationInfo = rlNumericSpec([4 1]);
ObservationInfo.Name = 'CarDitch States';
ObservationInfo.Description = 'x, xdot, vrelA, slipcond';

% Define the actions
ActionInfo = rlNumericSpec([1 1]);
ActionInfo.LowerLimit = -params.dim.maxT;
ActionInfo.UpperLimit = params.dim.maxT;
ActionInfo.Name = 'CarDitch Action';

% Define relevant parameters for learning
params.init_state_mean = [0, 0]';
params.init_state_bound = [0.1, 0.05]';
params.goal = 9.5;
params.goaltolerance = 0.1;
params.veltolerance = 0.1;
params.dt = 0.1;
Time_to_run = 28; % seconds
steps_to_run = Time_to_run*(1/params.dt) + 1;
episodecount = 15000;
params.experience = steps_to_run*episodecount;

% Define the reset and step functions
ResetHandle = @() myRWDResetFunction(params);
StepHandle = @(Action,LoggedSignals) myRWDStepFunction(Action,LoggedSignals,params,curves);

env = rlFunctionEnv(ObservationInfo,ActionInfo,StepHandle,ResetHandle);
        
%% Set up the agent
agent = getDDPGAgent(env,params);

%% Set up the training options
% % uncomment this code to restart training
% resumeTraining = true;
% agentOpts.resetExperienceBufferBeforeTraining = ~(resumeTraining);
% if resumeTraining
%     load('saved_agent');
%     agent = saved_agent;
% end

trainOpts = rlTrainingOptions;
trainOpts.MaxEpisodes = episodecount;
trainOpts.MaxStepsPerEpisode = steps_to_run;
trainOpts.StopTrainingCriteria = "EpisodeCount";
trainOpts.StopTrainingValue = episodecount;
trainOpts.SaveAgentCriteria = "EpisodeReward";
trainOpts.SaveAgentValue = 300;
trainOpts.Verbose = true;
trainOpts.Plots = "training-progress";

%% Train the Agent
trainStats = train(agent,env,trainOpts);

%% Simulate the trained agent
% simOptions = rlSimulationOptions('MaxSteps',steps_to_run);
% experience = sim(env,saved_agent,simOptions);
% 
% % Unpack Data
% times = experience.Observation.CarDitchStates.Time;
% action = zeros(size(times));
% states = zeros(length(times),4);
% reward = experience.Reward.Data;
% reward = [0; reward];
% for k = 1:(length(experience.Observation.CarDitchStates.Data) - 1)
%     action(k+1,1) = experience.Action.CarDitchAction.Data(:,:,k);
% end
% 
% for k = 1:length(experience.Observation.CarDitchStates.Data)
%     states(k,:) = experience.Observation.CarDitchStates.Data(:,:,k)';
% end
% 
% figure; 
% subplot(3,1,1);
% plot(times, states(:,1)); hold on
% plot([0, times(end)], [params.goal, params.goal], '--'); title('pos');
% subplot(3,1,2);
% plot(times, action); title('action');
% subplot(3,1,3);
% plot(times,reward); title('reward');
% 
% save('AgentWeLike','times','states','action','reward','saved_agent','savedAgentResultStruct');
