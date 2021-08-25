%% pilcoRWD_learn.m
% Levi Manring, Duke University
% 2021
%
% This is the primary learn function used to implement PILCO
%
% This code is modified for our scenario from 
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
%%

% Initialization
clear all;
close all; clc;

pilcoRWD_settings; % load scenario-specific settings

% filename used for saving data
basename = 'PILCOData';

% Initial J random rollouts
for jj = 1:J
    [xx, yy, realCost{jj}, latent{jj}] = ...
        rollout(gaussian(mu0, S0), struct('maxU',policy.maxU), H, plant, cost);
    x = [x; xx]; y = [y; yy];       % augment training sets for dynamics model
end

mu0Sim(odei,:) = mu0; S0Sim(odei,odei) = S0;
mu0Sim = mu0Sim(dyno); S0Sim = S0Sim(dyno,dyno);

% Controlled learning (N iterations)
j = 1;
while j < N
    close all
    trainDynModel;   % train (GP) dynamics model
    learnPolicy;     % learn policy
    applyController; % apply controller to system
    disp(['controlled trial # ' num2str(j)]);
    
    data = latent{j};
    
    X1 = data(:,1);
    X2 = data(:,2);
    Ta = data(:,3);
    time = ((1:1:length(X1)) - 1)'*params.dt;
    reward = zeros(length(X1),1);
    Force_scaled = Ta./max(abs(Ta));
    
    for m = 1:length(X1)
        reward(m,1) = 1-cost.fcn(cost,[X1(m,1), 0]', zeros(2));
    end
    
    experience = j*T;
    trial_num = j;
    
    hfig = figure;
    subplot(3,1,1);
    plot(time,X1); hold on;
    plot([time(1), time(end)],[cost.target(1), cost.target(1)],'k--');
    subplot(3,1,2);
    plot(time,Force_scaled);
    subplot(3,1,3);
    plot(time,reward);
    
    j = j + 1;
 
end


