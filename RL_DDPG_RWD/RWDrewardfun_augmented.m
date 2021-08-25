function Reward = RWDrewardfun_augmented(x,xdot,vrelA,slipcond,params)
%% RWDrewardfun_augmented
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the total reward given some states for
% the AWD vehicle dynamics training problem
%
% Inputs:
%   x: 1x1 double denoting vehicle position
%   xdot: 1x1 double denoting vehicle velosity
%   vrelA: 1x1 double denoting relative velocity between wheel A and the
%           surface
%   slipcond: 1x1 double indicating whether or not wheel A is slipping
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%   Reward: 1x1 double indicating the achieve reward

goal = params.goal;

if x >= 0 && x <= goal
    distance_reward = sinh(x/pi);
end
if x < 0
    distance_reward = 0;
end
if x > goal
    if x < 2*goal
        distance_reward = sinh(-(x - 2*goal)/pi);
    else
        distance_reward = 0;
    end
end

zero_v = 10;
if x > -goal/2 && x < goal && abs(xdot) < zero_v
    if xdot > 0
        velocity_add_reward = 0.25*abs(xdot);
    else
        velocity_add_reward = 0.1*abs(xdot);
    end
else
    velocity_add_reward = 0;
end

dist_locate_vel = -(x - goal)^2 + 1;
if dist_locate_vel < 0
    dist_locate_vel = 0;
end

if abs(xdot) > zero_v
    velocity_reward = 0;
else
    if xdot >= 0
        velocity_reward = sinh((-xdot + zero_v)*(1/3))*dist_locate_vel;
    else
        velocity_reward = sinh((xdot + zero_v)*(1/3))*dist_locate_vel;
    end
end

if abs(goal - x) < params.goaltolerance && abs(xdot) < params.veltolerance
    bonus_reward = 100;
else
    bonus_reward = 0;
end

vrelA_reward = -0.001*abs(vrelA);
if slipcond == 1
    slipcond_reward = -0.1;
else
    slipcond_reward = 0;
end

Reward = distance_reward + velocity_add_reward + velocity_reward + bonus_reward + vrelA_reward + slipcond_reward;