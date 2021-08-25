function Reward = AWDrewardfun_augmented(x,xdot,vrelA,vrelB,cond,params)
%% AWDrewardfun_augmented
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
%   vrelB: 1x1 double denoting relative velocity between wheel B and the
%           surface
%   cond: 1x1 double indicating the dynamic scenario:
%           1 -> noslip, 2 -> allslip, 3 -> rearslipfrontstick, 4 -> rearstickfrontslip
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


dist_locate_vel = -(x - goal)^2 + 1;
if dist_locate_vel < 0
    dist_locate_vel = 0;
end

zero_v = 2;
if abs(xdot) > zero_v
    velocity_reward = 0;
else
    if xdot >= 0
        velocity_reward = sinh((-xdot + zero_v)*(3/2))*dist_locate_vel;
    else
        velocity_reward = sinh((xdot + zero_v)*(3/2))*dist_locate_vel;
    end
end

if abs(goal - x) < params.goaltolerance && abs(xdot) < params.veltolerance
    bonus_reward = 100;
else
    bonus_reward = 0;
end

vrelA_reward = -0.001*abs(vrelA);
vrelB_reward = -0.001*abs(vrelB);

cond_reward = 0;
switch cond
    case 1
        cond_reward = 0;
    case 2
        cond_reward = -0.1;
    case 3
        cond_reward = -0.01;
    case 4
        cond_reward = -0.01;
end

Reward = distance_reward + velocity_reward + bonus_reward + vrelA_reward + vrelB_reward + cond_reward;