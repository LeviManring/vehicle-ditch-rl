function Reward = rewardfun(x,xdot,params)

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
    
% velocity_reward = increase_factor*exp(-(2*xdot/3)^2)*dist_locate_vel;

if abs(goal - x) < params.goaltolerance && abs(xdot) < params.veltolerance
    bonus_reward = 100;
else
    bonus_reward = 0;
end

% could also scale velocity reward between 0 and 1 and multiply by distance
% reward

Reward = distance_reward + velocity_reward + bonus_reward;