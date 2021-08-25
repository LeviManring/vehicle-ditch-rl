function slip_cond_start = curves_startslipcheck_fun(init_cond,init_Torque,params,curves)
%% curves_startslipcheck_fun
% Levi Manring, Duke University
% 2021
%
% This function is used to check if the vehicle begins slipping depending 
% on the torque applied to the wheel A at the start of the time step
%
% Inputs:
%   init_cond: 1xN array inticating the initial states of the dynamics
%           model at the beginning of a given timestep
%   init_Torque: 1x1 double inidicating the starting applied torque at
%           wheel A at the beginning of a given timestep
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs:
%   slip_cond_start: 1x1 boulean describing whether or not the vehicle is
%           slipping, 1 means it is slipping, 0 means it is not slipping

%%

% Calculate the Normal and Friction Forces
Fn_ns_init = curves_fnormal_ns(init_cond(1),init_cond(2),init_Torque,params,curves);
Ff_ns_init = curves_ffriction_ns(init_cond(1),init_cond(2),init_Torque,params,curves);

% Unpack the static friction coefficient
mu_s = params.dim.mu_s;

if abs(Ff_ns_init) > mu_s*abs(Fn_ns_init)
    slip_cond_start = 1;
else
    slip_cond_start = 0;
end