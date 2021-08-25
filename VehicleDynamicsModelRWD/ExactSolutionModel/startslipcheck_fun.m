function slip_cond_start = startslipcheck_fun(init_cond,init_Torque,guess_Dx,params)
%% startslipcheck_fun
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
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%   slip_cond_start: 1x1 boulean describing whether or not the vehicle is
%           slipping, 1 means it is slipping, 0 means it is not slipping

%%

% Calculate the initial normal and frictional forces at wheel A
Fn_ns_init = fnormal_ns(init_cond(1),init_cond(2),init_Torque,guess_Dx,params);
Ff_ns_init = ffriction_ns(init_cond(1),init_cond(2),init_Torque,guess_Dx,params);

% Unpack the static friction coefficient
mu_s = params.dim.mu_s;

% Check to see if the vehicle is slipping or not
if abs(Ff_ns_init) > mu_s*abs(Fn_ns_init)
    slip_cond_start = 1;
else
    slip_cond_start = 0;
end