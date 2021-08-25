function dy = curves_no_slipminimum_ode(t,y,torque,params,curves)
%% curves_no_slipminimum_ode
% Levi Manring, Duke University
% 2021
%
% This function is used by the ode45 dynamics solver to integrate the
% no-slip dynamics model. For this case we DO NOT use an augmented state-space.
%
% Inputs:
%   t: 1x1 double inidicating the timestep of ode integration
%   y: 1xN array inticating the states of the dynamics model from the ode.
%           In particular y(1) = x, y(2) = x', y(3) = theta_A, y(4) = theta_A' 
%   torque: 1x1 double indicating the torque applied to wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs: 
%   dy: 1xN array describing the outputs from the dynamic model statespace derivatives

%%

% Get the no-slip dynamics model parameters
H = curves.H(y(1));
J = curves.J(y(1));
Ta = curves.TA(y(1));

% Calculate the dynamic model states
% states: x, xdot
dy = zeros(2,1); 

dy(1) = y(2);

dy(2) = -(H*(y(2)^2) + J*params.dim.g + Ta*torque);

