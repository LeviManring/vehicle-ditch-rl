function dy = curves_slip_ode(t,y,torque,params,mu_s,order,curves)
%% curves_slip_ode
% Levi Manring, Duke University
% 2021
%
% This function is used by the ode45 dynamics solver to integrate the
% slip dynamics model.
%
% Inputs:
%   t: 1x1 double inidicating the timestep of ode integration
%   y: 1xN array inticating the states of the dynamics model from the ode.
%           In particular y(1) = x, y(2) = x', y(3) = theta_A, y(4) = theta_A' 
%   torque: 1x1 double indicating the torque applied to wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   mu_s: 1x1 double indicating the maximum friction coefficient for wheel A
%   order: 1x1 double indicating a selector for potentially different
%           friction models, simply set at 2 for now
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs: 
%   dy: 1xN array describing the outputs from the dynamic model statespace derivatives

%%

% Unpack the parameters structure
R = params.dim.R;
I_A = params.dim.I_A;

% Calculate the relative velocity
v_rel = v_rel_fcn(y,params);

% Calculate the directional friction coefficient
mu = mu_fcn(v_rel,mu_s,order);

% Calculate the dynamics model parameters
H_s = curves.H_s(y(1),mu);
J_s = curves.J_s(y(1),mu);
H_N_s = curves.H_N_s(y(1),mu);
J_N_s = curves.J_N_s(y(1),mu);
H_theta = (-R/I_A)*H_N_s;
J_theta = (-R/I_A)*J_N_s;
T_theta = (1/I_A);

% Calculate the dynamic model states
% states: x, xdot, theta_A, theta_Adot
dy = zeros(4,1);

dy(1) = y(2);

dy(2) = -(H_s*y(2)^2 + J_s*params.dim.g);

dy(3) = y(4);

dy(4) = mu*H_theta*y(2)^2 + mu*J_theta*params.dim.g + T_theta*torque;



