function dy = slip_ode(t,y,torque,guess_Dx,params,mu_s,order)
%% slip_ode
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
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   mu_s: 1x1 double indicating the maximum friction coefficient for wheel A
%   order: 1x1 double indicating a selector for potentially different
%           friction models, simply set at 2 for now
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
[H_s, J_s, H_N_s, J_N_s, ~, ~] = car_eom_slip(y(1),guess_Dx,mu,params);
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



