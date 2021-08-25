function dy = no_slip_ode(t,y,torque,guess_Dx,params)
%% no_slip_ode
% Levi Manring, Duke University
% 2021
%
% This function is used by the ode45 dynamics solver to integrate the
% no-slip dynamics model. For this case we use an augmented state-space by
% artificially forcing an extra DOF with theta_A and theta_Adot. This is to
% maintain the same state-space size when the vehicle slips.
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
%
% Outputs: 
%   dy: 1xN array describing the outputs from the dynamic model statespace derivatives

%%

% Calculate the no-slip dynamic model parameters
[H, J, TA, ~, ~, ~, ~, ~] = car_eom_no_slip(y(1),guess_Dx,params);
[phiA, phiAx] = phi_fun(y(1),'A',guess_Dx,params);

% Calculate the dynamic model states
% states: x, xdot, theta_A, theta_Adot
dy = zeros(4,1); 

dy(1) = y(2);

dy(2) = -(H*(y(2)^2) + J*params.dim.g + TA*torque);

dy(3) = y(4);

dy(4) = -phiA*(H*(y(2)^2) + J*params.dim.g + TA*torque) + (y(2)^2)*phiAx;