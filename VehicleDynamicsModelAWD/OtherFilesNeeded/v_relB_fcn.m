function v_relB = v_relB_fcn(y,guess_Dx,params)
%% v_relB_fcn
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the relative velocity between wheel B
% and the surface profile using the dynamic states of the system
%
% Inputs:
%   y: 1xN array inticating the states of the dynamics model from the ode.
%           In particular:
%               y(1) = x
%               y(2) = x'
%               y(3) = theta_A
%               y(4) = theta_A'
%               y(5) = theta_B
%               y(6) = theta_B'  
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs: 
%   v_relB: 1x1 double describing the relative velocity between wheel B and
%           the surface profile

%%

% Get psiB
psiB = psiBfun(y(1),guess_Dx,params);

% Calculate relative velocity
v_relB = params.dim.R*y(6) - y(2)*psiB;