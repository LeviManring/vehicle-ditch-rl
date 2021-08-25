function v_relA = v_relA_fcn(y,params)
%% v_relA_fcn
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the relative velocity between wheel A
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
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs: 
%   v_relA: 1x1 double describing the relative velocity between wheel B and
%           the surface profile

%%

% Get surface profile derivative
[~, yAx, ~, ~, ~] = y_eval_fcn(y(1),params);

% Calculate relative velocity
v_relA = params.dim.R*y(4) - y(2)*sqrt(1 + yAx^2);