function dvrel = curves_dv_rel_fcn(y,torque,params,curves)
%% curves_dv_rel_fcn
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the derivative of the relative
% velocity between wheel A and the surface profile
%
% Inputs:
%   y: 1xN array inticating the states of the dynamics model from the ode.
%           In particular y(1) = x, y(2) = x', y(3) = theta_A, y(4) = theta_A' 
%   torque: 1x1 double indicating the torque applied to wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs:
%   dFnormal: 1x1 double indicating the value for the derivative of the 
%           normal force at wheel A while wheel A is not slipping

%%

% Unpack the parameters structure
R = params.dim.R;
I_A = params.dim.I_A;

% Calculate the relative velocity
v_rel = v_rel_fcn(y,params);

% Calculate the directional friction coefficient
mu = mu_fcn(v_rel,params.dim.mu_s,params.settings.friction_model);

% Calculate the dynamics model parameters
H_s = curves.H_s(y(1),mu);
J_s = curves.J_s(y(1),mu);
H_N_s = curves.H_N_s(y(1),mu);
J_N_s = curves.J_N_s(y(1),mu);
H_theta = (-R/I_A)*H_N_s;
J_theta = (-R/I_A)*J_N_s;
T_theta = (1/I_A);

% Get the surface profile derivatives and parameters
[~, yAx, yAxx, ~, ~] = y_eval_fcn(y(1),params);
phiA = sqrt(yAx^2+1)/R;
phiAx = yAx*yAxx/(R*sqrt(yAx^2+1));

% Calculate the derivative of relative velocity
ddx = -(H_s*y(2)^2 + J_s*params.dim.g);
ddthetaA = mu*H_theta*y(2)^2 + mu*J_theta*params.dim.g + T_theta*torque;
dvrel = R*ddthetaA- ddx*R*phiA - phiAx*R*y(2)^2;