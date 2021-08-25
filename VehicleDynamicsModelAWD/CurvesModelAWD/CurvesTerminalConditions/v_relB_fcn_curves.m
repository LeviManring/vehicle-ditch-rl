function v_relB = v_relB_fcn_curves(y,params,curves)
%% v_relB_fcn_curves
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
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs: 
%   v_relB: 1x1 double describing the relative velocity between wheel B and
%           the surface profile

%%

% Get psiB
psiB = curves.psiB(y(1));

% Calculate relative velocity
v_relB = params.dim.R*y(6) - y(2)*psiB;