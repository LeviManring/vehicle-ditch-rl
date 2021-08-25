function Fn_s = curves_fnormal_s(xA,xdotA,params,curves)
%% curves_fnormal_s
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the normal force 
% at wheel A while the vehicle is slipping.
%
% Inputs:
%   xA: 1x1 double indicating the x-position of wheel A
%   xdotA: 1x1 double indicating the x-velocity of wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs:
%   Fn_s: 1x1 double indicating the value for the normal force at  
%           wheel A while wheel A is slipping

%%

% Get dynamic model parameters
H_N_s = curves.H_N_s(xA);
J_N_s = curves.J_N_s(xA);

% Calculate the normal force at wheel A
Fn_s = H_N_s*xdotA^2 + J_N_s*params.dim.g;