function Fn_ns = curves_fnormal_ns(xA,xdotA,torque,params,curves)
%% curves_fnormal_ns
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the normal force 
% at wheel A while the vehicle is not slipping.
%
% Inputs:
%   xA: 1x1 double indicating the x-position of wheel A
%   xdotA: 1x1 double indicating the x-velocity of wheel A
%   torque: 1x1 double indicating the torque applied to wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs:
%   Fn_ns: 1x1 double indicating the value for the normal force at  
%           wheel A while wheel A is not slipping

%%

% Get dynamic model parameters
H_N_ns = curves.H_N_ns(xA);
J_N_ns = curves.J_N_ns(xA);
TA_N_ns = curves.TA_N_ns(xA);

% Calculate the normal force at wheel A
Fn_ns = H_N_ns*xdotA^2 + J_N_ns*params.dim.g + TA_N_ns*torque;