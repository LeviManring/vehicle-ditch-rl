function dFnormal = curves_dfnormal_ns(xA,xdotA,torque,params,curves)
%% curves_dfnormal_ns
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the derivative of the normal force 
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
%   dFnormal: 1x1 double indicating the value for the derivative of the 
%           normal force at wheel A while wheel A is not slipping

%%

% Get dynamic model parameters
FN_H_x_ns = curves.FN_H_x_ns(xA);
FN_J_x_ns = curves.FN_J_x_ns(xA);
FN_TA_x_ns = curves.FN_TA_x_ns(xA);

% Calculate derivative of normal force
dFnormal = FN_H_x_ns*xdotA^3 + FN_J_x_ns*xdotA*params.dim.g + FN_TA_x_ns*xdotA*torque;