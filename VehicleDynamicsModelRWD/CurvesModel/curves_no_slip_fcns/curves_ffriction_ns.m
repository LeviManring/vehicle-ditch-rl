function Ff_ns = curves_ffriction_ns(xA,xdotA,torque,params,curves)
%% curves_ffriction_ns
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the friction force 
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
%   Ff_ns: 1x1 double indicating the value for the friction force at  
%           wheel A while wheel A is not slipping

%%

% Get dynamic model parameters
H = curves.H(xA);
J = curves.J(xA);
TA = curves.TA(xA);

% Get surface profile derivatives and parameters
[~, yAx, yAxx, ~, ~] = y_eval_fcn(xA,params);
phiA = sqrt(yAx^2+1)/params.dim.R;
phiAx = yAx*yAxx/(params.dim.R*sqrt(yAx^2+1));

% Calculate the friction force at wheel A
Ff_ns = (1/params.dim.R)*(torque + params.dim.I_A*(phiA*(H*xdotA^2 + J*params.dim.g + TA*torque) - phiAx*xdotA^2));