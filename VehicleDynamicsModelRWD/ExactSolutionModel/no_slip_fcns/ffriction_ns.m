function Ff_ns = ffriction_ns(xA,xdotA,torque,guess_Dx,params)
%% ffriction_ns
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the directional friction force at wheel
% A while the vehicle is not slipping.
%
% Inputs:
%   xA: 1x1 double indicating the x-position of wheel A
%   xdotA: 1x1 double indicating the x-velocity of wheel A
%   torque: 1x1 double indicating the torque applied to wheel A
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%   Ff_ns: 1x1 double indicating the value for the directional friction
%   force at wheel A while wheel A is not slipping

%%

% Calculate certain parameters
[H, J, TA, ~, ~, ~, ~, ~] = car_eom_no_slip(xA,guess_Dx,params);
[phiA, phiAx] = phi_fun(xA,'A',guess_Dx,params);

% Calculate the friction force
Ff_ns = (1/params.dim.R)*(torque + params.dim.I_A*(phiA*(H*xdotA^2 + J*params.dim.g + TA*torque) - phiAx*xdotA^2));
