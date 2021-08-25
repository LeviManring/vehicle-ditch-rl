function Fn_s = fnormal_s(xA,xdotA,guess_Dx,mu,params)
%% fnormal_s
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the normal force at wheel A while the
% vehicle is slipping.
%
% Inputs:
%   xA: 1x1 double indicating the x-position of wheel A
%   xdotA: 1x1 double indicating the x-velocity of wheel A
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   mu: 1x1 double indicating the directional friction coefficient for wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%   Fn_s: 1x1 double indicating the value for the normal force at wheel A
%           while wheel A is slipping

%%

[~, ~, H_N_s, J_N_s, ~, ~] = car_eom_slip(xA,guess_Dx,mu,params);

Fn_s = H_N_s*xdotA^2 + J_N_s*params.dim.g;