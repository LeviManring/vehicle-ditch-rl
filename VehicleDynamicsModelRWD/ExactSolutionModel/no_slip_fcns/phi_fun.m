function [phi, phi_x] = phi_fun(xA,str,guess_Dx,params)
%% phi_fun
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate phi and phi_x or (phi'(x)). These
% parameters are commonly needed in our dynamic model. We can choose to
% calculate these values either at x-position of wheel A or x-position of
% wheel B
%
% Inputs:
%   xA: 1x1 double indicating the x-position of wheel A
%   str: a string, either 'A' or 'B' indicating which wheel we want to
%           calculate phi and phi_x with respect to
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%   phi: 1x1 double
%   phi_x: 1x1 double

%%

% unpack parameters structure
R = params.dim.R;
l = params.dim.l;

% grab the function derivatives and features of the ditch profile at wheel A
[~, yAx, yAxx, ~, flag] = y_eval_fcn(xA,params);

% Numerically solve for Dx
x_k = xA;

if flag(1) == 0
    Dx_solve = l;
end

if flag(2) == 0 && flag(1) ~= 0
    Dx_solve = fzero(@(Dx) sqrt(Dx^2 + (params.fun.y{1}(x_k + Dx) - params.fun.y{1}(x_k))^2) - l, guess_Dx);
end

if flag(2) == 1
    Dx_solve = fzero(@(Dx) sqrt((Dx + R*((params.fun.y{2}(x_k)/sqrt(1 + params.fun.y{2}(x_k)^2)) - ...
        (params.fun.y{2}(x_k + Dx)/sqrt(1 + params.fun.y{2}(x_k + Dx)^2))))^2 + ...
        (params.fun.y{1}(x_k + Dx) - params.fun.y{1}(x_k) - R*((1/sqrt(1 + params.fun.y{2}(x_k)^2)) - ...
        (1/sqrt(1 + params.fun.y{2}(x_k + Dx)^2))))^2) - l, guess_Dx);
end

dx = Dx_solve;
xB = xA + dx; % determine where wheel B is located

% grab the function derivatives and features of the ditch profile at wheel B
[~, yBx, yBxx, ~, ~] = y_eval_fcn(xB,params);

% Choose phi at either contact for wheel A or contact for wheel B
if str == 'A'
    
    phi = sqrt(yAx^2+1)/R;
    phi_x = yAx*yAxx/(R*sqrt(yAx^2+1));
    
elseif str == 'B'
    
    phi = sqrt(yBx^2+1)/R;
    phi_x = yBx*yBxx/(R*sqrt(yBx^2+1));
    
end







