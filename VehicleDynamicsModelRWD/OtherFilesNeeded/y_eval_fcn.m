function [yK, yKx, yKxx, yKxxx, flag] = y_eval_fcn(xK,params)
%% y_eval_fcn
% Levi Manring, Duke University
% 2021
%
% This function is used to evaluate the surface profile function for the
% vehicle as well as it's 1st, 2nd, and 3rd derivatives with respect
% to x. Additionally, it determines a flag indicating the nature of the
% function.
%
% Inputs:
%   xK: 1x1 double indicating the location x we want to evaluate the
%           surface profile at
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs: 
%   yK: 1x1 double, surface profile function evaluated at xK
%   yKx: 1x1 double, 1st derivative with respect to x of the surface
%       profile function evaluated at xK
%   yKxx: 1x1 double, 2st derivative with respect to x of the surface
%       profile function evaluated at xK
%   yKxxx: 1x1 double, 3st derivative with respect to x of the surface
%       profile function evaluated at xK
%   flag: 1x1 double indicating how the params.fun.y{k} should be evaluated

%%

% Initialization
flag = zeros(length(params.fun.y),1);
y_p = zeros(length(params.fun.y),1);

% Evaluate the function and it's derivatives by using the params structure
for k = 1:length(params.fun.y)
    if isnumeric(params.fun.y{k})
        flag(k,1) = 0;
        y_p(k,1) = params.fun.y{k};
    else
        flag(k,1) = 1;
        y_p(k,1) = params.fun.y{k}(xK);
    end
end

% Output
yK      = y_p(1,1);
yKx     = y_p(2,1);
yKxx    = y_p(3,1);
yKxxx   = y_p(4,1);