function params = f_params(f,x)
%% f_params
% Levi Manring, Duke University
% 2021
%
% This function is used to create the params structure needed to contain
% the masses, moments of intertia, and dimensions of the planar vehicle
% model as well as the first and second derivatives of the surface profile
% the vehicle moves on.
%
% Inputs:
%   f: an anonymous function describing the shape of the surface profile
%           the vehicle moves on
%   x: 1x1 symbolic variable
%
% Outputs:
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model

%%

% f is the function that describes the surface the vehicle moves on
sim.y(1)   = f;
sim.y(2)   = diff(sim.y(1),x);
sim.y(3)   = diff(sim.y(2),x);
sim.y(4)   = diff(sim.y(3),x);

% convert to matlab Functions; this will be what is passed to this Function
if ~isnumeric(sim.y) % check to see if all any functions and derivatives are not just numbers/zeros
    for k = 1:length(sim.y)
        if isnumeric(eval(sim.y(k)))
            params.fun.y{k} = eval(sim.y(k));
        else
            params.fun.y{k} = matlabFunction(sim.y(k));
        end
    end
else % if they are, set the function value equal to the symbolic value
    for k = 1:length(sim.y)
        params.fun.y{k} = sim.y(k);
    end
end

% Set Dimensional parameters for the vehicle
params.dim.R    = 0.3675;           % meters
params.dim.mA   = 25.5*2;             % kg
params.dim.mB   = params.dim.mA;    % meters
params.dim.m    = 2039.25;          % kg
params.dim.I_A  = 1.6851*2;           % kg-m^2
params.dim.I_B  = params.dim.I_A;   % kg-m^2
params.dim.I_M  = 5091;             % kg-m^2
params.dim.g    = 9.81;             % m/s^2
params.dim.l    = 3.517;            % m
params.dim.h    = 1.4925;           % m            
params.dim.xc   = 2.052;            % m
params.dim.yc   = 0.3335;           % m
params.dim.mu_s = 0.25;             % non-dimensional
params.dim.maxT = 700;              % N-m
params.dt       = 0.1;              % s
params.settings.friction_model = 2;
