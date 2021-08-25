function event = velocityterminalcond_curves(y,wheel,params,curves)
%% velocityterminalcond_exactsoln
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the relative velocity at either wheel
% A or wheel B. This is a terminal condition used to switch dynamic models,
% hence why the output is called 'event'.
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
%   wheel: a string, selecting which wheel to output dynamic parameters
%           for. 'A' indicates wheel A, 'B' indicates wheel B, and 'AB
%           indicates both
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs: 
%   event: 1x1 array describing the value of the relative velocity at
%           either wheel A or wheel B

%%

if strcmp(wheel,'A')
    
    vevent = v_relA_fcn(y,params);
    
elseif strcmp(wheel,'B')
    
    vevent = v_relB_fcn_curves(y,params,curves);
    
end

event = vevent;
