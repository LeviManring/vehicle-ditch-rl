function [value,isterminal,direction] = curves_slip_event(t,y,torque,mu_s,params,curves)
%% curves_slip_event
% Levi Manring, Duke University
% 2021
%
% This function is used by the curves_no_slip_ode to check and see if the vehicle
% starts slipping.
%
% Inputs:
%   t: 1x1 double inidicating the timestep of ode integration
%   y: 1xN array inticating the states of the dynamics model from the ode
%   torque: 1x1 double indicating the torque applied to wheel A at the
%           timestep in question
%   mu_s: 1x1 double indicating the maximum friction coefficient for wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs: see https://www.mathworks.com/help/matlab/math/ode-event-location.html
%   value: 1x1 double describing the value of the event function
%   isterminal: 1x1 boulean determining termination of integration (1
%           terminates, 0 continues)
%   direction: 1x1 double, 0 if all events are to be found, +1 if only
%           events while event function is increasing, -1 if only events while
%           event function is decreasing

%%

% Calculate the normal and friction forces at wheel A
Fn_ns = curves_fnormal_ns(y(1),y(2),torque,params,curves);
Ff_ns = curves_ffriction_ns(y(1),y(2),torque,params,curves);

% Determine when the normal force Ff_ns > mu_s*Fn_ns
value = mu_s*abs(Fn_ns) - abs(Ff_ns);

% Choose to output certain parameters to the command window so as to
% visually determine the progress
if abs(value) < 100
    fprintf('%4.10f y, %4.10f ydot, %4.12f tau, %3.8f value \n',y(1), y(2),torque,value);
end

isterminal = 1; % stop the integration once the event is detected
direction = 0;

end
