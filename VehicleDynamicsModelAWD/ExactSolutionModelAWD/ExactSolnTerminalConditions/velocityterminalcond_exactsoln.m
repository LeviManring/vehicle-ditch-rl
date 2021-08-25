function event = velocityterminalcond_exactsoln(y,wheel,guess_Dx,params)
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
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs: 
%   event: 1x1 array describing the value of the relative velocity at
%           either wheel A or wheel B

%%

if strcmp(wheel,'A')
    
    vevent = v_relA_fcn(y,params);
    
elseif strcmp(wheel,'B')
    
    vevent = v_relB_fcn(y,guess_Dx,params);
    
end

event = vevent;
