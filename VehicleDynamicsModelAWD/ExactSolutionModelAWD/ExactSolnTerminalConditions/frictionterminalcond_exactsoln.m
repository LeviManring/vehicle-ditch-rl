function event = frictionterminalcond_exactsoln(y,tauA,tauB,cond,wheel,guess_Dx,params)
%% frictionterminalcond_exactsoln
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the friction terminal condition at 
% either wheel A or wheel B.
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
%   tauA: 1x1 double indicating the torque applied to wheel A
%   tauB: 1x1 double indicating the torque applied to wheel B
%   cond: 1x1 double indicating the dynamic scenario: 
%           1 -> noslip, 2 -> allslip, 3 -> rearslipfrontstick, 4 -> rearstickfrontslip
%   wheel: a string, selecting which wheel to output dynamic parameters
%           for. 'A' indicates wheel A, 'B' indicates wheel B, and 'AB
%           indicates both. 
%           NOTE: Can only input wheel 'A' or 'B', not 'AB' in THIS instance
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs: 
%   event: 1x1 array describing the value of the friction terminal condition at
%           either wheel A or wheel B

%%

% Unpack the parameters structure
g = params.dim.g;
R = params.dim.R;
I_A = params.dim.I_A;
I_B = params.dim.I_B;

% Calculate the dynamic parameters based on the scenario
[~, FN_coeff, theta_coeff] = exactsoln_dynamicsselectorfunction(y,cond,wheel,guess_Dx,params);

% Unpack the dynamic parameters
HN = FN_coeff.HN;
JN = FN_coeff.JN;
TAN = FN_coeff.TAN;
TBN = FN_coeff.TBN;
Htheta = theta_coeff.Htheta;
Jtheta = theta_coeff.Jtheta;
TAtheta = theta_coeff.TAtheta;
TBtheta = theta_coeff.TBtheta;

% calculate the Normal Force
Fnormal = HN*(y(2)^2) + JN*g + TAN*tauA + TBN*tauB;

% calculate Friction Force
ddtheta = Htheta*(y(2)^2) + Jtheta*g + TAtheta*tauA + TBtheta*tauB;
if strcmp(wheel,'A')
    
    Ffriction = (tauA - I_A*ddtheta)/R;
    
elseif strcmp(wheel,'B')
    
    Ffriction = (tauB - I_B*ddtheta)/R;
    
end

% Calculate the friction terminal event
Fevent = params.dim.mu_s*abs(Fnormal) - abs(Ffriction);
event = Fevent;