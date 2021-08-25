function eventdata = dvelocityterminalcond_curves(y,tauA,tauB,cond,wheel,params,curves)
%% dvelocityterminalcond_curves
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the relative velocity at either wheel
% A or wheel B and the derivative of the relative velocity at either wheel
% A or B. This is event information is used by the Newton-Raphson solver to
% determine the moment the AWD model switches between dynamic models.
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
%           indicates both
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs: 
%   eventdata: 1x2 array describing the value of the relative velocity at
%           either wheel A or wheel B and the derivative of the relative
%           velocity at either wheel A or wheel B

%%

% Unpack the parameters structure
g = params.dim.g;
R = params.dim.R;

% Calculate the dynamic parameters based on the scenario
[x_coeff, ~, theta_coeff] = curves_dynamicsselectorfunction(y,cond,wheel,params,curves);

% Unpack the dynamic parameters
H = x_coeff.Hx;
J = x_coeff.Jx;
TA = x_coeff.TAx;
TB = x_coeff.TBx;
Htheta = theta_coeff.Htheta;
Jtheta = theta_coeff.Jtheta;
TAtheta = theta_coeff.TAtheta;
TBtheta = theta_coeff.TBtheta;

% grab the function derivatives and features of the ditch profile at wheel A
[~, yAx, yAxx, ~, ~] = y_eval_fcn(y(1),params);
phiA = sqrt(yAx^2+1)/R;
phiAx = yAx*yAxx/(R*sqrt(yAx^2+1));

% Calculate the 2nd derivatives of x and theta
ddx = -(H*(y(2)^2) + J*g + TA*tauA + TB*tauB);
ddtheta = Htheta*(y(2)^2) + Jtheta*g + TAtheta*tauA + TBtheta*tauB;

if strcmp(wheel,'A')
    
    % Calculate the event and the event derivative
    vevent = v_relA_fcn(y,params);
    dvevent = R*ddtheta - ddx*phiA*R - phiAx*(y(2)^2)*R;
    
elseif strcmp(wheel,'B')
    
    % Calculate the event and the event derivative
    vevent = v_relB_fcn_curves(y,params,curves);
    dvevent = R*ddtheta - ddx*curves.psiB(y(1)) - curves.psiBx(y(1))*(y(2)^2);
    
end

eventdata = [vevent, dvevent];
