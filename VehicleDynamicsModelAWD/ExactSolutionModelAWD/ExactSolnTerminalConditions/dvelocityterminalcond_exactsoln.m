function eventdata = dvelocityterminalcond_exactsoln(y,tauA,tauB,cond,wheel,guess_Dx,params)
%% dvelocityterminalcond_exactsoln
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
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
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
[x_coeff, ~, theta_coeff] = exactsoln_dynamicsselectorfunction(y,cond,wheel,guess_Dx,params);

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
    vevent  = v_relA_fcn(y,params);
    dvevent = R*ddtheta - ddx*phiA*R - phiAx*(y(2)^2)*R;
    
elseif strcmp(wheel,'B')
    
    % take numerical derivative of psiB to get psiBx
    h = R/20;
    xA = y(1);
    Dvec = [xA-2*h, xA-h, xA, xA+h, xA+2*h]';
    psiBsub = zeros(size(Dvec));
    for m = 1:length(Dvec)
        xD = Dvec(m,1);
        psiBsub(m,1) = psiBfun(xD,guess_Dx,params);
    end
    m = 3;
    psiB = psiBsub(m,1);
    psiBx = (-psiBsub(m+2,1) + 8*psiBsub(m+1,1) - 8*psiBsub(m-1,1) + psiBsub(m-2,1))/(12*h);
    
    % Calculate the event and the event derivative
    vevent  = v_relB_fcn(y,guess_Dx,params);
    dvevent = R*ddtheta - ddx*psiB - psiBx*(y(2)^2);
    
end

eventdata = [vevent, dvevent];
