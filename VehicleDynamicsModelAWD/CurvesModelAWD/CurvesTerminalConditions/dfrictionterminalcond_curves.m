function eventdata = dfrictionterminalcond_curves(y,tauA,tauB,cond,wheel,params,curves)
%% dfrictionterminalcond_curves
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate the friction terminal condition at 
% either wheel A or wheel B and the derivative of the friction terminal
% condition at either wheel A or wheel B. This is event information is used
% by the Newton-Raphson solver to determine the moment the AWD model switches
% between dynamic models.
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
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs: 
%   eventdata: 1x2 array describing the value of the friction terminal condition at 
%           either wheel A or wheel B and the derivative of the friction terminal
%           condition at either wheel A or wheel B

%%

% Unpack the parameters structure
g = params.dim.g;
R = params.dim.R;
I_A = params.dim.I_A;
I_B = params.dim.I_B;

% Initialize numerical derivative stepsize
xA = y(1);
h = R/20;
Dvec = [xA-2*h, xA-h, xA, xA+h, xA+2*h]';

% Initialize numerical derivative computed components
HNk = zeros(size(Dvec));
JNk = zeros(size(Dvec));
TANk = zeros(size(Dvec));
TBNk = zeros(size(Dvec));
Hthetak = zeros(size(Dvec));
Jthetak = zeros(size(Dvec));
TAthetak = zeros(size(Dvec));
TBthetak = zeros(size(Dvec));

% Compute the numerical derivative of dynamic parameters
ymod = y;
for k = 1:length(Dvec)
    
    ymod(1) = Dvec(k,1);
    
    % Get dynamic model parameters
    [x_coeff, FN_coeff, theta_coeff] = curves_dynamicsselectorfunction(ymod,cond,wheel,params,curves);
    
    if k == 3
        H = x_coeff.Hx;
        J = x_coeff.Jx;
        TA = x_coeff.TAx;
        TB = x_coeff.TBx;
    end
    
    HNk(k,1) = FN_coeff.HN;
    JNk(k,1) = FN_coeff.JN;
    TANk(k,1) = FN_coeff.TAN;
    TBNk(k,1) = FN_coeff.TBN;
    
    Hthetak(k,1) = theta_coeff.Htheta;
    Jthetak(k,1) = theta_coeff.Jtheta;
    TAthetak(k,1) = theta_coeff.TAtheta;
    TBthetak(k,1) = theta_coeff.TBtheta;

end

k = 3;
HN = HNk(k,1);
JN = JNk(k,1);
TAN = TANk(k,1);
TBN = TBNk(k,1);
Htheta = Hthetak(k,1);
Jtheta = Jthetak(k,1);
TAtheta = TAthetak(k,1);
TBtheta = TBthetak(k,1);

% Compute the finite difference derivative
HNx = (-HNk(k+2,1) + 8*HNk(k+1,1) - 8*HNk(k-1,1) + HNk(k-2,1))/(12*h);
JNx = (-JNk(k+2,1) + 8*JNk(k+1,1) - 8*JNk(k-1,1) + JNk(k-2,1))/(12*h);
TANx = (-TANk(k+2,1) + 8*TANk(k+1,1) - 8*TANk(k-1,1) + TANk(k-2,1))/(12*h);
TBNx = (-TBNk(k+2,1) + 8*TBNk(k+1,1) - 8*TBNk(k-1,1) + TBNk(k-2,1))/(12*h);
Hthetax = (-Hthetak(k+2,1) + 8*Hthetak(k+1,1) - 8*Hthetak(k-1,1) + Hthetak(k-2,1))/(12*h);
Jthetax = (-Jthetak(k+2,1) + 8*Jthetak(k+1,1) - 8*Jthetak(k-1,1) + Jthetak(k-2,1))/(12*h);
TAthetax = (-TAthetak(k+2,1) + 8*TAthetak(k+1,1) - 8*TAthetak(k-1,1) + TAthetak(k-2,1))/(12*h);
TBthetax = (-TBthetak(k+2,1) + 8*TBthetak(k+1,1) - 8*TBthetak(k-1,1) + TBthetak(k-2,1))/(12*h);

% calculate Fnormal and dFnormal
Fnormal = HN*(y(2)^2) + JN*g + TAN*tauA + TBN*tauB;
dFnormal = (HNx - 2*HN*H)*(y(2)^3) + (JNx - 2*HN*J)*y(2)*g + (TANx - 2*HN*TA)*y(2)*tauA + (TBNx - 2*HN*TB)*y(2)*tauB;

% calculate Ffrction and dFfriction
ddtheta = Htheta*(y(2)^2) + Jtheta*g + TAtheta*tauA + TBtheta*tauB;
if strcmp(wheel,'A')
    Ffriction = (tauA - I_A*ddtheta)/R;
    dFfriction = (-I_A/R)*((Hthetax - 2*Htheta*H)*(y(2)^3) + (Jthetax - 2*Htheta*J)*y(2)*g ...
        + (TAthetax - 2*Htheta*TA)*y(2)*tauA + (TBthetax - 2*Htheta*TB)*y(2)*tauB);
elseif strcmp(wheel,'B')
    Ffriction = (tauB - I_B*ddtheta)/R;
    dFfriction = (-I_B/R)*((Hthetax - 2*Htheta*H)*(y(2)^3) + (Jthetax - 2*Htheta*J)*y(2)*g ...
        + (TAthetax - 2*Htheta*TA)*y(2)*tauA + (TBthetax - 2*Htheta*TB)*y(2)*tauB);
end

% Compute the friction terminal event and its derivative
Fevent = params.dim.mu_s*abs(Fnormal) - abs(Ffriction);
dFevent = params.dim.mu_s*(Fnormal*dFnormal/abs(Fnormal)) - (Ffriction*dFfriction/abs(Ffriction));

eventdata = [Fevent, dFevent];