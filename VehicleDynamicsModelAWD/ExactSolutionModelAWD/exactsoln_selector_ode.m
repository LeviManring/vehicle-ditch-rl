function dy = exactsoln_selector_ode(t,y,tauA,tauB,cond,params)
%% exactsoln_selector_ode
% Levi Manring, Duke University
% 2021
%
% This function is used by the ode45 dynamics solver to integrate the
% AWD dynamics model. This function allows switching between different
% dynamic parameters depending on the 'cond' value (designating 1 of 4
% possible scenarios: noslip, allslip, rearslipfrontstick, and
% rearstickfrontslip). The governing dynamic equations are of this form:
%           x'' + Hx*x'^2 + Jx*g + TAx*tau_A + TBx*tau_B = 0
%           thetaA'' = HthetaA*x'^2 + JthetaA*g + TAthetaA*tau_A + TBthetaA*tau_B
%           thetaB'' = HthetaB*x'^2 + JthetaB*g + TAthetaB*tau_A + TBthetaB*tau_B
%
% Inputs:
%   t: 1x1 double inidicating the timestep of ode integration
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
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs: 
%   dy: 1xN array describing the outputs from the dynamic model statespace derivatives

%%

% Unpack the parameters structure
guess_Dx = params.dim.l;
g = params.dim.g;

% Calculate the dynamic parameters based on the scenario
[x_coeff, ~, theta_coeff] = exactsoln_dynamicsselectorfunction(y,cond,'AB',guess_Dx,params);

% Unpack the dynamic parameters
Hx = x_coeff.Hx;
Jx = x_coeff.Jx;
TAx = x_coeff.TAx;
TBx = x_coeff.TBx;
HthetaA = theta_coeff.HthetaA;
JthetaA = theta_coeff.JthetaA;
TAthetaA = theta_coeff.TAthetaA;
TBthetaA = theta_coeff.TBthetaA;
HthetaB = theta_coeff.HthetaB;
JthetaB = theta_coeff.JthetaB;
TAthetaB = theta_coeff.TAthetaB;
TBthetaB = theta_coeff.TBthetaB;

% Calculate the dynamic model states
% states: x, xdot, thetaA, thetaAdot, thetaB, thetaBdot
dy = zeros(6,1); 

dy(1) = y(2);

dy(2) = -(Hx*(y(2)^2) + Jx*g + TAx*tauA + TBx*tauB);

dy(3) = y(4);

dy(4) = HthetaA*(y(2)^2) + JthetaA*g + TAthetaA*tauA + TBthetaA*tauB;

dy(5) = y(6);

dy(6) = HthetaB*(y(2)^2) + JthetaB*g + TAthetaB*tauA + TBthetaB*tauB;