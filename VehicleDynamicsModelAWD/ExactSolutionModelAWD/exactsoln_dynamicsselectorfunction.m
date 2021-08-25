function [x_coeff, FN_coeff, theta_coeff] = exactsoln_dynamicsselectorfunction(y,cond,wheel,guess_Dx,params)
%% exactsoln_dynamicsselectorfunction
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate dynamic model parameters based on the 
% 'cond' value (designating 1 of 4 possible scenarios: noslip, allslip, 
% rearslipfrontstick, and rearstickfrontslip).
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
%       Outputs are elements necessary to solve the vehicle dynamics model: 
%                   x'' + Hx*x'^2 + Jx*g + TAx*tau_A + TBx*tau_B = 0
%       And solve the angular acceleration at both wheels as well: 
%                   thetaA'' = HthetaA*x'^2 + JthetaA*g + TAthetaA*tau_A + TBthetaA*tau_B
%                   thetaB'' = HthetaB*x'^2 + JthetaB*g + TAthetaB*tau_A + TBthetaB*tau_B
%       And solve the normal force at both wheels as well: 
%                   FNA = HNA*x'^2 + JNA*g + TANA*tau_A + TBNA*tau_B
%                   FNB = HNB*x'^2 + JNB*g + TANB*tau_A + TBNB*tau_B
%   x_coeff: structure containing dynamic parameters pertaining to x
%   theta_coeff: structure containing dynamic parameters pertaining to angular acceleration
%   FN_coeff: structure containing dynamic parameters pertaining to normal force

%%

% Unpack the parameters structure
model = params.settings.friction_model;
mu_s = params.dim.mu_s;

% Calculate the dynamic parameters based on the scenario indicated by 'cond'
switch cond
    
    case 1 % noslip
        
        [x_coeff, thetaA_coeff, thetaB_coeff, FNA_coeff, FNB_coeff, ~, ~] = noslip_eom(y(1),guess_Dx,params);
        
    case 2 % allslip
        
        % calculate relative velocity at both wheels, get muA, muB
        v_relA = v_relA_fcn(y,params);
        v_relB = v_relB_fcn(y,guess_Dx,params);
        muA = mu_fcn(v_relA,mu_s,model);
        muB = mu_fcn(v_relB,mu_s,model);
        
        [x_coeff, thetaA_coeff, thetaB_coeff, FNA_coeff, FNB_coeff, ~, ~] = allslip_eom(y(1),guess_Dx,muA,muB,params);
        
    case 3 % rearslipfrontstick
        
        % calculate relative velocity at wheel A, get muA
        v_relA = v_relA_fcn(y,params);
        muA = mu_fcn(v_relA,mu_s,model);
        
        [x_coeff, thetaA_coeff, thetaB_coeff, FNA_coeff, FNB_coeff, ~, ~] = rearslipfrontstick_eom(y(1),guess_Dx,muA,params);
        
    case 4 % rearstickfrontslip
        
        % calculate relative velocity at wheel B, get muB
        v_relB = v_relB_fcn(y,guess_Dx,params);
        muB = mu_fcn(v_relB,mu_s,model);
        
        [x_coeff, thetaA_coeff, thetaB_coeff, FNA_coeff, FNB_coeff, ~, ~] = rearstickfrontslip_eom(y(1),guess_Dx,muB,params);
        
end

% Decide which wheel (or both) to output normal force and angular
% acceleration parameters
if strcmp(wheel,'A') 
    
    FN_coeff.HN = FNA_coeff.HNA;
    FN_coeff.JN = FNA_coeff.JNA;
    FN_coeff.TAN = FNA_coeff.TANA;
    FN_coeff.TBN = FNA_coeff.TBNA;
    theta_coeff.Htheta = thetaA_coeff.HthetaA;
    theta_coeff.Jtheta = thetaA_coeff.JthetaA;
    theta_coeff.TAtheta = thetaA_coeff.TAthetaA;
    theta_coeff.TBtheta = thetaA_coeff.TBthetaA;
    
elseif strcmp(wheel,'B')
    
    FN_coeff.HN = FNB_coeff.HNB;
    FN_coeff.JN = FNB_coeff.JNB;
    FN_coeff.TAN = FNB_coeff.TANB;
    FN_coeff.TBN = FNB_coeff.TBNB;
    theta_coeff.Htheta = thetaB_coeff.HthetaB;
    theta_coeff.Jtheta = thetaB_coeff.JthetaB;
    theta_coeff.TAtheta = thetaB_coeff.TAthetaB;
    theta_coeff.TBtheta = thetaB_coeff.TBthetaB;
    
elseif strcmp(wheel,'AB')
    
    FN_coeff.HNA = FNA_coeff.HNA;
    FN_coeff.JNA = FNA_coeff.JNA;
    FN_coeff.TANA = FNA_coeff.TANA;
    FN_coeff.TBNA = FNA_coeff.TBNA;
    FN_coeff.HNB = FNB_coeff.HNB;
    FN_coeff.JNB = FNB_coeff.JNB;
    FN_coeff.TANB = FNB_coeff.TANB;
    FN_coeff.TBNB = FNB_coeff.TBNB;
    
    theta_coeff.HthetaA = thetaA_coeff.HthetaA;
    theta_coeff.JthetaA = thetaA_coeff.JthetaA;
    theta_coeff.TAthetaA = thetaA_coeff.TAthetaA;
    theta_coeff.TBthetaA = thetaA_coeff.TBthetaA;
    theta_coeff.HthetaB = thetaB_coeff.HthetaB;
    theta_coeff.JthetaB = thetaB_coeff.JthetaB;
    theta_coeff.TAthetaB = thetaB_coeff.TAthetaB;
    theta_coeff.TBthetaB = thetaB_coeff.TBthetaB;
    
end
