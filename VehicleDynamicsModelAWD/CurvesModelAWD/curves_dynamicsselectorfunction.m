function [x_coeff, FN_coeff, theta_coeff] = curves_dynamicsselectorfunction(y,cond,wheel,params,curves)
%% curves_dynamicsselectorfunction
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
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
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
        
        Hx = curves.noslip.x_coeff.Hx(y(1));
        Jx = curves.noslip.x_coeff.Jx(y(1));
        TAx = curves.noslip.x_coeff.TAx(y(1));
        TBx = curves.noslip.x_coeff.TBx(y(1));
        if strcmp(wheel,'A')
            HN = curves.noslip.FNA_coeff.HNA(y(1));
            JN = curves.noslip.FNA_coeff.JNA(y(1));
            TAN = curves.noslip.FNA_coeff.TANA(y(1));
            TBN = curves.noslip.FNA_coeff.TBNA(y(1));
            Htheta = curves.noslip.thetaA_coeff.HthetaA(y(1));
            Jtheta = curves.noslip.thetaA_coeff.JthetaA(y(1));
            TAtheta = curves.noslip.thetaA_coeff.TAthetaA(y(1));
            TBtheta = curves.noslip.thetaA_coeff.TBthetaA(y(1));
        elseif strcmp(wheel,'B')
            HN = curves.noslip.FNB_coeff.HNB(y(1));
            JN = curves.noslip.FNB_coeff.JNB(y(1));
            TAN = curves.noslip.FNB_coeff.TANB(y(1));
            TBN = curves.noslip.FNB_coeff.TBNB(y(1));
            Htheta = curves.noslip.thetaB_coeff.HthetaB(y(1));
            Jtheta = curves.noslip.thetaB_coeff.JthetaB(y(1));
            TAtheta = curves.noslip.thetaB_coeff.TAthetaB(y(1));
            TBtheta = curves.noslip.thetaB_coeff.TBthetaB(y(1));
        elseif strcmp(wheel,'AB')
            HNA = curves.noslip.FNA_coeff.HNA(y(1));
            JNA = curves.noslip.FNA_coeff.JNA(y(1));
            TANA = curves.noslip.FNA_coeff.TANA(y(1));
            TBNA = curves.noslip.FNA_coeff.TBNA(y(1));
            HthetaA = curves.noslip.thetaA_coeff.HthetaA(y(1));
            JthetaA = curves.noslip.thetaA_coeff.JthetaA(y(1));
            TAthetaA = curves.noslip.thetaA_coeff.TAthetaA(y(1));
            TBthetaA = curves.noslip.thetaA_coeff.TBthetaA(y(1));
            HNB = curves.noslip.FNB_coeff.HNB(y(1));
            JNB = curves.noslip.FNB_coeff.JNB(y(1));
            TANB = curves.noslip.FNB_coeff.TANB(y(1));
            TBNB = curves.noslip.FNB_coeff.TBNB(y(1));
            HthetaB = curves.noslip.thetaB_coeff.HthetaB(y(1));
            JthetaB = curves.noslip.thetaB_coeff.JthetaB(y(1));
            TAthetaB = curves.noslip.thetaB_coeff.TAthetaB(y(1));
            TBthetaB = curves.noslip.thetaB_coeff.TBthetaB(y(1));
        end
        
    case 2 % allslip
        
        % calculate relative velocity, get muA, muB
        v_relA = v_relA_fcn(y,params);
        v_relB = v_relB_fcn_curves(y,params,curves);
        muA = mu_fcn(v_relA,mu_s,model);
        muB = mu_fcn(v_relB,mu_s,model);
        Hx = curves.allslip.x_coeff.Hx(y(1),muA,muB);
        Jx = curves.allslip.x_coeff.Jx(y(1),muA,muB);
        TAx = curves.allslip.x_coeff.TAx(y(1),muA,muB);
        TBx = curves.allslip.x_coeff.TBx(y(1),muA,muB);
        if strcmp(wheel,'A')
            HN = curves.allslip.FNA_coeff.HNA(y(1),muA,muB);
            JN = curves.allslip.FNA_coeff.JNA(y(1),muA,muB);
            TAN = curves.allslip.FNA_coeff.TANA(y(1),muA,muB);
            TBN = curves.allslip.FNA_coeff.TBNA(y(1),muA,muB);
            Htheta = curves.allslip.thetaA_coeff.HthetaA(y(1),muA,muB);
            Jtheta = curves.allslip.thetaA_coeff.JthetaA(y(1),muA,muB);
            TAtheta = curves.allslip.thetaA_coeff.TAthetaA(y(1),muA,muB);
            TBtheta = curves.allslip.thetaA_coeff.TBthetaA(y(1),muA,muB);
        elseif strcmp(wheel,'B')
            HN = curves.allslip.FNB_coeff.HNB(y(1),muA,muB);
            JN = curves.allslip.FNB_coeff.JNB(y(1),muA,muB);
            TAN = curves.allslip.FNB_coeff.TANB(y(1),muA,muB);
            TBN = curves.allslip.FNB_coeff.TBNB(y(1),muA,muB);
            Htheta = curves.allslip.thetaB_coeff.HthetaB(y(1),muA,muB);
            Jtheta = curves.allslip.thetaB_coeff.JthetaB(y(1),muA,muB);
            TAtheta = curves.allslip.thetaB_coeff.TAthetaB(y(1),muA,muB);
            TBtheta = curves.allslip.thetaB_coeff.TBthetaB(y(1),muA,muB);
        elseif strcmp(wheel,'AB')
            HNA = curves.allslip.FNA_coeff.HNA(y(1),muA,muB);
            JNA = curves.allslip.FNA_coeff.JNA(y(1),muA,muB);
            TANA = curves.allslip.FNA_coeff.TANA(y(1),muA,muB);
            TBNA = curves.allslip.FNA_coeff.TBNA(y(1),muA,muB);
            HthetaA = curves.allslip.thetaA_coeff.HthetaA(y(1),muA,muB);
            JthetaA = curves.allslip.thetaA_coeff.JthetaA(y(1),muA,muB);
            TAthetaA = curves.allslip.thetaA_coeff.TAthetaA(y(1),muA,muB);
            TBthetaA = curves.allslip.thetaA_coeff.TBthetaA(y(1),muA,muB);
            HNB = curves.allslip.FNB_coeff.HNB(y(1),muA,muB);
            JNB = curves.allslip.FNB_coeff.JNB(y(1),muA,muB);
            TANB = curves.allslip.FNB_coeff.TANB(y(1),muA,muB);
            TBNB = curves.allslip.FNB_coeff.TBNB(y(1),muA,muB);
            HthetaB = curves.allslip.thetaB_coeff.HthetaB(y(1),muA,muB);
            JthetaB = curves.allslip.thetaB_coeff.JthetaB(y(1),muA,muB);
            TAthetaB = curves.allslip.thetaB_coeff.TAthetaB(y(1),muA,muB);
            TBthetaB = curves.allslip.thetaB_coeff.TBthetaB(y(1),muA,muB);
        end
        
    case 3 % rearslipfrontstick
        
        % calculate relative velocity at wheel A, get muA
        v_relA = v_relA_fcn(y,params);
        muA = mu_fcn(v_relA,mu_s,model);
        Hx = curves.rearslipfrontstick.x_coeff.Hx(y(1),muA);
        Jx = curves.rearslipfrontstick.x_coeff.Jx(y(1),muA);
        TAx = curves.rearslipfrontstick.x_coeff.TAx(y(1),muA);
        TBx = curves.rearslipfrontstick.x_coeff.TBx(y(1),muA);
        if strcmp(wheel,'A')
            HN = curves.rearslipfrontstick.FNA_coeff.HNA(y(1),muA);
            JN = curves.rearslipfrontstick.FNA_coeff.JNA(y(1),muA);
            TAN = curves.rearslipfrontstick.FNA_coeff.TANA(y(1),muA);
            TBN = curves.rearslipfrontstick.FNA_coeff.TBNA(y(1),muA);
            Htheta = curves.rearslipfrontstick.thetaA_coeff.HthetaA(y(1),muA);
            Jtheta = curves.rearslipfrontstick.thetaA_coeff.JthetaA(y(1),muA);
            TAtheta = curves.rearslipfrontstick.thetaA_coeff.TAthetaA(y(1),muA);
            TBtheta = curves.rearslipfrontstick.thetaA_coeff.TBthetaA(y(1),muA);
        elseif strcmp(wheel,'B')
            HN = curves.rearslipfrontstick.FNB_coeff.HNB(y(1),muA);
            JN = curves.rearslipfrontstick.FNB_coeff.JNB(y(1),muA);
            TAN = curves.rearslipfrontstick.FNB_coeff.TANB(y(1),muA);
            TBN = curves.rearslipfrontstick.FNB_coeff.TBNB(y(1),muA);
            Htheta = curves.rearslipfrontstick.thetaB_coeff.HthetaB(y(1),muA);
            Jtheta = curves.rearslipfrontstick.thetaB_coeff.JthetaB(y(1),muA);
            TAtheta = curves.rearslipfrontstick.thetaB_coeff.TAthetaB(y(1),muA);
            TBtheta = curves.rearslipfrontstick.thetaB_coeff.TBthetaB(y(1),muA);
        elseif strcmp(wheel,'AB')
            HNA = curves.rearslipfrontstick.FNA_coeff.HNA(y(1),muA);
            JNA = curves.rearslipfrontstick.FNA_coeff.JNA(y(1),muA);
            TANA = curves.rearslipfrontstick.FNA_coeff.TANA(y(1),muA);
            TBNA = curves.rearslipfrontstick.FNA_coeff.TBNA(y(1),muA);
            HthetaA = curves.rearslipfrontstick.thetaA_coeff.HthetaA(y(1),muA);
            JthetaA = curves.rearslipfrontstick.thetaA_coeff.JthetaA(y(1),muA);
            TAthetaA = curves.rearslipfrontstick.thetaA_coeff.TAthetaA(y(1),muA);
            TBthetaA = curves.rearslipfrontstick.thetaA_coeff.TBthetaA(y(1),muA);
            HNB = curves.rearslipfrontstick.FNB_coeff.HNB(y(1),muA);
            JNB = curves.rearslipfrontstick.FNB_coeff.JNB(y(1),muA);
            TANB = curves.rearslipfrontstick.FNB_coeff.TANB(y(1),muA);
            TBNB = curves.rearslipfrontstick.FNB_coeff.TBNB(y(1),muA);
            HthetaB = curves.rearslipfrontstick.thetaB_coeff.HthetaB(y(1),muA);
            JthetaB = curves.rearslipfrontstick.thetaB_coeff.JthetaB(y(1),muA);
            TAthetaB = curves.rearslipfrontstick.thetaB_coeff.TAthetaB(y(1),muA);
            TBthetaB = curves.rearslipfrontstick.thetaB_coeff.TBthetaB(y(1),muA);
        end
        
    case 4 % rearstickfrontslip
        
        % calculate relative velocity at wheel B, get muB
        v_relB = v_relB_fcn_curves(y,params,curves);
        muB = mu_fcn(v_relB,mu_s,model);
        Hx = curves.rearstickfrontslip.x_coeff.Hx(y(1),muB);
        Jx = curves.rearstickfrontslip.x_coeff.Jx(y(1),muB);
        TAx = curves.rearstickfrontslip.x_coeff.TAx(y(1),muB);
        TBx = curves.rearstickfrontslip.x_coeff.TBx(y(1),muB);
        if strcmp(wheel,'A')
            HN = curves.rearstickfrontslip.FNA_coeff.HNA(y(1),muB);
            JN = curves.rearstickfrontslip.FNA_coeff.JNA(y(1),muB);
            TAN = curves.rearstickfrontslip.FNA_coeff.TANA(y(1),muB);
            TBN = curves.rearstickfrontslip.FNA_coeff.TBNA(y(1),muB);
            Htheta = curves.rearstickfrontslip.thetaA_coeff.HthetaA(y(1),muB);
            Jtheta = curves.rearstickfrontslip.thetaA_coeff.JthetaA(y(1),muB);
            TAtheta = curves.rearstickfrontslip.thetaA_coeff.TAthetaA(y(1),muB);
            TBtheta = curves.rearstickfrontslip.thetaA_coeff.TBthetaA(y(1),muB);
        elseif strcmp(wheel,'B')
            HN = curves.rearstickfrontslip.FNB_coeff.HNB(y(1),muB);
            JN = curves.rearstickfrontslip.FNB_coeff.JNB(y(1),muB);
            TAN = curves.rearstickfrontslip.FNB_coeff.TANB(y(1),muB);
            TBN = curves.rearstickfrontslip.FNB_coeff.TBNB(y(1),muB);
            Htheta = curves.rearstickfrontslip.thetaB_coeff.HthetaB(y(1),muB);
            Jtheta = curves.rearstickfrontslip.thetaB_coeff.JthetaB(y(1),muB);
            TAtheta = curves.rearstickfrontslip.thetaB_coeff.TAthetaB(y(1),muB);
            TBtheta = curves.rearstickfrontslip.thetaB_coeff.TBthetaB(y(1),muB);
        elseif strcmp(wheel,'AB')
            HNA = curves.rearstickfrontslip.FNA_coeff.HNA(y(1),muB);
            JNA = curves.rearstickfrontslip.FNA_coeff.JNA(y(1),muB);
            TANA = curves.rearstickfrontslip.FNA_coeff.TANA(y(1),muB);
            TBNA = curves.rearstickfrontslip.FNA_coeff.TBNA(y(1),muB);
            HthetaA = curves.rearstickfrontslip.thetaA_coeff.HthetaA(y(1),muB);
            JthetaA = curves.rearstickfrontslip.thetaA_coeff.JthetaA(y(1),muB);
            TAthetaA = curves.rearstickfrontslip.thetaA_coeff.TAthetaA(y(1),muB);
            TBthetaA = curves.rearstickfrontslip.thetaA_coeff.TBthetaA(y(1),muB);
            HNB = curves.rearstickfrontslip.FNB_coeff.HNB(y(1),muB);
            JNB = curves.rearstickfrontslip.FNB_coeff.JNB(y(1),muB);
            TANB = curves.rearstickfrontslip.FNB_coeff.TANB(y(1),muB);
            TBNB = curves.rearstickfrontslip.FNB_coeff.TBNB(y(1),muB);
            HthetaB = curves.rearstickfrontslip.thetaB_coeff.HthetaB(y(1),muB);
            JthetaB = curves.rearstickfrontslip.thetaB_coeff.JthetaB(y(1),muB);
            TAthetaB = curves.rearstickfrontslip.thetaB_coeff.TAthetaB(y(1),muB);
            TBthetaB = curves.rearstickfrontslip.thetaB_coeff.TBthetaB(y(1),muB);
        end
        
end

x_coeff.Hx = Hx;
x_coeff.Jx = Jx;
x_coeff.TAx = TAx;
x_coeff.TBx = TBx;

% Decide which wheel (or both) to output normal force and angular
% acceleration parameters
if strcmp(wheel,'A') || strcmp(wheel,'B')
    FN_coeff.HN = HN;
    FN_coeff.JN = JN;
    FN_coeff.TAN = TAN;
    FN_coeff.TBN = TBN;
    
    theta_coeff.Htheta = Htheta;
    theta_coeff.Jtheta = Jtheta;
    theta_coeff.TAtheta = TAtheta;
    theta_coeff.TBtheta = TBtheta;
elseif strcmp(wheel,'AB')
    FN_coeff.HNA = HNA;
    FN_coeff.JNA = JNA;
    FN_coeff.TANA = TANA;
    FN_coeff.TBNA = TBNA;
    FN_coeff.HNB = HNB;
    FN_coeff.JNB = JNB;
    FN_coeff.TANB = TANB;
    FN_coeff.TBNB = TBNB;
    
    theta_coeff.HthetaA = HthetaA;
    theta_coeff.JthetaA = JthetaA;
    theta_coeff.TAthetaA = TAthetaA;
    theta_coeff.TBthetaA = TBthetaA;
    theta_coeff.HthetaB = HthetaB;
    theta_coeff.JthetaB = JthetaB;
    theta_coeff.TAthetaB = TAthetaB;
    theta_coeff.TBthetaB = TBthetaB;
end
