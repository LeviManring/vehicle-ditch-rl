function [x_coeff, thetaA_coeff, thetaB_coeff, FNA_coeff, FNB_coeff, tM, dx] = rearstickfrontslip_eom(xA,guess_Dx,mu,params)
%% rearstickfrontslip_eom
% Levi Manring, Duke University
% 2021
%
% This function calculates certain parameters needed for integration of the
% rearstickfrontslip dynamics model for the AWD complete model.
%
% Inputs:
%   xA: 1x1 double indicating the x-position of wheel A
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%       Outputs are elements necessary to solve the rearstickfrontslip vehicle dynamics
%       model: 
%                   x'' + Hx*x'^2 + Jx*g + TAx*tau_A + TBx*tau_B = 0
%       And solve the angular acceleration at both wheels as well: 
%                   thetaA'' = HthetaA*x'^2 + JthetaA*g + TAthetaA*tau_A + TBthetaA*tau_B
%                   thetaB'' = HthetaB*x'^2 + JthetaB*g + TAthetaB*tau_A + TBthetaB*tau_B
%       And solve the normal force at both wheels as well: 
%                   FNA = HNA*x'^2 + JNA*g + TANA*tau_A + TBNA*tau_B
%                   FNB = HNB*x'^2 + JNB*g + TANB*tau_A + TBNB*tau_B
%   x_coeff: structure containing dynamic parameters pertaining to x
%   thetaA_coeff: structure containing dynamic parameters pertaining to thetaA
%   thetaB_coeff: structure containing dynamic parameters pertaining to thetaB
%   FNA_coeff: structure containing dynamic parameters pertaining to normal force at wheel A
%   FNB_coeff: structure containing dynamic parameters pertaining to normal force at wheel B
%   tM: 1x1 double, the angle of the vehicle with respect to the horizontal
%   dx: 1x1 double, the accurately solved-for guess_Dx 

%%

% unpack the parameters structure
R = params.dim.R;
l = params.dim.l;
mM = params.dim.m;
mA = params.dim.mA;
mB = params.dim.mB;
I_B = params.dim.I_B;
I_A = params.dim.I_A;
I_M = params.dim.I_M;
xc = params.dim.xc;
yc = params.dim.yc;

% grab the function derivatives and features of the ditch profile at wheel A
[~, yAx, yAxx, yAxxx, flag] = y_eval_fcn(xA,params);

% determine the parameters for taking numerical derivative
h = R/20;
Dvec = [xA-2*h, xA-h, xA, xA+h, xA+2*h]';

% Numerically solve for Dx and tM
Dx_solve = zeros(length(Dvec),1);
f_k = zeros(length(Dvec),1);

for k = 1:length(Dvec)
    x_k = Dvec(k,1);
    
    if flag(1) == 0
        Dx_solve(k,1) = l;
        f_k(k,1)      = 0;
    end
    
    if flag(2) == 0 && flag(1) ~= 0
        Dx_solve(k,1) = fzero(@(Dx) sqrt(Dx^2 + (params.fun.y{1}(x_k + Dx) - params.fun.y{1}(x_k))^2) - l, guess_Dx);
        guess_Dx = Dx_solve(k,1);
        
        f_k(k,1) = atan2((params.fun.y{1}(x_k + Dx_solve(k,1)) - params.fun.y{1}(x_k)),Dx_solve(k,1));
    end
    
    if flag(2) == 1
        
        Dx_solve(k,1) = fzero(@(Dx) sqrt((Dx + R*((params.fun.y{2}(x_k)/sqrt(1 + params.fun.y{2}(x_k)^2)) - ...
            (params.fun.y{2}(x_k + Dx)/sqrt(1 + params.fun.y{2}(x_k + Dx)^2))))^2 + ...
            (params.fun.y{1}(x_k + Dx) - params.fun.y{1}(x_k) - R*((1/sqrt(1 + params.fun.y{2}(x_k)^2)) - ...
            (1/sqrt(1 + params.fun.y{2}(x_k + Dx)^2))))^2) - l, guess_Dx);
        guess_Dx = Dx_solve(k,1);
        
        f_k(k,1) = atan2((params.fun.y{1}(x_k + Dx_solve(k,1)) - params.fun.y{1}(x_k) - ...
            R*((1/sqrt(1 + params.fun.y{2}(x_k)^2)) - (1/sqrt(1 + params.fun.y{2}(x_k + Dx_solve(k,1))^2)))),...
            (Dx_solve(k,1) + R*((params.fun.y{2}(x_k)/sqrt(1 + params.fun.y{2}(x_k)^2)) - ...
            (params.fun.y{2}(x_k + Dx_solve(k,1))/sqrt(1 + params.fun.y{2}(x_k + Dx_solve(k,1))^2)))));
    end
    
end

% Calculate numerical derivatives of tM
k = 3;
tM = f_k(k,1);
tMx = (-f_k(k+2,1) + 8*f_k(k+1,1) - 8*f_k(k-1,1) + f_k(k-2,1))/(12*h);
tMxx = (-f_k(k+2,1) + 16*f_k(k+1,1) - 30*f_k(k,1) + 16*f_k(k-1,1) - f_k(k-2,1))/(12*h^2);

dx = Dx_solve(k,1);
xB = xA + dx; % determine where wheel B is located

% grab the function derivatives and features of the ditch profile at wheel B
[~, yBx, ~, ~, ~] = y_eval_fcn(xB,params);

% Calculate certain dynamic model parameters
PA = 1-R*yAxx/(yAx^2+1)^(3/2);
PAx = -R*(yAx^2*yAxxx-3*yAx*yAxx^2+yAxxx)/(yAx^2+1)^(5/2);

% PB = 1-R*yBxx/(yBx^2+1)^(3/2);
% PBx = -R*(yBx^2*yBxxx-3*yBx*yBxx^2+yBxxx)/(yBx^2+1)^(5/2);

phiA = sqrt(yAx^2+1)/R;
phiAx = yAx*yAxx/(R*sqrt(yAx^2+1));

phiB = sqrt(yBx^2+1)/R;
% phiBx = yBx*yBxx/(R*sqrt(yBx^2+1));


%% Calculating rearstickfrontslip dynamic model parameters
% Applying torque at wheel A and B and including moment of inertia effects for both wheels

% Coefficients for x

Hx = (l*(-mB*(((mu*tMxx+tMx^2)*yBx-tMx^2*mu+tMxx)*yAx+(-mu*tMx^2+tMxx)*yBx-tMx^2-tMxx*mu)*l+((((tMx^2*yc-tMxx*xc)*mu-tMx^2*xc-yc*tMxx)*yBx+(tMx^2*xc+tMxx*yc)*mu+yc*tMx^2-tMxx*xc)*yAx+((tMx^2*xc+tMxx*yc)*mu+yc*tMx^2-tMxx*xc)*yBx+(-tMx^2*yc+tMxx*xc)*mu+tMx^2*xc+yc*tMxx)*mM)*R^2*phiA*cos(tM)^2+l*(R^2*phiA*(mB*(((mu*tMx^2-tMxx)*yBx+tMx^2+tMxx*mu)*yAx+(mu*tMxx+tMx^2)*yBx-tMx^2*mu+tMxx)*l+((((tMx^2*xc+tMxx*yc)*mu+yc*tMx^2-tMxx*xc)*yBx+(-tMx^2*yc+tMxx*xc)*mu+tMx^2*xc+yc*tMxx)*yAx+((-tMx^2*yc+tMxx*xc)*mu+tMx^2*xc+yc*tMxx)*yBx+(-tMx^2*xc-tMxx*yc)*mu-yc*tMx^2+tMxx*xc)*mM)*sin(tM)-(PAx*(mA+mM)*(mu*yBx+1)*yAx^2+(yAxx*PA*(mu*yBx+1)*mM+(PA*mA*mu*yAxx+PAx*mB)*yBx-PAx*mB*mu+yAxx*PA*mA)*yAx+PAx*(mu*yBx+1)*mM+(PAx*(mA+mB)*mu+yAxx*PA*mB)*yBx-yAxx*PA*mB*mu+PAx*(mA+mB))*R^2*phiA-(yAx^2+1)*(mu*yBx+1)*I_A*phiAx)*cos(tM)-l*(R^2*(PAx*(mA+mB+mM)*(yBx-mu)*yAx^2+(yAxx*PA*(yBx-mu)*mM+(PAx*mB*mu+yAxx*PA*(mA+mB))*yBx-yAxx*PA*(mA+mB)*mu+PAx*mB)*yAx+PAx*(mA+mM)*(yBx-mu))*phiA+I_A*phiAx*(yAx^2+1)*(yBx-mu))*sin(tM)+(((mu*tMxx+tMx^2)*yBx-tMx^2*mu+tMxx)*mB*yAx*l^2-((-tMx^2*xc-tMxx*yc)*yAx+yc*tMx^2-tMxx*xc)*mM*(yBx-mu)*l+tMxx*(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*phiA)/((-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-PA*((mA+mM)*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB+mM)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((mA+mB+mM)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+(mA+mM)*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx)*phiA);

Jx = -R^2*((((mu*yBx+1)*(mA+mM)*yAx+mB*(yBx-mu))*l-xc*((mu*yBx+1)*yAx-yBx+mu)*mM)*cos(tM)+(yAx*(mA+mB+mM)*(yBx-mu)*l+yc*((mu*yBx+1)*yAx-yBx+mu)*mM)*sin(tM))/(-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-PA*((mA+mM)*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB+mM)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((mA+mB+mM)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+(mA+mM)*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx);

TAx = ((yBx-mu)*sin(tM)+(mu*yBx+1)*cos(tM))*(yAx^2+1)*l/(phiA*(-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-PA*((mA+mM)*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB+mM)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((mA+mB+mM)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+(mA+mM)*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx));

TBx = 0;

x_coeff.Hx = Hx;
x_coeff.Jx = Jx;
x_coeff.TAx = TAx;
x_coeff.TBx = TBx;

% Coefficients for thetaA

HthetaA = R^2*(-l*(-mB*((((phiA*tMxx-phiAx*tMx)*mu+tMx^2*phiA)*yBx-tMx^2*phiA*mu-tMx*phiAx+phiA*tMxx)*yAx+(-mu*phiA*tMx^2+phiA*tMxx-phiAx*tMx)*yBx+(-phiA*tMxx+phiAx*tMx)*mu-tMx^2*phiA)*l+((((phiA*tMx^2*yc-phiA*tMxx*xc+phiAx*tMx*xc)*mu-tMx^2*phiA*xc+yc*tMx*phiAx-yc*phiA*tMxx)*yBx+(phiA*tMx^2*xc+phiA*tMxx*yc-phiAx*tMx*yc)*mu+yc*tMx^2*phiA+tMx*phiAx*xc-phiA*tMxx*xc)*yAx+((phiA*tMx^2*xc+phiA*tMxx*yc-phiAx*tMx*yc)*mu+yc*tMx^2*phiA+tMx*phiAx*xc-phiA*tMxx*xc)*yBx+(-phiA*tMx^2*yc+phiA*tMxx*xc-phiAx*tMx*xc)*mu+tMx^2*phiA*xc-yc*tMx*phiAx+yc*phiA*tMxx)*mM)*cos(tM)^2-l*((mB*(((mu*phiA*tMx^2-phiA*tMxx+phiAx*tMx)*yBx+(phiA*tMxx-phiAx*tMx)*mu+tMx^2*phiA)*yAx+((phiA*tMxx-phiAx*tMx)*mu+tMx^2*phiA)*yBx-tMx^2*phiA*mu-tMx*phiAx+phiA*tMxx)*l+((((phiA*tMx^2*xc+phiA*tMxx*yc-phiAx*tMx*yc)*mu+yc*tMx^2*phiA+tMx*phiAx*xc-phiA*tMxx*xc)*yBx+(-phiA*tMx^2*yc+phiA*tMxx*xc-phiAx*tMx*xc)*mu+tMx^2*phiA*xc-yc*tMx*phiAx+yc*phiA*tMxx)*yAx+((-phiA*tMx^2*yc+phiA*tMxx*xc-phiAx*tMx*xc)*mu+tMx^2*phiA*xc-yc*tMx*phiAx+yc*phiA*tMxx)*yBx+(-phiA*tMx^2*xc-phiA*tMxx*yc+phiAx*tMx*yc)*mu-yc*tMx^2*phiA-tMx*phiAx*xc+phiA*tMxx*xc)*mM)*sin(tM)-(mA+mM)*(mu*yBx+1)*(-PA*phiAx+PAx*phiA)*yAx^2+(-phiA*yAxx*PA*(mu*yBx+1)*mM+(-phiA*yAxx*PA*mA*mu-mB*(-PA*phiAx+PAx*phiA))*yBx+mB*(-PA*phiAx+PAx*phiA)*mu-phiA*yAxx*PA*mA)*yAx-(mu*yBx+1)*(-PA*phiAx+PAx*phiA)*mM+(-(mA+mB)*(-PA*phiAx+PAx*phiA)*mu-phiA*yAxx*PA*mB)*yBx+phiA*yAxx*PA*mB*mu-(mA+mB)*(-PA*phiAx+PAx*phiA))*cos(tM)+l*((mA+mB+mM)*(yBx-mu)*(-PA*phiAx+PAx*phiA)*yAx^2+(phiA*yAxx*PA*(yBx-mu)*mM+(mB*(-PA*phiAx+PAx*phiA)*mu+phiA*yAxx*PA*(mA+mB))*yBx-phiA*yAxx*PA*(mA+mB)*mu+mB*(-PA*phiAx+PAx*phiA))*yAx+(mA+mM)*(yBx-mu)*(-PA*phiAx+PAx*phiA))*sin(tM)-(((phiA*tMxx-phiAx*tMx)*mu+tMx^2*phiA)*yBx-tMx^2*phiA*mu-tMx*phiAx+phiA*tMxx)*mB*yAx*l^2+mM*(yBx-mu)*((-phiA*tMx^2*xc-phiA*tMxx*yc+phiAx*tMx*yc)*yAx+yc*tMx^2*phiA+tMx*phiAx*xc-phiA*tMxx*xc)*l+(I_M+mM*(xc^2+yc^2))*(-phiA*tMxx+phiAx*tMx)*((mu*yBx+1)*yAx-yBx+mu))/(-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-(mu*yBx+1)*(PA*R^2*mA+PA*R^2*mM+I_A)*yAx^2-R^2*PA*mB*(yBx-mu)*yAx-(mu*yBx+1)*(PA*R^2*mA+PA*R^2*mB+PA*R^2*mM+I_A))*cos(tM)-l*((yBx-mu)*(PA*R^2*mA+PA*R^2*mB+PA*R^2*mM+I_A)*yAx^2+R^2*PA*mB*(mu*yBx+1)*yAx+(yBx-mu)*(PA*R^2*mA+PA*R^2*mM+I_A))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx);

JthetaA = phiA*R^2*((((mu*yBx+1)*(mA+mM)*yAx+mB*(yBx-mu))*l-xc*((mu*yBx+1)*yAx-yBx+mu)*mM)*cos(tM)+(yAx*(mA+mB+mM)*(yBx-mu)*l+yc*((mu*yBx+1)*yAx-yBx+mu)*mM)*sin(tM))/(-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-PA*((mA+mM)*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB+mM)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((mA+mB+mM)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+(mA+mM)*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx);

TAthetaA = -((yBx-mu)*sin(tM)+(mu*yBx+1)*cos(tM))*(yAx^2+1)*l/(-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-PA*((mA+mM)*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB+mM)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((mA+mB+mM)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+(mA+mM)*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx);

TBthetaA = 0;

thetaA_coeff.HthetaA = HthetaA;
thetaA_coeff.JthetaA = JthetaA;
thetaA_coeff.TAthetaA = TAthetaA;
thetaA_coeff.TBthetaA = TBthetaA;

% Coefficients for thetaB

HthetaB = (-2*l*(-(1/2)*(2*tMx^2*yAx*PA+(PA*yAx*yAxx+PAx*yAx^2-PAx)*tMx-yAx^2*tMxx*PA+tMxx*PA)*mB*l+(-(1/2)*PA*(yAx^2*yc+2*xc*yAx-yc)*tMx^2+(-(1/2)*PAx*yAx^2*xc+(-(1/2)*yAxx*PA*xc+yc*PAx)*yAx+(1/2)*yc*yAxx*PA+(1/2)*PAx*xc)*tMx-PA*tMxx*(-(1/2)*yAx^2*xc+yc*yAx+(1/2)*xc))*mM)*mB*R^2*phiA*cos(tM)^2+(-l*mB*((2*((-(1/2)*yAx^2*PA+(1/2)*PA)*tMx^2+(PAx*yAx+(1/2)*PA*yAxx)*tMx-yAx*tMxx*PA))*mB*l+(2*PA*(-(1/2)*yAx^2*xc+yc*yAx+(1/2)*xc)*tMx^2+(yc*PAx*yAx^2+(PA*yAxx*yc+2*PAx*xc)*yAx+yAxx*PA*xc-yc*PAx)*tMx-tMxx*PA*(yAx^2*yc+2*xc*yAx-yc))*mM)*R^2*phiA*sin(tM)+(tMx^3*l^3*mB^2+tMx^3*mB*mM*(yAx*yc+xc)*l^2+((I_M+mM*(xc^2+yc^2))*tMx^3+yAxx*PA^2*(mA+mB+mM))*mB*l+(I_M+mM*(xc^2+yc^2))*tMx^3*(yAx*yc+xc)*mM)*phiA*R^2+l*I_A*((PA*yAxx+PAx*yAx)*phiA-yAx*phiAx*PA)*mB*(yAx^2+1))*cos(tM)+(-(-tMx^3*yAx*l^3*mB^2+tMx^3*mB*mM*(-xc*yAx+yc)*l^2-((I_M+mM*(xc^2+yc^2))*tMx^3+yAxx*PA^2*(mA+mB+mM))*mB*yAx*l+(I_M+mM*(xc^2+yc^2))*tMx^3*(-xc*yAx+yc)*mM)*phiA*R^2-I_A*l*mB*(yAx^2+1)*(-PA*phiAx+PAx*phiA))*sin(tM)-((tMx^2*yAx*PA*mB+((PA*yAx*yAxx+PAx*yAx^2+PAx)*mM+PAx*(mA+mB)*yAx^2+yAxx*PA*(mA+mB)*yAx+mA*PAx)*tMx-((yAx^2+1)*mM+(mA+mB)*yAx^2+mA)*PA*tMxx)*mB*l^2-mB*(PA*(-xc*yAx+yc)*tMx^2+PAx*(yAx*yc+xc)*tMx-tMxx*PA*(yAx*yc+xc))*mM*l+(I_M+mM*(xc^2+yc^2))*(mA+mB+mM)*((PA*yAx*yAxx+PAx*yAx^2+PAx)*tMx-tMxx*PA*(yAx^2+1)))*phiA*R^2-(-phiA*tMxx+phiAx*tMx)*I_A*(yAx^2+1)*(l^2*mB+mM*(xc^2+yc^2)+I_M))*mu*R^2*phiB/(phiA*I_B*(-l*tMx*R^2*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*yBx+1)*xc+yc*(yBx-mu))*yAx+(yBx-mu)*xc+(-mu*yBx-1)*yc)*mM)*cos(tM)^2+l*(tMx*(-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu-yBx)*xc+yc*(mu*yBx+1))*yAx+(mu*yBx+1)*xc+yc*(yBx-mu))*mM)*R^2*sin(tM)-PA*((yAx^2+1)*(mu*yBx+1)*mM+mA*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((yAx^2+1)*(yBx-mu)*mM+(mA+mB)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+mA*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+tMx*(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu+yAx-yBx))*R^2));

JthetaB = -((-yAx*l^2*mB^2+2*mB*mM*(-xc*yAx+yc)*l+mM^2*((-xc^2+yc^2)*yAx+2*yc*xc))*R^2*tMx*cos(tM)^2+(-(-l^2*mB^2-2*mB*mM*(yAx*yc+xc)*l+mM^2*(-2*xc*yAx*yc-xc^2+yc^2))*R^2*tMx*sin(tM)-(l*mB+mM*xc*(yAx^2+1))*PA*(mA+mB+mM)*R^2-I_A*(yAx^2+1)*(l*mB+mM*xc))*cos(tM)+((-yAx*l*mB+yc*mM*(yAx^2+1))*PA*(mA+mB+mM)*R^2+yc*I_A*mM*(yAx^2+1))*sin(tM)+(yAx*mB*(mA+mB+mM)*l^2-yc*l*mB*mM+(xc^2*yAx-xc*yc)*mM^2+yAx*(mB*(xc^2+yc^2)+yc^2*mA+mA*xc^2+I_M)*mM+I_M*yAx*(mA+mB))*R^2*tMx)*mu*R^2*phiB/(I_B*(-l*tMx*R^2*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*mM)*cos(tM)^2+l*(tMx*(-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((yBx*yc+xc)*mu-xc*yBx+yc)*yAx+(xc*yBx-yc)*mu+yc*yBx+xc)*mM)*R^2*sin(tM)-PA*((yAx^2+1)*(mu*yBx+1)*mM+mA*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((yAx^2+1)*(yBx-mu)*mM+(mA+mB)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+mA*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+tMx*(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu+yAx-yBx))*R^2));

TAthetaB = phiB*mu*(yAx^2+1)*R^2*(cos(tM)*PA*l*mB*yAx-sin(tM)*PA*l*mB+(l^2*mB+mM*(xc^2+yc^2)+I_M)*tMx)/(I_B*(-l*tMx*R^2*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*mM)*cos(tM)^2+l*(tMx*(-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((yBx*yc+xc)*mu-xc*yBx+yc)*yAx+(xc*yBx-yc)*mu+yc*yBx+xc)*mM)*R^2*sin(tM)-PA*((mA+mM)*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB+mM)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((mA+mB+mM)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+(mA+mM)*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx)*phiA);

TBthetaB = 1/I_B;

thetaB_coeff.HthetaB = HthetaB;
thetaB_coeff.JthetaB = JthetaB;
thetaB_coeff.TAthetaB = TAthetaB;
thetaB_coeff.TBthetaB = TBthetaB;

% Coefficients for Normal force at Wheel A

HNA = ((((-mB*(PA*((mu*yAx+1)*yBx+yAx-mu)*tMx^2+((PA*yAxx-PAx*mu+PAx*yAx)*yBx+(-PA*yAxx-PAx*yAx)*mu-PAx)*tMx-tMxx*((yAx-mu)*yBx-mu*yAx-1)*PA)*l+(-PA*(((xc*yAx-yc)*mu+yc*yAx+xc)*yBx+(-yAx*yc-xc)*mu+yAx*xc-yc)*tMx^2+(((yc*yAxx*PA+PAx*(yAx*yc+xc))*mu-yAxx*PA*xc+PAx*(-xc*yAx+yc))*yBx+(yAxx*PA*xc-PAx*(-xc*yAx+yc))*mu+yc*yAxx*PA+PAx*(yAx*yc+xc))*tMx-tMxx*(((yAx*yc+xc)*mu-yAx*xc+yc)*yBx+(xc*yAx-yc)*mu+yc*yAx+xc)*PA)*mM)*(mA+mM)*R^2-I_A*(mB*((mu*yAx*yBx-mu+yAx+yBx)*tMx^2-tMxx*((yAx-mu)*yBx-mu*yAx-1))*l+((((xc*yAx-yc)*mu+yc*yAx+xc)*yBx+(-yAx*yc-xc)*mu+yAx*xc-yc)*tMx^2+tMxx*(((yAx*yc+xc)*mu-yAx*xc+yc)*yBx+(xc*yAx-yc)*mu+yc*yAx+xc))*mM))*phiA+I_A*tMx*phiAx*(-mB*((yAx-mu)*yBx-mu*yAx-1)*l+(((yAx*yc+xc)*mu-yAx*xc+yc)*yBx+(xc*yAx-yc)*mu+yc*yAx+xc)*mM))*l*cos(tM)^2+(l*(((mB*(-((yAx-mu)*yBx-mu*yAx-1)*PA*tMx^2+(((PA*yAxx+PAx*yAx)*mu+PAx)*yBx+PAx*yAx-PAx*mu+PA*yAxx)*tMx-tMxx*PA*((mu*yAx+1)*yBx+yAx-mu))*l+((((yAx*yc+xc)*mu-yAx*xc+yc)*yBx+(xc*yAx-yc)*mu+yc*yAx+xc)*PA*tMx^2+(((yAxx*PA*xc-PAx*(-xc*yAx+yc))*mu+yc*yAxx*PA+PAx*(yAx*yc+xc))*yBx+(-yc*yAxx*PA-PAx*(yAx*yc+xc))*mu+yAxx*PA*xc-PAx*(-xc*yAx+yc))*tMx-tMxx*PA*(((xc*yAx-yc)*mu+yc*yAx+xc)*yBx+(-yAx*yc-xc)*mu+yAx*xc-yc))*mM)*(mA+mM)*R^2+I_A*(-mB*(((yAx-mu)*yBx-mu*yAx-1)*tMx^2+tMxx*((mu*yAx+1)*yBx+yAx-mu))*l+((((yAx*yc+xc)*mu-yAx*xc+yc)*yBx+(xc*yAx-yc)*mu+yc*yAx+xc)*tMx^2-tMxx*(((xc*yAx-yc)*mu+yc*yAx+xc)*yBx+(-yAx*yc-xc)*mu+yAx*xc-yc))*mM))*phiA+I_A*tMx*phiAx*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((xc*yAx-yc)*mu+yc*yAx+xc)*yBx+(-yAx*yc-xc)*mu+yAx*xc-yc)*mM))*sin(tM)+((tMx^3*mB*mM*((-mu*xc+yc)*yBx-yc*mu-xc)*l^2+(mu*yBx+1)*(((-xc^2-yc^2)*tMx^3-yAxx*PA^2)*mM^2+(mB*(xc^2+yc^2)*tMx^3-2*PA^2*yAxx*(mA+(1/2)*mB))*mM+I_M*tMx^3*mB-yAxx*PA^2*mA*(mA+mB))*l+(I_M+mM*(xc^2+yc^2))*tMx^3*((mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2+l*I_A*(-yAxx*PA*(mu*yBx+1)*mM+((-PA*mA*yAxx+PAx*mB*yAx)*mu+yAx*mB*(PA*yAxx+PAx*yAx))*yBx-yAx*mB*(PA*yAxx+PAx*yAx)*mu+PAx*yAx*mB-yAxx*PA*mA))*phiA-l*I_A*phiAx*mB*PA*yAx*((yAx+mu)*yBx-mu*yAx+1))*cos(tM)+(((-tMx^3*mB*((mu*yc+xc)*yBx-mu*xc+yc)*mM*l^2+(yBx-mu)*(((-xc^2-yc^2)*tMx^3-yAxx*PA^2)*mM^2+(mB*(xc^2+yc^2)*tMx^3-2*PA^2*yAxx*(mA+(1/2)*mB))*mM+I_M*tMx^3*mB-yAxx*PA^2*mA*(mA+mB))*l-((mu*yc-xc)*yBx+mu*xc+yc)*(I_M+mM*(xc^2+yc^2))*tMx^3*mM)*R^2-l*I_A*(yAxx*PA*(yBx-mu)*mM+(PAx*mB*mu+yAxx*PA*(mA+mB)+PAx*yAx*mB)*yBx+(-yAxx*PA*(mA+mB)-PAx*yAx*mB)*mu+PAx*mB))*phiA+l*I_A*phiAx*mB*PA*((yAx+mu)*yBx-mu*yAx+1))*sin(tM)+((-mB*(-PA*(yBx-mu)*tMx^2+PAx*(mu*yBx+1)*tMx-PA*tMxx*(mu*yBx+1))*(mA+mM)*l^2-((-PA*(yAx*yc+xc)*tMx^2+(-yAxx*PA*xc+PAx*(-xc*yAx+yc))*tMx-PA*tMxx*(-xc*yAx+yc))*(yBx-mu)*mM-((yc*mB*mu+yc*(mA+mB)*yAx+mA*xc)*yBx+(-yc*(mA+mB)*yAx-mA*xc)*mu+yc*mB)*PA*tMx^2+((-PAx*mB*mu*xc-yAxx*xc*(mA+mB)*PA+(-xc*(mA+mB)*yAx+yc*mA)*PAx)*yBx+(yAxx*xc*(mA+mB)*PA-(-xc*(mA+mB)*yAx+yc*mA)*PAx)*mu-PAx*mB*xc)*tMx-tMxx*PA*((-mB*mu*xc-xc*(mA+mB)*yAx+yc*mA)*yBx+(xc*(mA+mB)*yAx-yc*mA)*mu-mB*xc))*mM*l-(((PA*yAxx+PAx*mu+PAx*yAx)*yBx+(-PA*yAxx-PAx*yAx)*mu+PAx)*tMx-tMxx*PA*((yAx+mu)*yBx-mu*yAx+1))*(I_M+mM*(xc^2+yc^2))*(mA+mB+mM))*R^2+I_A*(((yBx-mu)*tMx^2+tMxx*(mu*yBx+1))*mB*l^2+((yAx*yc+xc)*tMx^2+tMxx*(-xc*yAx+yc))*mM*(yBx-mu)*l+tMxx*(I_M+mM*(xc^2+yc^2))*((yAx+mu)*yBx-mu*yAx+1)))*phiA-I_A*tMx*phiAx*(mB*(mu*yBx+1)*l^2+mM*(yBx-mu)*(-xc*yAx+yc)*l+(I_M+mM*(xc^2+yc^2))*((yAx+mu)*yBx-mu*yAx+1)))*R/(-l*tMx*R^2*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((xc*yAx-yc)*mu+yc*yAx+xc)*yBx+(-yAx*yc-xc)*mu+yAx*xc-yc)*mM)*cos(tM)^2+l*(tMx*(-mB*((yAx-mu)*yBx-mu*yAx-1)*l+(((yAx*yc+xc)*mu-yAx*xc+yc)*yBx+(xc*yAx-yc)*mu+yc*yAx+xc)*mM)*R^2*sin(tM)-PA*((yAx^2+1)*(mu*yBx+1)*mM+((mA*yAx^2+mA+mB)*mu+yAx*mB)*yBx+yAx^2*mA-yAx*mB*mu+mA+mB)*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((yAx^2+1)*(yBx-mu)*mM+(yAx*mB*mu+(mA+mB)*yAx^2+mA)*yBx+((-mA-mB)*yAx^2-mA)*mu+yAx*mB)*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+tMx*(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu+yAx-yBx))*R^2);

JNA = -phiA*(R^2*((((-mu*yc+xc)*yBx-mu*xc-yc)*l+(2*mu*xc*yc-xc^2+yc^2)*yBx+(xc^2-yc^2)*mu+2*yc*xc)*mM^2-l*(-mB*(yBx-mu)*l+(mA-mB)*((mu*yc-xc)*yBx+mu*xc+yc))*mM+l^2*mA*mB*(yBx-mu))*tMx*cos(tM)^2+(-((((mu*xc+yc)*yBx-yc*mu+xc)*l+((-xc^2+yc^2)*mu-2*yc*xc)*yBx+2*yc*mu*xc+yc^2-xc^2)*mM^2+l*((mu*yBx+1)*mB*l+(mA-mB)*((mu*xc+yc)*yBx-yc*mu+xc))*mM+l^2*mA*mB*(mu*yBx+1))*R^2*tMx*sin(tM)-PA*(mA+mB+mM)*(((-mu*yBx-1)*l+xc*((yAx+mu)*yBx-mu*yAx+1))*mM-l*mA*(mu*yBx+1))*R^2-I_A*(((-mu*yBx-1)*l+xc*((yAx+mu)*yBx-mu*yAx+1))*mM+l*((-mA*mu+mB*yAx)*yBx-yAx*mB*mu-mA)))*cos(tM)+((((yBx-mu)*l+yc*((yAx+mu)*yBx-mu*yAx+1))*mM+l*mA*(yBx-mu))*PA*(mA+mB+mM)*R^2+I_A*(((yBx-mu)*l+yc*((yAx+mu)*yBx-mu*yAx+1))*mM+l*(mA+mB)*(yBx-mu)))*sin(tM)+(-xc*((yBx-mu)*l+(mu*yc-xc)*yBx+mu*xc+yc)*mM^2+(((-yc*mB*mu-xc*(mA+mB))*yBx+xc*(mA+mB)*mu-yc*mB)*l+(yBx-mu)*(mB*(xc^2+yc^2)+yc^2*mA+mA*xc^2+I_M))*mM+I_M*(mA+mB)*(yBx-mu))*R^2*tMx)*R/(-l*((((xc*yAx-yc)*mu+yc*yAx+xc)*yBx+(-yAx*yc-xc)*mu+yAx*xc-yc)*mM+l*mB*(mu*yAx*yBx-mu+yAx+yBx))*R^2*tMx*cos(tM)^2+l*(tMx*(-mB*((yAx-mu)*yBx-mu*yAx-1)*l+(((yAx*yc+xc)*mu-yAx*xc+yc)*yBx+(xc*yAx-yc)*mu+yc*yAx+xc)*mM)*R^2*sin(tM)-PA*((yAx^2+1)*(mu*yBx+1)*mM+((mA*yAx^2+mA+mB)*mu+yAx*mB)*yBx+yAx^2*mA-yAx*mB*mu+mA+mB)*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((yAx^2+1)*(yBx-mu)*mM+(mA*yAx^2+mB*mu*yAx+mB*yAx^2+mA)*yBx+(-mA*yAx^2-mB*yAx^2-mA)*mu+yAx*mB)*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(((yBx-mu)*(yAx*yc+xc)*l+((mu*yAx-1)*yBx+yAx+mu)*(xc^2+yc^2))*mM+mB*(mu*yBx+1)*yAx*l^2+I_M*((mu*yAx-1)*yBx+yAx+mu))*R^2*tMx);

TANA = R*(-l*(-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*tMx*cos(tM)^2-l*((mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*tMx*sin(tM)-yAx*mB*PA*((yBx-mu)*yAx+mu*yBx+1))*cos(tM)-l*mB*PA*((yBx-mu)*yAx+mu*yBx+1)*sin(tM)+(mB*(mu*yBx+1)*l^2+mM*(yBx-mu)*(-xc*yAx+yc)*l+(I_M+mM*(xc^2+yc^2))*((yBx-mu)*yAx+mu*yBx+1))*tMx)/(-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-PA*((mA+mM)*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB+mM)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((mA+mB+mM)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+(mA+mM)*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx);

TBNA = 0;

FNA_coeff.HNA = HNA;
FNA_coeff.JNA = JNA;
FNA_coeff.TANA = TANA;
FNA_coeff.TBNA = TBNA;

% Coefficients for Normal force at Wheel B

HNB = -phiB*(-2*l*(-(1/2)*(2*tMx^2*yAx*PA+(PA*yAx*yAxx+PAx*yAx^2-PAx)*tMx-yAx^2*tMxx*PA+tMxx*PA)*mB*l+(-(1/2)*PA*(yAx^2*yc+2*xc*yAx-yc)*tMx^2+(-(1/2)*PAx*yAx^2*xc+(-(1/2)*yAxx*PA*xc+yc*PAx)*yAx+(1/2)*yc*yAxx*PA+(1/2)*PAx*xc)*tMx-PA*tMxx*(-(1/2)*yAx^2*xc+yc*yAx+(1/2)*xc))*mM)*mB*R^2*phiA*cos(tM)^2+(-l*mB*((2*((-(1/2)*yAx^2*PA+(1/2)*PA)*tMx^2+(PAx*yAx+(1/2)*PA*yAxx)*tMx-yAx*tMxx*PA))*mB*l+(2*PA*(-(1/2)*yAx^2*xc+yc*yAx+(1/2)*xc)*tMx^2+(yc*PAx*yAx^2+(PA*yAxx*yc+2*PAx*xc)*yAx+yAxx*PA*xc-yc*PAx)*tMx-tMxx*PA*(yAx^2*yc+2*xc*yAx-yc))*mM)*R^2*phiA*sin(tM)+(tMx^3*l^3*mB^2+tMx^3*mB*mM*(yAx*yc+xc)*l^2+((I_M+mM*(xc^2+yc^2))*tMx^3+yAxx*PA^2*(mA+mB+mM))*mB*l+(I_M+mM*(xc^2+yc^2))*tMx^3*(yAx*yc+xc)*mM)*phiA*R^2+l*I_A*((PA*yAxx+PAx*yAx)*phiA-yAx*phiAx*PA)*mB*(yAx^2+1))*cos(tM)+(-(-tMx^3*yAx*l^3*mB^2+tMx^3*mB*mM*(-xc*yAx+yc)*l^2-((I_M+mM*(xc^2+yc^2))*tMx^3+yAxx*PA^2*(mA+mB+mM))*mB*yAx*l+(I_M+mM*(xc^2+yc^2))*tMx^3*(-xc*yAx+yc)*mM)*phiA*R^2-I_A*l*mB*(yAx^2+1)*(-PA*phiAx+PAx*phiA))*sin(tM)-((tMx^2*yAx*PA*mB+((PA*yAx*yAxx+PAx*yAx^2+PAx)*mM+PAx*(mA+mB)*yAx^2+yAxx*PA*(mA+mB)*yAx+mA*PAx)*tMx-((yAx^2+1)*mM+(mA+mB)*yAx^2+mA)*PA*tMxx)*mB*l^2-mB*(PA*(-xc*yAx+yc)*tMx^2+PAx*(yAx*yc+xc)*tMx-tMxx*PA*(yAx*yc+xc))*mM*l+(I_M+mM*(xc^2+yc^2))*(mA+mB+mM)*((PA*yAx*yAxx+PAx*yAx^2+PAx)*tMx-tMxx*PA*(yAx^2+1)))*phiA*R^2-(-phiA*tMxx+phiAx*tMx)*I_A*(yAx^2+1)*(l^2*mB+mM*(xc^2+yc^2)+I_M))*R/(phiA*(-l*tMx*R^2*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*yBx+1)*xc+yc*(yBx-mu))*yAx+(yBx-mu)*xc+(-mu*yBx-1)*yc)*mM)*cos(tM)^2+l*(tMx*(-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu-yBx)*xc+yc*(mu*yBx+1))*yAx+(mu*yBx+1)*xc+yc*(yBx-mu))*mM)*R^2*sin(tM)-PA*((yAx^2+1)*(mu*yBx+1)*mM+mA*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((yAx^2+1)*(yBx-mu)*mM+(mA+mB)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+mA*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+tMx*(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu+yAx-yBx))*R^2));

JNB = phiB*((-yAx*l^2*mB^2+2*mB*mM*(-xc*yAx+yc)*l+mM^2*((-xc^2+yc^2)*yAx+2*yc*xc))*R^2*tMx*cos(tM)^2+(-(-l^2*mB^2-2*mB*mM*(yAx*yc+xc)*l+mM^2*(-2*xc*yAx*yc-xc^2+yc^2))*R^2*tMx*sin(tM)-(l*mB+mM*xc*(yAx^2+1))*PA*(mA+mB+mM)*R^2-I_A*(yAx^2+1)*(l*mB+mM*xc))*cos(tM)+((-yAx*l*mB+yc*mM*(yAx^2+1))*PA*(mA+mB+mM)*R^2+yc*I_A*mM*(yAx^2+1))*sin(tM)+(yAx*mB*(mA+mB+mM)*l^2-yc*l*mB*mM+(xc^2*yAx-xc*yc)*mM^2+yAx*(mB*(xc^2+yc^2)+yc^2*mA+mA*xc^2+I_M)*mM+I_M*yAx*(mA+mB))*R^2*tMx)*R/(-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-PA*((yAx^2+1)*(mu*yBx+1)*mM+mA*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((yAx^2+1)*(yBx-mu)*mM+(mA+mB)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+mA*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+tMx*(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu+yAx-yBx))*R^2);

TANB = -(cos(tM)*PA*l*mB*yAx-sin(tM)*PA*l*mB+(l^2*mB+mM*(xc^2+yc^2)+I_M)*tMx)*R*phiB*(yAx^2+1)/(phiA*(-l*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-mu*xc-yc)*mM)*R^2*tMx*cos(tM)^2+l*((-((yBx-mu)*yAx-mu*yBx-1)*mB*l+(((mu*yc-xc)*yBx+mu*xc+yc)*yAx+(mu*xc+yc)*yBx-yc*mu+xc)*mM)*R^2*tMx*sin(tM)-PA*((mA+mM)*(mu*yBx+1)*yAx^2+mB*(yBx-mu)*yAx+(mA+mB+mM)*(mu*yBx+1))*R^2-(yAx^2+1)*(mu*yBx+1)*I_A)*cos(tM)-l*(PA*((mA+mB+mM)*(yBx-mu)*yAx^2+mB*(mu*yBx+1)*yAx+(mA+mM)*(yBx-mu))*R^2+I_A*(yAx^2+1)*(yBx-mu))*sin(tM)+(mB*(mu*yBx+1)*yAx*l^2+mM*(yBx-mu)*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*((mu*yBx+1)*yAx-yBx+mu))*R^2*tMx));

TBNB = 0;

FNB_coeff.HNB = HNB;
FNB_coeff.JNB = JNB;
FNB_coeff.TANB = TANB;
FNB_coeff.TBNB = TBNB;


