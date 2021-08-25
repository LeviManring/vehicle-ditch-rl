function [x_coeff, thetaA_coeff, thetaB_coeff, FNA_coeff, FNB_coeff, tM, dx] = allslip_eom(xA,guess_Dx,muA,muB,params)
%% allslip_eom
% Levi Manring, Duke University
% 2021
%
% This function calculates certain parameters needed for integration of the
% allslip dynamics model for the AWD complete model.
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
%       Outputs are elements necessary to solve the all-slip vehicle dynamics
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
% phiAx = yAx*yAxx/(R*sqrt(yAx^2+1));

phiB = sqrt(yBx^2+1)/R;
% phiBx = yBx*yBxx/(R*sqrt(yBx^2+1));


%% Calculating allslip dynamic model parameters
% Applying torque at wheel A and B and including moment of inertia effects for both wheels

% Coefficients for x

Hx = (l*(mB*((((muA*tMx^2-tMxx)*yBx+muA*tMxx+tMx^2)*muB+(-muA*tMxx-tMx^2)*yBx+muA*tMx^2-tMxx)*yAx+((muA*tMxx+tMx^2)*yBx-muA*tMx^2+tMxx)*muB+(muA*tMx^2-tMxx)*yBx+muA*tMxx+tMx^2)*l+(((((muA*xc+yc)*tMx^2+tMxx*(muA*yc-xc))*yBx+(-muA*yc+xc)*tMx^2+tMxx*(muA*xc+yc))*muB+((muA*yc-xc)*tMx^2-tMxx*(muA*xc+yc))*yBx+(muA*xc+yc)*tMx^2+tMxx*(muA*yc-xc))*yAx+(((-muA*yc+xc)*tMx^2+tMxx*(muA*xc+yc))*yBx+(-muA*xc-yc)*tMx^2-tMxx*(muA*yc-xc))*muB+((muA*xc+yc)*tMx^2+tMxx*(muA*yc-xc))*yBx+(-muA*yc+xc)*tMx^2+tMxx*(muA*xc+yc))*mM)*cos(tM)^2-l*((((((-muA*tMxx-tMx^2)*yBx+muA*tMx^2-tMxx)*muB+(-muA*tMx^2+tMxx)*yBx-muA*tMxx-tMx^2)*yAx+((muA*tMx^2-tMxx)*yBx+muA*tMxx+tMx^2)*muB+(-muA*tMxx-tMx^2)*yBx+muA*tMx^2-tMxx)*mB*l+(((((muA*yc-xc)*tMx^2-tMxx*(muA*xc+yc))*yBx+(muA*xc+yc)*tMx^2+tMxx*(muA*yc-xc))*muB+((-muA*xc-yc)*tMx^2-tMxx*(muA*yc-xc))*yBx+(muA*yc-xc)*tMx^2-tMxx*(muA*xc+yc))*yAx+(((muA*xc+yc)*tMx^2+tMxx*(muA*yc-xc))*yBx+(-muA*yc+xc)*tMx^2+tMxx*(muA*xc+yc))*muB+((muA*yc-xc)*tMx^2-tMxx*(muA*xc+yc))*yBx+(muA*xc+yc)*tMx^2+tMxx*(muA*yc-xc))*mM)*sin(tM)+(muB*yBx+1)*(-PA*muA*yAxx+PA*yAx*yAxx+PAx*yAx^2+PAx)*mM-PAx*((-mA*yBx+mB*muA)*muB-muA*yBx*mB-mA)*yAx^2+(((PA*mA*yAxx+PAx*mB*muA)*yBx-mB*(PA*muA*yAxx+PAx))*muB+mB*(PA*muA*yAxx+PAx)*yBx+PAx*muA*mB+yAxx*PA*mA)*yAx+((-muA*yAxx*PA*mA+PAx*(mA+mB))*yBx-yAxx*PA*mB)*muB+yBx*yAxx*PA*mB-muA*yAxx*PA*mA+PAx*(mA+mB))*cos(tM)-l*(-(muB-yBx)*(-PA*muA*yAxx+PA*yAx*yAxx+PAx*yAx^2+PAx)*mM-PAx*(mA+mB)*(muB-yBx)*yAx^2+((PAx*yBx*mB+PAx*muA*mB-yAxx*PA*(mA+mB))*muB+(-PAx*muA*mB+yAxx*PA*(mA+mB))*yBx+PAx*mB)*yAx+(-PAx*muA*yBx*mB+yAxx*PA*(mA+mB)*muA-mA*PAx)*muB+(-yAxx*PA*(mA+mB)*muA+mA*PAx)*yBx-PAx*muA*mB)*sin(tM)+(muA-yAx)*mB*((tMx^2-tMxx*yBx)*muB-tMx^2*yBx-tMxx)*l^2+(muB-yBx)*(((muA*yc-xc)*tMx^2-tMxx*(muA*xc+yc))*yAx+(muA*xc+yc)*tMx^2+tMxx*(muA*yc-xc))*mM*l+tMxx*(I_M+mM*(xc^2+yc^2))*(((muA+yBx)*muB-muA*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx))/(l*tMx*((((muA-yBx)*muB-muA*yBx-1)*yAx+muA*muB*yBx+muA+muB-yBx)*mB*l+((((muA*yc-xc)*yBx+muA*xc+yc)*muB+(-muA*xc-yc)*yBx+muA*yc-xc)*yAx+((muA*xc+yc)*yBx-muA*yc+xc)*muB+(muA*yc-xc)*yBx+muA*xc+yc)*mM)*cos(tM)^2+l*((mB*((muA*muB*yBx+muA+muB-yBx)*yAx+(-muA+yBx)*muB+muA*yBx+1)*l+((((muA*xc+yc)*yBx-muA*yc+xc)*muB+(muA*yc-xc)*yBx+muA*xc+yc)*yAx+((-muA*yc+xc)*yBx-muA*xc-yc)*muB+(muA*xc+yc)*yBx-muA*yc+xc)*mM)*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+((-mA*yBx+mB*muA)*muB-muA*yBx*mB-mA)*yAx^2-mB*((muA*yBx-1)*muB+muA+yBx)*yAx-(mA+mB)*(muB*yBx+1))*PA)*cos(tM)-(-(yAx^2+1)*(muB-yBx)*mM-(mA+mB)*(muB-yBx)*yAx^2+mB*((muA+yBx)*muB-muA*yBx+1)*yAx+(-mB*muA*yBx-mA)*muB-muA*mB+yBx*mA)*l*PA*sin(tM)-tMx*(mB*(muB*yBx+1)*(muA-yAx)*l^2+(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*mM*l-(I_M+mM*(xc^2+yc^2))*(((muA+yBx)*muB-muA*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx)));

Jx = ((((muB*yBx+1)*(muA-yAx)*mM+((-mA*yBx+mB*muA)*muB-muA*yBx*mB-mA)*yAx+(mA*muA*yBx+mB)*muB+muA*mA-mB*yBx)*l+xc*(((muA+yBx)*muB-muA*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx)*mM)*cos(tM)-((mA+mB+mM)*(muB-yBx)*(muA-yAx)*l+yc*(((muA+yBx)*muB-muA*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx)*mM)*sin(tM))/(l*tMx*((((muA-yBx)*muB-muA*yBx-1)*yAx+muA*muB*yBx+muA+muB-yBx)*mB*l+((((muA*yc-xc)*yBx+muA*xc+yc)*muB+(-muA*xc-yc)*yBx+muA*yc-xc)*yAx+((muA*xc+yc)*yBx-muA*yc+xc)*muB+(muA*yc-xc)*yBx+muA*xc+yc)*mM)*cos(tM)^2+l*((mB*((muA*muB*yBx+muA+muB-yBx)*yAx+(-muA+yBx)*muB+muA*yBx+1)*l+((((muA*xc+yc)*yBx-muA*yc+xc)*muB+(muA*yc-xc)*yBx+muA*xc+yc)*yAx+((-muA*yc+xc)*yBx-muA*xc-yc)*muB+(muA*xc+yc)*yBx-muA*yc+xc)*mM)*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+((-mA*yBx+mB*muA)*muB-muA*yBx*mB-mA)*yAx^2-mB*((muA*yBx-1)*muB+muA+yBx)*yAx-(mA+mB)*(muB*yBx+1))*PA)*cos(tM)-(-(yAx^2+1)*(muB-yBx)*mM-(mA+mB)*(muB-yBx)*yAx^2+mB*((muA+yBx)*muB-muA*yBx+1)*yAx+(-mB*muA*yBx-mA)*muB-muA*mB+yBx*mA)*l*PA*sin(tM)-tMx*(mB*(muB*yBx+1)*(muA-yAx)*l^2+(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*mM*l-(I_M+mM*(xc^2+yc^2))*(((muA+yBx)*muB-muA*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx)));

TAx = 0;

TBx = 0;

x_coeff.Hx = Hx;
x_coeff.Jx = Jx;
x_coeff.TAx = TAx;
x_coeff.TBx = TBx;

% Coefficients for thetaA

HthetaA = -R^2*muA*(((-(((xc*yAx-yc)*yBx-yAx*yc-xc)*muB+(yAx*yc+xc)*yBx+yAx*xc-yc)*PA*tMx^2+(yAxx*((yBx*yc+xc)*muB-xc*yBx+yc)*PA+(((yAx*yc+xc)*yBx+yAx*xc-yc)*muB+(-xc*yAx+yc)*yBx+yAx*yc+xc)*PAx)*tMx-(((yAx*yc+xc)*yBx+yAx*xc-yc)*muB+(-xc*yAx+yc)*yBx+yAx*yc+xc)*PA*tMxx)*mM+l*(-PA*(muB*yAx*yBx-muB+yAx+yBx)*tMx^2+(yAxx*(muB-yBx)*PA+((yAx+yBx)*muB-yAx*yBx+1)*PAx)*tMx-((yAx+yBx)*muB-yAx*yBx+1)*PA*tMxx)*mB)*l*(mA+mM)*cos(tM)^2+((((((yAx*yc+xc)*yBx+yAx*xc-yc)*muB+(-xc*yAx+yc)*yBx+yAx*yc+xc)*PA*tMx^2+(yAxx*((xc*yBx-yc)*muB+yc*yBx+xc)*PA+PAx*(((xc*yAx-yc)*yBx-yAx*yc-xc)*muB+(yAx*yc+xc)*yBx+yAx*xc-yc))*tMx-(((xc*yAx-yc)*yBx-yAx*yc-xc)*muB+(yAx*yc+xc)*yBx+yAx*xc-yc)*PA*tMxx)*mM+l*(((yAx+yBx)*muB-yAx*yBx+1)*PA*tMx^2+(yAxx*(muB*yBx+1)*PA+PAx*((yAx*yBx-1)*muB+yAx+yBx))*tMx-tMxx*PA*(muB*yAx*yBx-muB+yAx+yBx))*mB)*l*(mA+mM)*sin(tM)+(-((xc^2+yc^2)*tMx^3+yAxx*PA^2)*(muB*yBx+1)*l+(xc^2+yc^2)*((xc*yBx-yc)*muB+yc*yBx+xc)*tMx^3)*mM^2+(-((xc*yBx+yc)*muB-yc*yBx+xc)*tMx^3*mB*l^2+(mB*(xc^2+yc^2)*tMx^3-2*yAxx*(mA+(1/2)*mB)*PA^2)*(muB*yBx+1)*l+I_M*((xc*yBx-yc)*muB+yc*yBx+xc)*tMx^3)*mM+(I_M*tMx^3*mB-yAxx*PA^2*mA*(mA+mB))*l*(muB*yBx+1))*cos(tM)+(((muB-yBx)*((xc^2+yc^2)*tMx^3+yAxx*PA^2)*l-(xc^2+yc^2)*((yBx*yc+xc)*muB-xc*yBx+yc)*tMx^3)*mM^2+(-((yBx*yc-xc)*muB+xc*yBx+yc)*tMx^3*mB*l^2-(mB*(xc^2+yc^2)*tMx^3-2*yAxx*(mA+(1/2)*mB)*PA^2)*(muB-yBx)*l-I_M*((yBx*yc+xc)*muB-xc*yBx+yc)*tMx^3)*mM-(I_M*tMx^3*mB-yAxx*PA^2*mA*(mA+mB))*(muB-yBx)*l)*sin(tM)+(-(muB-yBx)*(PA*(yAx*yc+xc)*tMx^2+(yAxx*PA*xc+PAx*(xc*yAx-yc))*tMx-tMxx*PA*(xc*yAx-yc))*l+(xc^2+yc^2)*((yAxx*(muB-yBx)*PA+PAx*((yAx-yBx)*muB-yAx*yBx-1))*tMx-((yAx-yBx)*muB-yAx*yBx-1)*PA*tMxx))*mM^2+(-(tMx^2*PA*(muB-yBx)+PAx*(muB*yBx+1)*tMx-tMxx*PA*(muB*yBx+1))*mB*l^2+(-((-yBx*mB*yc+yc*(mA+mB)*yAx+mA*xc)*muB+(-yc*(mA+mB)*yAx-mA*xc)*yBx-mB*yc)*PA*tMx^2+(-yAxx*xc*(mA+mB)*(muB-yBx)*PA-((-yBx*mB*xc+xc*(mA+mB)*yAx-mA*yc)*muB+(-xc*(mA+mB)*yAx+mA*yc)*yBx-mB*xc)*PAx)*tMx+((-yBx*mB*xc+xc*(mA+mB)*yAx-mA*yc)*muB+(-xc*(mA+mB)*yAx+mA*yc)*yBx-mB*xc)*PA*tMxx)*l+(mB*(xc^2+yc^2)+mA*xc^2+mA*yc^2+I_M)*((yAxx*(muB-yBx)*PA+PAx*((yAx-yBx)*muB-yAx*yBx-1))*tMx-((yAx-yBx)*muB-yAx*yBx-1)*PA*tMxx))*mM-mA*(tMx^2*PA*(muB-yBx)+PAx*(muB*yBx+1)*tMx-tMxx*PA*(muB*yBx+1))*mB*l^2+(mA+mB)*I_M*((yAxx*(muB-yBx)*PA+PAx*((yAx-yBx)*muB-yAx*yBx-1))*tMx-((yAx-yBx)*muB-yAx*yBx-1)*PA*tMxx))*phiA/(I_A*((((((muA*yc-xc)*yAx+muA*xc+yc)*yBx+(muA*xc+yc)*yAx-muA*yc+xc)*muB+((-muA*xc-yc)*yAx+muA*yc-xc)*yBx+(muA*yc-xc)*yAx+muA*xc+yc)*mM+l*(((muA-yAx)*yBx+muA*yAx+1)*muB-muA*yAx*yBx+muA-yAx-yBx)*mB)*l*tMx*cos(tM)^2+l*((((((muA*xc+yc)*yAx-muA*yc+xc)*yBx+(-muA*yc+xc)*yAx-muA*xc-yc)*muB+((muA*yc-xc)*yAx+muA*xc+yc)*yBx+(muA*xc+yc)*yAx-muA*yc+xc)*mM+((muA*yAx*yBx-muA+yAx+yBx)*muB+(muA-yAx)*yBx+muA*yAx+1)*l*mB)*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+((-mA*yAx^2-mB*muA*yAx-mA-mB)*yBx+yAx*mB*(muA*yAx+1))*muB-yAx*mB*(muA*yAx+1)*yBx-muA*yAx*mB-yAx^2*mA-mA-mB)*PA)*cos(tM)-l*(-(yAx^2+1)*(muB-yBx)*mM+(-mB*(muA-yAx)*yBx+(-mA-mB)*yAx^2+muA*yAx*mB-mA)*muB+((mA+mB)*yAx^2-muA*yAx*mB+mA)*yBx-mB*(muA-yAx))*PA*sin(tM)+((-(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*l+(xc^2+yc^2)*(((yAx-muA)*yBx+muA*yAx+1)*muB-muA*yAx*yBx-muA+yAx-yBx))*mM-mB*(muB*yBx+1)*(muA-yAx)*l^2+I_M*(((yAx-muA)*yBx+muA*yAx+1)*muB-muA*yAx*yBx-muA+yAx-yBx))*tMx));

JthetaA = -R^2*phiA*muA*(((((yBx*yc+xc)*muB-xc*yBx+yc)*l+(-2*xc*yBx*yc-xc^2+yc^2)*muB+(xc^2-yc^2)*yBx-2*xc*yc)*mM^2+(mB*(muB-yBx)*l+((yBx*yc+xc)*muB-xc*yBx+yc)*(mA-mB))*l*mM+l^2*mA*mB*(muB-yBx))*tMx*cos(tM)^2+(tMx*((((xc*yBx-yc)*muB+yc*yBx+xc)*l+((-xc^2+yc^2)*yBx+2*xc*yc)*muB-2*yBx*xc*yc-xc^2+yc^2)*mM^2+(mB*(muB*yBx+1)*l+(mA-mB)*((xc*yBx-yc)*muB+yc*yBx+xc))*l*mM+l^2*mA*mB*(muB*yBx+1))*sin(tM)-(mA+mB+mM)*PA*(((muB*yBx+1)*l+xc*((yAx-yBx)*muB-yAx*yBx-1))*mM+l*mA*(muB*yBx+1)))*cos(tM)+(mA+mB+mM)*PA*(((muB-yBx)*l+((yAx-yBx)*muB-yAx*yBx-1)*yc)*mM+l*mA*(muB-yBx))*sin(tM)+(xc*((-muB+yBx)*l+(yBx*yc+xc)*muB-xc*yBx+yc)*mM^2+(((yBx*mB*yc-xc*(mA+mB))*muB+xc*(mA+mB)*yBx+mB*yc)*l+(mB*(xc^2+yc^2)+mA*xc^2+mA*yc^2+I_M)*(muB-yBx))*mM+I_M*(mA+mB)*(muB-yBx))*tMx)/(I_A*((((((muA*yc-xc)*yAx+muA*xc+yc)*yBx+(muA*xc+yc)*yAx-muA*yc+xc)*muB+((-muA*xc-yc)*yAx+muA*yc-xc)*yBx+(muA*yc-xc)*yAx+muA*xc+yc)*mM+l*(((muA-yAx)*yBx+muA*yAx+1)*muB-muA*yAx*yBx+muA-yAx-yBx)*mB)*l*tMx*cos(tM)^2+l*((((((muA*xc+yc)*yAx-muA*yc+xc)*yBx+(-muA*yc+xc)*yAx-muA*xc-yc)*muB+((muA*yc-xc)*yAx+muA*xc+yc)*yBx+(muA*xc+yc)*yAx-muA*yc+xc)*mM+((muA*yAx*yBx-muA+yAx+yBx)*muB+(muA-yAx)*yBx+muA*yAx+1)*l*mB)*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+((-mA*yAx^2-mB*muA*yAx-mA-mB)*yBx+yAx*mB*(muA*yAx+1))*muB-yAx*mB*(muA*yAx+1)*yBx-muA*yAx*mB-yAx^2*mA-mA-mB)*PA)*cos(tM)-l*(-(yAx^2+1)*(muB-yBx)*mM+(-mB*(muA-yAx)*yBx+yAx*(muA-yAx)*mB-yAx^2*mA-mA)*muB+(-yAx*(muA-yAx)*mB+mA*(yAx^2+1))*yBx-mB*(muA-yAx))*PA*sin(tM)+((-(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*l+(xc^2+yc^2)*(((yAx-muA)*yBx+muA*yAx+1)*muB-muA*yAx*yBx-muA+yAx-yBx))*mM-mB*(muB*yBx+1)*(muA-yAx)*l^2+I_M*(((yAx-muA)*yBx+muA*yAx+1)*muB-muA*yAx*yBx-muA+yAx-yBx))*tMx));

TAthetaA = 1/I_A;

TBthetaA = 0;

thetaA_coeff.HthetaA = HthetaA;
thetaA_coeff.JthetaA = JthetaA;
thetaA_coeff.TAthetaA = TAthetaA;
thetaA_coeff.TBthetaA = TBthetaA;

% Coefficients for thetaB

HthetaB = -R^2*((-(mB*(muA*yAx^2-muA+2*yAx)*l+mM*((muA*xc+yc)*yAx^2+(-2*muA*yc+2*xc)*yAx-muA*xc-yc))*PA*tMx^2+((2*(-(1/2)*PAx*yAx^2+(PAx*muA-(1/2)*PA*yAxx)*yAx+(1/2)*muA*yAxx*PA+(1/2)*PAx))*mB*l+mM*(PAx*(muA*yc-xc)*yAx^2+(yAxx*(muA*yc-xc)*PA+2*PAx*(muA*xc+yc))*yAx+yAxx*(muA*xc+yc)*PA-PAx*(muA*yc-xc)))*tMx-(mB*(2*muA*yAx-yAx^2+1)*l+mM*((muA*yc-xc)*yAx^2+(2*muA*xc+2*yc)*yAx-muA*yc+xc))*PA*tMxx)*l*mB*cos(tM)^2+(((mB*(2*muA*yAx-yAx^2+1)*l+mM*((muA*yc-xc)*yAx^2+(2*muA*xc+2*yc)*yAx-muA*yc+xc))*PA*tMx^2+((PAx*muA*yAx^2+(PA*muA*yAxx+2*PAx)*yAx-PAx*muA+PA*yAxx)*mB*l+mM*(PAx*(muA*xc+yc)*yAx^2+(yAxx*(muA*xc+yc)*PA-2*PAx*(muA*yc-xc))*yAx-yAxx*(muA*yc-xc)*PA-PAx*(muA*xc+yc)))*tMx-(mB*(muA*yAx^2-muA+2*yAx)*l+mM*((muA*xc+yc)*yAx^2+(-2*muA*yc+2*xc)*yAx-muA*xc-yc))*PA*tMxx)*l*mB*sin(tM)-(l*mB*(muA*yAx+1)+mM*((muA*xc+yc)*yAx-muA*yc+xc))*(mB*l^2+mM*(xc^2+yc^2)+I_M)*tMx^3-yAxx*l*PA^2*mB*(mA+mB+mM)*(muA*yAx+1))*cos(tM)+((l*mB*(muA-yAx)+((muA*yc-xc)*yAx+muA*xc+yc)*mM)*(mB*l^2+mM*(xc^2+yc^2)+I_M)*tMx^3+yAxx*l*PA^2*mB*(mA+mB+mM)*(muA-yAx))*sin(tM)-(l*mB*(muA-yAx)+((muA*yc-xc)*yAx+muA*xc+yc)*mM)*l*PA*mB*tMx^2+(-((PA*muA*yAxx-PA*yAx*yAxx-PAx*yAx^2-PAx)*mM+(muA-yAx)*(PA*yAxx+PAx*yAx)*mB-mA*(-PA*muA*yAxx+PA*yAx*yAxx+PAx*yAx^2+PAx))*mB*l^2-PAx*mM*mB*((muA*xc+yc)*yAx-muA*yc+xc)*l+(mA+mB+mM)*(I_M+mM*(xc^2+yc^2))*(-PA*muA*yAxx+PA*yAx*yAxx+PAx*yAx^2+PAx))*tMx-(-((-yAx^2-1)*mM+yAx*(muA-yAx)*mB-yAx^2*mA-mA)*mB*l^2-mM*mB*((muA*xc+yc)*yAx-muA*yc+xc)*l+(mA+mB+mM)*(I_M+mM*(xc^2+yc^2))*(yAx^2+1))*PA*tMxx)*phiB*muB/(I_B*(((((muB-yBx)*muA-muB*yBx-1)*yAx+muA*muB*yBx+muA+muB-yBx)*mB*l+mM*((((muB-yBx)*xc+yc*(muB*yBx+1))*muA+(-muB*yBx-1)*xc+yc*(muB-yBx))*yAx+((muB*yBx+1)*xc-yc*(muB-yBx))*muA+(muB-yBx)*xc+yc*(muB*yBx+1)))*l*tMx*cos(tM)^2+l*((((muA*muB*yBx+muA+muB-yBx)*yAx+(-muB+yBx)*muA+muB*yBx+1)*mB*l+mM*((((muB*yBx+1)*xc-yc*(muB-yBx))*muA+(muB-yBx)*xc+yc*(muB*yBx+1))*yAx+((-muB+yBx)*xc+(-muB*yBx-1)*yc)*muA+(muB*yBx+1)*xc-yc*(muB-yBx)))*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+((muB-yBx)*yAx-muB*yBx-1)*(muA*yAx+1)*mB-mA*(yAx^2+1)*(muB*yBx+1))*PA)*cos(tM)-l*(-(yAx^2+1)*(muB-yBx)*mM+((muB-yBx)*yAx-muB*yBx-1)*(muA-yAx)*mB-mA*(yAx^2+1)*(muB-yBx))*PA*sin(tM)+(-mB*(muB*yBx+1)*(muA-yAx)*l^2-(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*mM*l+(I_M+mM*(xc^2+yc^2))*(((muB-yBx)*muA+muB*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx))*tMx));

JthetaB = R^2*(-2*tMx*((1/2)*mB^2*(muA-yAx)*l^2+mM*((muA*yc-xc)*yAx+muA*xc+yc)*mB*l+((muA*xc*yc-(1/2)*xc^2+(1/2)*yc^2)*yAx+((1/2)*xc^2-(1/2)*yc^2)*muA+xc*yc)*mM^2)*cos(tM)^2+(-(mB^2*(muA*yAx+1)*l^2+2*mM*mB*((muA*xc+yc)*yAx-muA*yc+xc)*l+mM^2*(((xc^2-yc^2)*muA+2*xc*yc)*yAx-2*muA*xc*yc+xc^2-yc^2))*tMx*sin(tM)+(mA+mB+mM)*PA*(l*mB*(muA*yAx+1)+mM*xc*(yAx^2+1)))*cos(tM)-(l*mB*(muA-yAx)+mM*yc*(yAx^2+1))*(mA+mB+mM)*PA*sin(tM)+(mB*(mA+mB+mM)*(muA-yAx)*l^2+mB*mM*yc*(muA*yAx+1)*l+xc*((muA*yc-xc)*yAx+muA*xc+yc)*mM^2+(mB*(xc^2+yc^2)+mA*xc^2+mA*yc^2+I_M)*(muA-yAx)*mM+I_M*(mA+mB)*(muA-yAx))*tMx)*phiB*muB/(I_B*(tMx*l*((((muB-yBx)*muA-muB*yBx-1)*yAx+muA*muB*yBx+muA+muB-yBx)*mB*l+((((yBx*yc+xc)*muB-xc*yBx+yc)*muA+(-xc*yBx+yc)*muB-yc*yBx-xc)*yAx+((xc*yBx-yc)*muB+yc*yBx+xc)*muA+(yBx*yc+xc)*muB-xc*yBx+yc)*mM)*cos(tM)^2+((((muA*muB*yBx+muA+muB-yBx)*yAx+(-muB+yBx)*muA+muB*yBx+1)*mB*l+((((xc*yBx-yc)*muB+yc*yBx+xc)*muA+(yBx*yc+xc)*muB-xc*yBx+yc)*yAx+((-yBx*yc-xc)*muB+xc*yBx-yc)*muA+(xc*yBx-yc)*muB+yc*yBx+xc)*mM)*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+(mB*(muB-yBx)*muA-muB*yBx*mA-mA)*yAx^2-((muB*yBx+1)*muA-muB+yBx)*mB*yAx-(mA+mB)*(muB*yBx+1))*PA)*l*cos(tM)-l*(-(yAx^2+1)*(muB-yBx)*mM-(mA+mB)*(muB-yBx)*yAx^2+((muB-yBx)*muA+muB*yBx+1)*mB*yAx+(-muB*yBx-1)*mB*muA-mA*(muB-yBx))*PA*sin(tM)+(-mB*(muB*yBx+1)*(muA-yAx)*l^2-(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*mM*l+(I_M+mM*(xc^2+yc^2))*(((muB-yBx)*muA+muB*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx))*tMx));

TAthetaB = 0;

TBthetaB = 1/I_B;

thetaB_coeff.HthetaB = HthetaB;
thetaB_coeff.JthetaB = JthetaB;
thetaB_coeff.TAthetaB = TAthetaB;
thetaB_coeff.TBthetaB = TBthetaB;

% Coefficients for Normal force at Wheel A

HNA = R*phiA*(((-(((xc*yAx-yc)*yBx-yAx*yc-xc)*muB+(yAx*yc+xc)*yBx+yAx*xc-yc)*PA*tMx^2+(yAxx*((yBx*yc+xc)*muB-xc*yBx+yc)*PA+(((yAx*yc+xc)*yBx+yAx*xc-yc)*muB+(-xc*yAx+yc)*yBx+yAx*yc+xc)*PAx)*tMx-(((yAx*yc+xc)*yBx+yAx*xc-yc)*muB+(-xc*yAx+yc)*yBx+yAx*yc+xc)*PA*tMxx)*mM+l*(-PA*(muB*yAx*yBx-muB+yAx+yBx)*tMx^2+(yAxx*(muB-yBx)*PA+((yAx+yBx)*muB-yAx*yBx+1)*PAx)*tMx-((yAx+yBx)*muB-yAx*yBx+1)*PA*tMxx)*mB)*l*(mA+mM)*cos(tM)^2+((((((yAx*yc+xc)*yBx+yAx*xc-yc)*muB+(-xc*yAx+yc)*yBx+yAx*yc+xc)*PA*tMx^2+(yAxx*((xc*yBx-yc)*muB+yc*yBx+xc)*PA+PAx*(((xc*yAx-yc)*yBx-yAx*yc-xc)*muB+(yAx*yc+xc)*yBx+yAx*xc-yc))*tMx-(((xc*yAx-yc)*yBx-yAx*yc-xc)*muB+(yAx*yc+xc)*yBx+yAx*xc-yc)*PA*tMxx)*mM+l*(((yAx+yBx)*muB-yAx*yBx+1)*PA*tMx^2+(yAxx*(muB*yBx+1)*PA+PAx*((yAx*yBx-1)*muB+yAx+yBx))*tMx-tMxx*PA*(muB*yAx*yBx-muB+yAx+yBx))*mB)*l*(mA+mM)*sin(tM)+(-((xc^2+yc^2)*tMx^3+yAxx*PA^2)*(muB*yBx+1)*l+(xc^2+yc^2)*((xc*yBx-yc)*muB+yc*yBx+xc)*tMx^3)*mM^2+(-((xc*yBx+yc)*muB-yc*yBx+xc)*tMx^3*mB*l^2+(mB*(xc^2+yc^2)*tMx^3-2*yAxx*(mA+(1/2)*mB)*PA^2)*(muB*yBx+1)*l+I_M*((xc*yBx-yc)*muB+yc*yBx+xc)*tMx^3)*mM+(I_M*tMx^3*mB-yAxx*PA^2*mA*(mA+mB))*l*(muB*yBx+1))*cos(tM)+(((muB-yBx)*((xc^2+yc^2)*tMx^3+yAxx*PA^2)*l-(xc^2+yc^2)*((yBx*yc+xc)*muB-xc*yBx+yc)*tMx^3)*mM^2+(-((yBx*yc-xc)*muB+xc*yBx+yc)*tMx^3*mB*l^2-(mB*(xc^2+yc^2)*tMx^3-2*yAxx*(mA+(1/2)*mB)*PA^2)*(muB-yBx)*l-I_M*((yBx*yc+xc)*muB-xc*yBx+yc)*tMx^3)*mM-(I_M*tMx^3*mB-yAxx*PA^2*mA*(mA+mB))*(muB-yBx)*l)*sin(tM)+(-(muB-yBx)*(PA*(yAx*yc+xc)*tMx^2+(yAxx*PA*xc+PAx*(xc*yAx-yc))*tMx-tMxx*PA*(xc*yAx-yc))*l+(xc^2+yc^2)*((yAxx*(muB-yBx)*PA+PAx*((yAx-yBx)*muB-yAx*yBx-1))*tMx-((yAx-yBx)*muB-yAx*yBx-1)*PA*tMxx))*mM^2+(-(tMx^2*PA*(muB-yBx)+PAx*(muB*yBx+1)*tMx-tMxx*PA*(muB*yBx+1))*mB*l^2+(-((-yBx*mB*yc+yc*(mA+mB)*yAx+mA*xc)*muB+(-yc*(mA+mB)*yAx-mA*xc)*yBx-mB*yc)*PA*tMx^2+(-yAxx*xc*(mA+mB)*(muB-yBx)*PA-((-yBx*mB*xc+xc*(mA+mB)*yAx-mA*yc)*muB+(-xc*(mA+mB)*yAx+mA*yc)*yBx-mB*xc)*PAx)*tMx+((-yBx*mB*xc+xc*(mA+mB)*yAx-mA*yc)*muB+(-xc*(mA+mB)*yAx+mA*yc)*yBx-mB*xc)*PA*tMxx)*l+(mB*(xc^2+yc^2)+mA*xc^2+mA*yc^2+I_M)*((yAxx*(muB-yBx)*PA+PAx*((yAx-yBx)*muB-yAx*yBx-1))*tMx-((yAx-yBx)*muB-yAx*yBx-1)*PA*tMxx))*mM-mA*(tMx^2*PA*(muB-yBx)+PAx*(muB*yBx+1)*tMx-tMxx*PA*(muB*yBx+1))*mB*l^2+(mA+mB)*I_M*((yAxx*(muB-yBx)*PA+PAx*((yAx-yBx)*muB-yAx*yBx-1))*tMx-((yAx-yBx)*muB-yAx*yBx-1)*PA*tMxx))/((((((muA*yc-xc)*yAx+muA*xc+yc)*yBx+(muA*xc+yc)*yAx-muA*yc+xc)*muB+((-muA*xc-yc)*yAx+muA*yc-xc)*yBx+(muA*yc-xc)*yAx+muA*xc+yc)*mM+l*(((muA-yAx)*yBx+muA*yAx+1)*muB-muA*yAx*yBx+muA-yAx-yBx)*mB)*l*tMx*cos(tM)^2+l*((((((muA*xc+yc)*yAx-muA*yc+xc)*yBx+(-muA*yc+xc)*yAx-muA*xc-yc)*muB+((muA*yc-xc)*yAx+muA*xc+yc)*yBx+(muA*xc+yc)*yAx-muA*yc+xc)*mM+((muA*yAx*yBx-muA+yAx+yBx)*muB+(muA-yAx)*yBx+muA*yAx+1)*l*mB)*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+((-mA*yAx^2-mB*muA*yAx-mA-mB)*yBx+yAx*mB*(muA*yAx+1))*muB-yAx*mB*(muA*yAx+1)*yBx-muA*yAx*mB-yAx^2*mA-mA-mB)*PA)*cos(tM)-l*(-(yAx^2+1)*(muB-yBx)*mM+(-mB*(muA-yAx)*yBx+(-mA-mB)*yAx^2+muA*yAx*mB-mA)*muB+((mA+mB)*yAx^2-muA*yAx*mB+mA)*yBx-mB*(muA-yAx))*PA*sin(tM)+((-(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*l+(xc^2+yc^2)*(((yAx-muA)*yBx+muA*yAx+1)*muB-muA*yAx*yBx-muA+yAx-yBx))*mM-mB*(muB*yBx+1)*(muA-yAx)*l^2+I_M*(((yAx-muA)*yBx+muA*yAx+1)*muB-muA*yAx*yBx-muA+yAx-yBx))*tMx);

JNA = R*(((((yBx*yc+xc)*muB-xc*yBx+yc)*l+(-2*xc*yBx*yc-xc^2+yc^2)*muB+(xc^2-yc^2)*yBx-2*xc*yc)*mM^2+(mB*(muB-yBx)*l+((yBx*yc+xc)*muB-xc*yBx+yc)*(mA-mB))*l*mM+l^2*mA*mB*(muB-yBx))*tMx*cos(tM)^2+(tMx*((((xc*yBx-yc)*muB+yc*yBx+xc)*l+((-xc^2+yc^2)*yBx+2*xc*yc)*muB-2*yBx*xc*yc-xc^2+yc^2)*mM^2+(mB*(muB*yBx+1)*l+(mA-mB)*((xc*yBx-yc)*muB+yc*yBx+xc))*l*mM+l^2*mA*mB*(muB*yBx+1))*sin(tM)-(mA+mB+mM)*PA*(((muB*yBx+1)*l+xc*((yAx-yBx)*muB-yAx*yBx-1))*mM+l*mA*(muB*yBx+1)))*cos(tM)+(mA+mB+mM)*PA*(((muB-yBx)*l+((yAx-yBx)*muB-yAx*yBx-1)*yc)*mM+l*mA*(muB-yBx))*sin(tM)+(xc*((-muB+yBx)*l+(yBx*yc+xc)*muB-xc*yBx+yc)*mM^2+(((yBx*mB*yc-xc*(mA+mB))*muB+xc*(mA+mB)*yBx+mB*yc)*l+(mB*(xc^2+yc^2)+mA*xc^2+mA*yc^2+I_M)*(muB-yBx))*mM+I_M*(mA+mB)*(muB-yBx))*tMx)*phiA/((((((muA*yc-xc)*yAx+muA*xc+yc)*yBx+(muA*xc+yc)*yAx-muA*yc+xc)*muB+((-muA*xc-yc)*yAx+muA*yc-xc)*yBx+(muA*yc-xc)*yAx+muA*xc+yc)*mM+l*(((muA-yAx)*yBx+muA*yAx+1)*muB-muA*yAx*yBx+muA-yAx-yBx)*mB)*l*tMx*cos(tM)^2+l*((((((muA*xc+yc)*yAx-muA*yc+xc)*yBx+(-muA*yc+xc)*yAx-muA*xc-yc)*muB+((muA*yc-xc)*yAx+muA*xc+yc)*yBx+(muA*xc+yc)*yAx-muA*yc+xc)*mM+((muA*yAx*yBx-muA+yAx+yBx)*muB+(muA-yAx)*yBx+muA*yAx+1)*l*mB)*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+((-mA*yAx^2-mB*muA*yAx-mA-mB)*yBx+yAx*mB*(muA*yAx+1))*muB-yAx*mB*(muA*yAx+1)*yBx-muA*yAx*mB-yAx^2*mA-mA-mB)*PA)*cos(tM)-l*(-(yAx^2+1)*(muB-yBx)*mM+(-mB*(muA-yAx)*yBx+yAx*(muA-yAx)*mB-yAx^2*mA-mA)*muB+(-yAx*(muA-yAx)*mB+mA*(yAx^2+1))*yBx-mB*(muA-yAx))*PA*sin(tM)+((-(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*l+(xc^2+yc^2)*(((yAx-muA)*yBx+muA*yAx+1)*muB-muA*yAx*yBx-muA+yAx-yBx))*mM-mB*(muB*yBx+1)*(muA-yAx)*l^2+I_M*(((yAx-muA)*yBx+muA*yAx+1)*muB-muA*yAx*yBx-muA+yAx-yBx))*tMx);

TANA = 0;

TBNA = 0;

FNA_coeff.HNA = HNA;
FNA_coeff.JNA = JNA;
FNA_coeff.TANA = TANA;
FNA_coeff.TBNA = TBNA;

% Coefficients for Normal force at Wheel B

HNB = R*phiB*((-(mB*(muA*yAx^2-muA+2*yAx)*l+mM*((muA*xc+yc)*yAx^2+(-2*muA*yc+2*xc)*yAx-muA*xc-yc))*PA*tMx^2+((2*(-(1/2)*PAx*yAx^2+(PAx*muA-(1/2)*PA*yAxx)*yAx+(1/2)*muA*yAxx*PA+(1/2)*PAx))*mB*l+mM*(PAx*(muA*yc-xc)*yAx^2+(yAxx*(muA*yc-xc)*PA+2*PAx*(muA*xc+yc))*yAx+yAxx*(muA*xc+yc)*PA-PAx*(muA*yc-xc)))*tMx-(mB*(2*muA*yAx-yAx^2+1)*l+mM*((muA*yc-xc)*yAx^2+(2*muA*xc+2*yc)*yAx-muA*yc+xc))*PA*tMxx)*l*mB*cos(tM)^2+(((mB*(2*muA*yAx-yAx^2+1)*l+mM*((muA*yc-xc)*yAx^2+(2*muA*xc+2*yc)*yAx-muA*yc+xc))*PA*tMx^2+((PAx*muA*yAx^2+(PA*muA*yAxx+2*PAx)*yAx-PAx*muA+PA*yAxx)*mB*l+mM*(PAx*(muA*xc+yc)*yAx^2+(yAxx*(muA*xc+yc)*PA-2*PAx*(muA*yc-xc))*yAx-yAxx*(muA*yc-xc)*PA-PAx*(muA*xc+yc)))*tMx-(mB*(muA*yAx^2-muA+2*yAx)*l+mM*((muA*xc+yc)*yAx^2+(-2*muA*yc+2*xc)*yAx-muA*xc-yc))*PA*tMxx)*l*mB*sin(tM)-(l*mB*(muA*yAx+1)+mM*((muA*xc+yc)*yAx-muA*yc+xc))*(mB*l^2+mM*(xc^2+yc^2)+I_M)*tMx^3-yAxx*l*PA^2*mB*(mA+mB+mM)*(muA*yAx+1))*cos(tM)+((l*mB*(muA-yAx)+((muA*yc-xc)*yAx+muA*xc+yc)*mM)*(mB*l^2+mM*(xc^2+yc^2)+I_M)*tMx^3+yAxx*l*PA^2*mB*(mA+mB+mM)*(muA-yAx))*sin(tM)-(l*mB*(muA-yAx)+((muA*yc-xc)*yAx+muA*xc+yc)*mM)*l*PA*mB*tMx^2+(-((PA*muA*yAxx-PA*yAx*yAxx-PAx*yAx^2-PAx)*mM+(muA-yAx)*(PA*yAxx+PAx*yAx)*mB-mA*(-PA*muA*yAxx+PA*yAx*yAxx+PAx*yAx^2+PAx))*mB*l^2-PAx*mM*mB*((muA*xc+yc)*yAx-muA*yc+xc)*l+(mA+mB+mM)*(I_M+mM*(xc^2+yc^2))*(-PA*muA*yAxx+PA*yAx*yAxx+PAx*yAx^2+PAx))*tMx-(-((-yAx^2-1)*mM+yAx*(muA-yAx)*mB-yAx^2*mA-mA)*mB*l^2-mM*mB*((muA*xc+yc)*yAx-muA*yc+xc)*l+(mA+mB+mM)*(I_M+mM*(xc^2+yc^2))*(yAx^2+1))*PA*tMxx)/(((((muB-yBx)*muA-muB*yBx-1)*yAx+muA*muB*yBx+muA+muB-yBx)*mB*l+mM*((((muB-yBx)*xc+yc*(muB*yBx+1))*muA+(-muB*yBx-1)*xc+yc*(muB-yBx))*yAx+((muB*yBx+1)*xc-yc*(muB-yBx))*muA+(muB-yBx)*xc+yc*(muB*yBx+1)))*l*tMx*cos(tM)^2+l*((((muA*muB*yBx+muA+muB-yBx)*yAx+(-muB+yBx)*muA+muB*yBx+1)*mB*l+mM*((((muB*yBx+1)*xc-yc*(muB-yBx))*muA+(muB-yBx)*xc+yc*(muB*yBx+1))*yAx+((-muB+yBx)*xc+(-muB*yBx-1)*yc)*muA+(muB*yBx+1)*xc-yc*(muB-yBx)))*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+((muB-yBx)*yAx-muB*yBx-1)*(muA*yAx+1)*mB-mA*(yAx^2+1)*(muB*yBx+1))*PA)*cos(tM)-l*(-(yAx^2+1)*(muB-yBx)*mM+((muB-yBx)*yAx-muB*yBx-1)*(muA-yAx)*mB-mA*(yAx^2+1)*(muB-yBx))*PA*sin(tM)+(-mB*(muB*yBx+1)*(muA-yAx)*l^2-(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*mM*l+(I_M+mM*(xc^2+yc^2))*(((muB-yBx)*muA+muB*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx))*tMx);

JNB = -R*(-2*tMx*((1/2)*mB^2*(muA-yAx)*l^2+mM*((muA*yc-xc)*yAx+muA*xc+yc)*mB*l+((muA*xc*yc-(1/2)*xc^2+(1/2)*yc^2)*yAx+((1/2)*xc^2-(1/2)*yc^2)*muA+xc*yc)*mM^2)*cos(tM)^2+(-(mB^2*(muA*yAx+1)*l^2+2*mM*mB*((muA*xc+yc)*yAx-muA*yc+xc)*l+mM^2*(((xc^2-yc^2)*muA+2*xc*yc)*yAx-2*muA*xc*yc+xc^2-yc^2))*tMx*sin(tM)+(mA+mB+mM)*PA*(l*mB*(muA*yAx+1)+mM*xc*(yAx^2+1)))*cos(tM)-(l*mB*(muA-yAx)+mM*yc*(yAx^2+1))*(mA+mB+mM)*PA*sin(tM)+(mB*(mA+mB+mM)*(muA-yAx)*l^2+mB*mM*yc*(muA*yAx+1)*l+xc*((muA*yc-xc)*yAx+muA*xc+yc)*mM^2+(mB*(xc^2+yc^2)+mA*xc^2+mA*yc^2+I_M)*(muA-yAx)*mM+I_M*(mA+mB)*(muA-yAx))*tMx)*phiB/(tMx*l*((((muB-yBx)*muA-muB*yBx-1)*yAx+muA*muB*yBx+muA+muB-yBx)*mB*l+((((yBx*yc+xc)*muB-xc*yBx+yc)*muA+(-xc*yBx+yc)*muB-yc*yBx-xc)*yAx+((xc*yBx-yc)*muB+yc*yBx+xc)*muA+(yBx*yc+xc)*muB-xc*yBx+yc)*mM)*cos(tM)^2+((((muA*muB*yBx+muA+muB-yBx)*yAx+(-muB+yBx)*muA+muB*yBx+1)*mB*l+((((xc*yBx-yc)*muB+yc*yBx+xc)*muA+(yBx*yc+xc)*muB-xc*yBx+yc)*yAx+((-yBx*yc-xc)*muB+xc*yBx-yc)*muA+(xc*yBx-yc)*muB+yc*yBx+xc)*mM)*tMx*sin(tM)+(-(yAx^2+1)*(muB*yBx+1)*mM+(mB*(muB-yBx)*muA-muB*yBx*mA-mA)*yAx^2-((muB*yBx+1)*muA-muB+yBx)*mB*yAx-(mA+mB)*(muB*yBx+1))*PA)*l*cos(tM)-l*(-(yAx^2+1)*(muB-yBx)*mM-(mA+mB)*(muB-yBx)*yAx^2+((muB-yBx)*muA+muB*yBx+1)*mB*yAx+(-muB*yBx-1)*mB*muA-mA*(muB-yBx))*PA*sin(tM)+(-mB*(muB*yBx+1)*(muA-yAx)*l^2-(muB-yBx)*((muA*xc+yc)*yAx-muA*yc+xc)*mM*l+(I_M+mM*(xc^2+yc^2))*(((muB-yBx)*muA+muB*yBx+1)*yAx-muA*muB*yBx-muA+muB-yBx))*tMx);

TANB = 0;

TBNB = 0;

FNB_coeff.HNB = HNB;
FNB_coeff.JNB = JNB;
FNB_coeff.TANB = TANB;
FNB_coeff.TBNB = TBNB;




