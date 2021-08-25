function [x_coeff, thetaA_coeff, thetaB_coeff, FNA_coeff, FNB_coeff, tM, dx] = rearslipfrontstick_eom(xA,guess_Dx,mu,params)
%% rearslipfrontstick_eom
% Levi Manring, Duke University
% 2021
%
% This function calculates certain parameters needed for integration of the
% rearslipfrontstick dynamics model for the AWD complete model.
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
%       Outputs are elements necessary to solve the rearslipfrontstick vehicle dynamics
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
[~, yBx, yBxx, yBxxx, ~] = y_eval_fcn(xB,params);

% Calculate certain dynamic model parameters
PA = 1-R*yAxx/(yAx^2+1)^(3/2);
PAx = -R*(yAx^2*yAxxx-3*yAx*yAxx^2+yAxxx)/(yAx^2+1)^(5/2);

PB = 1-R*yBxx/(yBx^2+1)^(3/2);
PBx = -R*(yBx^2*yBxxx-3*yBx*yBxx^2+yBxxx)/(yBx^2+1)^(5/2);

phiA = sqrt(yAx^2+1)/R;
% phiAx = yAx*yAxx/(R*sqrt(yAx^2+1));

phiB = sqrt(yBx^2+1)/R;
phiBx = yBx*yBxx/(R*sqrt(yBx^2+1));


%% Calculating rearslipfrontstick dynamic model parameters
% Applying torque at wheel A and B and including moment of inertia effects for both wheels

% Coefficients for x

Hx = (((mB*R^2*phiB*(((mu*tMxx+tMx^2)*yAx-tMx^2*mu+tMxx)*yBx+(-mu*tMx^2+tMxx)*yAx-tMx^2-tMxx*mu)*PB^2-(((mu*tMx^2-tMxx)*yAx+tMx^2+tMxx*mu)*phiB-tMx*phiBx*(yAx-mu))*(yBx^2+1)*I_B*PB-I_B*PBx*tMx*phiB*(yBx^2+1)*(yAx-mu))*l+PB^2*R^2*phiB*mM*((((-mu*yc+xc)*tMx^2+tMxx*(mu*xc+yc))*yAx+(-mu*xc-yc)*tMx^2+tMxx*(-mu*yc+xc))*yBx+((-mu*xc-yc)*tMx^2+tMxx*(-mu*yc+xc))*yAx+(mu*yc-xc)*tMx^2-tMxx*(mu*xc+yc)))*l*cos(tM)^2-(((mB*R^2*(((mu*tMx^2-tMxx)*yAx+tMx^2+tMxx*mu)*yBx+(mu*tMxx+tMx^2)*yAx-tMx^2*mu+tMxx)*phiB*PB^2+(yBx^2+1)*I_B*(((mu*tMxx+tMx^2)*yAx-tMx^2*mu+tMxx)*phiB+tMx*phiBx*(mu*yAx+1))*PB-I_B*PBx*tMx*phiB*(yBx^2+1)*(mu*yAx+1))*l+PB^2*((((mu*xc+yc)*tMx^2-tMxx*(-mu*yc+xc))*yAx+(-mu*yc+xc)*tMx^2+tMxx*(mu*xc+yc))*yBx+((-mu*yc+xc)*tMx^2+tMxx*(mu*xc+yc))*yAx+(-mu*xc-yc)*tMx^2+tMxx*(-mu*yc+xc))*R^2*phiB*mM)*sin(tM)-(mB*(mu*yAx+1)*(PA*yAxx+PAx*yAx)*yBx+PAx*(mA+mM)*yAx^2+(PAx*mB*mu+yAxx*PA*(mA+mM))*yAx-yAxx*PA*(mA+mM)*mu+PAx*(mA+mB+mM))*R^2*phiB*PB^2-I_B*(yBx^2+1)*(mu*yAx+1)*(PA*phiBx+PAx*phiB)*PB+I_B*PBx*phiB*PA*(yBx^2+1)*(mu*yAx+1))*l*cos(tM)+(R^2*phiB*((PAx*(mA+mB+mM)*yAx^2+(-PAx*mB*mu+yAxx*PA*(mA+mB+mM))*yAx-yAxx*PA*(mA+mB+mM)*mu+PAx*(mA+mM))*yBx+PAx*mB*(yAx-mu))*PB^2+I_B*(yBx^2+1)*(yAx-mu)*(PA*phiBx+PAx*phiB)*PB-I_B*PBx*phiB*PA*(yBx^2+1)*(yAx-mu))*l*sin(tM)-(yAx-mu)*(R^2*phiB*mB*(tMx^2*yBx+tMxx)*PB^2+I_B*(yBx^2+1)*(phiB*tMxx+phiBx*tMx)*PB-I_B*PBx*tMx*phiB*(yBx^2+1))*l^2-PB^2*R^2*phiB*mM*(((-mu*yc+xc)*tMx^2+tMxx*(mu*xc+yc))*yAx+(-mu*xc-yc)*tMx^2+tMxx*(-mu*yc+xc))*yBx*l+PB^2*tMxx*R^2*phiB*(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx))/(PB*phiB*(((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*tMx*l*cos(tM)^2+(((mB*R^2*((yAx-mu)*yBx-yAx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*R^2*mM)*tMx*sin(tM)+(R^2*(yAx*mB*(mu*yAx+1)*yBx+(mA+mM)*yAx^2+yAx*mB*mu+mA+mB+mM)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+(R^2*(((mA+mB+mM)*yAx^2-yAx*mB*mu+mA+mM)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*PA*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx));

Jx = PB*R^2*(((mB*(mu*yAx+1)*yBx+(mA+mM)*(yAx-mu))*l+xc*mM*((mu*yAx+1)*yBx-yAx+mu))*cos(tM)-(-yBx*(mA+mB+mM)*(yAx-mu)*l+mM*((mu*yAx+1)*yBx-yAx+mu)*yc)*sin(tM))/(((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*tMx*l*cos(tM)^2+(((mB*R^2*((yAx-mu)*yBx-yAx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*R^2*mM)*tMx*sin(tM)+(R^2*(yAx*mB*(mu*yAx+1)*yBx+(mA+mM)*yAx^2+yAx*mB*mu+mA+mB+mM)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+(R^2*(((mA+mB+mM)*yAx^2-yAx*mB*mu+mA+mM)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*PA*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx);

TAx = 0;

TBx = -l*(yBx^2+1)*((yAx-mu)*sin(tM)+(mu*yAx+1)*cos(tM))*PB/(phiB*(((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*mM*(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-xc*mu-yc))*tMx*l*cos(tM)^2+(((mB*R^2*((yBx-mu)*yAx-yBx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-mu*yc+xc)*yBx-xc*mu-yc)*yAx+(-mu*xc-yc)*yBx+yc*mu-xc)*R^2*mM)*tMx*sin(tM)+(R^2*((mB*mu*yBx+mA+mM)*yAx^2+mB*(yBx+mu)*yAx+mA+mB+mM)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+PA*(R^2*(yBx*(mA+mB+mM)*yAx^2+(-mB*mu*yBx+mB)*yAx+(mA+mM)*yBx-mB*mu)*PB+I_B*(yBx^2+1)*(yAx-mu))*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx));

x_coeff.Hx = Hx;
x_coeff.Jx = Jx;
x_coeff.TAx = TAx;
x_coeff.TBx = TBx;

% Coefficients for thetaA

HthetaA = phiA*R^2*mu*(-(mA+mM)*(phiB*R^2*(mB*(PA*(yAx+yBx)*tMx^2+(yBx*yAxx*PA+PAx*(yAx*yBx-1))*tMx-tMxx*PA*(yAx*yBx-1))*l+(PA*((yAx*yc+xc)*yBx+xc*yAx-yc)*tMx^2+(yAxx*(xc*yBx-yc)*PA+PAx*((xc*yAx-yc)*yBx-yc*yAx-xc))*tMx-PA*tMxx*((xc*yAx-yc)*yBx-yc*yAx-xc))*mM)*PB^2-((-PA*tMx^2*yAx-PA*tMxx+PAx*tMx)*phiB-tMx*phiBx*PA)*(yBx^2+1)*I_B*l*PB-I_B*PBx*tMx*l*phiB*PA*(yBx^2+1))*l*cos(tM)^2+((mA+mM)*(phiB*(mB*((-yAx*yBx+1)*PA*tMx^2+(PA*yAxx+PAx*(yAx+yBx))*tMx-tMxx*PA*(yAx+yBx))*l+(-PA*((xc*yAx-yc)*yBx-yc*yAx-xc)*tMx^2+(yAxx*(yBx*yc+xc)*PA+((yAx*yc+xc)*yBx+xc*yAx-yc)*PAx)*tMx-PA*((yAx*yc+xc)*yBx+xc*yAx-yc)*tMxx)*mM)*R^2*PB^2+(yBx^2+1)*I_B*((tMx^2*PA+(PA*yAxx+PAx*yAx)*tMx-yAx*tMxx*PA)*phiB-tMx*yAx*phiBx*PA)*l*PB+I_B*PBx*tMx*yAx*l*phiB*PA*(yBx^2+1))*l*sin(tM)+phiB*R^2*(-tMx^3*mB*mM*(-yBx*yc+xc)*l^2+(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*l+(I_M+mM*(xc^2+yc^2))*mM*(yBx*yc+xc)*tMx^3)*PB^2+(-xc*tMx^3*l*phiB*mM+((I_M+mM*(xc^2+yc^2))*tMx^3-yAxx*PA^2*(mA+mM))*phiB+yAx*phiBx*PA^2*(mA+mM))*(yBx^2+1)*I_B*l*PB-I_B*PBx*yAx*l*phiB*PA^2*(yBx^2+1)*(mA+mM))*cos(tM)+((-tMx^3*mB*mM*(xc*yBx+yc)*l^2+(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*yBx*l+(I_M+mM*(xc^2+yc^2))*mM*(xc*yBx-yc)*tMx^3)*phiB*R^2*PB^2+(yBx^2+1)*(-tMx^2*mM*(phiB*tMx*yc+phiBx*xc)*l+phiBx*((I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM)))*I_B*l*PB-phiB*PBx*(yBx^2+1)*(-xc*tMx^2*l*mM+(I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM))*I_B*l)*sin(tM)-phiB*R^2*(mB*(mA+mM)*(-PA*tMx^2*yBx-PA*tMxx+PAx*tMx)*l^2-(PA*((yAx*yc+xc)*yBx*mM+(yc*(mA+mB)*yAx+xc*mA)*yBx+yc*mB)*tMx^2+((xc*yAxx*PA+PAx*(xc*yAx-yc))*yBx*mM+xc*yBx*yAxx*(mA+mB)*PA+((xc*(mA+mB)*yAx-yc*mA)*yBx+xc*mB)*PAx)*tMx-((xc*yAx-yc)*yBx*mM+(xc*(mA+mB)*yAx-yc*mA)*yBx+xc*mB)*PA*tMxx)*mM*l+((yBx*yAxx*PA+PAx*(yAx*yBx+1))*tMx-tMxx*PA*(yAx*yBx+1))*(I_M+mM*(xc^2+yc^2))*(mA+mB+mM))*PB^2-(yBx^2+1)*(((-PA*tMxx+PAx*tMx)*phiB-tMx*phiBx*PA)*(mA+mM)*l^2-mM*((PA*tMx^2*yc-PA*tMxx*xc+PAx*tMx*xc)*phiB+xc*tMx*phiBx*PA)*l+((-PA*tMxx+PAx*tMx)*phiB+tMx*phiBx*PA)*(I_M+mM*(xc^2+yc^2)))*I_B*PB+phiB*((-mA-mM)*l^2-xc*l*mM+mM*(xc^2+yc^2)+I_M)*PBx*(yBx^2+1)*PA*I_B*tMx)/(PB*phiB*((R^2*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*PB+I_B*l*(yBx^2+1)*(yAx-mu))*tMx*l*cos(tM)^2+(((mB*((yAx-mu)*yBx-yAx*mu-1)*l+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*mM)*R^2*PB-I_B*l*(yBx^2+1)*(mu*yAx+1))*tMx*sin(tM)+(R^2*((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+(R^2*(yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*PA*l*sin(tM)+(R^2*(-mB*(yAx-mu)*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*mM*yBx*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx))*PB-I_B*l^2*(yBx^2+1)*(yAx-mu))*tMx)*I_A);

JthetaA = R^2*mu*((PB*(-yBx*mB*(mA+mM)*l^2-mM*(mA-mB+mM)*(xc*yBx-yc)*l+((xc^2-yc^2)*yBx-2*xc*yc)*mM^2)*R^2-yc*I_B*l*mM*(yBx^2+1))*tMx*cos(tM)^2+(-tMx*(PB*(-mB*(mA+mM)*l^2-mM*(mA-mB+mM)*(yBx*yc+xc)*l+mM^2*(2*xc*yBx*yc+xc^2-yc^2))*R^2+(yBx^2+1)*((-mA-mM)*l+mM*xc)*I_B*l)*sin(tM)+PA*(PB*((-mA-mM)*l+xc*mM*(yAx*yBx+1))*(mA+mB+mM)*R^2+(yBx^2+1)*((-mA-mM)*l+mM*xc)*I_B))*cos(tM)-PA*((yBx*(mA+mM)*l+yc*mM*(yAx*yBx+1))*PB*(mA+mB+mM)*R^2+yc*I_B*mM*(yBx^2+1))*sin(tM)-tMx*(PB*(-(xc*yBx*mM+xc*(mA+mB)*yBx+yc*mB)*mM*l+xc*(xc*yBx-yc)*mM^2+((mA+mB)*xc^2+(mA+mB)*yc^2+I_M)*yBx*mM+I_M*yBx*(mA+mB))*R^2-yc*I_B*l*mM*(yBx^2+1)))*phiA/(I_A*((R^2*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*PB+I_B*l*(yBx^2+1)*(yAx-mu))*tMx*l*cos(tM)^2+(((mB*((yAx-mu)*yBx-yAx*mu-1)*l+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*mM)*R^2*PB-I_B*l*(yBx^2+1)*(mu*yAx+1))*tMx*sin(tM)+(R^2*((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+(R^2*(yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*PA*l*sin(tM)+(R^2*(-mB*(yAx-mu)*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*mM*yBx*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx))*PB-I_B*l^2*(yBx^2+1)*(yAx-mu))*tMx));

TAthetaA = 1/I_A;

TBthetaA = PB*(-yAx*l*PA*(mA+mM)*cos(tM)+l*PA*(mA+mM)*sin(tM)+(-xc*l*mM+mM*(xc^2+yc^2)+I_M)*tMx)*R^2*mu*(yBx^2+1)*phiA/((((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*mM*(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-xc*mu-yc))*tMx*l*cos(tM)^2+(((mB*R^2*((yBx-mu)*yAx-yBx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-mu*yc+xc)*yBx-xc*mu-yc)*yAx+(-mu*xc-yc)*yBx+yc*mu-xc)*R^2*mM)*tMx*sin(tM)+(R^2*((mB*mu*yBx+mA+mM)*yAx^2+mB*(yBx+mu)*yAx+mA+mB+mM)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+PA*(R^2*(yBx*(mA+mB+mM)*yAx^2+(-mB*mu*yBx+mB)*yAx+(mA+mM)*yBx-mB*mu)*PB+I_B*(yBx^2+1)*(yAx-mu))*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx)*phiB*I_A);

thetaA_coeff.HthetaA = HthetaA;
thetaA_coeff.JthetaA = JthetaA;
thetaA_coeff.TAthetaA = TAthetaA;
thetaA_coeff.TBthetaA = TBthetaA;

% Coefficients for thetaB

HthetaB = R^2*(-(mB*((yAx-mu)*yBx-yAx*mu-1)*l+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*mM)*(-PB*phiBx+PBx*phiB)*tMx^2*l^2*cos(tM)^3+((-PB*phiBx+PBx*phiB)*tMx^2*l*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*sin(tM)+(-PB*(mB*(mu*yAx^2-mu+2*yAx)*yBx+(yAx^2+1)*(mA+mM))*PA*phiB*tMx^2+((-((yAxx*(yAx-mu)*PA+PAx*(yAx^2+1))*mM+yAxx*(mA+mB)*(yAx-mu)*PA+((mA+mB)*yAx^2-2*yAx*mB*mu-mB+mA)*PAx)*yBx*phiB+((yAx^2+1)*(mA+mB+mM)*yBx+2*mB*(yAx-mu))*phiBx*PA)*PB-((yAx^2+1)*(mA+mB+mM)*yBx+2*mB*(yAx-mu))*PBx*PA*phiB)*tMx+PB*tMxx*PA*phiB*((yAx^2+1)*mM+(mA+mB)*yAx^2-2*yAx*mB*mu-mB+mA)*yBx)*l+mM*(-PB*(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*PA*phiB*tMx^2+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*((PA*phiBx+PAx*phiB)*PB-PA*PBx*phiB)*tMx-PB*tMxx*PA*phiB*(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)))*l*cos(tM)^2-(((PB*PA*phiB*((yAx^2+1)*mM+(mA+mB)*yAx^2-2*yAx*mB*mu-mB+mA)*yBx*tMx^2+(((-mB*(yAxx*(mu*yAx+1)*PA+PAx*(mu*yAx^2-mu+2*yAx))*yBx-(yAxx*(yAx-mu)*PA+PAx*(yAx^2+1))*(mA+mM))*phiB+(mB*mu*(yAx^2+1)*yBx+(yAx^2+1)*mM+yAx^2*mA+2*yAx*mB*mu+mA+2*mB)*phiBx*PA)*PB-(mB*mu*(yAx^2+1)*yBx+(yAx^2+1)*mM+yAx^2*mA+2*yAx*mB*mu+mA+2*mB)*PBx*PA*phiB)*tMx+PB*tMxx*(mB*(mu*yAx^2-mu+2*yAx)*yBx+(yAx^2+1)*(mA+mM))*PA*phiB)*l-mM*(PB*PA*phiB*(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*tMx^2+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*((PA*phiBx+PAx*phiB)*PB-PA*PBx*phiB)*tMx-PB*tMxx*(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*PA*phiB))*sin(tM)-mB*tMx^2*(-yBx*phiB*PB*(mu*yAx+1)*tMx+((yAx-mu)*yBx-yAx*mu-1)*(-PB*phiBx+PBx*phiB))*l^2-(-PB*phiB*((-mu*yc+xc)*yAx-xc*mu-yc)*tMx+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*(-PB*phiBx+PBx*phiB))*mM*tMx^2*l+PB*(mu*yAx*yBx+mu-yAx+yBx)*phiB*(I_M+mM*(xc^2+yc^2))*tMx^3+PA^2*((yAxx*(mB*(mu*yAx+1)*yBx+(mA+mM)*(yAx-mu))*phiB-phiBx*((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB))*PB+PBx*phiB*((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB)))*l*cos(tM)+(-tMx^2*mB*(yAx-mu)*(PB*phiB*tMx*yBx-PB*phiBx+PBx*phiB)*l^2-mM*(PB*phiB*((-mu*yc+xc)*yAx-xc*mu-yc)*tMx+((mu*xc+yc)*yAx-yc*mu+xc)*(-PB*phiBx+PBx*phiB))*tMx^2*yBx*l+(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2))*(-PB*phiBx+PBx*phiB)*tMx^2-PA^2*((yBx*yAxx*(mA+mB+mM)*(yAx-mu)*phiB-phiBx*(((yAx^2+1)*mM+(mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu)))*PB+PBx*phiB*(((yAx^2+1)*mM+(mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu))))*l*sin(tM)+(yBx*phiB*PA*PB*mB*(yAx-mu)*tMx^2+((((yAxx*(yAx-mu)*PA+PAx*(yAx^2+1))*mM+yAxx*(mA+mB)*(yAx-mu)*PA+((mA+mB)*yAx^2-yAx*mB*mu+mA)*PAx)*yBx*phiB-phiBx*PA*(((yAx^2+1)*mM+(mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+2*mB*(yAx-mu)))*PB+PBx*PA*phiB*(((yAx^2+1)*mM+(mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+2*mB*(yAx-mu)))*tMx-PB*tMxx*PA*phiB*((yAx^2+1)*mM+(mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx)*l^2-(-PB*PA*phiB*((-mu*yc+xc)*yAx-xc*mu-yc)*tMx^2+((mu*xc+yc)*yAx-yc*mu+xc)*((PA*phiBx+PAx*phiB)*PB-PA*PBx*phiB)*tMx-PB*tMxx*((mu*xc+yc)*yAx-yc*mu+xc)*PA*phiB)*mM*yBx*l+(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2))*(((PA*phiBx+PAx*phiB)*PB-PA*PBx*phiB)*tMx-phiB*tMxx*PA*PB))/(PB*(tMx*l*((mB*R^2*((mu*yAx+1)*yBx+yAx-mu)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*cos(tM)^2+(((mB*R^2*((yAx-mu)*yBx-yAx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*R^2*mM)*tMx*sin(tM)+(R^2*((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+PA*(R^2*(((yAx^2+1)*mM+(mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx));

JthetaB = R^2*phiB*((-yBx*(mA+mB+mM)*(yAx-mu)*l+mM*((mu*yAx+1)*yBx-yAx+mu)*yc)*tMx*l*cos(tM)^2+((mB*(mu*yAx+1)*yBx+(mA+mM)*(yAx-mu))*l+xc*mM*((mu*yAx+1)*yBx-yAx+mu))*(l*tMx*sin(tM)-PA)*cos(tM)-(l*tMx-sin(tM)*PA)*(-yBx*(mA+mB+mM)*(yAx-mu)*l+mM*((mu*yAx+1)*yBx-yAx+mu)*yc))/(((I_B*(yAx-mu)*yBx^2+R^2*PB*mB*(mu*yAx+1)*yBx+(yAx-mu)*(PB*R^2*mB+I_B))*l+PB*R^2*(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*tMx*l*cos(tM)^2+(tMx*(((-I_B*mu*yAx-I_B)*yBx^2+R^2*PB*mB*(yAx-mu)*yBx-(mu*yAx+1)*(PB*R^2*mB+I_B))*l+PB*(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*R^2*mM)*sin(tM)+(I_B*(mu*yAx+1)*yBx^2+R^2*yAx*PB*mB*(mu*yAx+1)*yBx+R^2*PB*(mA+mM)*yAx^2+mu*(PB*R^2*mB+I_B)*yAx+PB*R^2*mM+R^2*(mA+mB)*PB+I_B)*PA)*l*cos(tM)+PA*(I_B*(yAx-mu)*yBx^2+PB*R^2*((mA+mB+mM)*yAx^2-yAx*mB*mu+mA+mM)*yBx+(yAx-mu)*(PB*R^2*mB+I_B))*l*sin(tM)+(-(yAx-mu)*(PB*R^2*mB+I_B*yBx^2+I_B)*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*((mu*yAx+1)*yBx-yAx+mu)*(I_M+mM*(xc^2+yc^2)))*tMx);

TAthetaB = 0;

TBthetaB = (tMx*l*(yAx-mu)*cos(tM)^2-(mu*yAx+1)*(l*tMx*sin(tM)-PA)*cos(tM)-(yAx-mu)*(l*tMx-sin(tM)*PA))*(yBx^2+1)*l/(((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*mM*(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-xc*mu-yc))*tMx*l*cos(tM)^2+(((mB*R^2*((yBx-mu)*yAx-yBx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-mu*yc+xc)*yBx-xc*mu-yc)*yAx+(-mu*xc-yc)*yBx+yc*mu-xc)*R^2*mM)*tMx*sin(tM)+(R^2*((mB*mu*yBx+mA+mM)*yAx^2+mB*(yBx+mu)*yAx+mA+mB+mM)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+PA*(R^2*(yBx*(mA+mB+mM)*yAx^2+(-mB*mu*yBx+mB)*yAx+(mA+mM)*yBx-mB*mu)*PB+I_B*(yBx^2+1)*(yAx-mu))*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx);

thetaB_coeff.HthetaB = HthetaB;
thetaB_coeff.JthetaB = JthetaB;
thetaB_coeff.TAthetaB = TAthetaB;
thetaB_coeff.TBthetaB = TBthetaB;

% Coefficients for Normal force at Wheel A

HNA = -R*phiA*(-(mA+mM)*(phiB*R^2*(mB*(PA*(yAx+yBx)*tMx^2+(yBx*yAxx*PA+PAx*(yAx*yBx-1))*tMx-tMxx*PA*(yAx*yBx-1))*l+(PA*((yAx*yc+xc)*yBx+xc*yAx-yc)*tMx^2+(yAxx*(xc*yBx-yc)*PA+PAx*((xc*yAx-yc)*yBx-yc*yAx-xc))*tMx-PA*tMxx*((xc*yAx-yc)*yBx-yc*yAx-xc))*mM)*PB^2-((-PA*tMx^2*yAx-PA*tMxx+PAx*tMx)*phiB-tMx*phiBx*PA)*(yBx^2+1)*I_B*l*PB-I_B*PBx*tMx*l*phiB*PA*(yBx^2+1))*l*cos(tM)^2+((mA+mM)*(phiB*(mB*((-yAx*yBx+1)*PA*tMx^2+(PA*yAxx+PAx*(yAx+yBx))*tMx-tMxx*PA*(yAx+yBx))*l+(-PA*((xc*yAx-yc)*yBx-yc*yAx-xc)*tMx^2+(yAxx*(yBx*yc+xc)*PA+((yAx*yc+xc)*yBx+xc*yAx-yc)*PAx)*tMx-PA*((yAx*yc+xc)*yBx+xc*yAx-yc)*tMxx)*mM)*R^2*PB^2+(yBx^2+1)*I_B*((tMx^2*PA+(PA*yAxx+PAx*yAx)*tMx-yAx*tMxx*PA)*phiB-tMx*yAx*phiBx*PA)*l*PB+I_B*PBx*tMx*yAx*l*phiB*PA*(yBx^2+1))*l*sin(tM)+phiB*R^2*(-tMx^3*mB*mM*(-yBx*yc+xc)*l^2+(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*l+(I_M+mM*(xc^2+yc^2))*mM*(yBx*yc+xc)*tMx^3)*PB^2+(-xc*tMx^3*l*phiB*mM+((I_M+mM*(xc^2+yc^2))*tMx^3-yAxx*PA^2*(mA+mM))*phiB+yAx*phiBx*PA^2*(mA+mM))*(yBx^2+1)*I_B*l*PB-I_B*PBx*yAx*l*phiB*PA^2*(yBx^2+1)*(mA+mM))*cos(tM)+((-tMx^3*mB*mM*(xc*yBx+yc)*l^2+(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*yBx*l+(I_M+mM*(xc^2+yc^2))*mM*(xc*yBx-yc)*tMx^3)*phiB*R^2*PB^2+(yBx^2+1)*(-tMx^2*mM*(phiB*tMx*yc+phiBx*xc)*l+phiBx*((I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM)))*I_B*l*PB-phiB*PBx*(yBx^2+1)*(-xc*tMx^2*l*mM+(I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM))*I_B*l)*sin(tM)-phiB*R^2*(mB*(mA+mM)*(-PA*tMx^2*yBx-PA*tMxx+PAx*tMx)*l^2-(PA*((yAx*yc+xc)*yBx*mM+(yc*(mA+mB)*yAx+xc*mA)*yBx+yc*mB)*tMx^2+((xc*yAxx*PA+PAx*(xc*yAx-yc))*yBx*mM+xc*yBx*yAxx*(mA+mB)*PA+((xc*(mA+mB)*yAx-yc*mA)*yBx+xc*mB)*PAx)*tMx-((xc*yAx-yc)*yBx*mM+(xc*(mA+mB)*yAx-yc*mA)*yBx+xc*mB)*PA*tMxx)*mM*l+((yBx*yAxx*PA+PAx*(yAx*yBx+1))*tMx-tMxx*PA*(yAx*yBx+1))*(I_M+mM*(xc^2+yc^2))*(mA+mB+mM))*PB^2-(yBx^2+1)*(((-PA*tMxx+PAx*tMx)*phiB-tMx*phiBx*PA)*(mA+mM)*l^2-mM*((PA*tMx^2*yc-PA*tMxx*xc+PAx*tMx*xc)*phiB+xc*tMx*phiBx*PA)*l+((-PA*tMxx+PAx*tMx)*phiB+tMx*phiBx*PA)*(I_M+mM*(xc^2+yc^2)))*I_B*PB+phiB*((-mA-mM)*l^2-xc*l*mM+mM*(xc^2+yc^2)+I_M)*PBx*(yBx^2+1)*PA*I_B*tMx)/(PB*phiB*((R^2*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*PB+I_B*l*(yBx^2+1)*(yAx-mu))*tMx*l*cos(tM)^2+(((mB*((yAx-mu)*yBx-yAx*mu-1)*l+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*mM)*R^2*PB-I_B*l*(yBx^2+1)*(mu*yAx+1))*tMx*sin(tM)+(R^2*((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+(R^2*(yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*PA*l*sin(tM)+(R^2*(-mB*(yAx-mu)*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*mM*yBx*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx))*PB-I_B*l^2*(yBx^2+1)*(yAx-mu))*tMx));

JNA = -((PB*(-yBx*mB*(mA+mM)*l^2-mM*(mA-mB+mM)*(xc*yBx-yc)*l+((xc^2-yc^2)*yBx-2*xc*yc)*mM^2)*R^2-yc*I_B*l*mM*(yBx^2+1))*tMx*cos(tM)^2+(-tMx*(PB*(-mB*(mA+mM)*l^2-mM*(mA-mB+mM)*(yBx*yc+xc)*l+mM^2*(2*xc*yBx*yc+xc^2-yc^2))*R^2+(yBx^2+1)*((-mA-mM)*l+mM*xc)*I_B*l)*sin(tM)+PA*(PB*((-mA-mM)*l+xc*mM*(yAx*yBx+1))*(mA+mB+mM)*R^2+(yBx^2+1)*((-mA-mM)*l+mM*xc)*I_B))*cos(tM)-PA*((yBx*(mA+mM)*l+yc*mM*(yAx*yBx+1))*PB*(mA+mB+mM)*R^2+yc*I_B*mM*(yBx^2+1))*sin(tM)-tMx*(PB*(-(xc*yBx*mM+xc*(mA+mB)*yBx+yc*mB)*mM*l+xc*(xc*yBx-yc)*mM^2+((mA+mB)*xc^2+(mA+mB)*yc^2+I_M)*yBx*mM+I_M*yBx*(mA+mB))*R^2-yc*I_B*l*mM*(yBx^2+1)))*R*phiA/((R^2*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*PB+I_B*l*(yBx^2+1)*(yAx-mu))*tMx*l*cos(tM)^2+(((mB*((yAx-mu)*yBx-yAx*mu-1)*l+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*mM)*R^2*PB-I_B*l*(yBx^2+1)*(mu*yAx+1))*tMx*sin(tM)+(R^2*((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+(R^2*(yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*PA*l*sin(tM)+(R^2*(-mB*(yAx-mu)*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*mM*yBx*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx))*PB-I_B*l^2*(yBx^2+1)*(yAx-mu))*tMx);

TANA = 0;

TBNA = -PB*(yBx^2+1)*R*phiA*(-yAx*l*PA*(mA+mM)*cos(tM)+l*PA*(mA+mM)*sin(tM)+(-xc*l*mM+mM*(xc^2+yc^2)+I_M)*tMx)/((((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*mM*(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-xc*mu-yc))*tMx*l*cos(tM)^2+(((mB*R^2*((yBx-mu)*yAx-yBx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-mu*yc+xc)*yBx-xc*mu-yc)*yAx+(-mu*xc-yc)*yBx+yc*mu-xc)*R^2*mM)*tMx*sin(tM)+(R^2*((mB*mu*yBx+mA+mM)*yAx^2+mB*(yBx+mu)*yAx+mA+mB+mM)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+PA*(R^2*(yBx*(mA+mB+mM)*yAx^2+(-mB*mu*yBx+mB)*yAx+(mA+mM)*yBx-mB*mu)*PB+I_B*(yBx^2+1)*(yAx-mu))*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx)*phiB);

FNA_coeff.HNA = HNA;
FNA_coeff.JNA = JNA;
FNA_coeff.TANA = TANA;
FNA_coeff.TBNA = TBNA;

% Coefficients for Normal force at Wheel B

HNB = R*((-PB*phiBx+PBx*phiB)*I_B*(mB*((mu*yBx+1)*yAx+yBx-mu)*l+(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*mM)*tMx^2*l^2*cos(tM)^3+((mB*((yBx-mu)*yAx-yBx*mu-1)*l+(((-yBx*yc-xc)*mu+xc*yBx-yc)*yAx+(-xc*yBx+yc)*mu-yc*yBx-xc)*mM)*(-PB*phiBx+PBx*phiB)*I_B*tMx^2*l*sin(tM)+(mB^2*R^2*phiB*(PA*(mu*yAx^2-mu+2*yAx)*tMx^2+(yAxx*(yAx-mu)*PA+PAx*(-2*mu*yAx+yAx^2-1))*tMx-tMxx*PA*(-2*mu*yAx+yAx^2-1))*PB^2+I_B*(-PA*phiB*(yBx*(yAx^2+1)*mM+(mA*yBx-mB*mu)*yAx^2-2*yAx*mB+yBx*mA+mB*mu)*tMx^2+(((yAxx*(yAx-mu)*PA+PAx*(yAx^2+1))*mM+yAxx*(mA+mB)*(yAx-mu)*PA+((mA+mB)*yAx^2-2*yAx*mB*mu-mB+mA)*PAx)*phiB-phiBx*PA*((yAx^2+1)*mM+(mA+mB)*yAx^2-2*yAx*yBx*mB+2*yBx*mB*mu+mA+mB))*tMx-((yAx^2+1)*mM+(mA+mB)*yAx^2-2*yAx*mB*mu-mB+mA)*tMxx*PA*phiB)*PB+PBx*PA*phiB*I_B*((yAx^2+1)*mM+(mA+mB)*yAx^2-2*yAx*yBx*mB+2*yBx*mB*mu+mA+mB)*tMx)*l+mM*((PA*((mu*xc+yc)*yAx^2+(-2*mu*yc+2*xc)*yAx-xc*mu-yc)*tMx^2+(((-mu*yc+xc)*yAx-xc*mu-yc)*yAxx*PA+PAx*((-mu*yc+xc)*yAx^2+(-2*mu*xc-2*yc)*yAx+yc*mu-xc))*tMx-tMxx*PA*((-mu*yc+xc)*yAx^2+(-2*mu*xc-2*yc)*yAx+yc*mu-xc))*mB*R^2*phiB*PB^2+(PA*phiB*(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*tMx^2+(((-yBx*yc-xc)*mu+xc*yBx-yc)*yAx+(-xc*yBx+yc)*mu-yc*yBx-xc)*(PA*phiBx+PAx*phiB)*tMx-tMxx*(((-yBx*yc-xc)*mu+xc*yBx-yc)*yAx+(-xc*yBx+yc)*mu-yc*yBx-xc)*PA*phiB)*I_B*PB-(((-yBx*yc-xc)*mu+xc*yBx-yc)*yAx+(-xc*yBx+yc)*mu-yc*yBx-xc)*PBx*PA*phiB*I_B*tMx))*l*cos(tM)^2+(-((mB^2*R^2*phiB*(-PA*(-2*mu*yAx+yAx^2-1)*tMx^2+(yAxx*(mu*yAx+1)*PA+PAx*(mu*yAx^2-mu+2*yAx))*tMx-tMxx*PA*(mu*yAx^2-mu+2*yAx))*PB^2-I_B*(((yAx^2+1)*mM+(mA+mB)*yAx^2-2*yAx*mB*mu-mB+mA)*PA*phiB*tMx^2+(((yAxx*(yAx-mu)*PA+PAx*(yAx^2+1))*yBx*mM+((mA*yBx-mB*mu)*yAx-yBx*mA*mu-mB)*yAxx*PA+((mA*yBx-mB*mu)*yAx^2-2*yAx*mB+yBx*mA+mB*mu)*PAx)*phiB-phiBx*PA*(yBx*(yAx^2+1)*mM+(mA*yBx-mB*mu)*yAx^2+2*yAx*yBx*mB*mu-mB*mu+yBx*(mA+2*mB)))*tMx-tMxx*PA*phiB*(yBx*(yAx^2+1)*mM+(mA*yBx-mB*mu)*yAx^2-2*yAx*mB+yBx*mA+mB*mu))*PB-PBx*PA*phiB*(yBx*(yAx^2+1)*mM+(mA*yBx-mB*mu)*yAx^2+2*yAx*yBx*mB*mu-mB*mu+yBx*(mA+2*mB))*I_B*tMx)*l+(mB*R^2*phiB*(-PA*((-mu*yc+xc)*yAx^2+(-2*mu*xc-2*yc)*yAx+yc*mu-xc)*tMx^2+(((mu*xc+yc)*yAx-yc*mu+xc)*yAxx*PA+((mu*xc+yc)*yAx^2+(-2*mu*yc+2*xc)*yAx-xc*mu-yc)*PAx)*tMx-tMxx*PA*((mu*xc+yc)*yAx^2+(-2*mu*yc+2*xc)*yAx-xc*mu-yc))*PB^2+I_B*(-(((-yBx*yc-xc)*mu+xc*yBx-yc)*yAx+(-xc*yBx+yc)*mu-yc*yBx-xc)*PA*phiB*tMx^2+(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*(PA*phiBx+PAx*phiB)*tMx-tMxx*PA*phiB*(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc))*PB-PBx*PA*phiB*(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*I_B*tMx)*mM)*l*sin(tM)+mB*(R^2*tMx*phiB*mB*(mu*yAx+1)*PB^2+I_B*(phiB*(mu*yAx+1)*tMx+phiBx*((mu*yBx+1)*yAx+yBx-mu))*PB-I_B*PBx*phiB*(mu*yAx*yBx-mu+yAx+yBx))*tMx^2*l^3+(mB*((mu*xc+yc)*yAx-yc*mu+xc)*R^2*phiB*tMx*PB^2-(phiB*((-mu*yc+xc)*yAx-xc*mu-yc)*yBx*tMx-phiBx*(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc))*I_B*PB-PBx*phiB*(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*I_B)*mM*tMx^2*l^2+(mB*(mu*yAx+1)*R^2*phiB*((I_M+mM*(xc^2+yc^2))*tMx^3+yAxx*PA^2*(mA+mB+mM))*PB^2+(((yBx+mu)*yAx-yBx*mu+1)*phiB*(I_M+mM*(xc^2+yc^2))*tMx^3+(-yAxx*(yBx*(yAx-mu)*mM+(mA*yBx-mB*mu)*yAx-yBx*mA*mu-mB)*phiB+phiBx*(yBx*(yAx^2+1)*mM+(mA*yBx-mB*mu)*yAx^2+mB*(mu*yBx-1)*yAx+yBx*(mA+mB)))*PA^2)*I_B*PB-PBx*PA^2*phiB*(yBx*(yAx^2+1)*mM+(mA*yBx-mB*mu)*yAx^2+mB*(mu*yBx-1)*yAx+yBx*(mA+mB))*I_B)*l+PB^2*((mu*xc+yc)*yAx-yc*mu+xc)*R^2*phiB*(I_M+mM*(xc^2+yc^2))*mM*tMx^3)*cos(tM)+((yAx-mu)*mB*(R^2*tMx*phiB*PB^2*mB+I_B*(phiB*tMx+phiBx*yBx)*PB-I_B*PBx*yBx*phiB)*tMx^2*l^3+(mB*R^2*phiB*((-mu*yc+xc)*yAx-xc*mu-yc)*tMx*PB^2+(phiB*((-mu*yc+xc)*yAx-xc*mu-yc)*tMx-phiBx*((mu*xc+yc)*yAx-yc*mu+xc))*I_B*PB+PBx*((mu*xc+yc)*yAx-yc*mu+xc)*phiB*I_B)*mM*tMx^2*l^2+((yAx-mu)*mB*R^2*phiB*((I_M+mM*(xc^2+yc^2))*tMx^3+yAxx*PA^2*(mA+mB+mM))*PB^2+(phiBx*((yBx+mu)*yAx-yBx*mu+1)*(I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(-yAxx*(mA+mB+mM)*(yAx-mu)*phiB+((yAx^2+1)*mM+(mA+mB)*yAx^2-mB*(yBx+mu)*yAx+yBx*mB*mu+mA)*phiBx))*I_B*PB-PBx*(((yBx+mu)*yAx-yBx*mu+1)*(I_M+mM*(xc^2+yc^2))*tMx^2-((yAx^2+1)*mM+(mA+mB)*yAx^2-mB*(yBx+mu)*yAx+yBx*mB*mu+mA)*PA^2)*phiB*I_B)*l+PB^2*R^2*phiB*(I_M+mM*(xc^2+yc^2))*((-mu*yc+xc)*yAx-xc*mu-yc)*mM*tMx^3)*sin(tM)+(-mB*R^2*(PA*mB*(yAx-mu)*tMx^2+((yAxx*(yAx-mu)*PA+PAx*(yAx^2+1))*mM+yAxx*(mA+mB)*(yAx-mu)*PA+((mA+mB)*yAx^2-yAx*mB*mu+mA)*PAx)*tMx-tMxx*PA*((yAx^2+1)*mM+(mA+mB)*yAx^2-yAx*mB*mu+mA))*phiB*PB^2-I_B*(phiB*PA*mB*(yAx-mu)*tMx^2+(((yAxx*(yAx-mu)*PA+PAx*(yAx^2+1))*mM+yAxx*(mA+mB)*(yAx-mu)*PA+((mA+mB)*yAx^2-yAx*mB*mu+mA)*PAx)*phiB-phiBx*PA*((yAx^2+1)*mM+(mA+mB)*yAx^2-2*mB*(yBx+(1/2)*mu)*yAx+2*yBx*mB*mu+mA))*tMx-tMxx*PA*phiB*((yAx^2+1)*mM+(mA+mB)*yAx^2-yAx*mB*mu+mA))*PB-PBx*PA*phiB*((yAx^2+1)*mM+(mA+mB)*yAx^2-2*mB*(yBx+(1/2)*mu)*yAx+2*yBx*mB*mu+mA)*I_B*tMx)*l^2+(mB*R^2*(-PA*((-mu*yc+xc)*yAx-xc*mu-yc)*tMx^2+((mu*xc+yc)*yAx-yc*mu+xc)*PAx*tMx-tMxx*((mu*xc+yc)*yAx-yc*mu+xc)*PA)*phiB*PB^2+(-PA*phiB*((-mu*yc+xc)*yAx-xc*mu-yc)*tMx^2+((mu*xc+yc)*yAx-yc*mu+xc)*(PA*phiBx+PAx*phiB)*tMx-tMxx*((mu*xc+yc)*yAx-yc*mu+xc)*PA*phiB)*I_B*PB-PBx*((mu*xc+yc)*yAx-yc*mu+xc)*PA*phiB*I_B*tMx)*mM*l-(I_M+mM*(xc^2+yc^2))*(((yAxx*(yAx-mu)*PA+PAx*(yAx^2+1))*tMx-tMxx*PA*(yAx^2+1))*R^2*phiB*(mA+mB+mM)*PB^2+((yBx+mu)*yAx-yBx*mu+1)*I_B*((PA*phiBx+PAx*phiB)*tMx-phiB*tMxx*PA)*PB-PBx*((yBx+mu)*yAx-yBx*mu+1)*PA*phiB*I_B*tMx))/(PB*(((mB*R^2*((mu*yBx+1)*yAx+yBx-mu)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*mM)*tMx*l*cos(tM)^2+(((mB*R^2*((yBx-mu)*yAx-yBx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-yBx*yc-xc)*mu+xc*yBx-yc)*yAx+(-xc*yBx+yc)*mu-yc*yBx-xc)*R^2*mM)*tMx*sin(tM)+(R^2*((yAx^2+1)*mM+(mB*mu*yBx+mA)*yAx^2+mB*(yBx+mu)*yAx+mA+mB)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+(R^2*(yBx*(yAx^2+1)*mM+yBx*(mA+mB)*yAx^2+(-mB*mu*yBx+mB)*yAx+yBx*mA-mB*mu)*PB+I_B*(yBx^2+1)*(yAx-mu))*PA*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx));

JNB = R*phiB*(tMx*((R^2*PB*mB^2+I_B*(mA+mB+mM))*(yAx-mu)*l^2+(2*(PB*mB*((-mu*yc+xc)*yAx-xc*mu-yc)*R^2-(1/2)*((yBx+mu)*yAx-yBx*mu+1)*I_B*yc))*mM*l+PB*R^2*mM^2*((-2*mu*xc*yc+xc^2-yc^2)*yAx+(-xc^2+yc^2)*mu-2*xc*yc))*cos(tM)^2+(-tMx*((PB*mB^2*(mu*yAx+1)*R^2-(yBx*(yAx-mu)*mM+(mA*yBx-mB*mu)*yAx-yBx*mA*mu-mB)*I_B)*l^2+2*mM*(PB*mB*((mu*xc+yc)*yAx-yc*mu+xc)*R^2+(1/2)*((yBx+mu)*yAx-yBx*mu+1)*I_B*xc)*l+(((xc^2-yc^2)*mu+2*xc*yc)*yAx-2*xc*yc*mu+xc^2-yc^2)*PB*R^2*mM^2)*sin(tM)+((PB*mB*(mA+mB+mM)*(mu*yAx+1)*R^2-(yBx*(yAx-mu)*mM+(mA*yBx-mB*mu)*yAx-yBx*mA*mu-mB)*I_B)*l+(R^2*(yAx^2+1)*(mA+mB+mM)*PB+((yBx+mu)*yAx-yBx*mu+1)*I_B)*mM*xc)*PA)*cos(tM)-PA*(-(mA+mB+mM)*(yAx-mu)*(PB*R^2*mB+I_B)*l+(R^2*(yAx^2+1)*(mA+mB+mM)*PB+((yBx+mu)*yAx-yBx*mu+1)*I_B)*mM*yc)*sin(tM)-tMx*((mA+mB+mM)*(yAx-mu)*(PB*R^2*mB+I_B)*l^2-(PB*mB*(mu*yAx+1)*R^2+((yBx+mu)*yAx-yBx*mu+1)*I_B)*mM*yc*l+PB*R^2*(((-mu*yc+xc)*yAx-xc*mu-yc)*xc*mM^2+(yAx-mu)*((xc^2+yc^2)*mB+xc^2*mA+yc^2*mA+I_M)*mM+I_M*(mA+mB)*(yAx-mu))))/(((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+PB*R^2*(((xc*yBx-yc)*mu+yc*yBx+xc)*yAx+(-yBx*yc-xc)*mu+xc*yBx-yc)*mM)*tMx*l*cos(tM)^2+(((mB*R^2*((yBx-mu)*yAx-yBx*mu-1)*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+PB*(((-yBx*yc-xc)*mu+xc*yBx-yc)*yAx+(-xc*yBx+yc)*mu-yc*yBx-xc)*R^2*mM)*tMx*sin(tM)+(R^2*((yAx^2+1)*mM+(mB*mu*yBx+mA)*yAx^2+mB*(yBx+mu)*yAx+mA+mB)*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+(R^2*(yBx*(yAx^2+1)*mM+yBx*(mA+mB)*yAx^2+(-mB*mu*yBx+mB)*yAx+yBx*mA-mB*mu)*PB+I_B*(yBx^2+1)*(yAx-mu))*PA*l*sin(tM)+(-(yAx-mu)*(R^2*PB*mB+I_B*(yBx^2+1))*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*(mu*yAx*yBx+mu-yAx+yBx)*(I_M+mM*(xc^2+yc^2)))*tMx);

TANB = 0;

TBNB = PB*R*(-tMx*(mB*((yBx-mu)*yAx-yBx*mu-1)*l+mM*(((-mu*yc+xc)*yBx-xc*mu-yc)*yAx+(-mu*xc-yc)*yBx+yc*mu-xc))*l*cos(tM)^2+((mB*((mu*yBx+1)*yAx+yBx-mu)*l+mM*(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-xc*mu-yc))*tMx*sin(tM)-PA*(((mA+mM)*yBx-mB*mu)*yAx^2+mB*(mu*yBx-1)*yAx+yBx*(mA+mB+mM)))*l*cos(tM)+((mA+mB+mM)*yAx^2-mB*(yBx+mu)*yAx+yBx*mB*mu+mA+mM)*PA*l*sin(tM)+(yBx*mB*(yAx-mu)*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*mM*l+((yBx+mu)*yAx-yBx*mu+1)*(I_M+mM*(xc^2+yc^2)))*tMx)/((((PB*R^2*mB*mu*yBx+PB*R^2*mB+I_B*yBx^2+I_B)*yAx-I_B*mu*yBx^2+R^2*yBx*PB*mB-mu*(PB*R^2*mB+I_B))*l+PB*R^2*mM*(((mu*xc+yc)*yBx-yc*mu+xc)*yAx+(-mu*yc+xc)*yBx-xc*mu-yc))*tMx*l*cos(tM)^2+((((-I_B*mu*yBx^2+R^2*yBx*PB*mB-mu*(PB*R^2*mB+I_B))*yAx-R^2*yBx*PB*mB*mu-R^2*PB*mB-I_B*yBx^2-I_B)*l+PB*(((-mu*yc+xc)*yBx-xc*mu-yc)*yAx+(-mu*xc-yc)*yBx+yc*mu-xc)*R^2*mM)*tMx*sin(tM)+(R^2*PB*(mB*mu*yBx+mA+mM)*yAx^2+(I_B*mu*yBx^2+R^2*yBx*PB*mB+mu*(PB*R^2*mB+I_B))*yAx+I_B*yBx^2+PB*R^2*mM+R^2*(mA+mB)*PB+I_B)*PA)*l*cos(tM)+PA*(R^2*yBx*PB*(mA+mB+mM)*yAx^2+(-PB*R^2*mB*mu*yBx+PB*R^2*mB+I_B*yBx^2+I_B)*yAx-I_B*mu*yBx^2+R^2*PB*(mA+mM)*yBx-mu*(PB*R^2*mB+I_B))*l*sin(tM)+tMx*(-(yAx-mu)*(PB*R^2*mB+I_B*yBx^2+I_B)*l^2-((mu*xc+yc)*yAx-yc*mu+xc)*PB*R^2*mM*yBx*l+PB*R^2*((mu*yBx-1)*yAx+yBx+mu)*(I_M+mM*(xc^2+yc^2))));

FNB_coeff.HNB = HNB;
FNB_coeff.JNB = JNB;
FNB_coeff.TANB = TANB;
FNB_coeff.TBNB = TBNB;


