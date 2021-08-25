function [x_coeff, thetaA_coeff, thetaB_coeff, FNA_coeff, FNB_coeff, tM, dx] = noslip_eom(xA,guess_Dx,params)
%% noslip_eom
% Levi Manring, Duke University
% 2021
%
% This function calculates certain parameters needed for integration of the
% noslip dynamics model for the AWD complete model.
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
%       Outputs are elements necessary to solve the no-slip vehicle dynamics
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
phiAx = yAx*yAxx/(R*sqrt(yAx^2+1));

phiB = sqrt(yBx^2+1)/R;
phiBx = yBx*yBxx/(R*sqrt(yBx^2+1));


%% Calculating noslip dynamic model parameters
% Applying torque at wheel A and B and including moment of inertia effects for both wheels

% Coefficients for x

Hx = (phiA*((-R^2*phiB*mB*((tMx^2*yBx+tMxx)*yAx-tMx^2+yBx*tMxx)*PB^2+(yBx^2+1)*((tMx^2-tMxx*yAx)*phiB-tMx*yAx*phiBx)*I_B*PB+I_B*PBx*tMx*yAx*phiB*(yBx^2+1))*l+R^2*phiB*PB^2*(((-tMx^2*xc-tMxx*yc)*yBx+yc*tMx^2-tMxx*xc)*yAx+(tMx^2*yc-tMxx*xc)*yBx+tMx^2*xc+yc*tMxx)*mM)*l*cos(tM)^2+l*(((R^2*phiB*((tMx^2-tMxx*yBx)*yAx+tMx^2*yBx+tMxx)*mB*PB^2+(yBx^2+1)*((tMx^2*yAx+tMxx)*phiB+tMx*phiBx)*I_B*PB-I_B*PBx*tMx*phiB*(yBx^2+1))*l+R^2*(((tMx^2*yc-tMxx*xc)*yBx+tMx^2*xc+yc*tMxx)*yAx+(tMx^2*xc+tMxx*yc)*yBx-yc*tMx^2+tMxx*xc)*phiB*PB^2*mM)*phiA*sin(tM)-(R^2*(PAx*(mA+mM)*yAx^2+(PAx*yBx*mB+yAxx*PA*(mA+mM))*yAx+yBx*yAxx*PA*mB+PAx*(mA+mB+mM))*phiA+I_A*phiAx*(yAx^2+1))*phiB*PB^2-I_B*phiA*(yBx^2+1)*(PA*phiBx+PAx*phiB)*PB+I_B*PBx*phiA*phiB*PA*(yBx^2+1))*cos(tM)-(phiB*(R^2*(PAx*yBx*(mA+mB+mM)*yAx^2+(yAxx*PA*(mA+mB+mM)*yBx+PAx*mB)*yAx+PAx*yBx*(mA+mM))*phiA+I_A*yBx*phiAx*(yAx^2+1))*PB^2+I_B*yAx*phiA*(yBx^2+1)*(PA*phiBx+PAx*phiB)*PB-I_B*PBx*yAx*phiA*phiB*PA*(yBx^2+1))*l*sin(tM)+phiA*(yAx*(R^2*phiB*mB*(tMx^2*yBx+tMxx)*PB^2+I_B*(yBx^2+1)*(phiB*tMxx+phiBx*tMx)*PB-I_B*PBx*tMx*phiB*(yBx^2+1))*l^2-R^2*phiB*PB^2*yBx*((-tMx^2*xc-tMxx*yc)*yAx+yc*tMx^2-tMxx*xc)*mM*l+R^2*phiB*PB^2*(yAx-yBx)*tMxx*(I_M+mM*(xc^2+yc^2))))/(phiA*phiB*PB*(-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx));

Jx = -((((mA+mM)*yAx+mB*yBx)*l-mM*xc*(yAx-yBx))*cos(tM)+((mA+mB+mM)*yBx*yAx*l+mM*yc*(yAx-yBx))*sin(tM))*PB*R^2/(-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx);

TAx = PB*l*(yAx^2+1)*(yBx*sin(tM)+cos(tM))/(phiA*(-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx));

TBx = PB*l*(yBx^2+1)*(yAx*sin(tM)+cos(tM))/((-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*phiB);

x_coeff.Hx = Hx;
x_coeff.Jx = Jx;
x_coeff.TAx = TAx;
x_coeff.TBx = TBx;

% Coefficients for thetaA

HthetaA = (-l*((-mB*phiB*R^2*((phiA*tMx^2*yBx+phiA*tMxx-phiAx*tMx)*yAx+(phiA*tMxx-phiAx*tMx)*yBx-tMx^2*phiA)*PB^2+(((-phiA*tMxx+phiAx*tMx)*yAx+tMx^2*phiA)*phiB-tMx*yAx*phiA*phiBx)*(yBx^2+1)*I_B*PB+I_B*PBx*tMx*yAx*phiA*phiB*(yBx^2+1))*l+phiB*PB^2*mM*((((-tMx^2*xc-tMxx*yc)*phiA+yc*tMx*phiAx)*yBx+(tMx^2*yc-tMxx*xc)*phiA+tMx*phiAx*xc)*yAx+((tMx^2*yc-tMxx*xc)*phiA+tMx*phiAx*xc)*yBx+(tMx^2*xc+tMxx*yc)*phiA-yc*tMx*phiAx)*R^2)*cos(tM)^2-((((((-phiA*tMxx+phiAx*tMx)*yBx+tMx^2*phiA)*yAx+tMx^2*yBx*phiA-tMx*phiAx+phiA*tMxx)*mB*phiB*R^2*PB^2+((phiA*tMx^2*yAx+phiA*tMxx-phiAx*tMx)*phiB+tMx*phiA*phiBx)*(yBx^2+1)*I_B*PB-I_B*PBx*tMx*phiA*phiB*(yBx^2+1))*l+((((tMx^2*yc-tMxx*xc)*phiA+tMx*phiAx*xc)*yBx+(tMx^2*xc+tMxx*yc)*phiA-yc*tMx*phiAx)*yAx+((tMx^2*xc+tMxx*yc)*phiA-yc*tMx*phiAx)*yBx+(-tMx^2*yc+tMxx*xc)*phiA-tMx*phiAx*xc)*phiB*PB^2*mM*R^2)*sin(tM)-phiB*((mA+mM)*(-PA*phiAx+PAx*phiA)*yAx^2+(mB*(-PA*phiAx+PAx*phiA)*yBx+phiA*yAxx*PA*(mA+mM))*yAx+yBx*phiA*yAxx*PA*mB+(mA+mB+mM)*(-PA*phiAx+PAx*phiA))*R^2*PB^2-((-PA*phiAx+PAx*phiA)*phiB+phiA*phiBx*PA)*(yBx^2+1)*I_B*PB+I_B*PBx*phiA*phiB*PA*(yBx^2+1))*l*cos(tM)+(phiB*R^2*(yBx*(mA+mB+mM)*(-PA*phiAx+PAx*phiA)*yAx^2+(phiA*yAxx*PA*(mA+mB+mM)*yBx+mB*(-PA*phiAx+PAx*phiA))*yAx+yBx*(mA+mM)*(-PA*phiAx+PAx*phiA))*PB^2+((-PA*phiAx+PAx*phiA)*phiB+phiA*phiBx*PA)*(yBx^2+1)*I_B*yAx*PB-I_B*PBx*yAx*phiA*phiB*PA*(yBx^2+1))*l*sin(tM)-yAx*(mB*phiB*R^2*(phiA*tMx^2*yBx+phiA*tMxx-phiAx*tMx)*PB^2+((phiA*tMxx-phiAx*tMx)*phiB+tMx*phiA*phiBx)*(yBx^2+1)*I_B*PB-I_B*PBx*tMx*phiA*phiB*(yBx^2+1))*l^2+phiB*PB^2*yBx*mM*(((-tMx^2*xc-tMxx*yc)*phiA+yc*tMx*phiAx)*yAx+(tMx^2*yc-tMxx*xc)*phiA+tMx*phiAx*xc)*R^2*l+(I_M+mM*(xc^2+yc^2))*phiB*PB^2*(-phiA*tMxx+phiAx*tMx)*(yAx-yBx)*R^2)/(phiB*(-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*PB);

JthetaA = phiA*((((mA+mM)*yAx+mB*yBx)*l-mM*xc*(yAx-yBx))*cos(tM)+((mA+mB+mM)*yBx*yAx*l+mM*yc*(yAx-yBx))*sin(tM))*PB*R^2/(-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx);

TAthetaA = -PB*l*(yAx^2+1)*(yBx*sin(tM)+cos(tM))/(-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx);

TBthetaA = -phiA*PB*l*(yBx^2+1)*(yAx*sin(tM)+cos(tM))/((-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*phiB);

thetaA_coeff.HthetaA = HthetaA;
thetaA_coeff.JthetaA = JthetaA;
thetaA_coeff.TAthetaA = TAthetaA;
thetaA_coeff.TBthetaA = TBthetaA;

% Coefficients for thetaB

HthetaB = (-l^2*((-mB*yAx*yBx+mB)*l+mM*((-xc*yAx+yc)*yBx+yc*yAx+xc))*tMx^2*(-PB*phiBx+PBx*phiB)*R^2*phiA*cos(tM)^3-l*(l*tMx^2*(-PB*phiBx+PBx*phiB)*R^2*phiA*(mB*(yAx+yBx)*l+((yAx*yc+xc)*yBx+yAx*xc-yc)*mM)*sin(tM)+((((-phiB*(2*yAx*yBx*mB+(yAx^2+1)*(mA+mM))*PA*tMx^2+(-yBx*(PAx*(mA+mB+mM)*yAx^2+yAxx*PA*(mA+mB+mM)*yAx+PAx*(mA-mB+mM))*phiB+((yAx^2+1)*(mA+mB+mM)*yBx+2*yAx*mB)*phiBx*PA)*tMx+phiB*tMxx*yBx*((mA+mB+mM)*yAx^2+mA-mB+mM)*PA)*PB-PBx*phiB*tMx*((yAx^2+1)*(mA+mB+mM)*yBx+2*yAx*mB)*PA)*R^2-I_A*(yAx^2+1)*((phiB*tMx^2-phiB*tMxx*yBx-phiBx*tMx*yBx)*PB+PBx*tMx*yBx*phiB))*phiA-I_A*tMx*yBx*phiB*phiAx*PB*(yAx^2+1))*l+((phiB*((-xc*yAx+yc)*yBx+yc*yAx+xc)*PA*tMx^2+((yAx*yc+xc)*yBx+yAx*xc-yc)*(PA*phiBx+PAx*phiB)*tMx-((yAx*yc+xc)*yBx+yAx*xc-yc)*phiB*tMxx*PA)*PB-PBx*((yAx*yc+xc)*yBx+yAx*xc-yc)*phiB*tMx*PA)*mM*R^2*phiA)*cos(tM)^2-l*((((((-phiB*yBx*((mA+mB+mM)*yAx^2+mA-mB+mM)*PA*tMx^2+(((PA*mB*yAxx+2*PAx*mB*yAx)*yBx+(mA+mM)*(PA*yAx*yAxx+PAx*yAx^2+PAx))*phiB-((mA+mM)*yAx^2+mA+2*mB+mM)*phiBx*PA)*tMx-phiB*tMxx*(2*yAx*yBx*mB+(yAx^2+1)*(mA+mM))*PA)*PB+PBx*((mA+mM)*yAx^2+mA+2*mB+mM)*phiB*tMx*PA)*R^2+I_A*(yAx^2+1)*((-phiB*tMx^2*yBx-phiB*tMxx-phiBx*tMx)*PB+PBx*phiB*tMx))*phiA+I_A*tMx*phiB*phiAx*PB*(yAx^2+1))*l-((-((yAx*yc+xc)*yBx+yAx*xc-yc)*phiB*PA*tMx^2+(PA*phiBx+PAx*phiB)*((-xc*yAx+yc)*yBx+yc*yAx+xc)*tMx-phiB*tMxx*((-xc*yAx+yc)*yBx+yc*yAx+xc)*PA)*PB-PBx*phiB*tMx*((-xc*yAx+yc)*yBx+yc*yAx+xc)*PA)*mM*R^2*phiA)*sin(tM)+((-phiB*tMx*yBx-phiBx*yAx*yBx+phiBx)*PB+PBx*phiB*(yAx*yBx-1))*mB*tMx^2*R^2*phiA*l^2-mM*tMx^2*((-tMx*(-xc*yAx+yc)*phiB-phiBx*((-xc*yAx+yc)*yBx+yc*yAx+xc))*PB+PBx*phiB*((-xc*yAx+yc)*yBx+yc*yAx+xc))*R^2*phiA*l+((((I_M+mM*(xc^2+yc^2))*phiB*(yAx-yBx)*tMx^3+(-((mA+mM)*yAx+mB*yBx)*yAxx*phiB+((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*phiBx)*PA^2)*PB-PBx*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*phiB*PA^2)*R^2+((PA*phiBx+PAx*phiB)*PB-PBx*phiB*PA)*I_A*(yAx^2+1))*phiA-I_A*phiB*phiAx*PA*PB*(yAx^2+1))*cos(tM)+l*(((phiB*tMx*yBx-phiBx)*PB+PBx*phiB)*mB*tMx^2*R^2*phiA*yAx*l^2+yBx*mM*tMx^2*R^2*phiA*((-tMx*(-xc*yAx+yc)*phiB-phiBx*(yAx*yc+xc))*PB+PBx*phiB*(yAx*yc+xc))*l+(((-(I_M+mM*(xc^2+yc^2))*phiBx*(yAx-yBx)*tMx^2-(-yAx*yBx*yAxx*(mA+mB+mM)*phiB+(((mA+mB+mM)*yAx^2+mA+mM)*yBx+yAx*mB)*phiBx)*PA^2)*PB+PBx*((I_M+mM*(xc^2+yc^2))*(yAx-yBx)*tMx^2+(((mA+mB+mM)*yAx^2+mA+mM)*yBx+yAx*mB)*PA^2)*phiB)*R^2-((PA*phiBx+PAx*phiB)*PB-PBx*phiB*PA)*I_A*(yAx^2+1)*yBx)*phiA+I_A*yBx*phiB*phiAx*PA*PB*(yAx^2+1))*sin(tM)+((((-tMx^2*yAx*yBx*phiB*PA*mB+(-(PAx*(mA+mB+mM)*yAx^2+yAxx*PA*(mA+mB+mM)*yAx+PAx*(mA+mM))*yBx*phiB+phiBx*PA*(((mA+mB+mM)*yAx^2+mA+mM)*yBx+2*yAx*mB))*tMx+phiB*tMxx*yBx*PA*((mA+mB+mM)*yAx^2+mA+mM))*PB-PBx*phiB*tMx*PA*(((mA+mB+mM)*yAx^2+mA+mM)*yBx+2*yAx*mB))*R^2-((-phiB*tMxx-phiBx*tMx)*PB+PBx*phiB*tMx)*I_A*(yAx^2+1)*yBx)*phiA-I_A*tMx*yBx*phiB*phiAx*PB*(yAx^2+1))*l^2+yBx*mM*((phiB*PA*(-xc*yAx+yc)*tMx^2+(PA*phiBx+PAx*phiB)*(yAx*yc+xc)*tMx-phiB*tMxx*PA*(yAx*yc+xc))*PB-PBx*tMx*phiB*PA*(yAx*yc+xc))*R^2*phiA*l+(I_M+mM*(xc^2+yc^2))*(yAx-yBx)*R^2*phiA*(((PA*phiBx+PAx*phiB)*tMx-phiB*tMxx*PA)*PB-PBx*phiB*tMx*PA))/((-l*tMx*((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+((yAx*yc+xc)*yBx+yAx*xc-yc)*PB*mM*R^2)*cos(tM)^2+l*(tMx*((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+PB*mM*R^2*((-xc*yAx+yc)*yBx+yc*yAx+xc))*sin(tM)-((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*PB*PA*R^2+(-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*cos(tM)-l*((((mA+mB+mM)*yAx^2+mA+mM)*yBx+yAx*mB)*PB*PA*R^2+I_A*yBx*(yAx^2+1)*PB+I_B*yAx*PA*(yBx^2+1))*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*PB*phiA);

JthetaB = (l*((mA+mB+mM)*yBx*yAx*l+mM*yc*(yAx-yBx))*tMx*cos(tM)^2-(((mA+mM)*yAx+mB*yBx)*l-mM*xc*(yAx-yBx))*(l*tMx*sin(tM)-PA)*cos(tM)-((mA+mB+mM)*yBx*yAx*l+mM*yc*(yAx-yBx))*(l*tMx-sin(tM)*PA))*R^2*phiB/(-l*((yAx*(PB*R^2*mB+I_B*yBx^2+I_B)+R^2*yBx*PB*mB)*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*cos(tM)^2+l*(((-PB*R^2*mB*yAx*yBx+PB*R^2*mB+I_B*yBx^2+I_B)*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)-PB*(PA*(mA+mM)*R^2+I_A)*yAx^2-PA*PB*R^2*mB*yAx*yBx-I_B*PA*yBx^2+(-PA*(mA+mB+mM)*R^2-I_A)*PB-I_B*PA)*cos(tM)-l*((PA*(mA+mB+mM)*R^2+I_A)*PB*yBx*yAx^2+PA*(PB*R^2*mB+I_B*yBx^2+I_B)*yAx+PB*yBx*(PA*(mA+mM)*R^2+I_A))*sin(tM)+tMx*(yAx*(PB*R^2*mB+I_B*yBx^2+I_B)*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2))));

TAthetaB = -(cos(tM)^2*l*tMx*yBx+(PA-l*tMx*sin(tM))*cos(tM)-yBx*(l*tMx-sin(tM)*PA))*l*phiB*(yAx^2+1)/((-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*phiA);

TBthetaB = -(cos(tM)^2*l*tMx*yAx+(PA-l*tMx*sin(tM))*cos(tM)-yAx*(l*tMx-sin(tM)*PA))*l*(yBx^2+1)/(-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx);

thetaB_coeff.HthetaB = HthetaB;
thetaB_coeff.JthetaB = JthetaB;
thetaB_coeff.TAthetaB = TAthetaB;
thetaB_coeff.TBthetaB = TBthetaB;

% Coefficients for Normal force at Wheel A

HNA = (l*(phiB*((mA+mM)*(-mB*(PA*(yAx+yBx)*tMx^2+((PA*yAxx+PAx*yAx)*yBx-PAx)*tMx-PA*tMxx*(yAx*yBx-1))*l+mM*(-((yAx*yc+xc)*yBx+yAx*xc-yc)*PA*tMx^2+((-PA*yAxx*xc+PAx*(-xc*yAx+yc))*yBx+yc*yAxx*PA+PAx*(yAx*yc+xc))*tMx-tMxx*((-xc*yAx+yc)*yBx+yc*yAx+xc)*PA))*phiA*R^2-I_A*((mB*((yAx+yBx)*tMx^2-yAx*yBx*tMxx+tMxx)*l+(((yAx*yc+xc)*yBx+yAx*xc-yc)*tMx^2+tMxx*((-xc*yAx+yc)*yBx+yc*yAx+xc))*mM)*phiA-phiAx*((-mB*yAx*yBx+mB)*l+mM*((-xc*yAx+yc)*yBx+yc*yAx+xc))*tMx))*R^2*PB^2+((mA+mM)*phiA*((-PA*tMx^2*yAx-PA*tMxx+PAx*tMx)*phiB-tMx*phiBx*PA)*R^2-I_A*(((tMx^2*yAx+tMxx)*phiB+tMx*phiBx)*phiA-tMx*phiB*phiAx))*l*(yBx^2+1)*I_B*PB+PBx*l*phiB*tMx*(PA*(mA+mM)*R^2+I_A)*(yBx^2+1)*phiA*I_B)*cos(tM)^2+(l*(((mA+mM)*(mB*((-PA*yAx*yBx+PA)*tMx^2+(PA*yAxx+PAx*yAx+PAx*yBx)*tMx-PA*tMxx*(yAx+yBx))*l+mM*(((-xc*yAx+yc)*yBx+yc*yAx+xc)*PA*tMx^2+((yc*yAxx*PA+PAx*(yAx*yc+xc))*yBx+PA*yAxx*xc-PAx*(-xc*yAx+yc))*tMx-((yAx*yc+xc)*yBx+yAx*xc-yc)*tMxx*PA))*phiA*R^2+((-mB*((yAx*yBx-1)*tMx^2+tMxx*(yAx+yBx))*l+(((-xc*yAx+yc)*yBx+yc*yAx+xc)*tMx^2-((yAx*yc+xc)*yBx+yAx*xc-yc)*tMxx)*mM)*phiA+phiAx*tMx*(mB*(yAx+yBx)*l+((yAx*yc+xc)*yBx+yAx*xc-yc)*mM))*I_A)*phiB*R^2*PB^2+(((PA*tMx^2+(PA*yAxx+PAx*yAx)*tMx-PA*tMxx*yAx)*phiB-tMx*yAx*phiBx*PA)*(mA+mM)*phiA*R^2+I_A*(((tMx^2-tMxx*yAx)*phiB-tMx*yAx*phiBx)*phiA+tMx*yAx*phiB*phiAx))*l*(yBx^2+1)*I_B*PB+PBx*l*phiB*tMx*(PA*(mA+mM)*R^2+I_A)*(yBx^2+1)*phiA*I_B*yAx)*sin(tM)+phiB*R^2*((tMx^3*mB*mM*(yBx*yc-xc)*l^2+(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*l+(yBx*yc+xc)*(I_M+mM*(xc^2+yc^2))*mM*tMx^3)*phiA*R^2+((-PA*yAxx*mM+yAx*mB*(PA*yAxx+PAx*yAx)*yBx+PAx*yAx*mB-PA*yAxx*mA)*phiA-yAx*phiAx*PA*mB*(yAx*yBx+1))*I_A*l)*PB^2+l*(yBx^2+1)*I_B*(((-tMx^3*l*mM*xc+(I_M+mM*(xc^2+yc^2))*tMx^3-yAxx*PA^2*(mA+mM))*phiB+yAx*phiBx*PA^2*(mA+mM))*phiA*R^2+((PA*phiBx+PAx*phiB)*phiA-phiB*phiAx*PA)*I_A*yAx)*PB-PBx*l*phiB*(PA*(mA+mM)*R^2+I_A)*(yBx^2+1)*phiA*I_B*PA*yAx)*cos(tM)+(-((tMx^3*mB*mM*(xc*yBx+yc)*l^2-yBx*(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*l+(-xc*yBx+yc)*(I_M+mM*(xc^2+yc^2))*mM*tMx^3)*phiA*R^2+I_A*l*((yBx*yAxx*PA*mM+(yAxx*(mA+mB)*PA+PAx*yAx*mB)*yBx+PAx*mB)*phiA-phiAx*PA*mB*(yAx*yBx+1)))*phiB*R^2*PB^2+l*((-yc*tMx^3*l*phiB*mM+(-tMx^2*l*mM*xc+(I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM))*phiBx)*phiA*R^2-((PA*phiBx+PAx*phiB)*phiA-phiB*phiAx*PA)*I_A)*(yBx^2+1)*I_B*PB-PBx*l*phiB*((-tMx^2*l*mM*xc+(I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM))*R^2-I_A*PA)*(yBx^2+1)*phiA*I_B)*sin(tM)-((mB*(mA+mM)*(-PA*tMx^2*yBx-PA*tMxx+PAx*tMx)*l^2+mM*(-((yAx*yc+xc)*yBx*mM+(yc*(mA+mB)*yAx+mA*xc)*yBx+yc*mB)*PA*tMx^2+((-PA*yAxx*xc+PAx*(-xc*yAx+yc))*yBx*mM+(-yAxx*xc*(mA+mB)*PA+PAx*(-xc*(mA+mB)*yAx+yc*mA))*yBx-PAx*mB*xc)*tMx-tMxx*((-xc*yAx+yc)*yBx*mM+(-xc*(mA+mB)*yAx+yc*mA)*yBx-mB*xc)*PA)*l+(((PA*yAxx+PAx*yAx)*yBx+PAx)*tMx-PA*tMxx*(yAx*yBx+1))*(I_M+mM*(xc^2+yc^2))*(mA+mB+mM))*phiA*R^2+I_A*((-mB*(tMx^2*yBx+tMxx)*l^2-yBx*((yAx*yc+xc)*tMx^2+tMxx*(-xc*yAx+yc))*mM*l-(I_M+mM*(xc^2+yc^2))*(yAx*yBx+1)*tMxx)*phiA+phiAx*tMx*(l^2*mB+(-xc*yAx+yc)*yBx*mM*l+(I_M+mM*(xc^2+yc^2))*(yAx*yBx+1))))*phiB*R^2*PB^2-((((mA+mM)*(-PA*tMxx+PAx*tMx)*l^2-mM*(PA*tMx^2*yc-PA*tMxx*xc+PAx*tMx*xc)*l+(I_M+mM*(xc^2+yc^2))*(-PA*tMxx+PAx*tMx))*phiB+((-mA-mM)*l^2-l*mM*xc+mM*(xc^2+yc^2)+I_M)*tMx*phiBx*PA)*phiA*R^2-I_A*l^2*((phiB*tMxx+phiBx*tMx)*phiA-tMx*phiB*phiAx))*(yBx^2+1)*I_B*PB+PBx*phiB*tMx*(yBx^2+1)*(((-mA-mM)*l^2-l*mM*xc+mM*(xc^2+yc^2)+I_M)*PA*R^2-I_A*l^2)*phiA*I_B)/(phiB*(-l*tMx*(R^2*(mB*(yAx+yBx)*l+((yAx*yc+xc)*yBx+yAx*xc-yc)*mM)*PB+I_B*yAx*l*(yBx^2+1))*cos(tM)^2+l*(tMx*(((-mB*yAx*yBx+mB)*l+mM*((-xc*yAx+yc)*yBx+yc*yAx+xc))*R^2*PB+l*(yBx^2+1)*I_B)*sin(tM)+(-((yAx^2+1)*mM+yAx^2*mA+yAx*yBx*mB+mA+mB)*PA*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*cos(tM)-l*(((yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2+mA)*yBx+yAx*mB)*PA*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*sin(tM)+tMx*(R^2*(yAx*l^2*mB+(yAx*yc+xc)*yBx*mM*l+(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*PB+I_B*yAx*l^2*(yBx^2+1)))*PB*R);

JNA = -R*(tMx*((yBx*mB*(mA+mM)*l^2-mM*(mA-mB+mM)*(-xc*yBx+yc)*l+((-xc^2+yc^2)*yBx+2*yc*xc)*mM^2)*R^2*PB+yc*I_B*l*mM*(yBx^2+1))*cos(tM)^2+(-tMx*((mB*(mA+mM)*l^2+mM*(mA-mB+mM)*(yBx*yc+xc)*l+mM^2*(-2*xc*yBx*yc-xc^2+yc^2))*R^2*PB+l*((mA+mM)*l-mM*xc)*(yBx^2+1)*I_B)*sin(tM)+(-(mA+mB+mM)*PA*((-mA-mM)*l+mM*xc*(yAx*yBx+1))*R^2-I_A*((mB*yAx*yBx-mA-mM)*l+mM*xc*(yAx*yBx+1)))*PB+((mA+mM)*l-mM*xc)*(yBx^2+1)*I_B*PA)*cos(tM)+(((mA+mB+mM)*(yBx*(mA+mM)*l+yc*mM*(yAx*yBx+1))*PA*R^2+I_A*(l*(mA+mB+mM)*yBx+yc*mM*(yAx*yBx+1)))*PB+yc*I_B*PA*mM*(yBx^2+1))*sin(tM)+((-(yBx*mM*xc+xc*(mA+mB)*yBx+yc*mB)*mM*l+(xc^2*yBx-xc*yc)*mM^2+((mA+mB)*xc^2+(mA+mB)*yc^2+I_M)*yBx*mM+I_M*yBx*(mA+mB))*R^2*PB-yc*I_B*l*mM*(yBx^2+1))*tMx)*phiA/(-l*tMx*(R^2*(mB*(yAx+yBx)*l+((yAx*yc+xc)*yBx+yAx*xc-yc)*mM)*PB+I_B*yAx*l*(yBx^2+1))*cos(tM)^2+l*(tMx*(((-mB*yAx*yBx+mB)*l+mM*((-xc*yAx+yc)*yBx+yc*yAx+xc))*R^2*PB+I_B*l*(yBx^2+1))*sin(tM)+(-((yAx^2+1)*mM+yAx^2*mA+yAx*yBx*mB+mA+mB)*PA*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*cos(tM)-l*(((yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2+mA)*yBx+yAx*mB)*PA*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*sin(tM)+tMx*(R^2*(yAx*l^2*mB+(yAx*yc+xc)*yBx*mM*l+(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*PB+I_B*yAx*l^2*(yBx^2+1)));

TANA = (-((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*l*tMx*cos(tM)^2-(((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*sin(tM)-yAx*(R^2*mB*(yAx*yBx+1)*PB+I_B*(yBx^2+1))*PA)*l*cos(tM)-l*(R^2*mB*(yAx*yBx+1)*PB+I_B*(yBx^2+1))*PA*sin(tM)+((PB*R^2*mB+I_B*yBx^2+I_B)*l^2+R^2*yBx*PB*mM*(-xc*yAx+yc)*l+(I_M+mM*(xc^2+yc^2))*PB*R^2*(yAx*yBx+1))*tMx)/((-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*R);

TBNA = PB*(-yAx*l*(PA*(mA+mM)*R^2+I_A)*cos(tM)+l*(PA*(mA+mM)*R^2+I_A)*sin(tM)+tMx*R^2*(-l*mM*xc+mM*(xc^2+yc^2)+I_M))*(yBx^2+1)*phiA/((-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)-((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*PB*PA*R^2+(-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-l*(PB*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*PA*R^2+I_A*yBx*(yAx^2+1)*PB+I_B*yAx*PA*(yBx^2+1))*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*phiB*R);

FNA_coeff.HNA = HNA;
FNA_coeff.JNA = JNA;
FNA_coeff.TANA = TANA;
FNA_coeff.TBNA = TBNA;

% Coefficients for Normal force at Wheel B

HNB = (-l^2*tMx^2*(-PB*phiBx+PBx*phiB)*R^2*phiA*I_B*(mB*(yAx+yBx)*l+((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*cos(tM)^3+2*l*((1/2)*l*((-mB*yAx*yBx+mB)*l+mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx^2*(-PB*phiBx+PBx*phiB)*R^2*phiA*I_B*sin(tM)+(-(1/2)*(2*tMx^2*yAx*PA+(PA*yAx*yAxx+PAx*yAx^2-PAx)*tMx-yAx^2*tMxx*PA+tMxx*PA)*mB*l+mM*(-(1/2)*PA*(yAx^2*yc+2*xc*yAx-yc)*tMx^2+(-(1/2)*PAx*xc*yAx^2+(PAx*yc-(1/2)*PA*yAxx*xc)*yAx+(1/2)*yc*yAxx*PA+(1/2)*PAx*xc)*tMx-tMxx*PA*(-(1/2)*xc*yAx^2+yc*yAx+(1/2)*xc)))*mB*phiB*PB^2*phiA*R^4+(1/2)*(((((yBx*(mA+mM)*yAx^2-2*yAx*mB+yBx*(mA+mM))*PA*tMx^2+(-PAx*(mA+mB+mM)*yAx^2-PA*yAxx*(mA+mB+mM)*yAx-PAx*(mA-mB+mM))*tMx+tMxx*((mA+mB+mM)*yAx^2+mA-mB+mM)*PA)*phiB+((mA+mB+mM)*yAx^2-2*yAx*yBx*mB+mA+mB+mM)*tMx*phiBx*PA)*l+((-((yBx*yc+xc)*yAx+yBx*xc-yc)*PA*tMx^2+PAx*((-xc*yBx+yc)*yAx+yc*yBx+xc)*tMx-tMxx*((-xc*yBx+yc)*yAx+yc*yBx+xc)*PA)*phiB+tMx*phiBx*((-xc*yBx+yc)*yAx+yc*yBx+xc)*PA)*mM)*PB-PBx*phiB*tMx*(((mA+mB+mM)*yAx^2-2*yAx*yBx*mB+mA+mB+mM)*l+mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*PA)*phiA*I_B*R^2-(1/2)*I_A*l*(yAx^2+1)*((((-tMx^2*yBx-tMxx)*phiB-tMx*phiBx)*PB+PBx*phiB*tMx)*phiA+PB*phiAx*phiB*tMx)*I_B)*cos(tM)^2+(((2*mB*((-(1/2)*yAx^2*PA+(1/2)*PA)*tMx^2+((1/2)*PA*yAxx+PAx*yAx)*tMx-tMxx*PA*yAx)*l+mM*(2*PA*(-(1/2)*xc*yAx^2+yc*yAx+(1/2)*xc)*tMx^2+(PAx*yc*yAx^2+(PA*yAxx*yc+2*PAx*xc)*yAx+PA*yAxx*xc-PAx*yc)*tMx-tMxx*PA*(yAx^2*yc+2*xc*yAx-yc)))*mB*phiB*PB^2*phiA*R^4+((((-((mA+mB+mM)*yAx^2+mA-mB+mM)*PA*tMx^2+(-PAx*yBx*(mA+mM)*yAx^2+(-PA*mA*yAxx*yBx-PA*mM*yAxx*yBx+2*PAx*mB)*yAx-PAx*yBx*mA-PAx*yBx*mM+PA*yAxx*mB)*tMx+tMxx*(yBx*(mA+mM)*yAx^2-2*yAx*mB+yBx*(mA+mM))*PA)*phiB+((mA+mM)*yAx^2+mA+2*mB+mM)*yBx*tMx*phiBx*PA)*l+((((-xc*yBx+yc)*yAx+yc*yBx+xc)*PA*tMx^2+((yBx*yc+xc)*yAx+yBx*xc-yc)*PAx*tMx-((yBx*yc+xc)*yAx+yBx*xc-yc)*tMxx*PA)*phiB+((yBx*yc+xc)*yAx+yBx*xc-yc)*tMx*phiBx*PA)*mM)*PB-PBx*phiB*(((mA+mM)*yAx^2+mA+2*mB+mM)*yBx*l+((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*PA)*phiA*I_B*R^2-((((tMx^2-tMxx*yBx)*phiB-tMx*yBx*phiBx)*PB+PBx*tMx*yBx*phiB)*phiA+tMx*yBx*phiB*phiAx*PB)*I_A*l*(yAx^2+1)*I_B)*l*sin(tM)-phiB*PB^2*phiA*(tMx^3*l^3*mB^2+tMx^3*mB*mM*(yAx*yc+xc)*l^2+((I_M+mM*(xc^2+yc^2))*tMx^3+yAxx*PA^2*(mA+mB+mM))*mB*l+(I_M+mM*(xc^2+yc^2))*mM*tMx^3*(yAx*yc+xc))*R^4-((I_A*phiB*mB*(yAx^2+1)*(PA*yAxx+PAx*yAx)*PB^2+((tMx*phiB+phiBx*(yAx+yBx))*mB*tMx^2*l^2+(tMx*yBx*(-xc*yAx+yc)*phiB+((yBx*yc+xc)*yAx+yBx*xc-yc)*phiBx)*mM*tMx^2*l+((I_M+mM*(xc^2+yc^2))*(yAx*yBx+1)*tMx^3-(yBx*(mA+mM)*yAx-mB)*yAxx*PA^2)*phiB+phiBx*(yBx*(mA+mM)*yAx^2-yAx*mB+yBx*(mA+mB+mM))*PA^2)*I_B*PB-PBx*phiB*(tMx^2*mB*(yAx+yBx)*l^2+((yBx*yc+xc)*yAx+yBx*xc-yc)*mM*tMx^2*l+(yBx*(mA+mM)*yAx^2-yAx*mB+yBx*(mA+mB+mM))*PA^2)*I_B)*phiA-I_A*yAx*phiB*phiAx*PA*PB^2*mB*(yAx^2+1))*l*R^2-(((PA*phiBx+PAx*phiB)*PB-PBx*phiB*PA)*phiA-phiB*phiAx*PA*PB)*I_A*l*(yAx^2+1)*yBx*I_B)*cos(tM)+((-tMx^3*yAx*l^3*mB^2+tMx^3*mB*mM*(-xc*yAx+yc)*l^2-((I_M+mM*(xc^2+yc^2))*tMx^3+yAxx*PA^2*(mA+mB+mM))*mB*yAx*l+(I_M+mM*(xc^2+yc^2))*mM*(-xc*yAx+yc)*tMx^3)*phiB*PB^2*phiA*R^4+l*((I_A*PAx*phiB*mB*(yAx^2+1)*PB^2-(tMx^2*yAx*mB*(phiB*tMx+phiBx*yBx)*l^2-(tMx*(-xc*yAx+yc)*phiB+phiBx*(yAx*yc+xc))*mM*tMx^2*l+yAx*yAxx*PA^2*(mA+mB+mM)*phiB+((I_M+mM*(xc^2+yc^2))*(yAx*yBx+1)*tMx^2-((mA+mB+mM)*yAx^2-yAx*yBx*mB+mA+mM)*PA^2)*phiBx)*I_B*PB+PBx*(tMx^2*yAx*yBx*l^2*mB-tMx^2*mM*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*(yAx*yBx+1)*tMx^2-((mA+mB+mM)*yAx^2-yAx*yBx*mB+mA+mM)*PA^2)*phiB*I_B)*phiA-I_A*phiB*phiAx*PA*PB^2*mB*(yAx^2+1))*R^2+(((PA*phiBx+PAx*phiB)*PB-PBx*phiB*PA)*phiA-phiB*phiAx*PA*PB)*I_A*l*(yAx^2+1)*I_B)*sin(tM)+((PA*mB*tMx^2*yAx+(PAx*(mA+mB+mM)*yAx^2+PA*yAxx*(mA+mB+mM)*yAx+PAx*(mA+mM))*tMx-tMxx*((mA+mB+mM)*yAx^2+mA+mM)*PA)*mB*l^2-mB*(PA*(-xc*yAx+yc)*tMx^2+PAx*(yAx*yc+xc)*tMx-tMxx*PA*(yAx*yc+xc))*mM*l+(I_M+mM*(xc^2+yc^2))*(mA+mB+mM)*((PA*yAx*yAxx+PAx*yAx^2+PAx)*tMx-tMxx*PA*(yAx^2+1)))*phiB*PB^2*phiA*R^4+((-(l^2*mB+mM*(xc^2+yc^2)+I_M)*I_A*phiB*tMxx*(yAx^2+1)*PB^2+(((PA*mB*tMx^2*yAx+(PAx*(mA+mB+mM)*yAx^2+PA*yAxx*(mA+mB+mM)*yAx+PAx*(mA+mM))*tMx-tMxx*((mA+mB+mM)*yAx^2+mA+mM)*PA)*phiB-((mA+mB+mM)*yAx^2-2*yAx*yBx*mB+mA+mM)*tMx*phiBx*PA)*l^2-((PA*(-xc*yAx+yc)*tMx^2+PAx*(yAx*yc+xc)*tMx-tMxx*PA*(yAx*yc+xc))*phiB+tMx*phiBx*PA*(yAx*yc+xc))*mM*l+((-PA*tMxx+PAx*tMx)*phiB+tMx*phiBx*PA)*(I_M+mM*(xc^2+yc^2))*(yAx*yBx+1))*I_B*PB-PBx*(((-mA-mB-mM)*yAx^2+2*yAx*yBx*mB-mA-mM)*l^2-mM*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*(yAx*yBx+1))*phiB*tMx*I_B*PA)*phiA+(l^2*mB+mM*(xc^2+yc^2)+I_M)*I_A*phiAx*phiB*(yAx^2+1)*PB^2*tMx)*R^2+I_A*l^2*(yAx^2+1)*(((-phiB*tMxx-phiBx*tMx)*PB+PBx*phiB*tMx)*phiA+PB*phiAx*phiB*tMx)*I_B)/((-l*tMx*(PB*(mB*(yAx+yBx)*l+((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*R^2+I_B*yAx*l*(yBx^2+1))*cos(tM)^2+l*(tMx*(((-mB*yAx*yBx+mB)*l+mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*PB*R^2+I_B*l*(yBx^2+1))*sin(tM)-((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*PB*PA*R^2+(-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*cos(tM)-l*(PB*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*PA*R^2+I_A*yBx*(yAx^2+1)*PB+I_B*yAx*PA*(yBx^2+1))*sin(tM)+tMx*(R^2*(yAx*l^2*mB+(yAx*yc+xc)*yBx*mM*l+(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*PB+I_B*yAx*l^2*(yBx^2+1)))*PB*R*phiA);

JNB = R*(((-yAx*l^2*mB^2+2*mB*mM*(-xc*yAx+yc)*l+((-xc^2+yc^2)*yAx+2*yc*xc)*mM^2)*R^2*PB+l*(-yAx*(mA+mB+mM)*l+mM*yc*(yAx*yBx+1))*I_B)*tMx*cos(tM)^2+(-(R^2*(-l^2*mB^2-2*mB*mM*(yAx*yc+xc)*l+mM^2*(-2*xc*yAx*yc-xc^2+yc^2))*PB+((yBx*(mA+mM)*yAx-mB)*l-mM*xc*(yAx*yBx+1))*l*I_B)*tMx*sin(tM)+(-(PA*(mA+mB+mM)*R^2+I_A*(yAx^2+1))*mB*l-(PA*(mA+mB+mM)*R^2+I_A)*xc*mM*(yAx^2+1))*PB+((yBx*(mA+mM)*yAx-mB)*l-mM*xc*(yAx*yBx+1))*I_B*PA)*cos(tM)+((-R^2*yAx*PA*mB*(mA+mB+mM)*l+(PA*(mA+mB+mM)*R^2+I_A)*mM*(yAx^2+1)*yc)*PB+(-yAx*(mA+mB+mM)*l+mM*yc*(yAx*yBx+1))*I_B*PA)*sin(tM)+((yAx*mB*(mA+mB+mM)*l^2-yc*l*mB*mM+(mM^2*xc^2+((xc^2+yc^2)*mB+yc^2*mA+mA*xc^2+I_M)*mM+I_M*(mA+mB))*yAx-yc*mM^2*xc)*R^2*PB-l*(-yAx*(mA+mB+mM)*l+mM*yc*(yAx*yBx+1))*I_B)*tMx)*phiB/(-l*tMx*(PB*(mB*(yAx+yBx)*l+((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*R^2+I_B*yAx*l*(yBx^2+1))*cos(tM)^2+l*(tMx*(((-mB*yAx*yBx+mB)*l+mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*PB*R^2+I_B*l*(yBx^2+1))*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+tMx*(R^2*(yAx*l^2*mB+(yAx*yc+xc)*yBx*mM*l+(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*PB+I_B*yAx*l^2*(yBx^2+1)));

TANB = -(-I_B*cos(tM)^2*l^2*tMx+l*(-I_B*tMx*yBx*l*sin(tM)+PA*(PB*R^2*mB*yAx+I_B*yBx))*cos(tM)-l*PA*(PB*R^2*mB+I_B)*sin(tM)+((PB*R^2*mB+I_B)*l^2+(I_M+mM*(xc^2+yc^2))*PB*R^2)*tMx)*phiB*(yAx^2+1)/((-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*R*phiA);

TBNB = -PB*(l*tMx*R^2*((-mB*yAx*yBx+mB)*l+mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*cos(tM)^2+((mB*(yAx+yBx)*l+((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*R^2*sin(tM)-(yBx*(mA+mM)*yAx^2-yAx*mB+yBx*(mA+mB+mM))*PA*R^2-I_A*yBx*(yAx^2+1))*l*cos(tM)+(((mA+mB+mM)*yAx^2-yAx*yBx*mB+mA+mM)*PA*R^2+I_A*(yAx^2+1))*l*sin(tM)+tMx*(yAx*yBx*l^2*mB-mM*(yAx*yc+xc)*l+(I_M+mM*(xc^2+yc^2))*(yAx*yBx+1))*R^2)/((-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*l*tMx*cos(tM)^2+l*(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+R^2*PB*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)-((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*PB*PA*R^2-I_A*PB*yAx^2-I_B*PA*yBx^2-I_A*PB-I_B*PA)*cos(tM)-l*(PB*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*PA*R^2+I_A*PB*yAx^2*yBx+I_B*yAx*PA*(yBx^2+1)+I_A*PB*yBx)*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+(yAx-yBx)*R^2*PB*(I_M+mM*(xc^2+yc^2)))*tMx)*R);

FNB_coeff.HNB = HNB;
FNB_coeff.JNB = JNB;
FNB_coeff.TANB = TANB;
FNB_coeff.TBNB = TBNB;


