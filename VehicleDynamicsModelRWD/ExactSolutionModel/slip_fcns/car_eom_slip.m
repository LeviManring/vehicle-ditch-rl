function [H_s, J_s, H_N_s, J_N_s, tM, dx] = car_eom_slip(xA,guess_Dx,mu,params)
%% car_eom_slip
% Levi Manring, Duke University
% 2021
%
% This function calculates certain parameters needed for integration of the slip_ode.m
%
% Inputs:
%   xA: 1x1 double indicating the x-position of wheel A
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   mu: 1x1 double indicating the directional friction coefficient for wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%
% Outputs:
%       Outputs are elements necessary to solve the slip vehicle dynamics
%       model: x'' + H_s*x'^2 + J_s*g = 0
%       And solve the normal force as well: F_N = H_N_s*x'^2 + J_N_s*g
%   H_s: 1x1 double, ode parameter
%   J_s: 1x1 double, ode parameter
%   TA_s: 1x1 double, ode parameter
%   H_N_s: 1x1 double, normal force parameter
%   J_N_s: 1x1 double, normal force parameter
%   TA_N_s: 1x1 double, normal force parameter
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
% I_A = params.dim.I_A;
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

%% Calculating dynamic model parameters
% Applying torque at just back wheel A and including moment of inertia effects for both wheels

% These equations derivation can be seen in Maple files included (copied
% and pasted directly into Matlab)

H_s = (((mB*R^2*(((mu*tMxx+tMx^2)*yAx-tMx^2*mu+tMxx)*yBx+(-mu*tMx^2+tMxx)*yAx-tMx^2-tMxx*mu)*phiB*PB^2-(((mu*tMx^2-tMxx)*yAx+tMx^2+tMxx*mu)*phiB-tMx*phiBx*(yAx-mu))*I_B*(yBx^2+1)*PB-I_B*PBx*tMx*phiB*(yBx^2+1)*(yAx-mu))*l+R^2*mM*phiB*((((-mu*yc+xc)*tMx^2+tMxx*(mu*xc+yc))*yAx+(-mu*xc-yc)*tMx^2+tMxx*(-mu*yc+xc))*yBx+((-mu*xc-yc)*tMx^2+tMxx*(-mu*yc+xc))*yAx+(mu*yc-xc)*tMx^2-tMxx*(mu*xc+yc))*PB^2)*l*cos(tM)^2-l*(((mB*R^2*phiB*(((mu*tMx^2-tMxx)*yAx+tMx^2+tMxx*mu)*yBx+(mu*tMxx+tMx^2)*yAx-tMx^2*mu+tMxx)*PB^2+(((mu*tMxx+tMx^2)*yAx-tMx^2*mu+tMxx)*phiB+tMx*phiBx*(mu*yAx+1))*I_B*(yBx^2+1)*PB-I_B*PBx*tMx*phiB*(yBx^2+1)*(mu*yAx+1))*l+R^2*mM*phiB*((((mu*xc+yc)*tMx^2-tMxx*(-mu*yc+xc))*yAx+(-mu*yc+xc)*tMx^2+tMxx*(mu*xc+yc))*yBx+((-mu*yc+xc)*tMx^2+tMxx*(mu*xc+yc))*yAx+(-mu*xc-yc)*tMx^2+tMxx*(-mu*yc+xc))*PB^2)*sin(tM)-(mB*(mu*yAx+1)*(PA*yAxx+PAx*yAx)*yBx+PAx*(mA+mM)*yAx^2+(PAx*mB*mu+yAxx*PA*(mA+mM))*yAx-yAxx*PA*(mA+mM)*mu+PAx*(mA+mB+mM))*R^2*phiB*PB^2-I_B*(yBx^2+1)*(mu*yAx+1)*(PA*phiBx+PAx*phiB)*PB+I_B*PBx*phiB*PA*(yBx^2+1)*(mu*yAx+1))*cos(tM)+l*(R^2*((PAx*(mA+mB+mM)*yAx^2+(-PAx*mB*mu+yAxx*PA*(mA+mB+mM))*yAx-yAxx*PA*(mA+mB+mM)*mu+PAx*(mA+mM))*yBx+PAx*mB*(yAx-mu))*phiB*PB^2+I_B*(yBx^2+1)*(yAx-mu)*(PA*phiBx+PAx*phiB)*PB-I_B*PBx*phiB*PA*(yBx^2+1)*(yAx-mu))*sin(tM)-(R^2*phiB*mB*(tMx^2*yBx+tMxx)*PB^2+I_B*(yBx^2+1)*(phiB*tMxx+phiBx*tMx)*PB-I_B*PBx*tMx*phiB*(yBx^2+1))*(yAx-mu)*l^2-R^2*yBx*mM*phiB*PB^2*(((-mu*yc+xc)*tMx^2+tMxx*(mu*xc+yc))*yAx+(-mu*xc-yc)*tMx^2+tMxx*(-mu*yc+xc))*l+tMxx*R^2*(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx)*phiB*PB^2)/((tMx*((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*R^2*mM*PB)*l*cos(tM)^2+(tMx*((mB*((yAx-mu)*yBx-yAx*mu-1)*R^2*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+R^2*mM*PB*(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc))*sin(tM)+PA*(R^2*(yAx*mB*(mu*yAx+1)*yBx+(mA+mM)*yAx^2+yAx*mB*mu+mA+mB+mM)*PB+I_B*(yBx^2+1)*(mu*yAx+1)))*l*cos(tM)+PA*l*((((mA+mB+mM)*yAx^2-yAx*mB*mu+mA+mM)*yBx+mB*(yAx-mu))*R^2*PB+I_B*(yBx^2+1)*(yAx-mu))*sin(tM)+(-(R^2*PB*mB+I_B*(yBx^2+1))*(yAx-mu)*l^2-R^2*yBx*mM*PB*((mu*xc+yc)*yAx-yc*mu+xc)*l+R^2*(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx)*PB)*tMx)*phiB*PB);

J_s = R^2*(((mB*(mu*yAx+1)*yBx+(mA+mM)*(yAx-mu))*l+xc*mM*((mu*yAx+1)*yBx-yAx+mu))*cos(tM)-sin(tM)*(-yBx*(mA+mB+mM)*(yAx-mu)*l+mM*((mu*yAx+1)*yBx-yAx+mu)*yc))*PB/(tMx*((R^2*mB*(mu*yAx*yBx-mu+yAx+yBx)*PB+I_B*(yBx^2+1)*(yAx-mu))*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*R^2*mM*PB)*l*cos(tM)^2+(tMx*((mB*((yAx-mu)*yBx-yAx*mu-1)*R^2*PB-I_B*(yBx^2+1)*(mu*yAx+1))*l+R^2*mM*PB*(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc))*sin(tM)+PA*(R^2*(yAx*mB*(mu*yAx+1)*yBx+(mA+mM)*yAx^2+yAx*mB*mu+mA+mB+mM)*PB+I_B*(yBx^2+1)*(mu*yAx+1)))*l*cos(tM)+PA*l*((((mA+mB+mM)*yAx^2-yAx*mB*mu+mA+mM)*yBx+mB*(yAx-mu))*R^2*PB+I_B*(yBx^2+1)*(yAx-mu))*sin(tM)+(-(R^2*PB*mB+I_B*(yBx^2+1))*(yAx-mu)*l^2-R^2*yBx*mM*PB*((mu*xc+yc)*yAx-yc*mu+xc)*l+R^2*(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx)*PB)*tMx);

HGN = -phiA*R*(-l*(mA+mM)*(R^2*phiB*((PA*(yAx+yBx)*tMx^2+(yBx*yAxx*PA+PAx*(yAx*yBx-1))*tMx-tMxx*PA*(yAx*yBx-1))*mB*l+(PA*((yAx*yc+xc)*yBx+xc*yAx-yc)*tMx^2+(yAxx*(xc*yBx-yc)*PA+((xc*yAx-yc)*yBx-yc*yAx-xc)*PAx)*tMx-tMxx*PA*((xc*yAx-yc)*yBx-yc*yAx-xc))*mM)*PB^2-I_B*((-PA*tMx^2*yAx-PA*tMxx+PAx*tMx)*phiB-tMx*phiBx*PA)*l*(yBx^2+1)*PB-I_B*PBx*tMx*l*phiB*PA*(yBx^2+1))*cos(tM)^2+((R^2*phiB*(((-yAx*yBx+1)*PA*tMx^2+(PA*yAxx+PAx*(yAx+yBx))*tMx-tMxx*PA*(yAx+yBx))*mB*l+mM*(-PA*((xc*yAx-yc)*yBx-yc*yAx-xc)*tMx^2+(yAxx*(yBx*yc+xc)*PA+((yAx*yc+xc)*yBx+xc*yAx-yc)*PAx)*tMx-tMxx*PA*((yAx*yc+xc)*yBx+xc*yAx-yc)))*PB^2+I_B*((tMx^2*PA+(PA*yAxx+PAx*yAx)*tMx-yAx*tMxx*PA)*phiB-tMx*yAx*phiBx*PA)*l*(yBx^2+1)*PB+I_B*PBx*tMx*yAx*l*phiB*PA*(yBx^2+1))*l*(mA+mM)*sin(tM)+R^2*(-tMx^3*mB*mM*(-yBx*yc+xc)*l^2+(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*l+(yBx*yc+xc)*(I_M+mM*(xc^2+yc^2))*tMx^3*mM)*phiB*PB^2+I_B*(-xc*tMx^3*l*phiB*mM+((I_M+mM*(xc^2+yc^2))*tMx^3-yAxx*PA^2*(mA+mM))*phiB+yAx*phiBx*PA^2*(mA+mM))*l*(yBx^2+1)*PB-I_B*PBx*yAx*l*phiB*PA^2*(yBx^2+1)*(mA+mM))*cos(tM)+(R^2*phiB*(-tMx^3*mB*mM*(xc*yBx+yc)*l^2+yBx*(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*l+(I_M+mM*(xc^2+yc^2))*(xc*yBx-yc)*tMx^3*mM)*PB^2+I_B*(-tMx^2*mM*(phiB*tMx*yc+phiBx*xc)*l+((I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM))*phiBx)*l*(yBx^2+1)*PB-I_B*(-xc*tMx^2*l*mM+(I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM))*PBx*phiB*l*(yBx^2+1))*sin(tM)-(mB*(mA+mM)*(-PA*tMx^2*yBx-PA*tMxx+PAx*tMx)*l^2-(PA*((yAx*yc+xc)*yBx*mM+(yc*(mA+mB)*yAx+xc*mA)*yBx+yc*mB)*tMx^2+((xc*yAxx*PA+PAx*(xc*yAx-yc))*yBx*mM+xc*yBx*yAxx*(mA+mB)*PA+((xc*(mA+mB)*yAx-yc*mA)*yBx+xc*mB)*PAx)*tMx-tMxx*((xc*yAx-yc)*yBx*mM+(xc*(mA+mB)*yAx-yc*mA)*yBx+xc*mB)*PA)*mM*l+(I_M+mM*(xc^2+yc^2))*(mA+mB+mM)*((yBx*yAxx*PA+PAx*(yAx*yBx+1))*tMx-tMxx*PA*(yAx*yBx+1)))*R^2*phiB*PB^2-I_B*(((-PA*tMxx+PAx*tMx)*phiB-tMx*phiBx*PA)*(mA+mM)*l^2-((PA*tMx^2*yc-PA*tMxx*xc+PAx*tMx*xc)*phiB+xc*tMx*phiBx*PA)*mM*l+(I_M+mM*(xc^2+yc^2))*((-PA*tMxx+PAx*tMx)*phiB+tMx*phiBx*PA))*(yBx^2+1)*PB+I_B*PA*PBx*phiB*tMx*(yBx^2+1)*((-mA-mM)*l^2-xc*l*mM+mM*(xc^2+yc^2)+I_M))/(((R^2*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*PB+I_B*l*(yBx^2+1)*(yAx-mu))*l*tMx*cos(tM)^2+((R^2*(mB*((yAx-mu)*yBx-yAx*mu-1)*l+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*mM)*PB-I_B*l*(yBx^2+1)*(mu*yAx+1))*tMx*sin(tM)+(((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB)*R^2*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+PA*l*(R^2*(yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*sin(tM)+(R^2*(-mB*(yAx-mu)*l^2-yBx*mM*((mu*xc+yc)*yAx-yc*mu+xc)*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx))*PB-I_B*l^2*(yBx^2+1)*(yAx-mu))*tMx)*PB*phiB);

JGN = -R*phiA*(((-yBx*mB*(mA+mM)*l^2-mM*(mA-mB+mM)*(xc*yBx-yc)*l+((xc^2-yc^2)*yBx-2*xc*yc)*mM^2)*PB*R^2-yc*I_B*l*mM*(yBx^2+1))*tMx*cos(tM)^2+(-(PB*(-mB*(mA+mM)*l^2-mM*(mA-mB+mM)*(yBx*yc+xc)*l+mM^2*(2*xc*yBx*yc+xc^2-yc^2))*R^2+I_B*((-mA-mM)*l+xc*mM)*l*(yBx^2+1))*tMx*sin(tM)+PA*((mA+mB+mM)*PB*((-mA-mM)*l+xc*mM*(yAx*yBx+1))*R^2+I_B*((-mA-mM)*l+xc*mM)*(yBx^2+1)))*cos(tM)-PA*((mA+mB+mM)*(yBx*(mA+mM)*l+yc*mM*(yAx*yBx+1))*PB*R^2+yc*I_B*mM*(yBx^2+1))*sin(tM)-((-(xc*yBx*mM+xc*(mA+mB)*yBx+yc*mB)*mM*l+xc*(xc*yBx-yc)*mM^2+yBx*((mA+mB)*xc^2+(mA+mB)*yc^2+I_M)*mM+I_M*yBx*(mA+mB))*PB*R^2-yc*I_B*l*mM*(yBx^2+1))*tMx)/((R^2*(mB*((mu*yAx+1)*yBx+yAx-mu)*l+(((mu*xc+yc)*yAx-yc*mu+xc)*yBx+(-mu*yc+xc)*yAx-xc*mu-yc)*mM)*PB+I_B*l*(yBx^2+1)*(yAx-mu))*l*tMx*cos(tM)^2+((R^2*(mB*((yAx-mu)*yBx-yAx*mu-1)*l+(((-mu*yc+xc)*yAx-xc*mu-yc)*yBx+(-mu*xc-yc)*yAx+yc*mu-xc)*mM)*PB-I_B*l*(yBx^2+1)*(mu*yAx+1))*tMx*sin(tM)+(((yAx^2+1)*mM+yAx*mB*(mu*yAx+1)*yBx+yAx^2*mA+yAx*mB*mu+mA+mB)*R^2*PB+I_B*(yBx^2+1)*(mu*yAx+1))*PA)*l*cos(tM)+PA*l*(R^2*(yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2-yAx*mB*mu+mA)*yBx+mB*(yAx-mu))*PB+I_B*(yBx^2+1)*(yAx-mu))*sin(tM)+(R^2*(-mB*(yAx-mu)*l^2-yBx*mM*((mu*xc+yc)*yAx-yc*mu+xc)*l+(I_M+mM*(xc^2+yc^2))*(mu*yAx*yBx+mu-yAx+yBx))*PB-I_B*l^2*(yBx^2+1)*(yAx-mu))*tMx);

H_N_s = HGN;
J_N_s = JGN;


