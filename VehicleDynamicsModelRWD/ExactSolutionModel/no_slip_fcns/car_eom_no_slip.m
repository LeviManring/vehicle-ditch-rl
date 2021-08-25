function [H, J, TA, H_N_ns, J_N_ns, TA_N_ns, tM, dx] = car_eom_no_slip(xA,guess_Dx,params)
%% car_eom_no_slip
% Levi Manring, Duke University
% 2021
%
% This function calculates certain parameters needed for integration of the no_slip_ode.m
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
%       model: x'' + H*x'^2 + J*g + TA*tau_A = 0
%       And solve the normal force as well: F_N = H_N_ns*x'^2 + J_N_ns*g + TA_N_ns*tau_A
%   H: 1x1 double, ode parameter
%   J: 1x1 double, ode parameter
%   TA: 1x1 double, ode parameter
%   H_N_ns: 1x1 double, normal force parameter
%   J_N_ns: 1x1 double, normal force parameter
%   TA_N_ns: 1x1 double, normal force parameter
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


%% Calculating dynamic model parameters
% Applying torque at just back wheel A and including moment of inertia effects for both wheels

% These equations derivation can be seen in Maple files included (copied
% and pasted directly into Matlab)

H = (phiA*l*((-mB*R^2*((tMx^2*yBx+tMxx)*yAx-tMx^2+yBx*tMxx)*phiB*PB^2+I_B*((tMx^2-tMxx*yAx)*phiB-tMx*yAx*phiBx)*(yBx^2+1)*PB+I_B*PBx*tMx*yAx*phiB*(yBx^2+1))*l+(((-tMx^2*xc-tMxx*yc)*yBx+yc*tMx^2-tMxx*xc)*yAx+(tMx^2*yc-tMxx*xc)*yBx+tMx^2*xc+yc*tMxx)*R^2*phiB*PB^2*mM)*cos(tM)^2+l*(phiA*((mB*R^2*phiB*((tMx^2-tMxx*yBx)*yAx+tMx^2*yBx+tMxx)*PB^2+I_B*((tMx^2*yAx+tMxx)*phiB+tMx*phiBx)*(yBx^2+1)*PB-I_B*PBx*tMx*phiB*(yBx^2+1))*l+(((tMx^2*yc-tMxx*xc)*yBx+tMx^2*xc+yc*tMxx)*yAx+(tMx^2*xc+tMxx*yc)*yBx-yc*tMx^2+tMxx*xc)*R^2*phiB*PB^2*mM)*sin(tM)-phiB*(R^2*(PAx*(mA+mM)*yAx^2+(PAx*yBx*mB+yAxx*PA*(mA+mM))*yAx+yBx*yAxx*PA*mB+PAx*(mA+mB+mM))*phiA+I_A*phiAx*(yAx^2+1))*PB^2-I_B*phiA*(yBx^2+1)*(PA*phiBx+PAx*phiB)*PB+I_B*PBx*phiA*phiB*PA*(yBx^2+1))*cos(tM)-(((PAx*yBx*(mA+mB+mM)*yAx^2+(yAxx*PA*(mA+mB+mM)*yBx+PAx*mB)*yAx+PAx*yBx*(mA+mM))*R^2*phiA+I_A*yBx*phiAx*(yAx^2+1))*phiB*PB^2+I_B*yAx*phiA*(yBx^2+1)*(PA*phiBx+PAx*phiB)*PB-I_B*PBx*yAx*phiA*phiB*PA*(yBx^2+1))*l*sin(tM)+phiA*((R^2*phiB*mB*(tMx^2*yBx+tMxx)*PB^2+I_B*(yBx^2+1)*(phiB*tMxx+phiBx*tMx)*PB-I_B*PBx*tMx*phiB*(yBx^2+1))*yAx*l^2-R^2*phiB*PB^2*yBx*((-tMx^2*xc-tMxx*yc)*yAx+yc*tMx^2-tMxx*xc)*mM*l+(yAx-yBx)*R^2*phiB*PB^2*(I_M+mM*(xc^2+yc^2))*tMxx))/(phiA*(-l*((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+((-xc*yBx+yc)*yAx+yc*yBx+xc)*R^2*PB*mM)*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-l*(((yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*PA*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+R^2*PB*(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*tMx)*PB*phiB);

J = -R^2*PB*((((mA+mM)*yAx+mB*yBx)*l-mM*xc*(yAx-yBx))*cos(tM)+sin(tM)*((mA+mB+mM)*yBx*yAx*l+mM*yc*(yAx-yBx)))/(-l*((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+((-xc*yBx+yc)*yAx+yc*yBx+xc)*R^2*PB*mM)*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-l*(((yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*PA*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+R^2*PB*(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*tMx);

TA = PB*l*(yAx^2+1)*(yBx*sin(tM)+cos(tM))/(phiA*(-l*((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+R^2*PB*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+((-xc*yBx+yc)*yAx+yc*yBx+xc)*R^2*PB*mM)*tMx*sin(tM)+(-PA*((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-l*(((yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*PA*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*sin(tM)+(yAx*(R^2*PB*mB+I_B*(yBx^2+1))*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+R^2*PB*(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*tMx));

HGN = (l*(phiB*((-(PA*(yAx+yBx)*tMx^2+((PA*yAxx+PAx*yAx)*yBx-PAx)*tMx-PA*tMxx*(yAx*yBx-1))*mB*l+mM*(-PA*((yAx*yc+xc)*yBx+yAx*xc-yc)*tMx^2+((-yAxx*PA*xc+PAx*(-xc*yAx+yc))*yBx+yc*yAxx*PA+PAx*(yAx*yc+xc))*tMx-PA*tMxx*((-xc*yAx+yc)*yBx+yc*yAx+xc)))*(mA+mM)*phiA*R^2-I_A*((((yAx+yBx)*tMx^2-yAx*yBx*tMxx+tMxx)*mB*l+mM*(((yAx*yc+xc)*yBx+yAx*xc-yc)*tMx^2+tMxx*((-xc*yAx+yc)*yBx+yc*yAx+xc)))*phiA-tMx*phiAx*((-mB*yAx*yBx+mB)*l+mM*((-xc*yAx+yc)*yBx+yc*yAx+xc))))*R^2*PB^2+((mA+mM)*((-PA*tMx^2*yAx-PA*tMxx+PAx*tMx)*phiB-tMx*phiBx*PA)*phiA*R^2-(((tMx^2*yAx+tMxx)*phiB+tMx*phiBx)*phiA-tMx*phiB*phiAx)*I_A)*I_B*l*(yBx^2+1)*PB+phiB*tMx*(PA*(mA+mM)*R^2+I_A)*phiA*I_B*l*PBx*(yBx^2+1))*cos(tM)^2+(l*(phiB*R^2*((mA+mM)*phiA*(((-PA*yAx*yBx+PA)*tMx^2+(PA*yAxx+PAx*yAx+PAx*yBx)*tMx-PA*tMxx*(yAx+yBx))*mB*l+mM*(PA*((-xc*yAx+yc)*yBx+yc*yAx+xc)*tMx^2+((yc*yAxx*PA+PAx*(yAx*yc+xc))*yBx+yAxx*PA*xc-PAx*(-xc*yAx+yc))*tMx-PA*tMxx*((yAx*yc+xc)*yBx+yAx*xc-yc)))*R^2+I_A*((-((yAx*yBx-1)*tMx^2+tMxx*(yAx+yBx))*mB*l+mM*(((-xc*yAx+yc)*yBx+yc*yAx+xc)*tMx^2-tMxx*((yAx*yc+xc)*yBx+yAx*xc-yc)))*phiA+(mB*(yAx+yBx)*l+mM*((yAx*yc+xc)*yBx+yAx*xc-yc))*tMx*phiAx))*PB^2+((mA+mM)*phiA*((PA*tMx^2+(PA*yAxx+PAx*yAx)*tMx-PA*tMxx*yAx)*phiB-tMx*yAx*phiBx*PA)*R^2+I_A*(((tMx^2-tMxx*yAx)*phiB-tMx*yAx*phiBx)*phiA+tMx*yAx*phiB*phiAx))*I_B*l*(yBx^2+1)*PB+phiB*tMx*(PA*(mA+mM)*R^2+I_A)*phiA*I_B*yAx*l*PBx*(yBx^2+1))*sin(tM)+phiB*R^2*((tMx^3*mB*mM*(yBx*yc-xc)*l^2+(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*l+(I_M+mM*(xc^2+yc^2))*mM*tMx^3*(yBx*yc+xc))*phiA*R^2+I_A*l*((-yAxx*PA*mM+yAx*mB*(PA*yAxx+PAx*yAx)*yBx+PAx*yAx*mB-yAxx*PA*mA)*phiA-yAx*phiAx*PA*mB*(yAx*yBx+1)))*PB^2+(phiA*((-tMx^3*l*mM*xc+(I_M+mM*(xc^2+yc^2))*tMx^3-yAxx*PA^2*(mA+mM))*phiB+yAx*phiBx*PA^2*(mA+mM))*R^2+I_A*((PA*phiBx+PAx*phiB)*phiA-phiB*phiAx*PA)*yAx)*I_B*l*(yBx^2+1)*PB-phiB*PA*(PA*(mA+mM)*R^2+I_A)*phiA*I_B*l*yAx*PBx*(yBx^2+1))*cos(tM)+(-phiB*(phiA*(tMx^3*mB*mM*(xc*yBx+yc)*l^2-yBx*(((-xc^2-yc^2)*mM^2+mB*(xc^2+yc^2)*mM+I_M*mB)*tMx^3-yAxx*PA^2*(mA+mM)*(mA+mB+mM))*l+(I_M+mM*(xc^2+yc^2))*mM*tMx^3*(-xc*yBx+yc))*R^2+I_A*((yBx*yAxx*PA*mM+(yAxx*(mA+mB)*PA+PAx*yAx*mB)*yBx+PAx*mB)*phiA-phiAx*PA*mB*(yAx*yBx+1))*l)*R^2*PB^2+((-yc*tMx^3*l*phiB*mM+(-tMx^2*l*mM*xc+(I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM))*phiBx)*phiA*R^2-I_A*((PA*phiBx+PAx*phiB)*phiA-phiB*phiAx*PA))*I_B*l*(yBx^2+1)*PB-phiB*((-tMx^2*l*mM*xc+(I_M+mM*(xc^2+yc^2))*tMx^2-PA^2*(mA+mM))*R^2-I_A*PA)*phiA*I_B*l*PBx*(yBx^2+1))*sin(tM)-phiB*R^2*(phiA*(mB*(mA+mM)*(-PA*tMx^2*yBx-PA*tMxx+PAx*tMx)*l^2+mM*(-PA*((yAx*yc+xc)*yBx*mM+(yc*(mA+mB)*yAx+mA*xc)*yBx+yc*mB)*tMx^2+((-yAxx*PA*xc+PAx*(-xc*yAx+yc))*yBx*mM+(-yAxx*xc*(mA+mB)*PA+PAx*(-xc*(mA+mB)*yAx+yc*mA))*yBx-PAx*mB*xc)*tMx-PA*tMxx*((-xc*yAx+yc)*yBx*mM+(-xc*(mA+mB)*yAx+yc*mA)*yBx-mB*xc))*l+(I_M+mM*(xc^2+yc^2))*(mA+mB+mM)*(((PA*yAxx+PAx*yAx)*yBx+PAx)*tMx-PA*tMxx*(yAx*yBx+1)))*R^2+I_A*((-mB*(tMx^2*yBx+tMxx)*l^2-yBx*mM*((yAx*yc+xc)*tMx^2+tMxx*(-xc*yAx+yc))*l-(I_M+mM*(xc^2+yc^2))*tMxx*(yAx*yBx+1))*phiA+tMx*phiAx*(l^2*mB+(-xc*yAx+yc)*yBx*mM*l+(I_M+mM*(xc^2+yc^2))*(yAx*yBx+1))))*PB^2-I_B*(yBx^2+1)*(phiA*(((mA+mM)*(-PA*tMxx+PAx*tMx)*l^2-mM*(PA*tMx^2*yc-PA*tMxx*xc+PAx*tMx*xc)*l+(-PA*tMxx+PAx*tMx)*(I_M+mM*(xc^2+yc^2)))*phiB+phiBx*((-mA-mM)*l^2-l*mM*xc+mM*(xc^2+yc^2)+I_M)*PA*tMx)*R^2-I_A*((phiB*tMxx+phiBx*tMx)*phiA-tMx*phiB*phiAx)*l^2)*PB+phiB*(((-mA-mM)*l^2-l*mM*xc+mM*(xc^2+yc^2)+I_M)*PA*R^2-I_A*l^2)*tMx*phiA*I_B*PBx*(yBx^2+1))/(PB*phiB*R*(-tMx*l*(R^2*(mB*(yAx+yBx)*l+mM*((yAx*yc+xc)*yBx+yAx*xc-yc))*PB+I_B*yAx*l*(yBx^2+1))*cos(tM)^2+l*(tMx*(R^2*((-mB*yAx*yBx+mB)*l+mM*((-xc*yAx+yc)*yBx+yc*yAx+xc))*PB+I_B*l*(yBx^2+1))*sin(tM)+(-PA*((yAx^2+1)*mM+yAx^2*mA+yAx*yBx*mB+mA+mB)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*cos(tM)-((PA*(yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2+mA)*yBx+yAx*mB)*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+tMx*(R^2*(yAx*l^2*mB+(yAx*yc+xc)*yBx*mM*l+(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*PB+I_B*yAx*l^2*(yBx^2+1))));

JGN = -R*(tMx*(R^2*(yBx*mB*(mA+mM)*l^2-mM*(mA-mB+mM)*(-xc*yBx+yc)*l+((-xc^2+yc^2)*yBx+2*yc*xc)*mM^2)*PB+yc*I_B*l*mM*(yBx^2+1))*cos(tM)^2+(-tMx*(R^2*(mB*(mA+mM)*l^2+mM*(mA-mB+mM)*(yBx*yc+xc)*l+mM^2*(-2*xc*yBx*yc-xc^2+yc^2))*PB+I_B*(yBx^2+1)*l*((mA+mM)*l-mM*xc))*sin(tM)+(-(mA+mB+mM)*PA*((-mA-mM)*l+mM*xc*(yAx*yBx+1))*R^2-I_A*((mB*yAx*yBx-mA-mM)*l+mM*xc*(yAx*yBx+1)))*PB+PA*I_B*(yBx^2+1)*((mA+mM)*l-mM*xc))*cos(tM)+(((mA+mB+mM)*PA*(yBx*(mA+mM)*l+yc*mM*(yAx*yBx+1))*R^2+I_A*(l*(mA+mB+mM)*yBx+yc*mM*(yAx*yBx+1)))*PB+yc*I_B*PA*mM*(yBx^2+1))*sin(tM)+tMx*(R^2*(-mM*(yBx*mM*xc+xc*(mA+mB)*yBx+yc*mB)*l+(xc^2*yBx-xc*yc)*mM^2+yBx*((mA+mB)*xc^2+(mA+mB)*yc^2+I_M)*mM+I_M*yBx*(mA+mB))*PB-yc*I_B*l*mM*(yBx^2+1)))*phiA/(-tMx*l*(R^2*(mB*(yAx+yBx)*l+mM*((yAx*yc+xc)*yBx+yAx*xc-yc))*PB+I_B*yAx*l*(yBx^2+1))*cos(tM)^2+l*(tMx*(R^2*((-mB*yAx*yBx+mB)*l+mM*((-xc*yAx+yc)*yBx+yc*yAx+xc))*PB+I_B*l*(yBx^2+1))*sin(tM)+(-PA*((yAx^2+1)*mM+yAx^2*mA+yAx*yBx*mB+mA+mB)*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*cos(tM)-((PA*(yBx*(yAx^2+1)*mM+((mA+mB)*yAx^2+mA)*yBx+yAx*mB)*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+tMx*(R^2*(yAx*l^2*mB+(yAx*yc+xc)*yBx*mM*l+(I_M+mM*(xc^2+yc^2))*(yAx-yBx))*PB+I_B*yAx*l^2*(yBx^2+1)));

TAN = (-((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+PB*R^2*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*l*cos(tM)^2-(((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+PB*R^2*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*sin(tM)-PA*yAx*(R^2*mB*(yAx*yBx+1)*PB+I_B*(yBx^2+1)))*l*cos(tM)-PA*l*(R^2*mB*(yAx*yBx+1)*PB+I_B*(yBx^2+1))*sin(tM)+tMx*((PB*R^2*mB+I_B*yBx^2+I_B)*l^2+R^2*yBx*PB*mM*(-xc*yAx+yc)*l+PB*(I_M+mM*(xc^2+yc^2))*R^2*(yAx*yBx+1)))/(R*(-((R^2*mB*(yAx+yBx)*PB+I_B*yAx*(yBx^2+1))*l+PB*R^2*((yBx*yc+xc)*yAx+yBx*xc-yc)*mM)*tMx*l*cos(tM)^2+(((-R^2*mB*(yAx*yBx-1)*PB+I_B*(yBx^2+1))*l+PB*R^2*mM*((-xc*yBx+yc)*yAx+yc*yBx+xc))*tMx*sin(tM)+(-((mA+mM)*yAx^2+yAx*yBx*mB+mA+mB+mM)*PA*R^2-I_A*yAx^2-I_A)*PB-I_B*PA*(yBx^2+1))*l*cos(tM)-((PA*(yBx*(mA+mB+mM)*yAx^2+yAx*mB+yBx*(mA+mM))*R^2+I_A*yBx*(yAx^2+1))*PB+I_B*yAx*PA*(yBx^2+1))*l*sin(tM)+tMx*((R^2*PB*mB+I_B*(yBx^2+1))*yAx*l^2+R^2*yBx*PB*mM*(yAx*yc+xc)*l+PB*(I_M+mM*(xc^2+yc^2))*R^2*(yAx-yBx))));

H_N_ns = HGN;
J_N_ns = JGN;
TA_N_ns = TAN;




