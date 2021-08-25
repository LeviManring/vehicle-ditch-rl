%% GetCurves.m
% Levi Manring, Duke University
% 2021
%
% This m-file calculates the curefitted version of the ODEeventModel
% dynamics models. The result is a curves structure that includes 
% interpolated functions for the parameters needed for the dynamic curves model.
%
% Run this function to calculate the curves prior to running the curves
% models. This function generally can take a significant amount of time to
% run, depending on the refinement of the interpolation mesh.

%% Initialize
clear all; close all; clc;
% Generate the paths for the RWD ehicle dynamics model
% addpath(genpath('.../RL Car/VehicleDynamicsModelRWD'));

% load VehicleParameters if these are not loaded previously
load('VehicleParameters'); 

%% Get Dynamic Curve Fit Model Parameters

% Set the interpolation range for x and mu (note a finer mesh selection is
% chosen for the region of the curve that varies the most)
xlow = (-20:0.1:-7.1)';
xmid = (-7:0.001:7)';
xhigh = (7.1:0.1:20)';
xrange = [xlow; xmid; xhigh];

% xlow = (-8:0.1:-5.1)';
% xmid = (-5:0.001:5)';
% xhigh = (5.1:0.1:10)';
% xrange = [xlow; xmid; xhigh];

murange = (-0.3:0.01:0.3)';

Lx = length(xrange);
Lmu = length(murange);

% Initialize the curves with zeros
H = zeros(Lx,1);
J = zeros(Lx,1);
TA = zeros(Lx,1);
G_N_ns = zeros(Lx,1);
H_N_ns = zeros(Lx,1);
J_N_ns = zeros(Lx,1);
TA_N_ns = zeros(Lx,1);
tM = zeros(Lx,1);
H_s = zeros(Lx,Lmu);
J_s = zeros(Lx,Lmu);
G_N_s = zeros(Lx,Lmu);
H_N_s = zeros(Lx,Lmu);
J_N_s = zeros(Lx,Lmu);

guess_Dx = params.dim.l;
% Run through interpolated rates for x and mu for the noslip and slip
% dynamics models
for k = 1:length(xrange)
    x = xrange(k,1);
    
    % Calculate no-slip dynamic model parameters
    [H(k,1), J(k,1), TA(k,1), H_N_ns(k,1), J_N_ns(k,1), TA_N_ns(k,1), tM(k,1), dx] = car_eom_no_slip(x,guess_Dx,params);
    
    % Calculate slip dynamic model parameters
    for m = 1:length(murange)
        mu = murange(m,1);
        [H_s(k,m), J_s(k,m), H_N_s(k,m), J_N_s(k,m), ~, ~] = car_eom_slip(x,guess_Dx,mu,params);
    end
    guess_Dx = dx;
    
    % provide a status update
    fprintf('k %3.0f out of %3.0f \n',k,Lx);
end


% Create interpolated functions stored in curves structure for no-slip model
curves.H = fit(xrange,H,'linearinterp');
curves.J = fit(xrange,J,'linearinterp');
curves.TA = fit(xrange,TA,'linearinterp');
curves.H_N_ns = fit(xrange,H_N_ns,'linearinterp');
curves.J_N_ns = fit(xrange,J_N_ns,'linearinterp');
curves.TA_N_ns = fit(xrange,TA_N_ns,'linearinterp');
curves.tM = fit(xrange,tM,'linearinterp');

% Create interpolated functions stored in curves structure for slip model
[mumesh,xmesh] = meshgrid(murange,xrange);
curves.H_s = fit([xmesh(:), mumesh(:)],H_s(:),'linearinterp'); fprintf('1 done\n');
curves.J_s = fit([xmesh(:), mumesh(:)],J_s(:),'linearinterp'); fprintf('1 done\n');
curves.H_N_s = fit([xmesh(:), mumesh(:)],H_N_s(:),'linearinterp'); fprintf('1 done\n');
curves.J_N_s = fit([xmesh(:), mumesh(:)],J_N_s(:),'linearinterp'); fprintf('1 done\n');

%% Get Friction and Normal Force and their derivatives in an interpolated model for Newton-Raphson solver

% Set interpolation range
h = 0.001;
% xrange_2 = (-17:h:17)';
xrange_2 = (-10:h:10)';
Lx = length(xrange_2);
H = zeros(Lx,1);
J = zeros(Lx,1);
TA = zeros(Lx,1);
H_N_ns = zeros(Lx,1);
J_N_ns = zeros(Lx,1);
TA_N_ns = zeros(Lx,1);
FF_H_ns = zeros(Lx,1);
FF_J_ns = zeros(Lx,1);
FF_TA_ns = zeros(Lx,1);
for k = 1:length(xrange_2)
    x = xrange_2(k,1);
    
    % Calculate no-slip dynamic model parameters
    [H(k,1), J(k,1), TA(k,1), H_N_ns(k,1), J_N_ns(k,1), TA_N_ns(k,1), ~, dx] = car_eom_no_slip(x,guess_Dx,params);
    
    % Get surface profile derivatives
    [~, yAx, yAxx, ~, ~] = y_eval_fcn(x,params);
    phiA = sqrt(yAx^2+1)/params.dim.R;
    phiAx = yAx*yAxx/(params.dim.R*sqrt(yAx^2+1));
    
    % Calculate Friction Force parameters
    FF_H_ns(k,1) = (params.dim.I_A/params.dim.R)*(phiA*H(k,1) - phiAx);
    FF_J_ns(k,1) = (params.dim.I_A/params.dim.R)*(phiA*J(k,1));
    FF_TA_ns(k,1) = (1/params.dim.R + (params.dim.I_A/params.dim.R)*phiA*TA(k,1));
    
    guess_Dx = dx;
    
    % Provide status update
    fprintf('k %3.0f out of %3.0f \n',k,Lx);
end

FN_H_ns = H_N_ns;
FN_J_ns = J_N_ns;
FN_TA_ns = TA_N_ns;

% Initialize the derivatives
FN_H_x_ns = zeros(Lx-4,1);
FN_J_x_ns = zeros(Lx-4,1);
FN_TA_x_ns = zeros(Lx-4,1);
FF_H_x_ns = zeros(Lx-4,1);
FF_J_x_ns = zeros(Lx-4,1);
FF_TA_x_ns = zeros(Lx-4,1);

% Compute numerical derivatives of the Friction and Normal Forces
m = 1;
for k = 3:(Lx-2)
    FN_H_x = (-FN_H_ns(k+2,1) + 8*FN_H_ns(k+1,1) - 8*FN_H_ns(k-1,1) + FN_H_ns(k-2,1))/(12*h);
    FN_J_x = (-FN_J_ns(k+2,1) + 8*FN_J_ns(k+1,1) - 8*FN_J_ns(k-1,1) + FN_J_ns(k-2,1))/(12*h);
    FN_TA_x = (-FN_TA_ns(k+2,1) + 8*FN_TA_ns(k+1,1) - 8*FN_TA_ns(k-1,1) + FN_TA_ns(k-2,1))/(12*h);
    FF_H_x = (-FF_H_ns(k+2,1) + 8*FF_H_ns(k+1,1) - 8*FF_H_ns(k-1,1) + FF_H_ns(k-2,1))/(12*h);
    FF_J_x = (-FF_J_ns(k+2,1) + 8*FF_J_ns(k+1,1) - 8*FF_J_ns(k-1,1) + FF_J_ns(k-2,1))/(12*h);
    FF_TA_x = (-FF_TA_ns(k+2,1) + 8*FF_TA_ns(k+1,1) - 8*FF_TA_ns(k-1,1) + FF_TA_ns(k-2,1))/(12*h);
    
    FN_H_x_ns(m,1) = (FN_H_x - 2*FN_H_ns(k,1)*H(k,1));
    FN_J_x_ns(m,1) = (FN_J_x - 2*FN_H_ns(k,1)*J(k,1));
    FN_TA_x_ns(m,1) = (FN_TA_x - 2*FN_H_ns(k,1)*TA(k,1));
    
    FF_H_x_ns(m,1) = (FF_H_x - 2*FF_H_ns(k,1)*H(k,1));
    FF_J_x_ns(m,1) = (FF_J_x - 2*FF_H_ns(k,1)*J(k,1));
    FF_TA_x_ns(m,1) = (FF_TA_x - 2*FF_H_ns(k,1)*TA(k,1));
    
    m = m + 1;
end

% Create interpolated functions stored in a curves structure
xrange_3 = xrange_2(3:end-2);
curves.FN_H_x_ns = fit(xrange_3,FN_H_x_ns,'linearinterp');
curves.FN_J_x_ns = fit(xrange_3,FN_J_x_ns,'linearinterp');
curves.FN_TA_x_ns = fit(xrange_3,FN_TA_x_ns,'linearinterp');
curves.FF_H_x_ns = fit(xrange_3,FF_H_x_ns,'linearinterp');
curves.FF_J_x_ns = fit(xrange_3,FF_J_x_ns,'linearinterp');
curves.FF_TA_x_ns = fit(xrange_3,FF_TA_x_ns,'linearinterp');

% Save the curve structure, which stores all the information needed to run
% the curve-fit dynamics model
save('CurveFitModel','curves');


