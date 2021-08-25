function dy = pilcoRWD_ode(t,y,tauA,params,curves)
%% pilcoRWD_ode
% Levi Manring, Duke University
% 2021
%
% This function is used by PILCO to integrate the RWD vehicle dynamics.
%
% Inputs:
%   t: 1x1 double inidicating the timestep of ode integration
%   y: 1xN array inticating the states of the dynamics model from the ode.
%           In particular y(1) = x, y(2) = x'
%   tauA: 1x1 double indicating the torque applied to wheel A
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%
% Outputs: 
%   dy: 1xN array describing the outputs from the dynamic model statespace derivatives
% 
% This code is modified for our scenario from 
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
%%

% Calculate the dynamics model parameters
G = curves.G(y(1));
H = curves.H(y(1));
J = curves.J(y(1));
Ta = curves.TA(y(1));

% Calculate the dynamic model states
% states: x, xdot

dy = zeros(2,1); 

dy(1) = y(2);

dy(2) = (-1/G)*(H*(y(2)^2) + J*params.dim.g + Ta*tauA(t));
