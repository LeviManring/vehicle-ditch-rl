function mu = mu_fcn(v_rel,mu_s,model)
%% mu_fcn
% Levi Manring, Duke University
% 2021
%
% This function is used to calculate a value for the directional friction 
% coefficient based on a value for the static friction coefficient and the
% relative velocity between wheel A and the surface profile.
%
% Inputs:
%   v_rel: 1x1 double describing the relative velocity between wheel A and the surface profile
%   mu_s: 1x1 double describing the static friction coefficient
%   model: 1x1 double selecting which friction model to choose for the scenario
%
% Outputs:
%   mu: 1x1 double indicating the calculated friction coefficient

%%

% We only have one model here, others can be created and selected using the
% model number selector
if model == 2
    mu = tanh(40*v_rel)*(exp(-abs(v_rel)) + mu_s)/4;
end
