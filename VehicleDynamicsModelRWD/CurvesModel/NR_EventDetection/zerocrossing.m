function [zero_crossing, guess_ind] = zerocrossing(tcondition,threshhold)
%% zerocrossing
% Levi Manring, Duke University
% 2021
%
% This function is used to find a zero crossing in an event vector
%
% Inputs:
%   tcondition: Nx1 array of a function that we want to check for a zero crossing
%   threshhold: 1x1 double that prohibits a super small function change at
%               the first indice of tcondition (this can help avoid finding false
%               zerocrossings due to integration error)
%
% Outputs:
%   zero_crossing: 1x1 boulean that tells if a zero crossing occurred (1
%               means a zero crossing occurred, 0 means one did not occur)
%   guess_ind: 1x1 double that gives a guess index as to the closest index
%               in tcondition where the zero crossing occurred. In particular, in
%               particular if the zero crossing happened between indices 30 and 31 in
%               tcondition, guess_ind would be 31. If there is no zero
%               crossing in tcondition, guess_ind is NaN

%%

for k = 1:(length(tcondition) - 1)
    prior = tcondition(k,1);
    post = tcondition(k+1,1);
    if prior*post < 0 % check for a change in sign
        if k == 1 && abs(prior) < threshhold % include a condition if the first step condition is very small
            zero_crossing = 0;
        else
            zero_crossing = 1;
            guess_ind = k+1;
            break
        end
        
    else
        zero_crossing = 0;
        guess_ind = NaN;
    end
end