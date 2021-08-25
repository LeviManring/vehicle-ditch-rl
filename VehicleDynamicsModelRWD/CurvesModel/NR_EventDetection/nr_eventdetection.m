function [tevent, yevent, failed, nozero] = nr_eventdetection(odehandle, eventhandle, deventhandle, init_interval, init_conds, max_err, max_steps, t_final)
%% nr_eventdetection
% Levi Manring, Duke University
% 2021
%
% This function is used to find an event that occurs during a numerical
% integration of an ode model
%
% Inputs:
%   odehandle: function handle for the dynamics model to integrate
%   eventhandle: function handle for the event to detect
%   deventhandle: function handle for the derivative of the event to detect
%   init_interval: 1x2 array with the initial guess time interval for the Newton-Raphson solver
%   init_conds: 1xN array of the initial conditions needed to integrate the dynamics model
%   max_err: 1x1 double dictating the maximum percent error allowed for convergence
%   max_steps: 1x1 double that tells the Newton-Raphson solver how many
%           interations it is allowed to perform before a failure to converge
%           condition is determined
%
% Outputs:
%   tevent: 1x1 double that tells the exact time of the event if the NR
%           routine converged, or gives the initial time if the NR routine failed
%   yevent: 1xN array of the states of the dynamics model at the exact
%           event moment if the NR routine converged, or the init_conds if
%           the NR routine failed
%   failed: 1x1 boulean indicating 1 if the NR routine failed or 0 if the
%           NR routine converged
%   nozero: 1x1 double telling if the zero-crossing was found due to a
%           likely integration error

%%
t0_nr = init_interval(1); % this is our initial time for integration (does not change)
tk = init_interval(2); % this is our initial guess for the event location, this gets updated every NR iteration
err = 1;
n_count = 1;
nozero = 0;
failed = 0;
while err > max_err
    [~, y_nr] = ode45(odehandle,[t0_nr, tk], init_conds);
    fprintf('%3.0f nr-count \n',n_count);
    eventk = eventhandle(y_nr(end,:));
    deventk = deventhandle(y_nr(end,:));
    tk_new = tk - eventk/deventk;
    err = abs((tk_new - tk))/tk;
    tk = tk_new;
    
    if n_count > max_steps || isnan(tk) % check to see if our maximum allowed NR steps occurred or if there is an error in our tk
        fprintf('The Newton Raphson Routine has failed to converge \n');
        failed = 1;
        break
    else
        failed = 0;
    end
    
    if abs(t0_nr-tk) < 1e-6 || tk > t_final || tk <= t0_nr % prohibits an error that can kill the dynamics solver
        nozero = 1;
        break
    end
    
    n_count = n_count + 1;
    
end

if failed == 0 && nozero == 0 % if the NR routine was a success, output the event
    [~, y_nr] = ode45(odehandle,[t0_nr, tk], init_conds);
    tevent = tk;
    yevent = y_nr(end,:);
else % if the NR routine failed, output the initial inputs to nr_eventdetection
    tevent = t0_nr;
    yevent = init_conds;
end


