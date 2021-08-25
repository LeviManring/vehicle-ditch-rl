function [states, slip_cond_end, event_num, failed] = curves_stickslip_nr(init_conds,slip_cond_start,current_Action,t_int,params,curves,refine)
%% curves_stickslip_nr
% Levi Manring, Duke University
% 2021
%
% This function is used to integrate the discontinuous dynamics model
% combination of not-slipping/slipping.
%
% Inputs:
%   init_conds: 1xN array indicating the initial states of the dynamics
%           model at the beginning of a given timestep
%   slip_cond_start: 1x1 boulean describing whether or not the vehicle is
%           slipping at the start of the time step, 1 means it is slipping,
%           0 means it is not slipping
%   current_Action: 1x1 double inidicating the starting applied torque at
%           wheel A at the beginning of a given timestep
%   t_int: 2x1 array describing the time interval to intergrate over
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
%   curves: a curves structure including the interpolated dynamics model
%           parameters needed to define the planar vehicle model
%   refine: 1x1 double which is a refine parameter used as a setting in
%           Matlab's ODE solver, set to 20 typically for this situation
%
% Outputs:
%   states: 1x(N+1) array indicating the states of the dynamics model at
%           the end of a given timestep (x, x', theta_A, theta_A', v_relA)
%   slip_cond_end: 1x1 boulean describing whether or not the vehicle is
%           slipping, 1 means it is slipping, 0 means it is not slipping
%   event_num: 1x1 double describing the number of events that occurred
%           during a timestep
%   failed: 1x1 boulean describing whether or not the integration failed
%           based on a failed Newton-Raphson convergence

%%

% Initialization
y0 = init_conds;
t_init = t_int(1,1);
t_final = t_int(1,2);
event_num = 0;
event_num_threshhold = 4;
max_err = 1e-7;
v_err = 1e-7;
mu_err = 1e-7;
max_nr_steps = 10;
failed = 0;

% Unpack parameters
mu_s = params.dim.mu_s;
fmodel = params.settings.friction_model;

% create a function handle for the ode solvers and terminal conditions
slipode     = @(t,y) curves_slip_ode(t,y,current_Action,params,mu_s,fmodel,curves);
noslipode   = @(t,y) curves_no_slip_ode(t,y,current_Action,params,curves);
vkfun       = @(y) v_rel_fcn(y,params);
vdotkfun    = @(y) curves_dv_rel_fcn(y,current_Action,params,curves);
Fkfun      = @(y) mu_s*abs(curves_fnormal_ns(y(1),y(2),current_Action,params,curves)) - ...
    abs(curves_ffriction_ns(y(1),y(2),current_Action,params,curves));
Fdotkfun   = @(y) mu_s*curves_fnormal_ns(y(1),y(2),current_Action,params,curves)*...
    curves_dfnormal_ns(y(1),y(2),current_Action,params,curves)/abs(curves_fnormal_ns(y(1),y(2),current_Action,params,curves))- ...
    curves_ffriction_ns(y(1),y(2),current_Action,params,curves)*...
    curves_dffriction_ns(y(1),y(2),current_Action,params,curves)/abs(curves_ffriction_ns(y(1),y(2),current_Action,params,curves));

% Set the ode options
ode_options = odeset('Refine',refine);

% check to see if it is slipping at the beginning of the time step with the
% current action
slip_cond_end = slip_cond_start;
if slip_cond_start == 0
    slip_cond_start = curves_startslipcheck_fun(y0,current_Action,params,curves);
end

% Switch odehandles and eventhandles depending on what the starting
% condition is
switch slip_cond_start
    
    case 1 % Starts the timestep slipping
        
        firstode = slipode;
        secondode = noslipode;
        
        termcond1 = vkfun;
        dtermcond1 = vdotkfun;
        err1 = v_err;
        
        termcond2 = Fkfun;
        dtermcond2 = Fdotkfun;
        err2 = mu_err;
        
    case 0 % Starts the timestep not slipping
        
        firstode = noslipode;
        secondode = slipode;
        
        termcond1 = Fkfun;
        dtermcond1 = Fdotkfun;
        err1 = mu_err;
        
        termcond2 = vkfun;
        dtermcond2 = vdotkfun;
        err2 = v_err;
        
end

% Now simulate
while t_init < t_final
    % Simulate system over entire timestep
    [t, y] = ode45(firstode, [t_init, t_final], y0, ode_options);
    
    % calculate terminal condition
    term1 = zeros(length(y),1);
    for k = 1:length(y)
        term1(k,1) = termcond1(y(k,:));
    end
    
    % Check for zero-crossing of terminal condition
    [zero_crossing, guess_ind] = zerocrossing(term1,err1);
    
    % Check for a zero crossing, and solve for it
    if zero_crossing == 1
        ind = guess_ind;
        tk = t(ind,1); % initial guess for the zero crossing
        
        % altering the starting point for the NR integration to reduce integration time
        if ind > 2
            nr_ind = ind - 2;
        else
            nr_ind = 1;
        end
        y0_nr = y(nr_ind,:);
        t0_nr = t(nr_ind,1);
        
        % Newton Raphson Solver to find zero crossing
        [tevent, yevent, failed, nozero] = nr_eventdetection(firstode, termcond1, dtermcond1, [t0_nr, tk], y0_nr, max_err, max_nr_steps, t_final);
        
        % Check to see if Newton-Raphson solver failed to converge
        if failed == 1 || nozero == 1
            slip_cond_end = slip_cond_start;
            break
        end
        
        event_num = event_num + 1;
        if event_num > event_num_threshhold
            failed = 1;
            break
        end
        
        % Reset initial conditions with final value from NR routine
        t_init = tevent;
        y0 = yevent;
        if abs(t_final - t_init) < 1e-6 % prohibits too small timestep integration
            break
        end
        
    else % no zero crossing, the initial solution is fine, so terminate
        slip_cond_end = slip_cond_start;
        break
    end
    
    % Now we have switched our slipping condition, so integrate these dynamics
    
    [t, y] = ode45(secondode, [t_init, t_final], y0, ode_options);
    
    % Check the terminal conditions
    term2 = zeros(length(y),1);
    for k = 1:length(y)
        term2(k,1) = termcond2(y(k,:));
    end
    
    [zero_crossing, guess_ind] = zerocrossing(term2,err2);
    
    % Now check for a zero crossing and solve for it
    if zero_crossing == 1
        ind = guess_ind;
        tk = t(ind,1); % initial guess for the zero crossing
        
        % altering the starting point for the NR integration to reduce integration time
        if ind > 2
            nr_ind = ind - 2;
        else
            nr_ind = 1;
        end
        y0_nr = y(nr_ind,:);
        t0_nr = t(nr_ind,1);
        
        % Newton Raphson Solver to find zero crossing
        [tevent, yevent, failed, nozero] = nr_eventdetection(secondode, termcond2, dtermcond2, [t0_nr, tk], y0_nr, max_err, max_nr_steps, t_final);
        
        % Check to see if Newton-Raphson solver failed to converge
        if failed == 1 || nozero == 1
            slip_cond_end = ~slip_cond_start;
            break
        end
        
        event_num = event_num + 1;
        if event_num > event_num_threshhold
            failed = 1;
            break
        end
        
        y0 = yevent;
        t_init = tevent;
        if abs(t_final - t_init) < 1e-6 % prohibits too small timestep integration
            break
        end
        
    else % no zero crossing, the initial solution is fine, so terminate
        slip_cond_end = ~slip_cond_start;
        
        break
    end
    
    
    
end

% Calculate relative velocity
v_rel = v_rel_fcn(y(end,:),params);

% Update states
states = [y(end,:), v_rel];