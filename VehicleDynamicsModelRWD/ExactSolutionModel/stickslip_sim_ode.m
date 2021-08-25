function [states, slip_cond_end, event_num] = stickslip_sim_ode(init_conds,slip_cond_start,current_Action,t_int,guess_Dx,params,refine)
%% stickslip_sim_ode
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
%   guess_Dx: 1x1 double that provides a guess for the location of wheel B
%           relative to wheel A in the x-direction. On a flat surface this would be
%           equal to the length of the vehicle (params.l). This guess (if
%           close) can help the solvers converge faster.
%   params: a parameter structure including the masses, moments of
%           inertia, and dimensions needed to define the planar vehicle model
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

%%

% Initialization
y0 = init_conds;
t_init = t_int(1,1);
t_final = t_int(1,2);
event_num = 0;

% Unpack parameters
mu_s = params.dim.mu_s;
friction_model = params.settings.friction_model;


% establish the terminal condition for the ODE's
options_no_slip = odeset('Events',@(t,y) slip_event(t,y,current_Action,mu_s,guess_Dx,params),'OutputSel',1,...
    'Refine',refine);
options_slip = odeset('Events',@(t,y) no_slip_event(t,y,params),'OutputSel',1,...
    'Refine',refine);

% setup the function handles
slipode = @(t,y) slip_ode(t,y,current_Action,guess_Dx,params,mu_s,friction_model);
noslipode = @(t,y) no_slip_ode(t,y,current_Action,guess_Dx,params);

% check to see if it is slipping at the beginning of the time step with the
% current action
if slip_cond_start == 0
    slip_cond_start = startslipcheck_fun(init_conds,current_Action,guess_Dx,params);
end

% Switch odehandles and eventhandles depending on what the starting
% condition is
switch slip_cond_start
    
    case 1 % Starts the timestep slipping
        
        firstode = slipode;
        secondode = noslipode;
        
        firstoptions = options_slip;
        secondoptions = options_no_slip;
        
    case 0 % Starts the timestep not slipping
        
        firstode = noslipode;
        secondode = slipode;
        
        firstoptions = options_no_slip;
        secondoptions = options_slip;
end

% Now simulate
while t_init < t_final
    
    [t,y,~,~,~] = ode45(firstode, [t_init, t_final], y0, firstoptions);
    
    % It goes back to sticking after slip_ode
    
    event_num = event_num + 1;
    
    if t(end) == t_final
        slip_cond_end = 1;
        break
    else
        t_init = t(end);
    end
    y0 = y(end,:);
    
    % Now we have switched our slipping condition, so integrate these dynamics
    [t,y,~,~,~] = ode45(secondode, [t_init t_final], y0, secondoptions);
    
    if t(end) == t_final
        slip_cond_end = 0;
        break
    else
        t_init = t(end);
    end
    y0 = y(end,:);
    
    event_num = event_num + 1;
    
end

% Calculate relative velocity
v_rel = v_rel_fcn(y(end,:),params);

% Update states
states = [y(end,:), v_rel];