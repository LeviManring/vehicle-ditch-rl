function [states, cond_end, event_num, failed] = exactsoln_dynamicsAWD_nr(init_conds,cond_start,tauA,tauB,t_int,guess_Dx,params,refine)
%% exactsoln_dynamicsAWD_nr
% Levi Manring, Duke University
% 2021
%
% This function is used to integrate and switch between four different
% dynamic models of a vehicle moving on a user-defined surface profile.
% These four scenarios are representative of the four combinations of
% sticking/slipping that can occur for our planar vehicle model with rear
% wheel A and front wheel B. These four conditions are listed as:
% 1 -> NoSlip (neither wheel A or wheel B are slipping)
% 2 -> AllSlip (both wheel A and wheel B are slipping)
% 3 -> RearSlipFrontStick (wheel A is slipping and wheel B is not)
% 4 -> RearStickFrontSlip (wheel A is not slipping and wheel B is slipping)
%
% Inputs:
%   init_conds: 1xN array inticating the initial states of the dynamics model from the ode.
%           In particular:
%               init_conds(1) = x
%               init_conds(2) = x'
%               init_conds(3) = theta_A
%               init_conds(4) = theta_A'
%               init_conds(5) = theta_B
%               init_conds(6) = theta_B'
%   cond_start: 1x1 double indicating the dynamic scenario at the start of the timestep:
%           1 -> noslip, 2 -> allslip, 3 -> rearslipfrontstick, 4 -> rearstickfrontslip
%   tauA: 1x1 double indicating the torque applied to wheel A
%   tauB: 1x1 double indicating the torque applied to wheel B
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
%   states: 1x(N+2) array indicating the states of the dynamics model at
%           the end of a given timestep:
%               states(1) = x
%               states(2) = x'
%               states(3) = theta_A
%               states(4) = theta_A'
%               states(5) = theta_B
%               states(6) = theta_B'
%               states(7) = vrel_A
%               states(8) = vrel_B
%   cond_end: 1x1 double indicating the dynamic scenario at the end of the timestep:
%           1 -> noslip, 2 -> allslip, 3 -> rearslipfrontstick, 4 -> rearstickfrontslip
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
% v_err = 1e-7;
% f_err = 1e-7;
max_nr_steps = 15;
failed = 0;
ode_options = odeset('Refine',refine);

%% Determine the cond_start based on checking the friction condition with the applied torques
% Even though we are told the cond_start via function input, with the
% instantanenously applied torques at wheel A and wheel B at the beginning of the
% timestep, it could change the slip condition at that moment.
switch cond_start
    case 1 % noslip
        % if noslipping (default for start of simulation), check to see if
        % either tauA or tauB or both make it slip (switches to case 2,3,
        % or 4)
        FeventA = frictionterminalcond_exactsoln(y0,tauA,tauB,cond_start,'A',guess_Dx,params);
        FeventB = frictionterminalcond_exactsoln(y0,tauA,tauB,cond_start,'B',guess_Dx,params);
        if FeventA < 0 % wheel A starts stlipping
            cond_start = 3;
        end
        if FeventB < 0 % wheel B starts slipping
            cond_start = 4;
        end
        if FeventA < 0 && FeventB < 0 % both wheels slip
            cond_start = 2;
        end
        % no case 2 checked here because both wheel are already slipping
    case 3 % rearslipfrontstick
        % check if tauB makes the front wheel slip (switches to case 2)
        FeventB = frictionterminalcond_exactsoln(y0,tauA,tauB,cond_start,'B',guess_Dx,params);
        if FeventB < 0
            cond_start = 2;
        end
    case 4 % rearstickfrontslip
        % check to see if tauA makes the rear wheel slip (switches to case 2)
        FeventA = frictionterminalcond_exactsoln(y0,tauA,tauB,cond_start,'A',guess_Dx,params);
        if FeventA < 0
            cond_start = 2;
        end
end
% now we have the accurate cond_start

%% Start simulation
% for each case, define terminal event functions and derivative of terminal
% event functions
% define ode_dynamics functions

% Steps:
% 1 simulate ode
% 2 calculate both event functions
% 3 check to see which one crossed first, if they are close, interpolate to
% approximate which one crossed first, if no event happened, finish loop
% 4 use Newton-Raphson routine to integrate up until the end of the true terminal event
% 5 switch case
% 6 Back to 1

while t_init < t_final
    % Simulate selected ode based on cond_start
    odestart = @(t,y) exactsoln_selector_ode(t,y,tauA,tauB,cond_start,params);
    [t, y] = ode45(odestart, [t_init, t_final], y0, ode_options);
    term1 = zeros(length(y),1);
    term2 = zeros(length(y),1);
    
    % Grab the correct set of terminal event functions and derivative of
    % terminal event functions based on the case
    switch cond_start
        case 1 % no slip
            % Terminal conditions:
            % Wheel A slips: mu*|FNA|-|FfA| = 0, we switch to case 3
            % Wheel B slips: mu*|FNB|-|FfB| = 0, we switch to case 4
            eventhandle1 = @(y) frictionterminalcond_exactsoln(y,tauA,tauB,cond_start,'A',guess_Dx,params);
            eventhandle2 = @(y) frictionterminalcond_exactsoln(y,tauA,tauB,cond_start,'B',guess_Dx,params);
            deventhandle1 = @(y) dfrictionterminalcond_exactsoln(y,tauA,tauB,cond_start,'A',guess_Dx,params);
            deventhandle2 = @(y) dfrictionterminalcond_exactsoln(y,tauA,tauB,cond_start,'B',guess_Dx,params);
%             eventerr1 = f_err;
%             eventerr2 = v_err;
            caseswitch1 = 3;
            caseswitch2 = 4;
        case 2 % allslip
            % Terminal conditions:
            % Wheel A sticks: vrel_A = 0, we switch to case 4
            % Wheel B sticks: vrel_B = 0, we switch to case 3
            eventhandle1 = @(y) velocityterminalcond_exactsoln(y,'A',guess_Dx,params);
            eventhandle2 = @(y) velocityterminalcond_exactsoln(y,'B',guess_Dx,params);
            deventhandle1 = @(y) dvelocityterminalcond_exactsoln(y,tauA,tauB,cond_start,'A',guess_Dx,params);
            deventhandle2 = @(y) dvelocityterminalcond_exactsoln(y,tauA,tauB,cond_start,'B',guess_Dx,params);
%             eventerr1 = v_err;
%             eventerr2 = v_err;
            caseswitch1 = 4;
            caseswitch2 = 3;
        case 3 % rearslipfrontstick
            % Terminal conditions:
            % Wheel A sticks: vrel_A = 0, we switch to case 1
            % Wheel B sliips: mu*|FNB| - |FfB| = 0, we switch to case 2
            eventhandle1 = @(y) velocityterminalcond_exactsoln(y,'A',guess_Dx,params);
            eventhandle2 = @(y) frictionterminalcond_exactsoln(y,tauA,tauB,cond_start,'B',guess_Dx,params);
            deventhandle1 = @(y) dvelocityterminalcond_exactsoln(y,tauA,tauB,cond_start,'A',guess_Dx,params);
            deventhandle2 = @(y) dfrictionterminalcond_exactsoln(y,tauA,tauB,cond_start,'B',guess_Dx,params);
%             eventerr1 = v_err;
%             eventerr2 = f_err;
            caseswitch1 = 1;
            caseswitch2 = 2;
        case 4 % rearstickfrontslip
            % Terminal conditions:
            % Wheel A slips: mu*|FNA|-|FfA| = 0, we switch to case 2
            % Wheel B sticks: vrel_B = 0, we switch to case 1
            eventhandle1 = @(y) frictionterminalcond_exactsoln(y,tauA,tauB,cond_start,'A',guess_Dx,params);
            eventhandle2 = @(y) velocityterminalcond_exactsoln(y,'B',guess_Dx,params);
            deventhandle1 = @(y) dfrictionterminalcond_exactsoln(y,tauA,tauB,cond_start,'A',guess_Dx,params);
            deventhandle2 = @(y) dvelocityterminalcond_exactsoln(y,tauA,tauB,cond_start,'B',guess_Dx,params);
%             eventerr1 = f_err;
%             eventerr2 = v_err;
            caseswitch1 = 2;
            caseswitch2 = 1;
    end
    
    % Calculate the terminal event function values
    for k = 1:length(y)
        term1(k,:) = eventhandle1(y(k,:));
        term2(k,:) = eventhandle2(y(k,:));
    end
    
    % Now check to see if we have zero crossings (did an event occur?)
    [zero_crossing1, guess_ind1] = zerocrossing(term1(:,1));
    [zero_crossing2, guess_ind2] = zerocrossing(term2(:,1));
    
    if ~zero_crossing1 && ~zero_crossing2 % if no zero crossings, terminate
        
        break
        
    elseif zero_crossing1 && zero_crossing2 % if both events happen, determine which one happened first
        
        if guess_ind1 ~= guess_ind2 % if one happened first, use that one
            
            [~, eventtype] = min([guess_ind1, guess_ind2]);
            
        else % otherwise if guess_ind is the same for both, use the slope at the event to choose one
            fprintf('Both events happened at the same guess index! \n'); % Choose to output this unusual scenario
            
            slope1 = abs(term1(guess_ind1 - 1,1) - term1(guess_ind1,1));
            slope2 = abs(term2(guess_ind2 - 1,1) - term2(guess_ind2,1));
            [~, eventtype] = max([slope1, slope2]);
            
        end
        
    elseif zero_crossing1 % just the first event happened
        
        eventtype = 1;
        
    elseif zero_crossing2 % just the second event happened
        
        eventtype = 2;
        
    end
    
    % Now we know what event happened (if any)
    % Grab the correct event handle to use for the Newton-Raphson solver
    switch eventtype
        case 1
            eventhandle = deventhandle1;
            ind = guess_ind1;
            caseswitch = caseswitch1;
        case 2
            eventhandle = deventhandle2;
            ind = guess_ind2;
            caseswitch = caseswitch2;
    end
    
    % Grab initial time guess for the Newton-Raphson starting point
    tk = t(ind,1);
    
    % initial starting index to use as the initial point for integration
    if ind > 2
        nr_ind = ind - 2;
    else
        nr_ind = 1;
    end
    y0_nr = y(nr_ind,:);
    t0_nr = t(nr_ind,1);
    
    % Solve for the event
    [tevent, yevent, failed, nozero] = nr_eventdetectionAWD(odestart, eventhandle, [t0_nr, tk], y0_nr, max_err, max_nr_steps, t_final);
    % Check to see if Newton-Raphson solver failed to converge
    if failed == 1 || nozero == 1
        break
    end
    
    % Update the initial conditions to being the final states when the
    % terminal event occurred
    y0 = yevent;
    t_init = tevent;
    if abs(t_final - t_init) < 1e-6 % prohibits too small timestep integration
        break
    end
    
    event_num = event_num + 1;
    if event_num > event_num_threshhold
        failed = 1;
        break
    end
    % Switch to the new case/scenario
    cond_start = caseswitch;
    
end

% Output the updated slipping condition
cond_end = cond_start;

% Calculate the relative velocities at wheels A and B
v_relA = v_relA_fcn(y(end,:),params);
v_relB = v_relB_fcn(y(end,:),guess_Dx,params);

% Output the augmented states
states = [y(end,:), v_relA, v_relB];
