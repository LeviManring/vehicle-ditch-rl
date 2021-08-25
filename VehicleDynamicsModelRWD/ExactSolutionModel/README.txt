ExactSolutionModel

Levi Manring
7/5/21

This folder contains no_slip_fcns and slip_fcns folders. Each of these folders contains the ode dynamics file for the slip/no_slip dynamics scenario, the parameters for the equation's of motion file (eom), and the event function for transitioning between slip/no_slip dynamics scenarios.

This directory contains the model information to numerically integrate the entire model at each time step. This is a slower, but more accurate, method of integrating the dynamic model for the vehicle.

In addition, this folder includes startslipcheck_fun.m which is used to see if the vehicle is slipping at the beginning of each time step.

stickslip_sim_ode.m allows integration of the discontinuous dynamics, but the terminal switching condition is determined using Matlab's event detection algorithm. This works well when not applying random actions to the system. Unfortunately, if random actions are applied the model can be thrown into an unsolvable event location, in which case Matlab's event detection fails and the code will never stop running (since it will just bounce back and forth across a zero-point without solving for the actual event).

% If applying more random, discontinuous events, use the CurvesModel to do so
