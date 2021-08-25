CurvesModel

Levi Manring
7/5/21

This folder contains curves_no_slip_fcns and curves_slip_fcns folders. Each of these folders contains the ode dynamics file for the slip/no_slip dynamics scenario and the event function for transitioning between slip/no_slip dynamics scenarios.

This directory contains a curvefit model instead of an exact solution model. This is a much faster method for calculating the dynamic response of the vehicle.

In addition, this folder includes curves_startslipcheck_fun.m which is used to see if the vehicle is slipping at the beginning of each time step.

curves_stickslip_ode.m allows integration of the discontinuous dynamics, but the terminal switching condition is determined using Matlab's event detection algorithm. This works well when not applying random actions to the system. Unfortunately, if random actions are applied the model can be thrown into an unsolvable event location, in which case Matlab's event detection fails and the code will never stop running (since it will just bounce back and forth across a zero-point without solving for the actual event).

curves_stickslip_nr.m allows integration of the discontinuous dynamics, but the terminal switching condition is determined by a user-coded Newton-Raphson routine. This allows a limit on the number of iterations used for finding a zero-point. This is the preferred method for accuracy/flexibility.

GetCurves.m is used to get the curves structure of interpolating functions needed to evaluate the curvefit model. The results of this m-file are stored in CurveFitModel.mat .