CurvesModelAWD

Levi Manring
7/5/21

This directory contains a curvefit model instead of an exact solution model. 

GetCurvesAWD.m is used to get the curves structure of interpolating functions needed to evaluate the curvefit model. The results of this m-file are stored in CurvesAllModels.mat under the folder CurvesDataModels.

curves_dynamicsselectorfunction.m grabs the appropriate set of dynamic parameters depending on which of the 4 slip condition possibilities. curves_selector_ode then integrates those dynamics.

curves_dynamicsAWD_nr.m steps through the logic process of integrating the dynamics, checking to see if a terminal condition has been reached, and then switching dynamic parameters.

CurvesTerminalConditions contains the velocity and friction terminal conditions used to determine when the vehicle transitions between the four dynamic scenarios.