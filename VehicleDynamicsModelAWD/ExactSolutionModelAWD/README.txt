ExactSolutionModelAWD

Levi Manring
7/5/21

This directory contains an exact solution model instead of a curvefit model. 

exactsoln_dynamicsselectorfunction.m grabs the appropriate set of dynamic parameters depending on which of the 4 slip condition possibilities. exactsoln_selector_ode then integrates those dynamics.

exactsoln_dynamicsAWD_nr.m steps through the logic process of integrating the dynamics, checking to see if a terminal condition has been reached, and then switching dynamic parameters.

ExactSolnTerminalConditions contains the velocity and friction terminal conditions used to determine when the vehicle transitions between the four dynamic scenarios.

ExactSolnParameters includes the m-files needed to calculate the dynamic parameters for each of the four dynamic scenarios. This is essentially a replacement for the curvefit model used by the CurvesModelAWD.
