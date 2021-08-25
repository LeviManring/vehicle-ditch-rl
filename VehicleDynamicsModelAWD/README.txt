VehicleDynamicsModelAWD

Levi Manring
7/5/21

This folder contains the packages necessary to simulate a vehicle moving on a user-defined surface profile in planar motion. This is the AWD (All-Wheel-Drive) version, thus torque is applied to the rear wheel and the front wheel of the planar vehicle model. Thus both wheels are assumed to be able to slip. So, there are four dynamic regions modeled in this package: 
	1. No wheels slipping (nonslip)
	2. All wheels slipping (all slip)
	3. Rear wheel slipping, front wheel sticking (rearslipfrontstick)
	4. Rear wheel sticking, front wheel slipping (rearstickfrontslip)
In addition, we have developed methods for simulating such a discontinuous system transitioning instantaneously between these four different sets of dynamics. 

This folder contains two different dynamic model model methods for simulating this system. One is a curvefit model, located in CurvesModelAWD. This system uses curvefit interpolations of some of the dynamic parameters needed to simulate the system. This model is more memory intensive, and perhaps slightly faster computationally in operation.
The other is an exact solution model, located in ExactSolutionModelAWD. This model calculates all of the necessary dynamic parameters at each time step. This should be slower than the curves model, but due to the large memory requirements of the curvefit interpolated models, it actually ends up being somewhat similar in computational intensity. We chose to use the exact solution model for reinforcement learning training.

OtherFilesNeeded includes a video plotting function, models for friction coefficient, calculation of relative velocity, etc.

AWDExample contains data from results from these two models that can be used to initially check the functionality of the package.

VehicleAWDTheoreticalModels_Maple contains the equations used to derive the dynamics for all four regions of dynamic behavior. If modifying certain aspects of the dynamics, it is recommended to inspect these files first.


TO GET STARTED: Run TestingAWDEnvironments.m . An inspection of each of the scenarios presented in that m-file will help get a grip on what is contained in this package. If running TestingAWDEnvironments for the first time, either simulate just the ExactSoln method on its own first, or load the recorded data from AWDExample to cut down the amount of time needed to examine the functionality of these packages.



