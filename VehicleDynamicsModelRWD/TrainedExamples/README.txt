VehicleDynamicsModelRWD

Levi Manring
7/5/21

This folder contains the packages necessary to simulate a vehicle moving on a user-defined surface profile in planar motion. This is the RWD (Rear-Wheel-Drive) version, thus torque is only applied to the rear wheel of the planar vehicle model. Thus only one wheel is assumed to be able to slip. So, there are two dynamic regions modeled in this package: rear wheel slipping and rear wheel not slipping. In addition, we have developed methods for simulating such a discontinuous system transitioning instantaneously between two different sets of dynamics. 

This folder contains two different dynamic model model methods for simulating this system. One is a curvefit model, located in CurvesModel. This system uses curvefit interpolations of some of the dynamic parameters needed to simulate the system. This model is much faster computationally.
The other is an exact solution model, located in ODEeventModel. This model calculates all of the necessary dynamic parameters at each time step. This is significantly slow, and the curvefit model is sufficiently accurate to be useful for a typical user's purposes.

OtherFilesNeeded includes a video plotting function, models for friction coefficient, calculation of relative velocity, etc.

TrainedExamples contains two .mat files that have raw data that can be used to initially check the functionality of the package.

VehicleRWDTheoreticalModals_Maple contains the equations used to derive the dynamics for both the slip region and the no-slip region. If modifying certain aspects of the dynamics, it is recommended to inspect these files first.


TO GET STARTED: Run TestingRWDEnvironments.m . An inspection of each of the scenarios presented in that m-file will help get a grip on what is contained in this package.



