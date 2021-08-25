RL_DDPG_AWD

Levi Manring
8/25/21

This folder contains the RL framework using DDPG and assuming an AWD dynamic model with front and rear wheel drive and wheel slippage. It includes the: 

ResetFunction (which defines the initial conditions for each training episode), 

StepFunction (which defines the plant or dynamics of the environment), 

DDPGAgent (which defines the neural-network parameters for the agent), 

Rewardfun (which defines the achieved reward given states and actions),

Learn (which is the file used to train the agent).

To get STARTED: open AWDlearn.m and run.