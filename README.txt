Repository

Levi Manring
8/25/21

This repository contains the modeling, simulation tools, reinforcement learning (RL) framework and results for multiple control simulations of a vehicle stuck in a ditch. The goal of this work is to provide some insight into vehicle automation and training a vehicle to intelligently escape from a ditch without human intervention.

BestResults -> contains the data files and plotting tools to quickly reproduce our results

RL_DDPG_AWD -> contains the RL framework using DDPG and assuming an AWD dynamic model with front and rear wheel drive and wheel slippage

RL_DDPG_RWD -> contains the RL framework using DDPG and assuming a RWD dynamic model with rear wheel drive and wheel slippage

RL_DDPG_RWD_NOSLIP -> contains the RL framework using DDPG and assuming a RWD dynamic model with rear wheel drive and assuming no wheel slippage

RL_PILCO_RWD -> contains the RL framework using PILCO and assuming a RWD dynamic model with rear wheel drive and assuming no wheel slippage

VehicleDynamicsModelAWD -> contains the AWD dynamics model

VehicleDynamicsModelRWD -> contains the RWD dynamics model

Each one of these files has its own README files for further explanation.