VehicleRWDTheoreticalModels_Maple

Levi Manring
7/5/21

This folder four necessary Maple files needed to derive the dynamic model for the vehicle moving on an user-defined surface profile. These files derive two dynamic models: one when the rear wheel (A) is slipping and the other when the rear wheel is not slipping.

'CarRWD Noslip EOM Solution.mw' is used to combine the 9 equations derived using Newton's Method for the system to solve for all the free variables for the case when the vehicle is not slipping. 'CarRWD Noslip Model Parameters.mw' uses the results of the EOM solution to get the necessary dynamic model parameters and the normal and friction force parameters.

'CarRWD Slip EOM Solution.mw' is used to combine the 9 equations derived using Newton's Method for the system to solve for all the free variables for the case when the vehicle is slipping. 'CarRWD Slip Model Parameters.mw' uses the results of the EOM solution to get the necessary dynamic model parameters and the normal and friction force parameters.

The folder 'FBDs' contains eps and png files that illustrate the vehicle case.



