VehicleAWDTheoreticalModels_Maple

Levi Manring
7/5/21

This folder four necessary folders needed to derive the dynamic model for a vehicle moving on a user-defined surface profile in planar motion. This is the AWD (All-Wheel-Drive) version, thus torque is applied to the rear wheel and the front wheel of the planar vehicle model. Thus both wheels are assumed to be able to slip. So, there are four dynamic regions modeled in Maple: 
	1. No wheels slipping (nonslip)
	2. All wheels slipping (all slip)
	3. Rear wheel slipping, front wheel sticking (rearslipfrontstick)
	4. Rear wheel sticking, front wheel slipping (rearstickfrontslip)

'... EOM Solution.mw' is used to combine the 9 equations derived using Newton's Method for the system to solve for all the free variables for the case when the given scenario (1-4). '... Model Parameters.mw' uses the results of the EOM solution to get the necessary dynamic model parameters and the normal and friction force parameters.

The folder 'AWDFBDs' contains eps files that illustrate the vehicle case.



