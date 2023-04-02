# 3DOF Hand model

This repo presents a dynamic model for 3DOF and 2DOF handsystem. The 2DOF system is essentially, two links acting as forearm and hand. The 3DOF system is the same 2DOF system, but a vertical mass, spring and damper is attached to it. Note that in this model elbow and hand can only move in a planar movement (Flexion-Extension). This model is used in my Masters degree which is in it, we have used a single SMA spring instead of a spring-damper system for absorbing the vibrations. The constitutive model for the spring is in the repository below so we won't explain how the SMA spring works here. Please refere to the following repository for more information about the SMA spring model.
https://github.com/amirabbasja/SMA-Constitutive-Model-Brinson

## How to use the model

The codes are commented in detail so you can use and change the codes as you wish. But there are some notes that you should remember:
1)The files "dynamicEqn_2DOF_DEAD_MASS" and "dynamicEqn_3DOF_AbsorberOnElbow" contain the dynamic equations of the 2DOF and 3DOF systems respectively. The dynamic equations are written in the form of a state space model. The state space model is a very useful model for simulating the dynamic systems. The state space model is a set of differential equations that can be solved using numerical methods. The mechanism of acquiring the stiffness of the joints are commented out in these files and we have only used its output to avoid redundant calculations. because we have movement only in the FE direction. if you chnage the system equations so it supports RUD or PS as well, you need to uncomment the stiffness calculation for better accuracy. Also the method of exciting the elbow and hand are implemented in this file.
2)At first, i implemented RK4 to solve the dynamic equations. It worked fine for simple systems but as i noticed after alot of trial and errors, Adding SMA's to the dynamic model, stiffens the differential equations so we cant use RK4. To solve these equations, I have used matlab's ode15s solver but some changes needed been done to it (I hope they dont mind it! If there are any problems please contact me for changes or removals). Because we used SMA in the model, To get its excerted force for the current iteration, we needed to have it's state at some iterations prior to the current iteration. ode15s offered no functionality in this regard; Also stopping the simulation when we have reached the steady state, mostly is implemented with an event function but again, This function only lets us access the parameters in the current iteration. The modified solver for ode15s is located in the "solver" folder for two matlab versions.
3)Two files "stopSim" and "StopSolverEvent" stop the simulation before reaching the END time (When we have reached a steady state). This helps us avoid redundant calculations.
4)Two files "SpringClass" and "SysParams" are used to define the spring and system parameters. You can change the parameters in these files.
5)The file "SMA_Spring_Excite" is in charge of exciting the SMA spring. A spring is passed into it as an object with a certain displacement. And it returns updated spring (Martensite ratio in cross section may change) and the force it exerts on the system.
6)The files "HandSystem_2DOF_DEAD_MASS_CODE_SINGLE" and "HandSystem_3DOF_AbsorberOnElbow_CODE_SINGLE" run the simulation for once. But the files "HandSystem_2DOF_DEAD_MASS_CODE_MULTIPLE" and "HandSystem_3DOF_AbsorberOnElbow_CODE_MULTIPLE" run the simulation for a range if frequencies. This will help us get the FRF of the system. you can change the frequency range in the "sysParams" file. These two file have a "PARALLEL" version as well which will get the system's response with help of parallel processing.
7)It's important to note that i have removed some parameters from the sysParams which have been acquired with a heavy optimization process. So to run the codes you need to update the "sysParams" file with your own parameters.

I really hope The readers can use this little project in what they are doing. Thanks for your attention!