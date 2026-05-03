# Lap Time Simulation and Sensitivity Analysis

As a motorsports enthusiast and someone who is interested in vehicle modeling and controls, this is a project I am working on to try and learn. The framework is built upon established literature and references as mentioned below - 
References - 

**[1]** Giacomo Perantoni and David J.N. Limebeer, "Optimal Control of a Formula One Car on Three-Dimensional Track Part 1: Track Modelling and Identification".
**[2]** Giacomo Perantoni and David J.N. Limebeer, "Optimal Control of a Formula One Car on Three-Dimensional Track Part 2: Optimal Control", Journal of Dynamic Systems, Measurement and Control, 2015.
**[3]** William F. Milliken & Douglas L. Milliken, Racecar Vehicle Dynamics, PA: SAE International, 1995.
**[4]** Jorge Segers, Analysis Techniques for Racecar Data Aquisition, 2nd ed. Warrendale, PA: SAE International, 2014.

The lap time simulation framework includes a dynamic simulation of a vehicle around a given racetrack and includes four main components:

1. <u>**Track Model**</u> - A discrete ribbon model of a closed-loop track in Frenet coordinates, which provides the simulated vehicle with a path to follow. The model processes GPS coordinates (Latitude, Longitude, Height) of both right and left track boundaries to define the 3-dimensional ribbon.

2. <u>**Tire Model**</u> - A mathematical model of tire forces and moments that serves as the critical link between vehicle performance and the track surface.

3. <u>**Vehicle Model**</u> - A mathematical description of the vehicle dynamics. Inputs include vehicle data and a vehicle model of varying complexity.

4. <u>**Lap Time Optimizer**</u> - Combines the track, tire and vehicle model together along with an optimization strategy which could be online (in real time) which could be useful for autonomous driving considering a part of the track at once or offline which optimizes performance across entire track at once, useful for vehicle setup optimization. Strategies could include convex/Nonlinear Programming Numerical Optimization using gradient based solvers and/or interior point methods. Some of the readily available solvers include IPOPT (Interior Point OPTimization), SNOPT (Sparse Nonlinear OPTimizer) etc.

<img width="877" height="303" alt="image" src="https://github.com/user-attachments/assets/8b473e1e-1d27-4476-8f5c-3ce3dee69428" />
**Source: [4] Jorge Segers, Analysis Techniques for Racecar Data Acquisition, 2nd ed., 2014.**

