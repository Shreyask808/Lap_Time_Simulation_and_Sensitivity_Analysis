# Track Model Optimization

# Overview <br>
The track model provides the simulated vehicle with a path to follow, defined by the track curvature as a function of arc length distance. This module is implemented in Python using CasADi and is based on the methodology described in: <br>

[1] Giacomo Perantoni and David J.N. Limebeer, "Optimal Control of a Formula One Car on Three-Dimensional Track — Part 1: Track Modelling and Identification" <br>

This module takes raw GPS coordinates (latitude, longitude, altitude) of points along the left and right track boundaries (recorded in the direction of travel) and passes them through an optimizer that simultaneously converts them into a Frenet-coordinate ribbon model and smooths the track geometry, producing continuous curves that represent the track boundaries with minimal boundary errors. The resulting model provides well-defined track curvatures (Normal, torsion, and geodesic), heading angles, track widths as a function of centerline arc length, suitable for use with vehicle and tire models in simulations.

<p align="center">
  <img src="Iamges/Barcelona.png" alt="Circuit de Barcelona-Catalunya 3D Track Model" width="800"/>
  <br>
  <em>3D ribbon model of Circuit de Barcelona-Catalunya — raw GPS boundary data (dashed)</em>
</p>

**Key Capabilities of the Module**
1. Converts raw GPS track boundary data into srutcutred Frenet coordinate track model
2. Compatible with both **clockwise and anti-clockwise** track layouts
3. Accepts left and right boundary data with **unequal number of points**, no resampling or preprocessing required
4. Simultaneously optimizes and smoothens track geometry, producing continuous track curvatures
5. Outputs include relative track curvatures (torsion, normal, geodesic), and track width as a function of arc length
6. Compatible with vehicle dynamics and tire models for simulations
