# Track Model Optimization

# Overview <br>
The track model provides the simulated vehicle with a path to follow, defined by the track curvature as a function of arc length distance. This module is implemented in Python using CasADi and is based on the methodology described in: <br>

[1] Giacomo Perantoni and David J.N. Limebeer, "Optimal Control of a Formula One Car on Three-Dimensional Track — Part 1: Track Modelling and Identification" <br>

This module takes raw GPS coordinates (latitude, longitude, altitude) of points along the left and right track boundaries (recorded in the direction of travel) and passes them through an optimizer that simultaneously converts them into a Frenet-coordinate ribbon model and smooths the track geometry, producing continuous curves that represent the track boundaries with minimal boundary errors. The resulting model provides well-defined track curvatures (Normal, torsion, and geodesic), heading angles, track widths as a function of centerline arc length, suitable for use with vehicle and tire models in simulations.

<p align="center">
  <img src="Images/Barcelona.png" alt="Circuit de Barcelona-Catalunya 3D Track Model" width="800"/>
  <br>
  <em>3D ribbon model of Circuit de Barcelona-Catalunya — raw GPS boundary data (dashed) vs optimized Frenet-frame model (solid)</em>
</p>

**Key Capabilities of the Module**
1. Converts raw GPS track boundary data into srutcutred and closed Frenet coordinate track model
2. Compatible with both **clockwise and anti-clockwise** track layouts
3. Can be used to generate both, 2D and 3D track models
4. Accepts left and right boundary data with **unequal number of points**, no resampling or preprocessing required
5. Simultaneously optimizes and smoothens track geometry, producing continuous track curvatures
6. Outputs include relative track curvatures (torsion, normal, geodesic), and track width as a function of arc length
7. Compatible with vehicle dynamics and tire models for simulations

# How it Works <br>

Ribbons are centeral to road models and can be studied in terms of surfaces which in turn can be studied in terms of spines (curves) which may represent the centerline of the ribbon. A curve can be represented as : <br>
<br>
*C* = {**x**(s) = [x(s), y(s), z(s)]<sup>T</sup>} ∈ **R**<sup>3</sup>: s ∈ [s<sub>0</sub>,s<sub>f</sub>], with s being the arc length of *C*<br>
<br>
A moving coordinate frame is defined along the curve *C*, consisting of three orthonormal vectors that form a right-handed coordinate system. The tangent vector $\mathbf{t}(s) = \dot{\mathbf{x}}(s)$ points along the direction of travel, the principal normal $\mathbf{p}(s) = \dot{\mathbf{t}}/|\dot{\mathbf{t}}|$ points toward the center of curvature, and the binormal vector $\mathbf{b}(s) = \mathbf{t} \times \mathbf{p}$ is perpendicular to both. Together, these three vectors define the local Frenet frame at each point along the track centerline.


# Installation & Dependencies <br>

# Usage <br>

# Outputs <br>
