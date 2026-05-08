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

## Ribbon Modeling

Ribbons are centeral to road models and can be studied in terms of surfaces which in turn can be studied in terms of spines (curves) which may represent the centerline of the ribbon. A curve can be represented as : <br>
<br>
*C* = {**x**(s) = [x(s), y(s), z(s)]<sup>T</sup>} ∈ **R**<sup>3</sup>: s ∈ [s<sub>0</sub>, s<sub>f</sub>], with s being the arc length of *C*<br>
<br>
A moving coordinate frame is defined along the curve *C*, consisting of three orthonormal vectors that form a right-handed coordinate system. The tangent vector $\mathbf{t}(s) = \dot{\mathbf{x}}(s)$ points along the direction of travel, the principal normal $\mathbf{p}(s) = \dot{\mathbf{t}}/|\dot{\mathbf{t}}|$ points toward the center of curvature, and the binormal vector $\mathbf{b}(s) = \mathbf{t} \times \mathbf{p}$ is perpendicular to both. Together, these three vectors define the local Frenet frame at each point along the track centerline.

A ribbon *R* is constructed along the curve *C* by introducing a unit camber vector $\mathbf{n}(s)$ that lies in the plane of the ribbon and is perpendicular to the tangent vector $\mathbf{t}(s)$. This defines a surface with both width and twist, represented as: <br>
<br>
$$R = \{\mathbf{r}(s,n) = \mathbf{x}(s) + n.\mathbf{n}(s)\} \in \mathbb{R}^3, \quad s \in [s_0, s_f], \quad n \in [n_l(s), n_r(s)]$$ <br>
<br>
where $s$ is the arc length along *C* and $n$ is the lateral offset within the ribbon plane. The limits $n_l(s)$ and $n_r(s)$ define the left and right track boundaries respectively at each point along the centerline. \mathbf{n} = 0 represents the curve *C*. The unit camber vector \mathbf{n}(s) can be described in terms of vector \mathbf{p}(s) and \mathbf{b}(s) and a twist angle γ(s) as: <br>
<br>
$$\mathbf{n}(s) = \mathbf{p}(s)\cos(\boldsymbol{\gamma}(s)) - \mathbf{b}(s)\sin(\boldsymbol{\gamma}(s))$$ <br>
<br>
A unit normal vector $\mathbf{m}(s) = \mathbf{t}(s) \times \mathbf{n}(s)$ is defined perpendicular to the ribbon surface, completing the right-handed coordinate system for the ribbon frame. The curvatures of the ribbon can be derived as: <br>
<br>
1. Relative Torsion: $\omega_x = \tau + \dot{\gamma}$
2. Normal Curvature: $\omega_y = \kappa \sin(\gamma)$
3. Geodesic Curvature: $\omega_z = \kappa \cos(\gamma)$
<br>
where κ and τ are the curvature and torsion of the curve C respectively.
For the complete mathematical derivation, refer to [1].

## Optimizer

The track optimization is formulated as a Nonlinear ProBLEM (NLP) and solved using **IPOPT** (Interior Point OPTimizer) through the CasADi library. The continuous dynamics are discretized via **Direct Multiple Shooting** with 4th Order Runge-Kutta (RK4) integration across N uniformly spaced intervals along the track arc length, where N is user-defined. <br>
**A higher N produces a smoother, more accurate track model at the cost of increased computation time.**

### Dynamics Definition

The state vector (X) in terms of centerline arc length (s) for the track model can be defined as: <br>
<br>
X(s) = [x(s), y(s), z(s), $\theta(s)$, $\mu(s)$, $\phi(s)$, $\dot{\theta}(s)$, $\dot{\mu}(s)$, $\dot{\phi}(s)$, n<sub>l</sub>(s), n<sub>r</sub>(s)]<sup>T</sup> <br>
<br>

where $\theta$, $\mu$ and $\phi$ represent the euler angles describing elementary rotations about z-axis, y-axis and x-axis respectively. n<sub>l</sub>(s) and n<sub>r</sub>(s) represent the left- and right-hand track widths along s. <br>
The ribbon dynamics is given by: <br>
<br>
$\ddot{X}$ = f(s,X,U) = $[\cos(\theta).\cos(\mu), \sin(\theta).\cos(\mu), -\sin(\mu), \dot{\theta}, \dot{\mu}, \dot{\phi}, \ddot{\theta}, \ddot{\mu}, \ddot{\phi}, \dot{n_l}, \dot{n_r}]^T$ <br>
<br>
where U = $[\ddot{\theta}, \ddot{\mu}, \ddot{\phi}, \dot{n_l}, \dot{n_r}]^T$ is the **control vector**. **The track dynamics $\ddot{X}$ = f(s,X,U) serves as an equality constraint for optimizer**.

### Cost Function Definition

The cost function **J(X,U)** is defined as: J(s,X,U) = $$\int_{s_0}^{s_f} l(X(s),U(s).ds + l_f(X(s_f))$$. As the track model is continuous in s (centerline arc length), a Hamiltonian for the NLP can be defined as: <br>
<br>
$H(s,X,U) = J(s,X,U) + \lambda^T.f(s,X,U)$<br>
<br>
with $\lambda$ are the lagrange multiplers. According to the Pontryagin's Minimum Principle, the optimal solution (X*, U*) for this track model NLP satisfy the following conditions: <br>
<br>
1. $\dot{X*} = \frac{\partial \mathcal{H}}{\partial \lambda}$
2. $\dot{\lambda} = -\frac{\partial \mathcal{H}}{\partial X}$
3. $\frac{\partial \mathcal{H}}{\partial U*} = 0$

It is worth noting that Pontryagin's Minimum Principle provides the first-order necessary conditions for a continuous-time optimal control problem, which correspond to the KKT necessary conditions enforced by IPOPT in the discrete domain. When controls enter the Hamiltonian linearly, the optimal solution $U*$ has a possibility of exhibiting bang-bang controls or rapid oscilations between extreme control values. To correct this, a cost function which is quadratic in control inputs is defined. This method is called **Regularization**. <br>
<br>
The cost function consists of three terms, a tracking error term $(l_e)$, curvature rate term $(l_c)$ and a width rate term $(l_w)$. They are defined as: <br>
<br>
1. $l_e = w_c \|\mathbf{x} - \bar{\mathbf{c}}\|^2 + w_l \|\mathbf{b}_l - \bar{\mathbf{b}}_l\|^2 + w_r \|\mathbf{b}_r -\bar{\mathbf{b}}_r\|^2$
2. $l_c = w_\theta \ddot{\theta}^2 + w_\mu \ddot{\mu}^2 + w_\phi \ddot{\phi}^2$
3. $l_w = w_{nl} \dot{n}_l^2 + w_{nr} \dot{n}_r^2$
4. $J(s,X,U) = l_e + l_c + l_w$

where $\bar{\mathbf{b}}_l$, $\bar{\mathbf{b}}_r$ and $\bar{\mathbf{c}}$ are the left boundary, right boundary and centerline coordinates interpolated from raw GPS data. Weights $w_c, w_l, w_r, w_\theta, w_\mu, w_\phi, w_{nl}$ and $w_{nr}$ can be tuned to minimize boundary error while ensuring the track curvature outputs remain free of high frequency oscillations.

# Installation & Dependencies <br>

# Usage <br>

# Outputs <br>
