# MPC-Kalman-Optimization-for-Induction-Motor-Control

Project Overview
This project designs a Model Predictive Control (MPC) system for an induction motor, leveraging a Kalman Filter to estimate unmeasurable rotor fluxes and track dynamic sinusoidal references. The solution integrates:

-Discrete-time state-space modeling (10 kHz sampling).

-Steady-state Kalman Filter for robust state estimation under process and measurement noise.

-Constraint-aware MPC with quadratic programming (QP) to optimize control inputs while adhering to polytopic state and hexagonal input constraints.

Key Components
1. System Discretization & Kalman Filter Design
-Continuous-to-Discrete Conversion: System matrices Ad and Bd derived using matrix exponential methods.

-Noise Modeling:

-Kalman Gain Calculation: Solved via Discrete Algebraic Riccati Equation (DARE), yielding steady-state gain:

Result: Near-perfect state estimation (Figs. 1–2).

2. MPC with Dynamic Reference Tracking
-Cost Function: Minimizes state deviation, input changes, and terminal error using weights Q,R,S:


-QP Formulation

3. Performance & Constraints
-Input Constraints: Hexagonal voltage limits (peak: 272 V, spike: 235 V) (Fig. 3).

-State Tracking:

-Currents (x₁, x₂): Fast, oscillation-free tracking (Fig. 4).

-Rotor Fluxes (x₃, x₄): Slower convergence (~0.5 s) due to unmeasured states (Fig. 5).

4. Repository Structure
-matlab code

-figures:

-Report detailing theory, tuning, and results.

5. Dependencies
-MATLAB (Control System Toolbox, Optimization Toolbox).

-QP Solver (e.g., quadprog).

6. References
Lecture notes on Kalman Filters and QP optimization.
