clear all
close all 
clc

%% Define the Induction Motor for constant speed
p = 1;                      %[-] pole pair number
w_el = 50  * 2 * pi;        %[1/s] electrical frequency
w_me = pi*2850/30;          %[1/s] rotor frequency

L_m = 464.704e-3;                      %[Vs/A] mutual inductance
L_sigmar = 13.777e-3;                  %[Vs/A] rotor leakage inductance
L_sigmas = 9.868e-3;                 %[Vs/A] stator leakage inductance
R_r = 2.6;                      %[V/A] rotor winding resistance
R_s = 2.7;                      %[V/A] stator winding inductance 

L_r = L_m + L_sigmar;       %[Vs/A] rotor self-inductance
L_s = L_m + L_sigmas;       %[Vs/A] stator self-inductance

tau_r = L_r / R_r;          %[1/s] rotor time constant
sigma = (L_r * L_s - L_m^2) / (L_r * L_s); %[-] total leakage factor
tau_sigma = (sigma * L_s) / (R_s + R_r * L_m^2 / L_r^2); %[-] leakage time constant

A = [-1/tau_sigma            0     (R_r * L_m)/(sigma * L_r^2 * L_s) (p * w_me * L_m)/(sigma * L_r * L_s);
                0 -1/tau_sigma -(p * w_me * L_m)/(sigma * L_r * L_s)    (R_r * L_m)/(sigma * L_r^2 * L_s);
        L_m/tau_r            0                              -1/tau_r                              -p*w_me;
                0    L_m/tau_r                                p*w_me                             -1/tau_r;];
   
B = [1/(sigma * L_s)               0;
                   0 1/(sigma * L_s);
                   0               0;
                   0               0;];
 
C = [1 0 0 0;
     0 1 0 0];
 
x0 = zeros(4,1); %[A A Vs Vs] initial state

%% Discretize the system 
T_s = 100e-6; %[s] sampling time

A_d = expm(A*T_s);
B_d = inv(A)*(A_d-eye(4))*B; 
C_d = C;

%% Design the steady-state Kalman filter
M = 0.0001* [1 0 0 0;
    0 1 0 0
    0 0 1 0
    0 0 0 1];   % system noise covariance matrix

N = 10* [1 0;
    0 1];         % measurement noise covariance matrix

[P,~,~] = dare(A_d',C_d',M,N); % Estimation Covariance matrix (We use dual property because dare command is used for controller and not for estimators).
K = P*C_d'*inv(C_d*P*C_d'+N); % Kalman Gain

%% Design of the state-tracking MPC with preview 
kf = 2;          % prediction horizon

% Definition of the state constraints (approximation of a circular
% current constraint as regular polytope with 24 vertices)
W_x = [-0.1459   -0.0604         0         0;
       -0.1253   -0.0961         0         0;
       -0.0961   -0.1253         0         0;
       -0.0604   -0.1459         0         0;
       -0.0206   -0.1566         0         0;
        0.0206   -0.1566         0         0;
        0.0604   -0.1459         0         0;
        0.0961   -0.1253         0         0;
        0.1253   -0.0961         0         0;
        0.1459   -0.0604         0         0;
        0.1566   -0.0206         0         0;
        0.1566    0.0206         0         0;
        0.1459    0.0604         0         0;
        0.1253    0.0961         0         0;
        0.0961    0.1253         0         0;
        0.0604    0.1459         0         0;
        0.0206    0.1566         0         0;
       -0.0206    0.1566         0         0;
       -0.0604    0.1459         0         0;
       -0.0961    0.1253         0         0;
       -0.1253    0.0961         0         0;
       -0.1459    0.0604         0         0;
       -0.1566    0.0206         0         0;
       -0.1566   -0.0206         0         0];
   
omega_x = ones(length(W_x),1);

% Definition of the input constraints (hexagon shape due to the three-phase 
% two-level inverter)
W_u = [-0.002307684114745  -0.001332342044853;
       -0.000000000000000  -0.002664684089706;
        0.002307684114745  -0.001332342044853;
        0.002307684114745   0.001332342044853;
        0.000000000000000   0.002664684089706;
       -0.002307684114745   0.001332342044853];
   
omega_u = ones(length(W_u),1);
       
% State constraints of the terminal state
W_xf = W_x;
omega_xf = omega_x;
   
% Weighting matrices
Q = [100 0 0 0;
    0 150 0 0
    0 0 100 0
    0 0 0 100];
S = Q;
R = [0.001 0;0 0.05];
       
% Compute matrices for the condense system representation (hint: exercise, 
% do not change the condensed matrix names)
[Qkf, Rkf] = Compute_Qkf_Rkf(kf,Q,S,R);
[W_cal_x, Omega_x, W_cal_u, Omega_u] = Compute_constraints(kf,W_x,W_xf,W_u,omega_x,omega_xf,omega_u);
[Akf, Bkf] = Compute_Akf_Bkf(kf,A_d,B_d);
C_new = eye (4); % Kalman filter can be interpreted as a sensor that provides the information of all states, therefore, the C-matrix to calculate Ckf is adapted.
[D_cal, E_cal, Ckf] = Compute_D_cal_E_cal_Ckf(kf, B_d, C_new);


% From here on, the reference quantities for the states are calculated. 
% Since this is an underactuated system, the reference variables cannot be 
% selected independently of each other. For this reason, these reference 
% variables must be calculated by a model-based coupling. For the solution 
% of the task, however, this approach does not have to be understood, 
% since expert electrical drives engineering knowledge is required here.

%% Calculation of state reference phasors
% Here a state-space representation in a rotor-flux orientated coordinate 
% system is defined. 
% defined. 

psi_dq_ss=[0.6; 0]; %[Vs] rotorflux 

% state-space representation in rotoflux orientated dq-coordinate system
A_ss =[-1/tau_sigma          w_el      (R_r * L_m)/(sigma * L_r^2 * L_s) (p * w_me * L_m)/(sigma * L_r * L_s);
              -w_el  -1/tau_sigma  -(p * w_me * L_m)/(sigma * L_r * L_s)    (R_r * L_m)/(sigma * L_r^2 * L_s);
          L_m/tau_r             0                               -1/tau_r                         -p*w_me+w_el;
                  0     L_m/tau_r                            p*w_me-w_el                             -1/tau_r;];
              
% stationary stator current for a given rotorflux
i_dq_ss=-inv(A_ss(3:4,1:2))*A_ss([3 4],[3 4])*psi_dq_ss;

% stator current in the alpha-beta-coordinate system at t=0
i_alphabeta_0=i_dq_ss;

% rotor flux in the alpha-beta-coordinate system at t=0
psi_alphabeta_0=psi_dq_ss;



 