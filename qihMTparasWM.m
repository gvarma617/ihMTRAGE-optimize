% Matlab R2018a code with global variables relating to white matter tissue
% based on the paper by Varma et al. Low duty-cycle pulsed irradiation
% reduces magnetization transfer and increases the inhomogeneous
% magnetization transfer effect. J Magn Reson. 2018 Nov;296:60-71.
% doi: 10.1016/j.jmr.2018.08.004.

clear global T_2B
clear global R_A
clear global M_0A
clear global R
clear global M_0B1
clear global M_0B2
clear global T_2A
clear global R_B
clear global T_D

global T_2B
global R_A
global M_0A
global R
global M_0B1
global M_0B2
global T_2A
global R_B
global T_D

T_2B = 9.0e-6;
R_A  = 0.92;
M_0A = 1.0;
R = 60;
f = 1.00;
M_0B1 = f*0.100;
M_0B2 = (1-f)*0.100;
T_2A = 69e-3;
R_B = 1.0;
T_D = 6.2e-3;