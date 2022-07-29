# Parameters for the global_seir model
import math

r_d = (1 / 3.5)
r_l = (1 / 5)
r_R = (1 / 15)
p_V_H = 0.0296
p_H = 0.296
T_D = 45                       #30
w = 0.9
lmbda = 0.6
N = 1000000000                #50000
L = 30                        #20
p = 0.01
k = math.log((1-p)/p)/T_D     #natural log, ln()
rho = .95
p_e = 0.6
p_r = 0.6
p_D = 0.1