# Parameters for the global_seir model
import math

r_d = (1 / 3.5)
r_l = (1 / 5)
r_R = (1 / 15)
p_V_H = 0.04473
p_H = 0.08914
beta = 0.5
T_D = 45
w = 0.9
lmbda = .6
N = 50000
L = 30
p = 0.01
k = math.log((1-p)/p)/T_D #natural log, ln()
rho = .95
delta = 0.1
p_e = 0.6
p_r = 0.6
