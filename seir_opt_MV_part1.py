# SEIR-OPT model: Overall structure and Optimize
# Does not include: reading inputs, Simulate, or outputs
# Some parts are not fully written or not correct syntax. Where I don't know the syntax, I generally use comments
# Includes calls to Gurobi


# Calling functions: All variables are treated as global, so there are no arguments. 
# In some cases, when they compute new values they are give nnew names: V is always the input to simulate(), V_new is the output.
# In some cases, the same varaible is updated: tn and m get updated in optimize_inner().

# Mike Veatch 
# Based on Abraham Holleran's seir_opt.py

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd
import time

# Read input files

from areas import * #areas, donor, m, n, N(a), rho_V(a), rho_I(a), delta_r(a), gamma(a), rho(a) 
from scenario import * #rI, r0, rR, pH, pHV, pD, a0, delta_a, pe, pr, L, TD, p, T, B_0, b(t)
from parameters import * #pk, lambda_0, phi, eps_0, delta_I, delta, beta, iter_lmt, iter_lmt_search # Changed name from lambda to lambda_0

# Computed constants

donor = "area1"
rd = {a: r0 + delta_r[a]) for a in areas}
B = {t: B0 * b(t) for T in [0,T-1]}  # If no b(t) in input file, use b(t) = 1
k = log((1-p)/p)/TD

# num_areas = len(areas)

# Initial States

E[a,0]  = {a: (1 - rho_V[a]) / (pr*rho_V[a] + 1 - rho_V[a]) /rI * rho_I[a] * N[a] for a in areas}
EV[a,0] = {a: (pr*rho_V[a]) / (pr*rho_V[a] + 1 - rho_V[a]) /rI * rho_I[a] * N[a] for a in areas}
I[a,0]  = {a: (1 - rho_V[a]) / (pr*rho_V[a] + 1 - rho_V[a]) /(r0 + delta_r[a]) * rho_I[a] * N[a] for a in areas}
IV[a,0] = {a: (pr*rho_V[a]) / (pr*rho_V[a] + 1 - rho_V[a]) /(r0 + delta_r[a]) * rho_I[a] * N[a] for a in areas}
SV[a,0] = {a: rho_V[a]*N[a] - EV[a,0] - IV[a,0] for a in areas}
S[a,0]  = {a: N[a] - E[a,0] - I[a,0] - EV[a,0] - IV[a,0] -SV[a,0] for a in areas}
H[a,0]  = {a: 0 for a in areas}
D[a,0]  = {a: 0 for a in areas}
R[a,0]  = {a: 0 for a in areas}

# Initialize first Simulate

S_dot = S.sum(a,0) - S[donor,0]     # Sum of nondonor S at time 0
V[a,t]  = {a,t: (1 - pk) * S[a,0]/Sdot * B[t] for a in areas for t in [0,T-1]}
V[donor,t]  = {t: pk * B[t] for t in [0,T-1]}

#Create a simulate object

def simulate():
    # key input: V[a,t]
    # key outputs:  V_plus[a,t] (actual vacc), I_hat[a,t], IV_hat[a,t], tn, m

#Create the inner loop of Optimize

def optimize_inner():
    # key inputs: I_hat[a,t], IV_hat[a,t], tn, m, lambda[i]
    # key outputs: V_new[a,t] (optimal vacc), tn, m, z[i] (donor deaths)

    j = 0
    eps = eps_0
    zLP = {j: 0 for j in [0,iter_lmt]} # Define vector zLP and set to 0
    
    while abs(zLP[j] - z[max(0,j-1)]) > delta_I or i < 2:
        j = j + 1
        solve_LP()
        V[a,t]  = {a,t: V_new for a in areas for t in [0,T-1]}  # Update V using optimal V_new
        simulate(V)                                             
        eps = beta*eps
    z[i] = D[donor,T]   #Donor deaths from current LP solution

def solve_LP()          
    # key inputs: I_hat, IV_hat, tn, m, lambda[i], eps
    # key outputs:  V_new[a,t] (actual vacc), zLP[j] (optimal objective value)
    # TBD

## Simulate

If optflag == 0:      # or FALSE
    simulate(V)                                                     
    # Print output
    # Write output file
    # Terminate

## Optimize

simulate(V)                                             

# Initialize Phase 1 of Optimize 
# Notation: f = z[i] = donor deaths, x = lambda[i]

phase = 1
lambda = {i: 0 for i in [0,iter_lmt_search]} # Define vector
lamba[1] = lambda_0
lambda[2] = lambda_0*phi
z = {i: 0 for i in [0,iter_lmt_search]} # Define vector and set z[0] = 0

# Iterations 1 and 2
i = 1
optimize_inner(I_hat, IV_hat, tn, m,  lambda[i])    #These are the inputs, but tn and m are updated in the function 

fL = z[i]           #f at left end pt
i = 2
optimize_inner(I_hat, IV_hat, tn, m,  lambda[i]) 
z[i] = D[donor,T]   
fR = z[i]           #f at right end pt
If fL <= fR:        #f is increasing, search to left
    mult = 1/phi	#x multiplier < 1 
    x1 = lambda[2]	#largest of 3 current x values
    x2 = lambda[1]	#middle of 3 current x values
    f1 = fR	        #f(x1)
    f2 = fL	        #f(x2)
If fL > fR: 	    #f is decreasing, search to right
    mult = phi	    #x multiplier > 1
    x1 = lambda[1]	#smallest of 3 current x values
    x2 = lambda[2]	#middle of 3 current x values
    f1 = fL         	
    f2 = fR

# Phase 1 Main loop

While (iter < iter_lmt_search and phase == 1)
    i = i + 1
    lambda[i] = mult*x2
    x3 = lambda[i]	        #3rd of 3 current x values
    optimize_inner(I_hat, IV_hat, tn, m,  lambda[i]) 
    z[i] = D[donor,T] 
    f3 = z[i]               #f(x3)
    If f3 > f2: phase = 2	#x3 is past the minimum
    x1 = x2                 #shift x’s for next Phase 1 iteration 
    x2 = x3
    f1 = f2
    f2 = f3
    
### Phase 2: Golden ratio search on interval [a, b] with check for unimin

# Initialize Phase 2 of Optimize

If x1 < x3: 
    a = x1 
    b = x3 
    fa = f1 
    fb = f3
If x1 > x3: 
    a = x3
    b = x1
    fa = f3
    fb = f1

# Two more iterations, at x and y
x = a + 0.618 * (b - a) 	#current larger x value
y = b - 0.618 * (b - a)		#current smaller x value
i = i + 1
lambda[i] = X
optimize_inner(I_hat, IV_hat, tn, m,  lambda[i]) 
z[i] = D[donor,T]
fx = z[i]                   #f(x)
i = i + 1
lambda[i] = y
optimize_inner(I_hat, IV_hat, tn, m,  lambda[i]) 
z[i] = D[donor,T]
fy = z[i]                   #f(y)

# Phase 2 main loop

while (abs(fx - fy) > delta and iter < iter_lmt_search)
    i = i + 1
    if fx > fy           		# minimum is in [a,x], so update y
        b, x, fx = (x, y, fy)
        y = b - 0.618 * (b - a)
        lambda[i] = y
        optimize_inner(I_hat, IV_hat, tn, m,  lambda[i]) 
        z[i] = D[donor,T]
        fy = z[i]
    else				        # minimum is in [y,b], so update x
        a, y, fy = (y, x, fx)
        x = a + 0.618 * (b - a)
        lambda[i] = x
        optimize_inner(I_hat, IV_hat, tn, m,  lambda[i]) 
        z[i] = D[donor,T]
        fx = z[i]
    If (fy < fx and fx > fb) print "Warning: f is not unimin" # print lambda[i-3],z[i-3],...,lambda[i],z[i]
    If (fy > fx and fa < fy) print "Warning: f is not unimin" # print lambda[i-3],z[i-3],...,lambda[i],z[i]

if (fy <= fx): 
    z_opt = fy 
    lambda_opt = y
else
    z_opt = fx 
    lambda_opt = X

# Write output file
# Print outputs, including z_opt = optimal donor deaths