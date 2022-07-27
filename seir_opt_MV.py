# SEIR-OPT model: Overall structure and Optimize by Mike Veatch, based on Abraham Holleran's seir_opt.py

# Need to add:
    # Reading inputs: I have listed the names of the inputs, which you should recognize but might differ slightly from yours.
    # Simulate: I call a function simulate(), which is intended to be what you wrote. Printing/writing output should NOT be in simulate().
      # Computing alpha[a,t] should be in simulate(). If not, add it to formulate_LP() 
    # Writing "simulation data" file: I refer to Abraham's code. It won't work, partly b/c we changed inputs.
    # Error checking, e.g., whether the two simulation limits were reached
    # Better stopping criteria for outer loop: |fx - fy| < delta does NOT imply fx or fy is within delta of min. 
# Some parts are not correct syntax. Where I don't know the syntax, I generally use comments.
 
# Gurobi (LP) variables are different than other variables. When both are needed and had the same name, LP variables have "1"" appended: V1, etc.
# The non-LP versions appear as constants in the LP. Adopting this convention led to these changes:
#   Initial conditions W[a,0] etc. are constants, not LP variables, so constraints at t = 0 are written separately, e.g., W1[a,1] = W[a,0] + ... 
#   V_new was renamed V (that is sufficient b/c it is different than V1) 
#   I_hat, IV_hat were renamed I, IV
# After solving, the optimal values of variables are reference, e.g., V1.x[a,t] for V1[a,t] 

# Calling functions: All variables are treated as global, so there are no arguments. 
# In some cases, when they compute new values they are give new names: V is the input to simulate(), V_act is the output.
# In some cases, the same variable is updated: tn and m get updated in optimize_inner().

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd
import time

# Read input files

from areas import * #areas, donor, m, n, N[a], rho_V[a], rho_I[a], delta_r[a], gamma[a], rho[a] 
from scenario import * #rI, r0, rR, pH, pHV, pD, a0, delta_a, pe, pr, L, TD, p, T, B_0, b[t]
from parameters import * #pk, lambda_0, phi, eps_0, delta_I, delta, beta, iter_lmt, iter_lmt_search # Changed name from lambda to lambda_0

# Computed constants

donor = "area1"
rd = {a: r0 + delta_r[a]) for a in areas}
B = {t: B0 * b(t) for t in range(T)}  # If no b(t) in input file, use b(t) = 1. # range(T): 0, ..., T-1
k = log((1-p)/p)/TD

# num_areas = len(areas)    # Not needed?

# Initial States. These are constants of LP. Variables of LP are named E1, etc.
W = E = EV = I = IV = SV = S = H = D = R = {(a,t): 0 for a in areas for t in range(T+1)}
E[a,0]  = {a: (1 - rho_V[a]) / (pr*rho_V[a] + 1 - rho_V[a]) /rI * rho_I[a] * N[a] for a in areas}
EV[a,0] = {a: (pr*rho_V[a]) / (pr*rho_V[a] + 1 - rho_V[a]) /rI * rho_I[a] * N[a] for a in areas}
I[a,0]  = {a: (1 - rho_V[a]) / (pr*rho_V[a] + 1 - rho_V[a]) /(r0 + delta_r[a]) * rho_I[a] * N[a] for a in areas}
IV[a,0] = {a: (pr*rho_V[a]) / (pr*rho_V[a] + 1 - rho_V[a]) /(r0 + delta_r[a]) * rho_I[a] * N[a] for a in areas}
SV[a,0] = {a: rho_V[a]*N[a] - EV[a,0] - IV[a,0] for a in areas}
S[a,0]  = {a: N[a] - E[a,0] - I[a,0] - EV[a,0] - IV[a,0] -SV[a,0] for a in areas}
H[a,0]  = {a: 0 for a in areas}
D[a,0]  = {a: 0 for a in areas}
R[a,0]  = {a: 0 for a in areas}

W[a,0]  = {a: rho[a]*N[a] - SV[a,0] - EV[a,0] - IV[a,0] - rho[a]*E[a,0] - rho[a]*I[a,0] for a in areas} 

# Initialize first Simulate

S_dot = S.sum(a,0) - S[donor,0]     # Sum of nondonor S at time 0
V  = {(a,t): (1 - pk) * S[a,0]/Sdot * B[t] for a in areas for t in range(T)}
V[donor,t]  = {t: pk * B[t] for t in range(T)}

#Create a simulate object

def simulate():
    # key input: V[a,t]
    # key outputs:  V_act[a,t] (actual vacc, called V* in paper), I[a,t], IV[a,t] (simulated, called I_hat, IV_hat in paper), tn, m

# Create the inner loop of Optimize

def optimize_inner():
    # key inputs: I[a,t], IV[a,t]v(simulated, called I_hat, IV_hat in paper), tn, m, lambda[i]
    # key outputs: V1[a,t] (optimal vacc from LP, called V_new in paper), tn, m, z[i] (donor deaths)

    j = 0
    eps = eps_0
    zLP = {j: 0 for j in range(iter_lmt)} # Define vector zLP and set to 0
    
    while abs(zLP[j] - zLP[max(0,j-1)]) > delta_I or j < 2:
        j = j + 1
        solve_LP()
        V = {(a,t): V1.x[a,t] for a in areas for t in range(T)}  # Update V using optimal V1 from LP
        simulate()                                             
        eps = beta*eps
    z[i] = D1[donor,T].x   #Donor deaths from LP solution
                        #Could move outside of this function.

# Create object that formulates and solves LP

def formulate_LP()          
    # key inputs: (I, IV, alpha) for all areas and times, tn, m, lambda[i], eps
    # key outputs:  LP model
    
    t_int = ceil(tn)
    V_cal = {(a,t): I[a,t] + pe*IV[a,t] for a in areas for t in range(T) } # "vectors" (constants)

    # Compute alpha[a,t] here if not done in simulate()?? 
      # Use constant rate before tn or tn + L and variable rate after

    v = gp.Model("vaccine_opt")
    # v.setParam("NumericFocus",1)  # Increased protection from roundoff/instability. Shouldn't need.

    # Define LP variables, e.g., S1 for state var S. All are continuous and nonnegative by default. 

    S1 = v.addVars(areas, range(1,T+1), name="S1")  # range(1,T+1): 1, ..., T
    SV1 = v.addVars(areas, range(1,T+1), name="SV1")
    E1 = v.addVars(areas, range(1,T+1), name="E1")
    EV1 = v.addVars(areas, range(1,T+1), name="EV1")
    I1 = v.addVars(areas, range(1,T+1), name="I1")
    IV1 = v.addVars(areas, range(1,T+1), name="IV1")
    H1 = v.addVars(areas, range(1,T+1), name="H1")
    D1 = v.addVars(areas, range(1,T+1), name="D1")
    R1 = v.addVars(areas, range(1,T+1), name="R1")
    W1 = v.addVars(areas, range(1,T+1), name="W1")
    V1 = v.addVars(areas, T, name="V1")             # t = 0, ..., T-1


    if non_donor_deaths_flag == 0: #or FALSE
        v.setObjective(D1[donor, T] + lambda[i]*sum(I1[a,t] for a in areas for t in range(1, t_int + 1) ) \
            + 1/200000*D1.sum('*', T), GRB.MINIMIZE)  # include non-donor deaths w/ small weight
                                                    # Omit constant I[a,0] from objective, so t = 1, ..., t_int
    else:
            v.setObjective(D1[donor, T] + lambda[i]*sum(I1[a, t] for a in areas for t in range(1, t_int + 1) ), \
                 GRB.MINIMIZE)
                 
    v.addConstrs((V1[donor, t] <= pk*B[t] for t in range(T)), "Policy: donor limit") # constraint must be in () if name argument used
    v.addConstrs((V1.sum('*', t) <= B[t] for t in range(T)), "Vaccine budget")                  

    # Dynamics for start time t=1 to T-1
    v.addConstrs((W1[a,t+1] == W1[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*W1[a,t] - V1[a,t] \
                        for a in areas for t in range(1,T)), "Vaccine willingness")
    v.addConstrs((S1[a,t+1] == S1[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*S1[a,t] - V1[a,t] \
                        for a in areas for t in range(1,T)), "S")                    
    v.addConstrs((SV1[a,t+1] == SV1[a,t] - pr*alpha[a,t]*V_cal[a,t]/N[a]*SV1[a,t] + V1[a,t] \
                        for a in areas for t in range(1,T)), "SV")
    v.addConstrs((E1[a,t+1] == E1[a,t] + alpha[a,t]*V_cal[a,t]/N[a]*S[a,t] - rI*E1[a,t] \
                        for a in areas for t in range(1,T)), "E")  
    v.addConstrs((EV1[a,t+1] == EV1[a,t] + pr*alpha[a,t]*V_cal[a,t]/N[a]*SV1[a,t] - rI*EV1[a,t] \
                        for a in areas for t in range(1,T)), "EV")
    v.addConstrs((I1[a,t+1] == I1[a,t] + rI*E1[a,t] - rd*I1[a,t] \
                        for a in areas for t in range(1,T)), "I")
    v.addConstrs((IV1[a,t+1] == IV1[a,t] + rI*EV1[a,t] - rd*IV1[a,t] \
                        for a in areas for t in range(1,T)), "IV")
    v.addConstrs((H1[a,t+1] == H1[a,t] + rd*pH*I1[a,t] + rd*pHV*IV1[a,t] - rR*H1[a,t] \ 
                        for a in areas for t in range(1,T)), "H")
    v.addConstrs((D1[a,t+1] == D1[a,t] + rR*pD*H1[a,t] \
                        for a in areas for t in range(1,T)), "D")
    v.addConstrs((R1[a,t+1] == R1[a,t] + rR*(1 - pD)*H1[a,t] + rd*(1 - pH)*I1[a,t] + rd*(1 - pHV)*IV1[a,t] \
                        for a in areas for t in range(1,T)), "R")

    # Dynamics for start time t=0. Use constant W[a,0} etc. (but variable V1). Use same constraint names??
    v.addConstrs((W1[a,1] == W[a,0] - alpha[a,0]*V_cal[a,0]/N[a]*W[a,0] - V1[a,0] \
                        for a in areas), "Vaccine willingness t=0")
    v.addConstrs((S1[a,1] == S[a,0] - alpha[a,0]*V_cal[a,0]/N[a]*S[a,0] - V1[a,0] \
                        for a in areas), "S t=0")                   
    v.addConstrs((SV1[a,1] == SV[a,0] - pr*alpha[a,0]*V_cal[a,0]/N[a]*SV[a,0] + V1[a,0] \
                        for a in areas), "SV t=0")
    v.addConstrs((E1[a,1] == E[a,0] + alpha[a,t]*V_cal[a,0]/N[a]*S[a,0] - rI*E[a,0] \
                        for a in areas), "E t=0")  
    v.addConstrs((EV1[a,1] == EV[a,0] + pr*alpha[a,0]*V_cal[a,0]/N[a]*SV[a,0] - rI*EV[a,0] \
                        for a in areas), "EV t=0")
    v.addConstrs((I1[a,1] == I[a,0] + rI*E[a,0] - rd*I[a,0] \
                        for a in areas), "I t=0")
    v.addConstrs((IV1[a,1] == IV[a,0] + rI*EV[a,0] - rd*IV[a,0] \
                        for a in areas), "IV t=0")
    v.addConstrs((H1[a,1] == H[a,0] + rd*pH*I[a,0] + rd*pHV*IV[a,0] - rR*H[a,0] \ 
                        for a in areas), "H t=0")
    v.addConstrs((D1[a,1] == D[a,0] + rR*pD*H[a,0] \
                        for a in areas), "D t=0")
    v.addConstrs((R1[a,1] == R[a,0] + rR*(1 - pD)*H[a,0] + rd*(1 - pH)*I[a,0] + rd*(1 - pHV)*IV[a,0] \
                        for a in areas), "R t=0")

    # Regularization constraints on I, IV
    v.addConstrs((I1[a,t] - I[a,t] <= eps for a in areas for t in range(1,T+1)), "I upper bd")
    v.addConstrs((I1[a,t] - I[a,t] >= -eps for a in areas for t in range(1,T+1)), "I lower bd")
    v.addConstrs((IV1[a,t] - IV[a,t] <= eps for a in areas for t in range(1,T+1)), "IV upper bd")
    v.addConstrs((IV1[a,t] - IV[a,t] >= -eps for a in areas for t in range(1,T+1)), "IV lower bd")

 def solve_LP()          
    # key inputs: I, IV, tn, m, lambda[i], eps
    # key outputs:  V[a,t] (actual vacc), zLP[j] (optimal objective value)

    t_int = ceil(tn)
    V_cal = {(a,t): I[a,t] + pe*IV[a,t] for a in areas for t in range(T) } # "vectors" (constants)

    # Compute alpha[a,t] here if not done in simulate()?? 
      # Needs to include constant rate before tn or tn + L and variable rate after

    # Objective changes with lambda (outer loop) and tn (inner loop). 

    if non_donor_deaths_flag == 0: #or FALSE
        v.setObjective(D1[donor, T] + lambda[i]*sum(I1[a,t] for a in areas for t in range(1, t_int + 1) ) \
            + 1/200000*D1.sum('*', T), GRB.MINIMIZE)  # include non-donor deaths w/ small weight
                                                    # Omit constant I[a,0] from objective, so t = 1, ..., t_int
    else:
            v.setObjective(D1[donor, T] + lambda[i]*sum(I1[a, t] for a in areas for t in range(1, t_int + 1) ), \
                 GRB.MINIMIZE)

    # Constraints change with alpha, V_cal (inner loop)
    # Assumes that addConstrs REPLACES constraints with the same name.

    # Start time t=1 to T-1
    v.addConstrs((W1[a,t+1] == W1[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*W1[a,t] - V1[a,t] \
                        for a in areas for t in range(1,T)), "Vaccine willingness")
    v.addConstrs((S1[a,t+1] == S1[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*S1[a,t] - V1[a,t] \
                        for a in areas for t in range(1,T)), "S")                    
    v.addConstrs((SV1[a,t+1] == SV1[a,t] - pr*alpha[a,t]*V_cal[a,t]/N[a]*SV1[a,t] + V1[a,t] \
                        for a in areas for t in range(1,T)), "SV")
    v.addConstrs((E1[a,t+1] == E1[a,t] + alpha[a,t]*V_cal[a,t]/N[a]*S[a,t] - rI*E1[a,t] \
                        for a in areas for t in range(1,T)), "E")  
    v.addConstrs((EV1[a,t+1] == EV1[a,t] + pr*alpha[a,t]*V_cal[a,t]/N[a]*SV1[a,t] - rI*EV1[a,t] \
                        for a in areas for t in range(1,T)), "EV")

    # Start time t=0
    v.addConstrs((W1[a,1] == W[a,0] - alpha[a,0]*V_cal[a,0]/N[a]*W[a,0] - V1[a,0] \
                        for a in areas), "Vaccine willingness t=0")
    v.addConstrs((S1[a,1] == S[a,0] - alpha[a,0]*V_cal[a,0]/N[a]*S[a,0] - V1[a,0] \
                        for a in areas), "S t=0")                   
    v.addConstrs((SV1[a,1] == SV[a,0] - pr*alpha[a,0]*V_cal[a,0]/N[a]*SV[a,0] + V1[a,0] \
                        for a in areas), "SV t=0")
    v.addConstrs((E1[a,1] == E[a,0] + alpha[a,t]*V_cal[a,0]/N[a]*S[a,0] - rI*E[a,0] \
                        for a in areas), "E t=0")  
    v.addConstrs((EV1[a,1] == EV[a,0] + pr*alpha[a,0]*V_cal[a,0]/N[a]*SV[a,0] - rI*EV[a,0] \
                        for a in areas), "EV t=0")

    # Constraints change with I, IV, eps (inner loop)
      # Regularization constraints on I, IV
    v.addConstrs((I1[a,t] - I[a,t] <= eps for a in areas for t in range(1,T+1)), "I upper bd")
    v.addConstrs((I1[a,t] - I[a,t] >= -eps for a in areas for t in range(1,T+1)), "I lower bd")
    v.addConstrs((IV1[a,t] - IV[a,t] <= eps for a in areas for t in range(1,T+1)), "IV upper bd")
    v.addConstrs((IV1[a,t] - IV[a,t] >= -eps for a in areas for t in range(1,T+1)), "IV lower bd")

    try:
        v.optimize()
    except gp.GurobiError:
        print("Optimize failed at lambda =", lambda[i]. "LP iteration j =", j, "exploration tol eps =" eps)
    zLP[j] = v.objVal  # Store optimal objective value. Could move outside of this function. 

## Simulate

If optflag == 0:      # or FALSE
    simulate()                                                     
    # Print simulate output: donor deaths D[donor,T], deaths D.sum['*',T], day of variant tn, variant area m, 
        # vaccinations by area  V_act.sum[area,'*'] (called V* in paper), total vaccinations V_act.sum['*','*']
    # Write dynamics file containing all state variables, W, and V_act for a in areas for t in [0,T]. Abraham's code (more is at the end):
    filename = f"donor_{print_area}global_simulation_{ALPHA_VALUE}_alpha.csv"
    with open("simulation_data/" + filename, 'w') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerows(zip(var_names, var_values))
        wr.writerows(zip(vax_var_names, vax_var_values))
    print(f"State Variables written to {filename}.")
     
    # Terminate. All that follows is an "else" but I didn't indent it. 

## Optimize

simulate()                                             

### Initialize Phase 1 of Optimize 
# Notation: f = z[i] = donor deaths, x = lambda[i]

phase = 1
lambda = {i: 0 for i in range(1, iter_lmt_search + 1)} # Define vector
lamba[1] = lambda_0
lambda[2] = lambda_0*phi
z = {i: 0 for i in range(1, iter_lmt_search + 1)} # Define vector and set z[0] = 0

# Iterations 1 and 2
i = 1
optimize_inner()     

fL = z[i]           #f at left end pt
i = 2
optimize_inner() 
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

While (i < iter_lmt_search and phase == 1)
    i = i + 1
    lambda[i] = mult*x2
    x3 = lambda[i]	        #3rd of 3 current x values
    optimize_inner() 
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
lambda[i] = x
optimize_inner() 
fx = z[i]                   #f(x)
i = i + 1
lambda[i] = y
optimize_inner() 
fy = z[i]                   #f(y)

# Phase 2 main loop

while (abs(fx - fy) > delta and i < iter_lmt_search)
    i = i + 1
    if fx > fy           		# minimum is in [a,x], so update y
        b, x, fx = (x, y, fy)
        y = b - 0.618 * (b - a)
        lambda[i] = y
        optimize_inner() 
        fy = z[i]
    else				        # minimum is in [y,b], so update x
        a, y, fy = (y, x, fx)
        x = a + 0.618 * (b - a)
        lambda[i] = x
        optimize_inner() 
        fx = z[i]
    if (fy < fx and fx > fb) print "Warning: f is not unimin" # print lambda[i-3],z[i-3],...,lambda[i],z[i]
    if (fy > fx and fa < fy) print "Warning: f is not unimin" # print lambda[i-3],z[i-3],...,lambda[i],z[i]

# Save the optimal objective value 
if (fy <= fx): 
    z_opt = fy 
    lambda_opt = y
else
    z_opt = fx 
    lambda_opt = x

# v.write("seir_opt.lp")  # This writes a file in LP format: equations are written out, with values of constants. Omit?
        # If kept, name the file using names of input files. Doing it here writes the last LP solved.

# Print z_opt = optimal donor deaths and lambda_opt. These could be from the last or next-to-last LP, but z_opt should only differ by delta
# Print using solution of last LP: 
    # donor deaths = D1.x[donor,T]
    # deaths = sum a in areas of D1.x[a,T]
    # day of variant tn
    # variant area m, 
    # vaccinations by area  = sum t in [0,T-1] of V1.x[a,t]
    # total vaccinations = sum a in areas sum t in [0,T-1] of V1.x[a,t]

# If verbosity == 2 (more verbose than 1), print
#   table of optimal values for each outer loop: 
#       z[i], z[i] - z_opt, and lambda[i]. Could be in order of i, but increasing order of lambda[i] would be better.
#   table of optimal LP values for each inner loop (for the last outer loop):
#       j, zLP[j], zLP[j] - z_opt

# Write dynamics file containing state variables S1.x[a,t], etc., W1.x[a,t], and V1.x[a,t] for a in areas for t in [1,T]
#   This is the same range of t as Abraham. Could insert constants S[a,0], etc. at beginning, but not necessary.


### Abraham's old code:  
    
    for var in opt_Mod.getVars():

    obj_function_values[iteration_number] = opt_Mod.getObjective().getValue()
    vaccine_allocation_totals[iteration_number] = sum(vax_dictionary.values())
    obj_difference = abs(obj_function_values[iteration_number] - obj_function_values[iteration_number - 1])
    infection_difference = sum([abs(i) for i in [
        infections[iteration_number, ar, "I", i] - infections[iteration_number - 1, ar, "I", i] for ar in areas for i in
        range(T)]])
    if iteration_number > 6: #check for loops
        #rounds the obj dictionary to be at the term tolerance
        tol = termination_tolerance/2
        rounded_obj = [round(i / tol) * tol for i in list(obj_function_values.values())[2:]]
        number_of_repeated_obj_values = len(rounded_obj)-len(set(rounded_obj))
        if number_of_repeated_obj_values > 2:
            print("Warning: The algorithm is looping. Terminating due to lack of convergence.")
            no_loop = False
    print(f"The objective difference is {obj_difference}, and the infection difference is {infection_difference}.")
    exploration_tolerance *= gamma
    #if iteration_number == 3:
    #    no_loop = False
    print(f"The exploration tolerance is {exploration_tolerance}")


print_solution = {ar: [] for ar in areas}
var_names = []
var_values = []
vax_var_names = []
vax_var_values = []
total_deaths_opt = 0
for var in opt_Mod.getVars():
    var_name_list = var.VarName.split(",")
    if "vax_variables" in var_name_list[0]:
        # format is ['vax_variables[ISO', '79]']
        day = int(var_name_list[-1][:-1])
        if not four_binary:
            area_name = var_name_list[0][-3:]
            print_solution[area_name].append(float(var.x))
        else:

            area_name = var_name_list[0][-5:]
            print_solution[area_name].append(float(var.x))
        if var.x >= 0:
            vax_var_names.append(str(var.varName))
            vax_var_values.append(var.x)
    if str(var.varName)[:5] == "state" and var.x >= 0:
        if "D," + str(T-1) in str(var.varName):
            total_deaths_opt += var.x
        var_names.append(str(var.varName))
        var_values.append(var.X)
        #if "S" in var.VarName:
        #    print(var.VarName, var.x)
vax_decision_variables = []
print("In each area, the vaccines given in each day (rounded) is as follows.")
for ar in areas:
    print(ar, [round(j, 1) for j in print_solution[ar]])
    vax_decision_variables.append(print_solution[ar])
print("Total vaccines given each day (opt)", np.array(vax_decision_variables).sum(axis=0))
print("Total vaccines given over the time horizon (opt)", np.array(vax_decision_variables).sum(axis=0).sum(), "\n")
print(f"The initial feasible solution had objective value (sim) {first_obj}.\n")
print(f"After {iteration_number} iterations, we found the solution (opt after time 1) {obj_function_values[iteration_number]}.\n")
print(f"There were {total_deaths_opt} total deaths (opt).")
print(f"The objective values were (opt) {obj_function_values}.\n")
print(f"The vaccine allocation totals (optimization) were {opt_vax_totals}.\n"
      f"This is (probably) lower than the simulated total, which was {vaccine_allocation_totals[iteration_number]}")
print(f"The total infected people entering the donor area was (sim) {sum(entering_donor.values())}.")
print(f"The total deaths in the donor area was (opt) {state_variables[donor, 'D', T-1]}")
print(f"The variant occured on day {t_N} (sim).")

if not four_binary:
    if policy:
        if small_run:
            print_area = f"some_countries_{p_k}_"
        else:
            print_area = f"all_countries_{p_k}_"
    else:
        if small_run:
            print_area = f"some_countries_"
        else:
            print_area = f"all_countries_"
else:
    print_area = "four_countries_"
filename = f"donor_{print_area}global_simulation_{ALPHA_VALUE}_alpha.csv"
with open("simulation_data/" + filename, 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerows(zip(var_names, var_values))
    wr.writerows(zip(vax_var_names, vax_var_values))
print(f"State Variables written to {filename}.")
opt_Mod.write("solution_sol.sol")
if policy:
    print(f"Strict policy = {strict_policy}, p^k = {p_k}")
elif vaccine_incentive:
    print("We used a vaccine incentive in the objective function, and no policy.")
else:
    print("No policy was used.")
print(f"The python file ran for a duration of {time.strftime('%H:%M:%S', time.gmtime(time.time() - start_time))}.")
