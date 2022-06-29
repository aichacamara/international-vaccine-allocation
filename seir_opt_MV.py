# SEIR-OPT model: Overall structure and Optimize
# Need to add:
    # Reading inputs: I have listed the names of the inputs, which you should recognize but might differ slightly from yours.
    # Simulate: I call a function simulate(), which is intended to be what you wrote. Printing/writing output should NOT be in simulate().
      # Computing alpha[a,t] should be in simulate(). If not, add it to formulate_LP() 
    # Writing "simulation data" file: I refer to Abraham's code. It won't work, partly b/c we chaged from  I don't know how to modify it.
    # Error checking, e.g., whether the two simulation limits were reached
    # Better stopping criteria for outer loop: |fx - fy| < delta does NOT imply fx or fy is within delta of min. 
# Some parts are not fully written or not correct syntax. Where I don't know the syntax, I generally use comments. Need to fix:
    # Gurobi variable referencing. After solving, to get the optimal value of V[a,t], use V.x[a,t]. 
    #   This needs to be fixed when V, I, IV, D are used.
    # Variables and Gurobi variables with same name: V is used in many places. W and other state varaibles are used in simulate(). 
    # All these name conflicts need to be resolved, either by making variables local or changing names. 
    #   For example, Gurobi variables could have a one appended: V1, etc.
    # Initial conditions W[a,0] etc. are variables, not Gurobi variables. They will need a different name than the Gurobi variables 
    #   and constraints at t = 0 will have to be written separately, e.g., W1[a,1] = W[a,0] + ... (W1 is the Gurobi var)   

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

from areas import * #areas, donor, m, n, N[a], rho_V[a], rho_I[a], delta_r[a], gamma[a], rho[a] 
from scenario import * #rI, r0, rR, pH, pHV, pD, a0, delta_a, pe, pr, L, TD, p, T, B_0, b[t]
from parameters import * #pk, lambda_0, phi, eps_0, delta_I, delta, beta, iter_lmt, iter_lmt_search # Changed name from lambda to lambda_0

# Computed constants

donor = "area1"
rd = {a: r0 + delta_r[a]) for a in areas}
B = {t: B0 * b(t) for T in [0,T-1]}  # If no b(t) in input file, use b(t) = 1
k = log((1-p)/p)/TD

# num_areas = len(areas)    # Not needed?

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

W[a,0]  = {a: rho[a]*N[a] - SV[a,0] - EV[a,0] - IV[a,0] - rho[a]*E[a,0] - rho[a]*I[a,0] for a in areas} 

# Initialize first Simulate

S_dot = S.sum(a,0) - S[donor,0]     # Sum of nondonor S at time 0
V[a,t]  = {a,t: (1 - pk) * S[a,0]/Sdot * B[t] for a in areas for t in [0,T-1]}
V[donor,t]  = {t: pk * B[t] for t in [0,T-1]}

#Create a simulate object

def simulate():
    # key input: V[a,t]
    # key outputs:  V_act[a,t] (actual vacc, called V* in paper), I_hat[a,t], IV_hat[a,t], tn, m

# Create the inner loop of Optimize

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
        simulate()                                             
        eps = beta*eps
    z[i] = D[donor,T]   #Donor deaths from current LP solution ("D" is Gurobi variable, not simulate variable!)
                        #Could move outside of this function.

# Create object that formulates and solves LP

def formulate_LP()          
    # key inputs: (I_hat, IV_hat, alpha) for all areas and times, tn, m, lambda[i], eps
    # key outputs:  LP model
    
    t_int = ceil(tn)
    V_cal = {a, t: I_hat[a, t] + pe*IV_hat[a, t] for a in areas for t in range(T-1) } # "vectors" (constants)

    # Compute alpha[a,t] here if not done in simulate()?? 
      # Needs to include constant rate before tn or tn + L and variable rate after

    v = gp.Model("vaccine_opt")
    # v.setParam("NumericFocus",1)  # Increased protection from roundoff/instability. Shouldn't need.

    # Define Variables. All are continuous and nonnegative by default. State var's 

    S = v.AddVars(areas, range(1,T+1), name="S")
    SV = v.AddVars(areas, range(1,T+1), name="SV")
    E = v.AddVars(areas, range(1,T+1), name="E")
    EV = v.AddVars(areas, range(1,T+1), name="EV")
    I = v.AddVars(areas, range(1,T+1), name="I")
    IV = v.AddVars(areas, range(1,T+1), name="IV")
    H = v.AddVars(areas, range(1,T+1), name="H")
    D = v.AddVars(areas, range(1,T+1), name="D")
    R = v.AddVars(areas, range(1,T+1), name="R")
    W = v.AddVars(areas, range(1,T+1), name="W")
    V = v.AddVars(areas, T, name="V")

    if non_donor_deaths_flag == 0: #or FALSE
        m.setObjective(D[donor, T] + lambda[i]*sum(I[a, t] for a in areas for t in range(t_int) )
            + 1/200000*D.sum('*', T), GRB.MINIMIZE)  # include non-donor deaths w/ small weight
    else:
            m.setObjective(D[donor, T] + lambda[i]*sum(I[a, t] for a in areas for t in range(t_int) ),
            GRB.MINIMIZE)
   
    v.AddConstrs((V[donor, t] <= pk*B[t] for t in range(T-1)), "Policy: donor limit") # constraint must be in () if name argument used
    v.AddConstrs((V.sum('*', t) <= B[t] for t in range(T-1)), "Vaccine budget")                  

    v.AddConstrs((W[a,t+1] == W[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*W[a,t] - V[a,t]
                        for a in areas for t in range(T-1)), "Vaccine willingness")
    v.AddConstrs((S[a,t+1] == S[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*S[a,t] - V[a,t]
                        for a in areas for t in range(T-1)), "S")                    
    v.AddConstrs((SV[a,t+1] == SV[a,t] - pr*alpha[a,t]*V_cal[a,t]/N[a]*SV[a,t] + V[a,t]
                        for a in areas for t in range(T-1)), "SV")
    v.AddConstrs((E[a,t+1] == E[a,t] + alpha[a,t]*V_cal[a,t]/N[a]*S[a,t] - rI*E[a,t]
                        for a in areas for t in range(T-1)), "E")  
    v.AddConstrs((EV[a,t+1] == EV[a,t] + pr*alpha[a,t]*V_cal[a,t]/N[a]*SV[a,t] - rI*EV[a,t]
                        for a in areas for t in range(T-1)), "EV")
    v.AddConstrs((I[a,t+1] == I[a,t] + rI*E[a,t] - rd*I[a,t]
                        for a in areas for t in range(T-1)), "I")
    v.AddConstrs((IV[a,t+1] == IV[a,t] + rI*EV[a,t] - rd*IV[a,t]
                        for a in areas for t in range(T-1)), "IV")
    v.AddConstrs((H[a,t+1] == H[a,t] + rd*pH*I[a,t] + rd*pHV*IV[a,t] - rR*H[a,t] 
                        for a in areas for t in range(T-1)), "H")
    v.AddConstrs((D[a,t+1] == D[a,t] + rR*pD*H[a,t]
                        for a in areas for t in range(T-1)), "D")
    v.AddConstrs((R[a,t+1] == R[a,t] + rR*(1 - pD)*H[a,t] + rd*(1 - pH)*I[a,t] + rd*(1 - pHV)*IV[a,t]
                        for a in areas for t in range(T-1)), "R")

    v.AddConstrs((I[a,t] - I_hat[a,t] <= eps
                        for a in areas for t in range(1,T)), "I upper bd")
    v.AddConstrs((I[a,t] - I_hat[a,t] >= -eps
                        for a in areas for t in range(1,T)), "I lower bd")
    v.AddConstrs((IV[a,t] - IV_hat[a,t] <= eps
                        for a in areas for t in range(1,T)), "IV upper bd")
    v.AddConstrs((IV[a,t] - IV_hat[a,t] >= -eps
                        for a in areas for t in range(1,T)), "IV lower bd")

 def solve_LP()          
    # key inputs: I_hat, IV_hat, tn, m, lambda[i], eps
    # key outputs:  V_new[a,t] (actual vacc), zLP[j] (optimal objective value)

    t_int = ceil(tn)
    V_cal = {a, t: I_hat[a, t] + pe*IV_hat[a, t] for a in areas for t in range(T-1) } # "vectors" (constants)

    # Compute alpha[a,t] here if not done in simulate()?? 
      # Needs to include constant rate before tn or tn + L and variable rate after

    # Objective changes with lambda (outer loop) and tn (inner loop)

    if non_donor_deaths_flag == 0: #or FALSE
        m.setObjective(D[donor, T] + lambda[i]*sum(I[a, t] for a in areas for t in range(t_int) )
            + 1/200000*D.sum('*', T), GRB.MINIMIZE)  # include non-donor deaths w/ small weight
    else:
            m.setObjective(D[donor, T] + lambda[i]*sum(I[a, t] for a in areas for t in range(t_int) ),
            GRB.MINIMIZE)

    # Constraints change with alpha, V_cal (inner loop)
    # Assumes that AddConstrs REPLACES constraints with the same name.

    v.AddConstrs((W[a,t+1] == W[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*W[a,t] - V[a,t]
                        for a in areas for t in range(T-1)), "Vaccine willingness")
    v.AddConstrs((S[a,t+1] == S[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*S[a,t] - V[a,t]
                        for a in areas for t in range(T-1)), "S")                    
    v.AddConstrs((SV[a,t+1] == SV[a,t] - pr*alpha[a,t]*V_cal[a,t]/N[a]*SV[a,t] + V[a,t]
                        for a in areas for t in range(T-1)), "SV")
    v.AddConstrs((E[a,t+1] == E[a,t] + alpha[a,t]*V_cal[a,t]/N[a]*S[a,t] - rI*E[a,t]
                        for a in areas for t in range(T-1)), "E")  
    v.AddConstrs((EV[a,t+1] == EV[a,t] + pr*alpha[a,t]*V_cal[a,t]/N[a]*SV[a,t] - rI*EV[a,t]
                        for a in areas for t in range(T-1)), "EV")

    # Constraints change with I_hat, IV_hat, eps (inner loop)

    v.AddConstrs((I[a,t] - I_hat[a,t] <= eps
                        for a in areas for t in range(1,T)), "I upper bd")
    v.AddConstrs((I[a,t] - I_hat[a,t] >= -eps
                        for a in areas for t in range(1,T)), "I lower bd")
    v.AddConstrs((IV[a,t] - IV_hat[a,t] <= eps
                        for a in areas for t in range(1,T)), "IV upper bd")
    v.AddConstrs((IV[a,t] - IV_hat[a,t] >= -eps
                        for a in areas for t in range(1,T)), "IV lower bd")

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
lambda = {i: 0 for i in [0,iter_lmt_search]} # Define vector
lamba[1] = lambda_0
lambda[2] = lambda_0*phi
z = {i: 0 for i in [0,iter_lmt_search]} # Define vector and set z[0] = 0

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

While (iter < iter_lmt_search and phase == 1)
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
lambda[i] = X
optimize_inner() 
fx = z[i]                   #f(x)
i = i + 1
lambda[i] = y
optimize_inner() 
fy = z[i]                   #f(y)

# Phase 2 main loop

while (abs(fx - fy) > delta and iter < iter_lmt_search)
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
    lambda_opt = X

# v.write("seir_opt.lp")  # This writes a file in LP format: equations are written out, with values of constants. Omit?
        # If kept, name the file using names of input files. Doing it here writes the last LP solved.

# Print z_opt = optimal donor deaths and lambda_opt. These could be from the last or next-to-last LP, but z_opt should only differ by delta
# Print using solution of last LP: donor deaths D[donor,T], deaths D.sum['*',T], day of variant tn, variant area m, 
        # vaccinations by area  V.sum[area,'*'], total vaccinations V_act.sum['*','*']

# Write dynamics file containing all state variables, W, and V for a in areas for t in [0,T]
#   These are Gurobi variables, not from simulate()













    
    
    
    
    
    
    
    
    
    
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
