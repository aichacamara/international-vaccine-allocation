from io import TextIOWrapper
import os
import math
import numpy as np
import gurobipy as gp
from gurobipy import GRB
import xml.etree.ElementTree as ET
import argparse
import shutil
import csv
import sys

def simulate_switchover_policy():
    """Runs simulation for switchover
    """
    global new_priority
    if t_switch is None:
        # allocate vaccine to highest priority area
        for t in range(T - T0 + 1): # t=0,...,T-T0+1
            V[priority[0], t] = B[t]     # highest priority
    else:
        new_priority = priority # for consistency with simulate_only, where the priorities change
        t_prev = 0	# initialize previous switching time
        for z in range(n_a):      # area index 0, ..., n_a - 1
            if z < n_a - 1: 
                t_next =  t_switch[z]	            # set next switching time
                for t in range(t_prev, t_next):         # t=t_prev,..., t_next - 1
                    V[new_priority[z], t] = B[t] * (1 - split)	    # allocate to area specified in switching policy
                    for z1 in range(z + 1, n_a):
                        if z1 > z:          # divide proportion split b/t lower-priority areas
                            V[new_priority[z1], t] = B[t] * split / (n_a - 1 - z)
                t_prev = t_next                         # update for next area	 
            else: 
                t_next =  T	                            # for last area, next switching time is T
                for t in range(t_prev, t_next):
                    V[A[z], t] = B[t]                      # for last area, no splitting
          
def main():
    ## start_time = time.time() #work for SP but not MV
    global S0, SV0, E0, EV0, I0, IV0, W0, S1, SV1, E1, EV1, I1, IV1, D1, R1, W1, V1, v, z, i, phase, fn_base
    # global deaths, donor_deaths, tot_deaths, t_sim # from opt_inner
    global fn, csv_file #output files
    
    global new_priority

    # read input file, compute constants
    import_xml(xml_path=os.getcwd() + "/" + input_file)
    # initialize output filename and directory
    try:
        os.mkdir(os.getcwd() + "/" + "output")
    except:
        pass 
    fn_base = f'./output/{input_filename.split("/")[-1][0:-4]}_nu{nu:3.1f}' #_nu{nu:3.1f} output file name uses inputs
    ##sys.stdout = open(fn_base + "_con.out", "w") # console output redirected to file

    # Initialize state variables
    S0 = {a: 0 for a in A}
    SV0 = {a: 0 for a in A}
    E0 = {a: 0 for a in A}  
    EV0 = {a: 0 for a in A}
    I0 = {a: 0 for a in A}  
    IV0 = {a: 0 for a in A}
    W0 = {a: 0 for a in A}
    for a in A:
        E0[a] = ((1 - rho_V[a])/(p_r*rho_V[a] +
                      1 - rho_V[a]))*(rho_I_N[a]/r_I)
        EV0[a] = ((p_r*rho_V[a])/(p_r*rho_V[a] +
                        1 - rho_V[a]))*(rho_I_N[a]/r_I)
        I0[a] = ((1 - rho_V[a])/(p_r*rho_V[a] + 1 - rho_V[a]))* \
                rho_I_N[a]/r_d[a]
        IV0[a] = ((p_r*rho_V[a])/(p_r*rho_V[a] + 1 - rho_V[a]))* \
                rho_I_N[a]/r_d[a]
        SV0[a] = rho_V[a]*N[a] - EV0[a] - IV0[a]
        S0[a] = N[a] - E0[a] - EV0[a] - I0[a] - IV0[a] - SV0[a]
        W0[a] = rho[a]*N[a] - SV0[a] - EV0[a] - IV0[a] - rho[a]*E0[a] - rho[a]*I0[a]
        
    # Initialize V for initial sim using priority list. Use splitting between areas.
    global V
    V = {(a, t): 0 for a in A for t in range(T)}       # t=0,...,T-1
    simulate_switchover_policy()

    if not simulate_only:
        # Initialize LP Variables. LP will be updated (objective, constraints) in solve_LP
        v = gp.Model("vaccine_opt")
        v.Params.LogToConsole = 0 # 1: Gurobi output turned on
        v.Params.LPWarmStart = 2    # Allows presolve with PStart, which otherwise might prevent presolve. See PStart docum.
        #v.Params.Method = 2 # 2: barrier method, -1: default: automatic; for LP, usually concurrently runs all 3 methods 
        # LP variables are S1 for state var S, etc. All are continuous and nonnegative by default.
        S1 = v.addVars(A, range(1, T + 1), name="S1")          # t = 1, ..., T
        SV1 = v.addVars(A, range(1, T + 1), name="SV1")
        E1 = v.addVars(A, range(1, T + 1), name="E1")
        EV1 = v.addVars(A, range(1, T + 1), name="EV1")
        I1 = v.addVars(A, range(1, T + 1), name="I1")
        IV1 = v.addVars(A, range(1, T + 1), name="IV1")
        #H1 = v.addVars(A, range(1, T + 1), name="H1")
        D1 = v.addVars(A, range(1, T + 1), name="D1")
        R1 = v.addVars(A, range(1, T + 1), name="R1")
        W1 = v.addVars(A, range(1, T + 1), name="W1")
        V1 = v.addVars(A, range(T), name="V1")             # t = 0, ..., T-1

        # Iteration 0. Uses t_n, alpha, V_cal from first sim, lambda_0. Updates t_n, alpha, V_cal 
        # but z[0] is NOT used in search b/c initial V may give better z[0] than z[1], giving wrong search direction for lambda
        # Define z, l (lambda); initialize l[0], l[1], l[2]
        z = {i1: 0 for i1 in range(iter_lmt_search + 4)}  
        l = {i1: 0 for i1 in range(iter_lmt_search + 4)}
        l[0] = lambda_0
        l[1] = lambda_0
        l[2] = lambda_0*phi
        i = 0    # outer iteration
        phase = 1
        z[i] = optimize_inner(l[i], V)
        if iter_lmt_search > 1:
            # Iterations 1 and 2
            i = 1
            z[i] = optimize_inner(l[i], V)
            fL = z[i]  # f at left end pt

            i = 2
            z[i] = optimize_inner(l[i], V)
            fR = z[i]  # f at right end pt

            if fL <= fR:  # f is increasing, search to left
                mult = 1/phi  # x multiplier < 1
                x1 = l[2]  # largest of 3 current x values
                x2 = l[1]  # middle of 3 current x values
                f1 = fR  # f(x1)
                f2 = fL  # f(x2)
            if fL > fR:  # f is decreasing, search to right
                mult = phi  # x multiplier > 1
                x1 = l[1]  # smallest of 3 current x values
                x2 = l[2]  # middle of 3 current x values
                f1 = fL
                f2 = fR

            # Phase 1 loop
            while (i < max(iter_lmt_search, 3) and phase == 1): # At least 3 iter needed to execute this loop, set x3
                i = i + 1
                l[i] = mult*x2
                x3 = l[i]  # 3rd of 3 current x values
                z[i] = optimize_inner(l[i], V)
                f3 = z[i]  # f(x3)
                if f3 > f2 or (f3 == f2 and f3 < z[1]):
                    phase = 2  # x3 is past the minimum or (flat here but not initially) 
                x1 = x2  # shift x for next Phase 1 iteration
                x2 = x3
                f1 = f2
                f2 = f3
            if i < iter_lmt_search:
                # Phase 2: Golden ratio search on interval [a, b] with check for unimin
                # Initialize Phase 2
                if x1 <= x3:
                    a0 = x1
                    b0 = x3
                    fa = f1
                    fb = f3

                if x1 > x3:
                    a0 = x3
                    b0 = x1
                    fa = f3
                    fb = f1

                # Two more iterations, at x and y
                x = a0 + 0.618 * (b0 - a0)  # current larger x value
                y = b0 - 0.618 * (b0 - a0)  # current smaller x value
                i = i + 1
                l[i] = x
                z[i] = optimize_inner(l[i], V)

                fx = z[i]  # f(x)
                i = i + 1
                l[i] = y
                z[i] = optimize_inner(l[i], V)
                fy = z[i]  # f(y)

                # Phase 2 loop
                while (abs(fx - fy) > delta and i < iter_lmt_search):
                    i = i + 1
                    if fx > fy:        		    # minimum is in [a1,x], so update b
                        b0, fb, x, fx = (x, fx, y, fy)
                        y = b0 - 0.618 * (b0 - a0)
                        l[i] = y
                        z[i] = optimize_inner(l[i], V)
                        fy = z[i]
                    else:	                    # minimum is in [y,b], so update a1
                        a0, fa, y, fy = (y, fy, x, fx)
                        x = a0 + 0.618 * (b0 - a0)
                        l[i] = x
                        z[i] = optimize_inner(l[i], V)
                        fx = z[i]
                    if (fy < fx and fx > fb):
                        print("Warning: f is not unimin: y < fx and fx > fb")
                    if (fy > fx and fa < fy):
                        print("Warning: f is not unimin: fy > fx and fa < fy")
        elapsed_time = 0
        ## elapsed_time = time.time() - start_time # Worked for SM but not MV
        print("\nLP count: ", LP_count, "Infeas count: ", infeas_count, 
                "Time elapsed: ", elapsed_time, "s")
        print("Number of areas: ", len(A), " Iter_lmt: ", iter_lmt,
                " Iter_lmt_search: ", iter_lmt_search)

    # open output files
    fn = open(fn_base + ".out", "w")
    csv_file = open(fn_base + "_plot.csv", "w") 
    
    o_input_echo()
    
    if not simulate_only:     
        # Compute total vacc by area (sim, opt, and min)
        global V_tot_sim, V_tot_min, V_tot_opt
        V_tot_sim = {a: 0 for a in A}
        V_tot_min = {a: 0 for a in A}
        V_tot_opt = {a: 0 for a in A}
        for a in A:
            for t in range(T - T0 + 1): # t=0,...,T-T0
                V_tot_sim[a] += V_table[a, t, 0, 0]
                V_tot_opt[a] += V_table[a, t, i_opt, j_opt]
                V_tot_min[a] += V_table[a, t, i, j_min[i]]
        # first sim
        o_loop_report()

    else: # Simulate only
        fn.write("Simulate. Policy gives vaccine to one area, then reallocates using priorities  " + input_file + "\n\n") 
        fn.write("policy     wtd_deaths donor_deaths tot_deaths t_n   vacc by area\n")
        global V_sim
        V_sim = {(a1, a, t): 0 for a1 in A for a in A for t in range(T)}       # To store V by priority policy
        for a1 in A: 
            if t_switch is None:
                # Initialize V giving priority to area a1
                V = {(a, t): 0 for a in A for t in range(T)}       # t=0,...,T-1
                for t in range(T - T0 + 1): # t=0,...,T-T0+1
                    V[a1, t] = B[t]
                continue
             # Initialize V giving initial priority to area a1, then use priority list. Use splitting between areas.

            priority.remove(a1)
            priority.insert(0,a1)

            V = {(a, t): 0 for a in A for t in range(T)}    # t=0,...,T-1
            t_prev = 0	# initialize previous switching time
            for z in range(n_a):      # area index 0, ..., a_n - 1
                if z < n_a - 1: 
                    t_next =  t_switch[z]	                        # set next switching time
                    for t in range(t_prev, t_next):                 # t=t_prev,..., t_next - 1
                        V[new_priority[z], t] = B[t] * (1 - split)	# allocate to area specified in switching policy
                        for z1 in range(z + 1, n_a):
                            if z1 > z:          # divide proportion split b/t lower-priority areas
                                V[new_priority[z1], t] = B[t] * split / (n_a - 1 - z)
                    t_prev = t_next                                 # update for next area	 
                else: 
                    t_next =  T	                                    # for last area, next switching time is T
                    for t in range(t_prev, t_next):
                        V[A[z], t] = B[t]                           # for last area, no splitting
            # Simulate   
            t_sim, alpha, V_cal, V, D = simulate(V)
            deaths_sim_only = (1 - nu)*D[donor, T] 
            for a in A:
                deaths_sim_only += nu*D[a, T]
            donor_deaths_sim_only = D[donor, T]
            tot_deaths_sim_only = 0
            for a in A:
                tot_deaths_sim_only += D[a, T]
            for a in A:                     # Store V for this policy
                for t in range(T - T0 + 1): # t=0,...,T-T0+1
                    V_sim[a1, a, t] = V[a, t]

            # Compute total vacc by area
            V_tot_sim = {a: 0 for a in A} 
            for a in A:
                for t in range(T - T0 + 1):
                    V_tot_sim[a] += V[a, t]
                    
            o_policy_report(a1, deaths_sim_only, donor_deaths_sim_only, tot_deaths_sim_only, V_tot_sim)
        o_loop_report()

def optimize_inner(l, V): 
    global deaths, donor_deaths, tot_deaths, t_n, zNLP, V_table, \
        j_min, j_max, alpha_min, V_cal_min, V_min, \
        alpha_prev, V_cal_prev, V_prev, \
        deaths_opt, i_opt, j_opt, \
        j, LP_count, infeas_count, dVmax, dVmin

    if i == 0:
        LP_count = infeas_count = 0
        # Define arrays
        j_min = {i1: 0 for i1 in range(iter_lmt_search + 4)} 
            # iter j of inner loop that achieves best zNLP for each i 
        deaths = {(i1, j1): 0 for i1 in range(iter_lmt_search + 4) for j1 in range(iter_lmt + 1)} 
            # wtd deaths for each i, j
        donor_deaths = {(i1, j1): 0 for i1 in range(iter_lmt_search + 4) for j1 in range(iter_lmt + 1)} 
        tot_deaths = {(i1, j1): 0 for i1 in range(iter_lmt_search + 4) for j1 in range(iter_lmt + 1)}
        t_n = {(i1, j1): 0 for i1 in range(iter_lmt_search + 4) for j1 in range(iter_lmt + 1)} 
            # Time variant emerges
        zNLP = {(i1, j1): 0 for i1 in range(iter_lmt_search + 4) for j1 in range(iter_lmt + 1)} 
            # LP objective fcn value. Includes Lagrangain
        V_table = {(a, t, i1, j1): 0 for a in A for t in range(T+1) \
            for i1 in range(iter_lmt_search + 4) for j1 in range(iter_lmt + 1)} 
            # Vaccinations by area and day. t=0,...,T for use in CSV
        dVmax = {(i1, j1): 0 for i1 in range(iter_lmt_search + 4) for j1 in range(iter_lmt + 1)} 
        dVmin = {(i1, j1): 0 for i1 in range(iter_lmt_search + 4) for j1 in range(iter_lmt + 1)} 
            # Change in V_cal, sim - LP. Measures convergence.

        # First simulate. Record as i = 0, j = 0
        t_n[0, 0], alpha, V_cal, V, D = simulate(V)
        # Store sim V in V_table
        for t in range(T - T0 + 1): # t=0,...,T-T0+1
            for a in A:
                V_table[a, t, 0, 0] = V[a, t]
        # Compute deaths, zNLP from sim. This zNLP is used in output if it is best for i = 0 (min)
        deaths[0, 0] = (1 - nu)*D[donor, T]
        donor_deaths[0, 0] = D[donor, T]
        tot_deaths[0, 0] = 0
        for a in A:
            deaths[0, 0] += nu*D[a, T]
            tot_deaths[0, 0] += D[a, T]
        # Use fixed t_n in LP for all j, but include dT more days since t_n may increase
        t_LP = min(math.ceil(t_n[0, 0]) + dT, T) 
        zNLP[0, 0] = deaths[0, 0] + l*sum(I[a, t] for a in A_D for t in range(1, t_LP + 1)) \
          - (1e-9)*sum(V[a, t]*(T-t) for a in A for t in range(T)) # sum I over A_D: no donor in Lagr
        # Initialize min from this sim (for solve_LP)
        alpha_min = alpha
        V_cal_min = V_cal
        V_min = V
        # Initialize opt from this sim
        deaths_opt = deaths[0, 0]
        i_opt = i
        j_opt = 0 # iter j of inner loop that achieves deaths_opt

    elif i > 0 and phase == 1:  # During phase 1 use best solution from i = 0 as the init solution 
        # to make them more comparable. Init solution for i goes in j = 0. 
        deaths[i, 0] = deaths[0, j_min[0]]
        donor_deaths[i, 0] = donor_deaths[0, j_min[0]]
        tot_deaths[i, 0] = tot_deaths[0, j_min[0]]
        zNLP[i, 0] = zNLP[0, j_min[0]]
        t_n[i, 0] = t_n[0, j_min[0]]
        alpha = alpha_prev
        V_cal = V_cal_prev
        V = V_prev

    else:   # Phase 2: use best solution from previous i as the init solution
            # Init solution for i goes in j = 0. (alpha, V_cal needed for solve_lp)
        deaths[i, 0] = deaths[i-1, j_min[i-1]]
        donor_deaths[i, 0] = donor_deaths[i-1, j_min[i-1]]
        tot_deaths[i, 0] = tot_deaths[i-1, j_min[i-1]]
        zNLP[i, 0] = zNLP[i-1, j_min[i-1]]
        t_n[i, 0] = t_n[i-1, j_min[i-1]]
        alpha = alpha_min
        V_cal = V_cal_min
        V = V_min
        # Store V_min in V_table at j = 0 for output
        for t in range(T - T0 + 1): # t=0,...,T-T0+1
            for a in A:
                V_table[a, t, i, 0] = V_min[a, t]
 
    # Initialize for j loop
    # Use fixed t_n in LP for all j, but include dT more days since t_n may increase
    t_LP = min(math.ceil(t_n[i, 0]) + dT, T) # For i=0, t_LP was already computed. This doesn't change it.
    j = 0 # inner loop counter
    eps = eps_prev = epsilon_0 # eps_prev is eps at the previous best sim
    zLP = 1e14              # initialize value of LP. Used for stopping.
    zLP_prev = 2*zLP        # initialize previous value of zLP. They must differ. Used for stopping.

    while (abs(zLP - zLP_prev) >= delta_I or j < 2) and j < iter_lmt:
        j += 1
        LP_count += 1
        zLP_prev = zLP
        # solve LP
        if improving == 1:
            zLP = solve_LP(l, t_LP, alpha_min, V_cal_min, eps)  # Improving: alpha_min, V_cal_min are best for this i
        else:
            zLP = solve_LP(l, t_LP, alpha, V_cal, eps)   # Not improving: alpha, V_cal from sim of current LP
        if v.status == GRB.OPTIMAL: # If LP infeas, update eps and iterate. All remaining j likely infeas. 
            V = {(a, t): V1[a, t].x for a in A for t in range(T)} # update V using optimal V1 from LP
            t_n[i, j], alpha, V_cal, V, D = simulate(V)   # simulate LP solution 
            # compute deaths, NLP objective from sim
            deaths[i, j] = (1 - nu)*D[donor, T]
            donor_deaths[i, j] = D[donor, T]
            tot_deaths[i, j] = 0  
            for a in A:
                deaths[i, j] += nu*D[a, T]
                tot_deaths[i, j] += D[a, T]
            zNLP[i, j] = deaths[i, j] + l*sum(I1[a, t].x for a in A_D for t in range(1, t_LP + 1)) \
                - (1e-9)*sum(V1[a, t].x*(T-t) for a in A for t in range(T)) # A_D: no donor in Lagr
            # Store sim V in V_table
            for t in range(T - T0 + 1): # t=0,...,T-T0
                for a in A:
                    V_table[a, t, i, j] = V[a, t]
            # Compute difference in V_cal, simulate - LP
            # Ignore behavior G. Easier than using G from simulate but must recompute LP. 
            dVmax_0 = -10000
            dVmin_0 = 10000
            for t in range(1, T+1):
                for a in A:
                    dV = (I[a, t] + p_e*I_V[a, t]) - (I1[a, t].x + p_e*IV1[a, t].x)
                    dVmax_0 = max(dVmax_0, dV) # max dV so far
                    dVmin_0 = min(dVmin_0, dV) # min dV so far
            dVmax[i, j] = dVmax_0
            dVmin[i, j] = dVmin_0

            if zNLP[i, j] <= zNLP[i, j_min[i]]:  # update best j, alpha, V_cal, V for this i (min)
                j_min[i] = j
                alpha_min = alpha
                V_cal_min = V_cal
                V_min = V
                if i == 0:
                    alpha_prev = alpha # previous is best j for i = 0, used in Phase 1 as init
                    V_cal_prev = V_cal
                    V_prev = V 
                eps = eps_prev      # reset eps to its value at the previous best zNLP
                eps_prev *= beta    # reduce eps b/c better zNLP achieved

            if deaths[i, j] < deaths_opt:  # update best deaths (opt)
                i_opt = i
                j_opt = j
                deaths_opt = deaths[i, j]

        else:
            infeas_count += 1
        eps *= beta              # reduce current eps
    if i == 0: # Save iter 0 min solution in "prev" for use initializing all phase 1 iter 
        alpha_prev = alpha_min
        V_cal_prev = V_cal_min
        V_prev = V_min
    return deaths[i, j_min[i]]

def solve_LP(l, t_LP, alpha, V_cal, eps):
    global S1, SV1, E1, EV1, I1, IV1, D1, R1, W1, V1, v

    # Warm start using bases (can also use solution)
    if LP_count - infeas_count > 1: #if v.status == GRB.OPTIMAL: 
        vbas = {i: v.getVars()[i].VBasis for i in range(v.NumVars)} # Save basis for variables
        ##psol = {i: v.getVars()[i].x for i in range(v.NumVars)}
        cbas = {i: v.getConstrs()[i].CBasis for i in range(v.NumConstrs)} # Save basis for constraints (dual)
        ##dsol = {i: v.getConstrs()[i].Pi for i in range(v.NumConstrs)}

    v.setObjective((1 - nu)*D1[donor, T] + nu*D1.sum('*', T)+ l*sum(I1[a, t] for a in A_D for t in range(1, t_LP + 1))
        - (1e-9)*sum(V1[a, t]*(T - t) for a in A for t in range(T)), GRB.MINIMIZE)  
    
    # Some constraints change with alpha, V_cal (inner loop). Rewrite them all.
    if LP_count - infeas_count > 1:
        v.remove([constraint for constraint in v.getConstrs()]) # Remove all constraints

    v.addConstrs((V1[donor, t] <= p_k*B[t]
                 for t in range(T - T0 + 1)), "Policy_donor_limit")
    v.addConstrs((V1.sum('*', t) <= B[t] for t in range(T - T0 + 1)), "Vaccine_budget")
    v.addConstrs((V1.sum('*', t) == 0 for t in range(T - T0 + 1, T)), "No_vaccine")

    # Dynamics for start time t=1 to T-1
    v.addConstrs((W1[a, t+1] == W1[a, t] - alpha[a, t]*V_cal[a, t]/N[a]*W1[a, t] - V1[a, t]
                  for a in A for t in range(1, T)), "Vaccine_willingness")
    v.addConstrs((S1[a, t+1] == S1[a, t] - alpha[a, t]*V_cal[a, t]/N[a]*S1[a, t] - V1[a, t]
                  for a in A for t in range(1, T)), "S")
    v.addConstrs((SV1[a, t+1] == SV1[a, t] - p_r*alpha[a, t]*V_cal[a, t]/N[a]*SV1[a, t] + V1[a, t]
                  for a in A for t in range(1, T)), "SV")
    v.addConstrs((E1[a, t+1] == E1[a, t] + alpha[a, t]*V_cal[a, t]/N[a]*S1[a, t] - r_I*E1[a, t]
                  for a in A for t in range(1, T)), "E")
    v.addConstrs((EV1[a, t+1] == EV1[a, t] + p_r*alpha[a, t]*V_cal[a, t]/N[a]*SV1[a, t] - r_I*EV1[a, t]
                  for a in A for t in range(1, T)), "EV")
    v.addConstrs((I1[a, t+1] == I1[a, t] + r_I*E1[a, t] - r_d[a]*I1[a, t]
                  for a in A for t in range(1, T)), "I")
    v.addConstrs((IV1[a, t+1] == IV1[a, t] + r_I*EV1[a, t] - r_d[a]*IV1[a, t]
                  for a in A for t in range(1, T)), "IV")
    v.addConstrs((D1[a, t+1] == D1[a, t] + r_d[a]*p_D*I1[a, t] + r_d[a]*p_V_D*IV1[a, t]
                  for a in A for t in range(1, T)), "D")
    v.addConstrs((R1[a, t+1] == R1[a, t] + r_d[a]*(1 - p_D)*I1[a, t] + r_d[a]*(1 - p_V_D)*IV1[a, t]
                  for a in A for t in range(1, T)), "R")

    # Dynamics for start time t=0. Use global constants S0[a],SV0[a],E0[a],EV0[a],I0[a],IV0[a]  
    # constant V_cal[a, 0] but variable V1. Assume D[a, 0] = R[a, 0] = 0 (they aren't defined)
    v.addConstrs((W1[a, 1] == W0[a] - alpha[a, 0]*V_cal[a, 0]/N[a]*W0[a] - V1[a, 0]
                  for a in A), "Vaccine_willingness_t=0")
    v.addConstrs((S1[a, 1] == S0[a] - alpha[a, 0]*V_cal[a, 0]/N[a]*S0[a] - V1[a, 0]
                  for a in A), "S_t=0")
    v.addConstrs((SV1[a, 1] == SV0[a] - p_r*alpha[a, 0]*V_cal[a, 0]/N[a]*SV0[a] + V1[a, 0]
                  for a in A), "SV_t=0")
    v.addConstrs((E1[a, 1] == E0[a] + alpha[a, 0]*V_cal[a, 0]/N[a]*S0[a] - r_I*E0[a]
                  for a in A), "E_t=0")
    v.addConstrs((EV1[a, 1] == EV0[a] + p_r*alpha[a, 0]*V_cal[a, 0]/N[a]*SV0[a] - r_I*EV0[a]
                  for a in A), "EV_t=0")
    v.addConstrs((I1[a, 1] == I0[a] + r_I*E0[a] - r_d[a]*I0[a]
                  for a in A), "I_t=0")
    v.addConstrs((IV1[a, 1] == IV0[a] + r_I*EV0[a] - r_d[a]*IV0[a]
                  for a in A), "IV_t=0")
    v.addConstrs((D1[a, 1] == r_d[a]*p_D*I0[a] + r_d[a]*p_V_D*IV0[a]
                  for a in A), "D_t=0")
    v.addConstrs((R1[a, 1] == r_d[a]*(1 - p_D)*I0[a] + r_d[a]*(1 - p_V_D)*IV0[a]
                  for a in A), "R_t=0")

    # Regularization constraints on V_cal: constant in t
    v.addConstrs((G[a, t]*I1[a, t] + G[a, t]*p_e*IV1[a, t] <= V_cal[a, t] + eps #*(t/T) ##*(t/T)^2 #Behavior
                 for a in A for t in range(1, T)), "V_cal_upper_bd")
    v.addConstrs((G[a, t]*I1[a, t] + G[a, t]*p_e*IV1[a, t] >=  max(V_cal[a, t] - eps, 0) ##V_cal[a, t] - eps #*(t/T) ##*(t/T)^2 #Behavior 
                 for a in A for t in range(1, T)), "V_cal_lower_bd")

    if LP_count -infeas_count > 1: # Warm start using VBasis/CBasis if v.status == GRB.OPTIMAL: #
        v.update() # Must update to refer to new constraints
        for i in range(v.NumVars):
            v.getVars()[i].VBasis = vbas[i]
            #v.getVars()[i].PStart = psol[i]
        for i in range(v.NumConstrs):
            v.getConstrs()[i].CBasis = cbas[i]
            #v.getConstrs()[i].DStart = dsol[i]   
    v.optimize()
    if v.status == GRB.OPTIMAL:
        return v.ObjVal
    elif v.status == GRB.INFEASIBLE:
        print("Infeasible LP at lambda =", l, "LP iteration j =", j, "exploration tol eps =", eps)
        return 100000000
    else:
        print("Optimize failed at lambda =", l, "LP iteration j =", j, "exploration tol eps =", eps)
        print("Status =", v.status)
        return 100000000

def simulate(V):
    global I, I_V, G
    # Define state variables, alpha, delta_E, V_star, V_minus
    S = {(a, t): 0 for a in A for t in range(T+1)}     # t=0,...,T b/c computed in diff eq
    S_V = {(a, t): 0 for a in A for t in range(T+1)}
    E = {(a, t): 0 for a in A for t in range(T+1)}   
    E_V = {(a, t): 0 for a in A for t in range(T+1)}
    I = {(a, t): 0 for a in A for t in range(T+1)}   
    I_V = {(a, t): 0 for a in A for t in range(T+1)}
    R = {(a, t): 0 for a in A for t in range(T+1)}
    #H = {(a, t): 0 for a in A for t in range(T+1)}
    D = {(a, t): 0 for a in A for t in range(T+1)}
    W = {(a, t): 0 for a in A for t in range(T+1)}
    V_cal = {(a, t): 0 for a in A for t in range(T+1)}  # t=0,...,T b/c used in output
    G = {(a, t): 0 for a in A for t in range(T+1)}      # Behavior
    alpha = {(a, t): 0 for a in A for t in range(T)}    #  t=0,...,T-1 b/c not computed in diff eq
    delta_E = {(a, t): 0 for a in A for t in range(T)}
    V_star = {(a, t): 0 for a in A for t in range(T)}
    V_minus = {a: 0 for a in A}                       # stores V_star before realloc

    # Assign initial states
    for a in A:
        S[a, 0] = S0[a]
        S_V[a, 0] = SV0[a]
        E[a, 0] = E0[a]
        E_V[a, 0] = EV0[a]
        I[a, 0] = I0[a]
        I_V[a, 0] = IV0[a]
        W[a, 0] = W0[a]

    # Compute alpha_0. Initialize t_sim, variant_emerge
    alpha_0 = {(a): a_0 * gamma[a] for a in A}
    t_sim = T
    variant_emerge = False

    for t in range(T):   # t = 0,...,T-1
        # compute alpha if variant found
        if variant_emerge:
            a_t = a_0 + delta_a/(1 + math.exp(-k*(t - (t_sim + T_D))))
            alpha[m, t] = a_t * gamma[m]
            for a in A:
                if a != m:
                    if t - L < 0:
                        alpha[a, t] = alpha_0[a]
                    else:
                        alpha[a, t] = a_0 + delta_a / \
                            (1 + math.exp(-k*(t - L - (t_sim + T_D))))
        else:                           # constant alpha
            for a in A:
                alpha[a, t] = alpha_0[a]

        # Compute V_cal, delta_E, and V_star (w/o realloc)
        for a in A:                                 
            V_cal[a, t] = I[a, t] + p_e*I_V[a, t]
            G[a, t] = 1                             # Behavior...
            if V_cal[a, t] > N[a]*v_l[a]:
               G[a, t] = 1 - (1 - g[a])/(v_u[a] - v_l[a])*(V_cal[a, t] - N[a]*v_l[a])/N[a]
            if V_cal[a, t] > N[a]*v_u[a]:
                G[a, t] = g[a]
            V_cal[a, t] *= G[a, t] 
            delta_E[a, t] = min(S[a, t], alpha[a, t]*S[a, t]*V_cal[a, t]/N[a])
            if S[a, t] < 0.0000001:
                V_star[a, t] = min(W[a, t], V[a, t])
            else:
                V_star[a, t] = min(W[a, t] - W[a, t]*delta_E[a, t]/S[a, t], V[a, t])

        # Realloc: any number of areas. Uses priority.
        realloc = 0 # amount avail to realloc
        for a in A:
            realloc += V[a, t] - V_star[a, t] # unused vacc
            V_minus[a] = V_star[a, t] # Store vacc before realloc
        for a in priority: # loop over areas in priority order
            if S[a, t] < 0.0000001:
                Wnew = W[a, t] # limit on V_star
            else:
                Wnew = W[a, t] - W[a, t]*delta_E[a, t]/S[a, t]
            V_star[a, t] = min(Wnew, V_minus[a] + realloc) # Realloc as much as allowed to this area
            realloc -= V_star[a, t] - V_minus[a] # subtract amount realloc to this area 

        # Difference equations including W
        for a in A:
            if S[a, t] < 0.0000001:
                W[a, t + 1] = W[a, t] - V_star[a, t]
            else:
                W[a, t + 1] = W[a, t] - W[a, t]*delta_E[a, t]/S[a, t] - V_star[a, t]
            S[a, t + 1] = S[a, t] - delta_E[a, t] - V_star[a, t]
            S_V[a, t + 1] = S_V[a, t] + V_star[a, t] - p_r * \
                alpha[a, t]*(S_V[a, t]/N[a])*V_cal[a, t]
            E[a, t + 1] = E[a, t] + delta_E[a, t] - r_I*E[a, t]
            E_V[a, t + 1] = E_V[a, t] + p_r*alpha[a, t] * \
                (S_V[a, t]/N[a]) * V_cal[a, t] - r_I*E_V[a, t]
            I[a, t + 1] = I[a, t] + r_I*E[a, t] - r_d[a]*I[a, t]
            I_V[a, t + 1] = I_V[a, t] + r_I * E_V[a, t] - r_d[a]*I_V[a, t]
            D[a, t + 1] = D[a, t] + r_d[a]*p_D*I[a, t] + r_d[a]*p_V_D*I_V[a, t]
            R[a, t + 1] = R[a, t] + r_d[a]*(1 - p_D)*I[a, t] + r_d[a]*(1 - p_V_D)*I_V[a, t]

        # check for the variant occuring, do not calculate if variant already emerged
        if not variant_emerge:
            I_sum = I_max = 0 # Isum = sum over a and t, Imax = max over a of sum over t
            for a in A:
                I_current = 0 # sum over t for current a
                if a != donor: ## donor doesn't contribute to variant
                    for t1 in range(t+1): # t1=0,...,t 
                        I_sum += I[a, t1]     # sum all infected up to current time
                        I_current += I[a, t1] # sum current area up to current time
                    if I_max < I_current:     # update area that has max infections
                        area_max = a
                        I_max = I_current
            # If cum infections > n, variant emerges. Compute t_sim and m = area where variant emerges
            if I_sum > n:
                variant_emerge = True
                m = area_max # variant emerges in area with max cum I
                I_tot = 0
                for a in A:
                    if a != donor:
                        I_tot += I[a, t]
                t_sim = t + 1 - (I_sum - n)/(I_tot)
    # Compute V_cal at T (used in LP and opt output)
    for a in A:
        V_cal[a, T] = I[a, T] + p_e*I_V[a, T]
        G[a, T] = 1                          # Behavior...
        if V_cal[a, T] > N[a]*v_l[a]:
            G[a, T] = 1 - (1 - g[a])/(v_u[a] - v_l[a])*(V_cal[a, T] - N[a]*v_l[a])/N[a]
        if V_cal[a, T] > N[a]*v_u[a]:
            G[a, T] = g[a]
        V_cal[a, T] *= G[a, T]               

    if simulate_only:
        # Write the csv (simulate)
        with open(fn_base + "_plot" + ".csv", "w") as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                ["area", "t", "S", "SV", "E", "EV", "I", "IV", "H", "D", "R", "W", "V", "t_n", "L"])
            csv_writer.writerow(
                [m, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_sim, L])
            for a in A:
                for t in range(T + 1):
                    if t != T:
                        csv_writer.writerow([a, t, S[a, t], S_V[a, t], E[a, t], E_V[a, t], I[a, t],
                                            I_V[a, t], 0, D[a, t], R[a, t], W[a, t], V_star[a, t], t_sim, L])
                    else:
                        csv_writer.writerow([a, t, S[a, t], S_V[a, t], E[a, t], E_V[a, t], I[a, t],
                                             I_V[a, t], 0, D[a, t], R[a, t], W[a, t], 0, t_sim, L])
    return t_sim, alpha, V_cal, V_star, D

######################################## INPUT HELPERS ########################################

def xml_text(element: ET.Element, name: str):
    try:
        return convert_num(element.find(name).text)
    except:
        return None

def import_xml(xml_path: str): # Read inputs from XML file. xml_path: path to the XML file
    root = ET.parse(xml_path).getroot()

    area_data = root.find("area_data")
    scenario_data = root.find("scenario_data")
    params = root.find("params")

    # read area data
    global A, A_D, gamma, rho_I_N, rho, rho_V, delta_r, N, priority, t_switch, split, donor, m, n, n_a
    A = [] # all areas
    A_D = [] # areas except donor
    gamma = {}
    rho_I_N = {}
    rho = {}
    rho_V = {}
    delta_r = {}
    N = {}

    priority = area_data.find("priority").text
    # if there is no priority, assign priority to empty array
    if priority == None:
        priority = []
    priority = priority.split(sep=",")

    donor = area_data.find("donor").text
    # m = area_data.find("m").text  # computing varaint area, so don't read it
    n = convert_num(area_data.find("n").text)
    for child in area_data.findall("area"):
        area = child.attrib["name"]
        A.append(area)
        if area != donor:
            A_D.append(area)
        gamma[area] = convert_num(child.find("gamma").text)
        rho_I_N[area] = convert_num(child.find("rho_I_N").text)
        rho[area] = convert_num(child.find("rho").text)
        rho_V[area] = convert_num(child.find("rho_V").text)
        delta_r[area] = convert_num(child.find("delta_r").text)
        N[area] = convert_num(child.find("N").text)
    n_a = len(A)
    
    # read scenario data
    global T, B_0, nu, p_k, r_I, r_0, p_D, p_V_D, a_0, delta_a, \
        p_e, p_r, L, T_D, p, b_arr, v_u, v_l, g
    
    # switchover policy
    global t_switch, split
    
    T = convert_num(scenario_data.find("T").text)
    B_0 = convert_num(scenario_data.find("B_0").text) 
    nu = convert_num(scenario_data.find("nu").text)
    p_k = convert_num(scenario_data.find("p_k").text)
    r_I = convert_num(scenario_data.find("r_I").text) 
    r_0 = convert_num(scenario_data.find("r_0").text)
    p_D = convert_num(scenario_data.find("p_D").text)
    p_V_D = convert_num(scenario_data.find("p_V_D").text)
    a_0 = convert_num(scenario_data.find("a_0").text)
    delta_a = convert_num(scenario_data.find("delta_a").text)
    p_e = convert_num(scenario_data.find("p_e").text)
    p_r = convert_num(scenario_data.find("p_r").text)
    L = convert_num(scenario_data.find("L").text)
    T_D = convert_num(scenario_data.find("T_D").text)
    p = convert_num(scenario_data.find("p").text)
    b_arr = scenario_data.find("b").text
    v_u = convert_num(scenario_data.find("v_u").text) 
    t_switch = scenario_data.find("t_switch")
    split = scenario_data.find("switch_split")
    
    if t_switch is not None:
        t_switch = t_switch.text.split(",")
        for i in range(len(t_switch)):
            t_switch[i] = convert_num(t_switch[i])
    
    if split is not None:
        split = convert_num(split.text)
    
    if not b_arr == None:
        b_arr = b_arr.split(sep=",")
        for i in range(len(b_arr)):
            b_arr[i] = convert_num(b_arr[i])
    else:
        b_arr = []
    # Behavior dynamics hardcoded, using input v_u
    v_input = v_u
    v_l = {} 
    v_u = {}
    g = {}
    if v_input: # linear behavior
        for a in A:
            v_l[a] = 0
            v_u[a] = v_input 
            g[a] = 0  
    else:           # no behavior
        for a in A:
            v_l[a] = 0
            v_u[a] = 1 
            g[a] = 1   

    # read params
    global simulate_only, lambda_0, phi, epsilon_0, delta_I, \
        delta, beta, iter_lmt, iter_lmt_search, dT, verbosity, T0,\
        improving, realloc_flag, t_priority, t_priority_vector  #realloc_flag not currently used
    simulate_only = bool(convert_num(params.find("simulate_only").text))
    lambda_0 = convert_num(params.find("lambda_0").text)
    phi = convert_num(params.find("phi").text)
    epsilon_0 = convert_num(params.find("epsilon_0").text)
    delta_I = convert_num(params.find("delta_I").text)
    delta = convert_num(params.find("delta").text)
    beta = convert_num(params.find("beta").text)
    iter_lmt = convert_num(params.find("iter_lmt").text)
    iter_lmt_search = convert_num(params.find("iter_lmt_search").text)
    dT = convert_num(params.find("dT").text)
    verbosity = convert_num(params.find("verbosity").text)
    

    # Hard-coded parameters
    T0 = 3
    improving = 0
    # Compute k, B, r_d
    global k, B, r_d
    k = math.log((1-p)/p)/T_D  # natural log, ln()
    B = {(t): B_0 * 1 for t in range(T)}
    for t in range(0,len(b_arr)):
        B[t] = B_0 * b_arr[t]
    r_d = {(a): r_0 + delta_r[a] for a in A}

def convert_num(num: str): #Converts an input string "num" to float or int
    if float(num) % 1 != 0:
        return float(num)
    return int(num)

######################################## OUTPUT HELPERS ########################################

def o_state_equations(fn: TextIOWrapper,
                           num_length: int,
                           lower_limit: int,
                           upper_limit: int,
                           state_name: str,
                           state: dict):
    """
    Helper function to generate file output listing state values over area and time in a formatted way

    Params:
        fn: the file to have the output written to
        num_length: the length of the spacing numbers are centered inside of
        lower_limit: the lower limit of the time bounds
        upper_limit: the upper limit of the time bounds
        state_name: the name of the outputted state values
        state: the object containing the values for area and time

    Returns:
        None
    """
    fn.write(state_name + "\n")
    fn.write(f'{"day": ^{num_length}}')
    for a in A:
        fn.write(f'{a: ^{num_length}}')
    fn.write("\n")
    for t in range(lower_limit, upper_limit):
        fn.write(f'{t: ^{num_length}}')
        for a in A:
            fn.write(
                f'{round(state[a, t], 2): ^{num_length}}')
        fn.write("\n")

def o_input_echo():
    """
    Generates input file echo
    If optimization scheme is in play, function will output to .csv
    
    Returns:
    (int) -> 1 on Success, 0 on Fail
    """
    try:
        if not simulate_only:
            #output to .csv
            global csv_file
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                ["area", "t", "S", "SV", "E", "EV", "I",
                    "IV", "H", "D", "R", "W", "V", "t_n", "L"]
            )
            csv_writer.writerow(
                [m, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_n[i_opt,j_opt], L]
            )
            for a in A:
                csv_writer.writerow(
                    [a, 0, S0[a], SV0[a], E0[a], EV0[a], I0[a],
                    IV0[a], 0, 0, 0, S0[a], V_table[a, 0, i_opt, j_opt], t_n[i_opt,j_opt], L] # was V_min[a, 0], t_min 
                )
                for t in range(1, T + 1):
                    csv_writer.writerow(
                        [a, t, S1[a, t].x, SV1[a, t].x, E1[a, t].x, EV1[a, t].x, I1[a, t].x,
                            IV1[a, t].x, 0, D1[a, t].x, R1[a, t].x, W1[a, t].x, 
                            V_table[a, 0, i_opt, j_opt], t_n[i_opt,j_opt], L] 
                        )       # was V_min[a, t], t_min, which may be from a diff LP than S1,...
        
    
        #output to .out
        if verbosity >= 1:
            # input echo
            fn.write("--------------------------------Area data--------------------------------" + "\n")
            fn.write("Area name:                        ")
            for name in A:
                fn.write(name + " ")
            fn.write("\n")
            fn.write("Population:                      ")
            for area in A:
                fn.write(str(N[area]) + " ")
            fn.write("\n")            
            fn.write("Vaccination Rate:                ")
            for area in A:
                fn.write(str(rho_V[area]) + " ")
            fn.write("\n")            
            fn.write("Initial cases per day = rho_I*N: ")
            for area in A:
                fn.write(str(rho_I_N[area]) + " ")
            fn.write("\n")
            fn.write("Rate out of state I:             ")
            for area in A:
                fn.write(str(r_d[area]) + " ")
            fn.write("\n")
            fn.write("Behavior infection multiplier:   ")
            for area in A:
                fn.write(str(gamma[area]) + " ")
            fn.write("\n")
            fn.write("Proportion willing to be vacc:   ")
            for area in A:
                fn.write(str(rho[area]) + " ")
            fn.write("\n")

            fn.write("n: " + str(n) + "\n")
            fn.write("------------------------------Scenario data------------------------------" + "\n")
            fn.write("Time horizon (days): " + str(T) + "\n")
            fn.write("Vaccine avaiable day 0: " + str(B_0) + "\n")
            fn.write("Weight for non-donor deaths in objective: " + str(nu) + "\n")
            fn.write("Upper limit on proportion infectious due to behavior: " + str(v_u) + "\n")
            fn.write("Max. prop. of vaccine allocated to donor areas: " + str(p_k) + "\n")
            fn.write("Rate out of state E into I:" + str(r_I) + "\n")
            fn.write("Rate out of state I w/o testing: " + str(r_0) + "\n")
            fn.write("P(death | infected, not vacc): " + str(p_D) + "\n")
            fn.write("P(death | infected, vacc): " + str(p_V_D) + "\n")
            fn.write("Initial infection rate: " + str(a_0) + "\n")
            fn.write("Change in infection rate for variant: " + str(delta_a) + "\n")
            fn.write("Prop. transmission from a vaccinated person: " + str(p_e) + "\n")
            fn.write("Prop. transmission to a vaccinated person: " + str(p_r) + "\n")
            fn.write("Lag for variant to reach other areas (days): " + str(L) + "\n")
            fn.write("Time for variant to dominate (days): " + str(T_D) + "\n")
            fn.write("Prop. of people in state I that have the new variant when introduced: " + str(p) + "\n")
    
            fn.write("-------------------------------Parameters--------------------------------" + "\n")
            fn.write("Simulate only: " + str(simulate_only) + "\n")
            fn.write("Priority (decreasing): ")
            for a1 in range(len(priority)): 
                fn.write(str(priority[a1]) + " ")
            fn.write("\n")
            fn.write("Lagrange multiplier for infection: " + str(lambda_0) + "\n")
            fn.write("Exploration multiplier for lambda: " + str(phi) + "\n")
            fn.write("Exploration tolerance for LP: " + str(epsilon_0) + "\n")
            fn.write("Termination tolerance for LP: " + str(delta_I) + "\n")
            fn.write("Termination tolerance for lambda: " + str(delta) + "\n")
            fn.write("Exploration convergence parameter for LP: " + str(beta) + "\n")
            fn.write("Iteration limit for LP: " + str(iter_lmt) + "\n")
            fn.write("Iteration limit for lambda: " + str(iter_lmt_search) + "\n")
            fn.write("Days after t_n[0] in Lagrangian: " + str(dT) + "\n")
            fn.write("Verbosity: " + str(verbosity) + "\n")
            fn.write("-------------------------------------------------------------------------" + "\n")
        return 1
    except:
        print("Failed to output Input Echo")
        return 0

def o_policy_report(a1, deaths_sim_only, donor_deaths_sim_only, tot_deaths_sim_only, V_tot_sim):
    """Outputs policy report

    Args:
        a1 : element of list A, current priority view
        deaths_sim_only (int): deaths by simulation
        donor_deaths_sim_only (int): donor deaths by simulation
        tot_deaths_sim_only (int): total deaths
        V_tot_sim (int): vaccination total by simulation

    Returns:
        int: 0 on fail, 1 on success
    """
    if simulate_only:
        try:
            fn.write(f'{a1: ^{9}}  {deaths_sim_only: 8.2f}  {donor_deaths_sim_only: 8.2f}  {tot_deaths_sim_only: 12.2f}  {t_sim: 6.2f} ')
            for a in A:
                fn.write(f'{V_tot_sim[a]: 5.0f} ')                    
            fn.write("\n")
            return 1
        except:
            return 0
    else:
        try:
            return 1
        except:
            return 0

def o_loop_report():
    """
    Outputs loop report based on verbosity and (global) simulate_only
    
    Returns:
    (int) -> 1 on Success, 0 on Fail
    """
    if simulate_only:
        try:
            # Verbosity 1
            if verbosity >= 1:
                for a1 in A:
                    fn.write("\nVaccinations, priority to " + a1)
                    fn.write("\n  day    V by area \n")
                    for t in range(T - T0 + 1):
                        fn.write(f'{t: ^{7}}')
                        for a in A:
                            fn.write("  " + str(V_sim[a1, a, t]) + "  ")
                        fn.write("\n")
            return 1
        except:
            return 0
    else:
        try:
            # Verbosity 0
            fn.write("Convergence: Min/Max change in V_cal, (sim - LP)\n\n")
               
            fn.write( "                           --------deaths-------_\n")
            fn.write( "i  j     lambda    zNLP    weighted donor   total    t_n    conv of V_cal     vacc by area\n")
            fn.write("first simulation\n")
            fn.write(f'0  0       0        0    {deaths[0,0]: 8.2f} {donor_deaths[0,0]: 8.2f} {tot_deaths[0,0]: 8.2f} {t_n[0,0]: 6.2f}                  ')
            for a in A:
                fn.write(f'{V_tot_sim[a]: 5.0f} ')                    
            fn.write("\n")

            fn.write(f'optimal (best deaths found)\n')

            fn.write(f'{i_opt: ^{2}} {j_opt: ^{2}} {l[i_opt]: 9.4f} {zNLP[i_opt,j_opt]: 8.2f} ')
            fn.write(f'{deaths[i_opt,j_opt]: 8.2f} {donor_deaths[i_opt,j_opt]: 8.2f} {tot_deaths[i_opt,j_opt]: 8.2f} ')
            fn.write(f'{t_n[i_opt,j_opt]: 6.2f} ({dVmin[i_opt,j_opt]: 7.1f},{dVmax[i_opt,j_opt]: 6.1f}) ')
            for a in A:
                fn.write(f'{V_tot_opt[a]: 5.0f} ')                    
            fn.write("\n") 

            fn.write(f'minimum (best zNLP w/ Lagrangian for last lambda) w/ convergence for last LP\n')

            fn.write(f'{i: ^{2}} {j_min[i]: ^{2}} {l[i]: 9.4f} {zNLP[i,j_min[i]]: 8.2f} ')
            fn.write(f'{deaths[i,j_min[i]]: 8.2f} {donor_deaths[i,j_min[i]]: 8.2f} {tot_deaths[i,j_min[i]]: 8.2f} ')
            fn.write(f'{t_n[i,j_min[i]]: 6.2f} ({dVmin[i,j_min[i]]: 7.1f},{dVmax[i,j_min[i]]: 6.1f}) ')
            for a in A:
                fn.write(f'{V_tot_min[a]: 5.0f} ')                    
            fn.write("\n\n") 
            
            # Verbosity 2
            if verbosity >= 2:
                fn.write("Outer Loop over lambda. j_min = iter of inner loop that achieves best wtd deaths\n")
                fn.write("iter  lambda j_min  zNLP  wtd_deaths  subopt  t_n   conv of V_cal\n")
                
                for i1 in range(i+1):
                    fn.write(f'{i1: ^{2}} {l[i1]: 9.4f}  {j_min[i1]: ^2} {zNLP[i1,j_min[i1]]: 8.2f} ')
                    fn.write(f'{z[i1]: 8.2f} {z[i1] - deaths_opt: 9.2f} {t_n[i1,j_min[i1]]: 6.2f} ')
                    fn.write(f'({dVmin[i1,j_min[i1]]: 7.1f},{dVmax[i1,j_min[i1]]: 6.1f})')
                    fn.write("\n")
                fn.write("\n")

                fn.write("Inner Loop at last i (last lambda)\n")
                fn.write("iter  zNLP subopt w/in this i wtd_deaths subopt  t_n   conv of V_cal\n")
                for j1 in range(j+1):
                    fn.write(f'{j1: ^{2}} {zNLP[i,j1]: 8.2f} {zNLP[i,j1] - zNLP[i,j_min[i]]: 9.2f}       ')
                    fn.write(f'{deaths[i,j1]: 8.2f} {deaths[i,j1] - deaths_opt: 9.2f}  {t_n[i,j_min[i]]: 6.2f} ')
                    fn.write(f'({dVmin[i,j1]: 7.1f},{dVmax[i,j1]: 6.1f})')
                    fn.write("\n")
                fn.write("\n")

            # Verbosity 1
            if verbosity >= 1:
                fn.write("Vaccinations: sim of best LP (best j), last lambda (min) \n"
                        "  day    Vacc by area \n")
                for t in range(T - T0 + 1):
                    fn.write(f'{t: ^{7}}')
                    for a in A:
                        fn.write("  " + str(V_table[a, t, i, j_min[i]]) + "  ")
                    fn.write("\n")
                fn.write("\nOptimal vaccinations \n"
                        "  day    Vacc by area \n")
                for t in range(T - T0 + 1):
                    fn.write(f'{t: ^{7}}')
                    for a in A:
                        fn.write("  " + str(V_table[a, t, i_opt, j_opt]) + "  ")
                    fn.write("\n")
            return 1
        except:
            return 0
    

########################################### Script Run ###########################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # First positional argument (this must be present)
    parser.add_argument('input', type=str, help='Directory of input xml file')

    # Parse the command line
    args = parser.parse_args()    
    
    input_dir = args.input
    
    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir,f))]
    
    for f in files:
        """
        @TODO
            Highly suggested to remove the global property of the 
            `input_file` 
            `input_filename` 
            as these files are only used in main() for printing and csv parsing,
            they can be passed in as parameters
        """
        global input_file
        global input_filename
        input_filename = f
        input_file = os.path.join(input_dir,input_filename)
        main()