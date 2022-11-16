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

def main():
    global S0, SV0, E0, EV0, I0, IV0, W0, S1, SV1, E1, EV1, I1, IV1, D1, R1, W1, V1, v, iter, V_opt, \
        z, donor_deaths, tot_deaths, t_sim, dVmax, dVmin, phase, fn_base #H1,

    # read input file, compute k and B
    import_xml(xml_path=os.getcwd() + "/" + input_file)
    # initialize output filename and directory
    try:
            os.mkdir(os.getcwd() + "/" + "output")
    except:
        pass 
    fn_base = f'./output/{input_filename.split("/")[-1][0:-4]}_tsw{t_switch:03d}' #_nu{nu:3.1f} output file name uses inputs
    sys.stdout = open(fn_base + "_con.out", "w") # console output redirected to file

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
        I0[a] = ((1 - rho_V[a])/(p_r*rho_V[a] + 1 -
                      rho_V[a]))*(rho_I_N[a]/(r_0 + delta_r[a]))
        IV0[a] = ((p_r*rho_V[a])/(p_r*rho_V[a] +
                        1 - rho_V[a]))*(rho_I_N[a]/(r_0 + delta_r[a]))
        SV0[a] = rho_V[a]*N[a] - EV0[a] - IV0[a]
        S0[a] = N[a] - E0[a] - EV0[a] - I0[a] - IV0[a] - SV0[a]
        W0[a] = rho[a]*N[a] - SV0[a] - EV0[a] - IV0[a] - rho[a]*E0[a] - rho[a]*I0[a]

    # Initialize V
    V = {(a, t): 0 for a in A for t in range(T)}       # t=0,...,T-1
    # 2 areas: priority to area1 for t < t_switch
    if n_a == 2:
        for t in range(min(t_switch,T - T0 + 1)): # t = 0,...,t_switch-1
            V[A[0], t] = B[t]
        for t in range(t_switch, T - T0 + 1): # t = t_switch,...,T-T0. Keep V=0 for t = T-T0+1,...,T
            V[A[1], t] = B[t]
    else:
    # More than 2 areas: Initialize V using p_k, compute S_dot = Sum of nondonor S at time 0
        S_dot = 0
        for a in A:
            if a != donor:
                S_dot += S0[a]
        for a in A:
            for t in range(T - T0 + 1): # Keep V=0 for t = T-T0+1,...,T
                if a != donor:
                    V[a, t] = (1 - p_k) * (S0[a]/S_dot) * B[t]
                else:
                    V[a, t] = p_k * B[t]

    if not simulate_only:
        # Initialize LP Variables. LP will be updated (objective, constraints) in solve_LP
        v = gp.Model("vaccine_opt")
        #v.Params.DualReductions = 0 ##test: to see if unbounded or infeasible
        v.Params.LogToConsole = 0 # 1: Gurobi output turned on
        v.Params.LPWarmStart = 2    # Allows presolve with PStart, which otherwise might prevent presolve. See PStart docum. 
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
        z = {i: 0 for i in range(iter_lmt_search + 1)}  
        l = {i: 0 for i in range(iter_lmt_search + 1)}
        l[0] = lambda_0
        l[1] = lambda_0
        l[2] = lambda_0*phi
        iter = i = 0    # i is reassigned writing outputs, so store number of outer iterations in global "iter"
        phase = 1
        dVmin = dVmax = 0 # In case init sim is optimal
        z[i] = optimize_inner(l[i], V)
        if iter_lmt_search > 1:
            # Iterations 1 and 2
            iter = i = 1
            z[i] = optimize_inner(l[i], V)
            fL = z[i]  # f at left end pt

            iter = i = 2
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
                iter = i = i + 1
                l[i] = mult*x2
                x3 = l[i]  # 3rd of 3 current x values
                z[i] = optimize_inner(l[i], V)
                f3 = z[i]  # f(x3)
                if f3 > f2 or (f3 == f2 and f3 < z[1]): ##tested: + delta
                    phase = 2  # x3 is past the minimum or (flat here but not initially) 
                x1 = x2  # shift x’s for next Phase 1 iteration
                x2 = x3
                f1 = f2
                f2 = f3

            # Phase 2: Golden ratio search on interval [a, b] with check for unimin
            # Initialize Phase 2
            if x1 <= x3:
                a = x1
                b = x3
                fa = f1
                fb = f3

            if x1 > x3:
                a = x3
                b = x1
                fa = f3
                fb = f1

            # Two more iterations, at x and y
            x = a + 0.618 * (b - a)  # current larger x value
            y = b - 0.618 * (b - a)  # current smaller x value
            iter = i = i + 1
            l[i] = x
            z[i] = optimize_inner(l[i], V)

            fx = z[i]  # f(x)
            iter = i = i + 1
            l[i] = y
            z[i] = optimize_inner(l[i], V)
            fy = z[i]  # f(y)

            # Phase 2 loop
            while (abs(fx - fy) > delta and i < iter_lmt_search):
                iter = i = i + 1
                if fx > fy:        		    # minimum is in [a,x], so update b
                    b, fb, x, fx = (x, fx, y, fy)
                    y = b - 0.618 * (b - a)
                    l[i] = y
                    z[i] = optimize_inner(l[i], V)
                    fy = z[i]
                else:	                    # minimum is in [y,b], so update a
                    a, fa, y, fy = (y, fy, x, fx)
                    x = a + 0.618 * (b - a)
                    l[i] = x
                    z[i] = optimize_inner(l[i], V)
                    fx = z[i]
                if (fy < fx and fx > fb):
                    print("Warning: f is not unimin: y < fx and fx > fb")
                if (fy > fx and fa < fy):
                    print("Warning: f is not unimin: fy > fx and fa < fy")
        print("LP count: ", LP_count, "Infeas count: ", infeas_count)

        # Write the csv (optimize)
        #with open("./output/plot_" + input_file.split("/")[-1][0:-4] + ".csv", "w") as csv_file: OLD: no param in fn; plot_ ...
        with open(fn_base + "_plot.csv", "w") as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                ["area", "t", "S", "SV", "E", "EV", "I",
                    "IV", "H", "D", "R", "W", "V", "t_n", "L"]
            )
            csv_writer.writerow(
                [m, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_opt, L]
            )
            for a in A:
                csv_writer.writerow(
                    [a, 0, S0[a], SV0[a], E0[a], EV0[a], I0[a],
                     IV0[a], 0, 0, 0, S0[a], V_min[a, 0], t_min, L] # was V_opt[a, 0], t_opt 
                )
                for t in range(1, T + 1):
                    csv_writer.writerow(
                        [a, t, S1[a, t].x, SV1[a, t].x, E1[a, t].x, EV1[a, t].x, I1[a, t].x,
                            IV1[a, t].x, 0, D1[a, t].x, R1[a, t].x, W1[a, t].x, V_min[a, t], t_min, L] # was V_opt[a, 0], t_opt
                        )           # V_min, t_min may be from a diff LP than S1,...
 
        # Write output file (optimize)
        #with open("./output/" +  input_file.split("/")[-1][0:-4] + ".log", "w") as fn: OLD: no param in fn
        with open(fn_base + ".out", "w") as fn:
            # input echo
            fn.write("Optimization " + input_file + "\n\n")
            fn.write("Time Horizon: " + str(T) +  "  nu: " + str(nu) + "\n")
            fn.write("Donor in t_n: no\n")
            fn.write("Improving search: " + str(improving) + "\n")
            if n_a == 2:
                fn.write("Initial policy priority to " + str(2 - np.sign(t_switch)) + "  t_switch: " + str(t_switch) + "\n")
            fn.write("variant emergence n: " + str(n) + "  L: " + str(L) + "  T_D: " + str(T_D) + "\n")
            fn.write("Mortality p_D: " + str(p_D) + "  p_V_D: " + str(p_V_D) \
                + "  Vacc effectiveness p_e: " + str(p_e) + "  p_r: " + str(p_r) +"\n")
            fn.write("Rates r_I: " + str(r_I) + "  r_0: " + str(r_0) + "  a_0: " + str(a_0) + "  delta_a: " + str(delta_a) + "\n")
            fn.write("gamma by area: ") 
            for a in A:
                fn.write(str(gamma[a]) + "  ")
            fn.write("\n")
            fn.write("Lambda Converg. lambda: " + str(lambda_0) + "  phi: " + str(phi) + "  dT: " + str(dT) + \
                    "  delta: " + str(delta) + "  iter_search: " + str(iter_lmt_search) + "\n")
            fn.write("LP Converg. delta_I: " + str(delta_I) + "  beta: " + str(beta) + "  iter: " + str(iter_lmt) + "\n\n")
            # Verbosity 0
            fn.write("Convergence: Min/Max change in V_cal, (sim - LP)\n\n")        

            # first sim
            fn.write( "i  j     lambda    zNLP    deaths tot_deaths t_n    conv of V_cal     vacc by area\n")
            fn.write("first simulation\n")
            fn.write(f'0  0       0        0    {deaths_sim: 8.2f} {tot_deaths_sim: 8.2f}  {t_sim: 6.2f}                   ')
            for a in A:
                fn.write(f'{V_tot_sim[a]: 5.0f} ')                    
            fn.write("\n")

            fn.write(f'optimal (best deaths found)\n')

            # Compute total vacc by area (opt)
            V_total = {a: 0 for a in A} 
            for a in A:
                for t in range(T):
                    V_total[a] += V_opt[a, t] 
            fn.write(f'{i_opt: ^{2}} {j_opt: ^{2}} {lambda_opt: 9.4f}          ')
            fn.write(f'{z_opt: 8.2f} {tot_deaths_opt: 8.2f}  {t_opt: 6.2f} ({dVmin: 7.1f},{dVmax: 6.1f}) ')
            for a in A:
                fn.write(f'{V_total[a]: 5.0f} ')                    
            fn.write("\n") 

            fn.write(f'minimum (best zNLP for last lambda) w/ convergence for last LP\n')
            # Compute difference in V_cal, last simulate - last LP
            dVmax = -10000
            dVmin = 10000
            for t in range(1, T+1):
                for a in A:
                    dV = (I[a, t] + p_e*I_V[a, t]) - (I1[a, t].x + p_e*IV1[a, t].x)  
                    dVmax = max(dVmax, dV) # max dV so far
                    dVmin = min(dVmin, dV) # min dV so far
            # Compute total vacc by area (min)
            V_total = {a: 0 for a in A} 
            for a in A:
                for t in range(T):
                    V_total[a] += V_min[a, t] 
            fn.write(f'{iter: ^{2}} {j_min[iter]: ^{2}} {l[iter]: 9.4f} {zNLP_min: 8.2f} ')
            fn.write(f'{deaths_min: 8.2f} {tot_deaths_min: 8.2f}  {t_min: 6.2f} ({dVmin: 7.1f},{dVmax: 6.1f}) ')
            for a in A:
                fn.write(f'{V_total[a]: 5.0f} ')                    
            fn.write("\n\n") 
     
            # Verbosity 2
            if verbosity >= 2:
                fn.write("Outer Loop at last lambda. j_min = iter of inner loop that achieves best zNLP\n")
                fn.write("iter   lambda j_min wtd deaths subopt   t_n\n")
                
                for i in range(iter+1):
                    fn.write(f'{i: ^{2}}  {l[i]: 9.4f}  {j_min[i]: ^2} {z[i]: 8.2f} {z[i] - z_opt: 9.2f} {t_n[i]: 6.2f} ')
                    fn.write("\n")
                fn.write("\n")

                fn.write("Inner Loop at last lambda\n")
                #fn.write("iter     zLP     sim zNLP subopt_zNLP deaths   t_n   vacc by area\n")
                fn.write("iter   zNLP     subopt   deaths donor_deaths tot_deaths t_n   vacc by area\n")
                for i in range(j+1):
                    fn.write(f'{i: ^{2}} {zNLP[i]: 10.2f} {zNLP[i] - zNLP_min: 9.2f} ')
                    fn.write(f'{deaths[i]: 8.2f}  {donor_deaths[i]: 8.2f}  {tot_deaths[i]: 8.2f}  {t_n[i]: 6.2f} ')
                    for a in A:
                        fn.write(f'{V_tot[a, i]: 5.0f} ')                    
                    fn.write("\n")
                fn.write("\n")

            # Verbosity 1
            if verbosity >= 1:
                fn.write("Vaccinations: sim of best LP (best j), last lambda (V_min) \n"
                        "  day    V_opt by area \n")
                for t in range(T - T0 + 1):
                    fn.write(f'{t: ^{7}}')
                    for a in A:
                        fn.write("  " + str(V_min[a, t]) + "  ")
                    fn.write("\n")
                fn.write("\nOptimal vaccinations (V_opt) \n"
                        "  day    V_opt by area \n")
                for t in range(T - T0 + 1):
                    fn.write(f'{t: ^{7}}')
                    for a in A:
                        fn.write("  " + str(V_opt[a, t]) + "  ")
                    fn.write("\n")

    else: # Simulate only
        t_sim, alpha, V_cal, V, D = simulate(V)
    """ Output not finished. wtd deaths, vacc by area from this sim 
        deaths = (1 - nu)*D[donor, T] 
        for a in A:
            deaths += nu*D[a, T]
        donor_deaths = D[donor, T]
        total_deaths = 0
        for a in A:
            total_deaths += D[a, T]
        V_opt = {(a, t): (V[a, t] if t < T else 0) for a in A for t in range(T+1)} """

def optimize_inner(l, V): 
    global deaths, donor_deaths, tot_deaths, t_n, zLP, zNLP, j_min, V_tot, \
        deaths_min, donor_deaths_min, tot_deaths_min, V_min, t_min, alpha_min, V_cal_min, zNLP_min, \
        deaths_prev, donor_deaths_prev, tot_deaths_prev, t_prev, alpha_prev, V_cal_prev, \
        z_opt, z_init_sim, donor_deaths_opt, tot_deaths_opt, V_opt, t_opt, V_cal_opt, lambda_opt, j_opt, i_opt,  \
        deaths_sim, donor_deaths_sim, tot_deaths_sim, t_sim, V_tot_sim, j, LP_count, infeas_count, \
        dVmax, dVmin

    if iter == 0:
        LP_count = infeas_count = 0
        # Define arrays
        j_min = {i: 0 for i in range(iter_lmt_search + 1)} # iter j of inner loop that achieves best zNLP 
        deaths = {j1: 0 for j1 in range(iter_lmt + 1)} # wtd deaths, iteration j of inner loop
        donor_deaths = {j1: 0 for j1 in range(iter_lmt + 1)} # donor deaths ...
        tot_deaths = {j1: 0 for j1 in range(iter_lmt + 1)} # total deaths
        t_n = {j1: 0 for j1 in range(iter_lmt + 1)} # t_n
        zLP = {j1: 0 for j1 in range(iter_lmt + 1)} # LP objective fcn value
        zNLP = {j1: 0 for j1 in range(iter_lmt + 1)} # NLP objective from sim
        V_tot = {(a, j1): 0 for a in A for j1 in range(iter_lmt + 1)} # total vacc by area, j

        # First simulate. Record as j = 0
        t_n[0], alpha, V_cal, V, D = simulate(V)
        # Compute wtd deaths, vacc by area from this sim 
        deaths[0] = (1 - nu)*D[donor, T]
        donor_deaths[0] = D[donor, T]
        tot_deaths[0] = 0
        for a in A:
            deaths[0] += nu*D[a, T]
            tot_deaths[0] += D[a, T]
            for t in range(T):
                V_tot[a, 0] += V[a, t] 
        # Record first sim for output
        V_tot_sim = {a: 0 for a in A}
        for a in A:
            V_tot_sim[a] = V_tot[a, 0]
        deaths_sim = deaths[0]
        donor_deaths_sim = donor_deaths[0]
        tot_deaths_sim = tot_deaths[0]
        t_sim = t_n[0]
        # Initialize min solution from this sim (needed for solve_LP)
        alpha_min = alpha
        V_cal_min = V_cal
        # Initialize opt from this sim
        z_init_sim = deaths[0]
        z_opt = deaths[0]
        donor_deaths_opt = donor_deaths[0]
        tot_deaths_opt = tot_deaths[0]
        t_opt = t_n[0]
        lambda_opt = 0 # there is no lambda for init sim
        i_opt = iter
        j_opt = -1 # iter j of inner loop that achieves z_opt
        V_opt = {(a, t): (V[a, t] if t < T else 0) for a in A for t in range(T+1)}
        V_cal_opt = {(a, t): V_cal[a, t] for a in A for t in range(T+1)}
                                            # t=0,...,T b/c used in output

    elif iter > 0 and phase == 1:  # Update j = 0 solution using best solution from iter 0 
        # All phase 1 iter use iter 0 as their init solution to make them more comparable 
        deaths[0] = deaths_prev
        donor_deaths[0] = donor_deaths_prev
        tot_deaths[0] = tot_deaths_prev
        zLP[0] = 0
        zNLP[0] = 0
        t_n[0] = t_prev
        alpha = alpha_prev
        V_cal = V_cal_prev
        for a in A:
            V_tot[a, 0] = V_tot[a, j_min[0]] # for output

    else:   # Phase 2: Update j = 0 solution using best solution from last lambda (alpha, V_cal needed for solve_lp)
        deaths[0] = deaths_min
        donor_deaths[0] = donor_deaths_min
        tot_deaths[0] = tot_deaths_min
        zLP[0] = 0
        zNLP[0] = 0
        t_n[0] = t_min
        alpha = alpha_min
        V_cal = V_cal_min
        for a in A:
            V_tot[a, 0] = V_tot[a, j_min[iter-1]] # for output
 
    # Initialize for j loop
    t_LP = min(math.ceil(t_n[0]) + dT, T) # Use fixed t_n in LP for all j, but include dT more days 
        # since t_n may increase. Use t_n from first sim or best from last iter
    j = 0 # inner loop counter
    eps = eps_prev = epsilon_0 # eps_prev is eps at the previous best sim
    zNLP_min = 1e14    # initialize min ZNLP
    zNLP_prev = 2*zNLP_min   # initialize previous value of min zNLP. They must differ.
 
#    while (abs(zNLP_prev - zNLP_min) >= delta_I or j < 2) and j < iter_lmt:
    while (abs(zLP[j] - zLP[max(0,j-1)]) >= delta_I or j < 2) and j < iter_lmt:
        j += 1
        LP_count += 1
        # solve LP
        if improving == 1:
            zLP[j] = solve_LP(l, t_LP, alpha_min, V_cal_min, eps)  # Improving: alpha_min, V_cal_min are best for this lambda
        else:
            zLP[j] = solve_LP(l, t_LP, alpha, V_cal, eps)   # Not improving: alpha, V_cal from sim of current LP
        if v.status == GRB.OPTIMAL: # If LP infeas, update eps and iterate. All remaining j likely infeas. 
            V = {(a, t): V1[a, t].x for a in A for t in range(T)} # update V using optimal V1 from LP
            t_n[j], alpha, V_cal, V, D = simulate(V)   # simulate LP solution 
            # compute deaths, total vacc by area, NLP objective from sim
            deaths[j] = (1 - nu)*D[donor, T]
            donor_deaths[j] = D[donor, T]
            tot_deaths[j] = 0  
            for a in A:
                deaths[j] += nu*D[a, T]
                tot_deaths[j] += D[a, T]
                V_tot[a, j] = 0
                for t in range(T):
                    V_tot[a, j] += V[a, t] 
            zNLP[j] = deaths[j] + l*sum(I1[a, t].x for a in A_D for t in range(1, t_LP + 1)) \
                - (1e-7)*sum(V1[a, t].x for a in A for t in range(T)) # A_D: no donor in Lagr

            if zNLP[j] <= zNLP_min:  # update best deaths, obj, vacc, alpha, t_n, V_cal for this lambda
                j_min[iter] = j
                deaths_min = deaths[j]
                donor_deaths_min = donor_deaths[j]
                tot_deaths_min = tot_deaths[j] 
                zNLP_prev = zNLP_min  # save previous value to use in stopping criterion
                zNLP_min = zNLP[j]
                V_min = {(a, t): (V[a, t] if t < T else 0) for a in A for t in range(T+1)}
                alpha_min = alpha
                V_cal_min = V_cal
                t_min = t_n[j]
                eps = eps_prev      # reset eps to its value at the previous best sim objective
                eps_prev *= beta    # reduce eps b/c better sim objective achieved

            if deaths[j] < z_opt:  # update best deaths, t_n, lambda, vacc, V_cal for all lambda
                i_opt = iter
                j_opt = j
                z_opt = deaths[j]
                donor_deaths_opt = donor_deaths[j]
                tot_deaths_opt = tot_deaths[j] 
                t_opt = t_n[j]
                lambda_opt = l
                V_opt = {(a, t): (V[a, t] if t < T else 0) for a in A for t in range(T+1)}
                V_cal_opt = V_cal
                # Compute difference in V_cal, simulate - LP (opt)
                dVmax = -10000
                dVmin = 10000
                for t in range(1, T+1):
                    for a in A:
                        dV = (I[a, t] + p_e*I_V[a, t]) - (I1[a, t].x + p_e*IV1[a, t].x)
                    dVmax = max(dVmax, dV) # max dV so far
                    dVmin = min(dVmin, dV) # min dV so far
        else:
            infeas_count += 1
        eps *= beta              # reduce current eps
    if iter == 0: # Save iter 0 min solution in "prev" for use initializing all phase 1 iter 
                 ## Save iter 1 min solution in "prev" for use initializing iter 3 
        deaths_prev = deaths_min
        donor_deaths_prev = donor_deaths_min
        tot_deaths_prev = tot_deaths_min 
        alpha_prev = alpha_min
        V_cal_prev = V_cal_min
        t_prev = t_min
    return deaths_min

def solve_LP(l, t_LP, alpha, V_cal, eps):
    global S1, SV1, E1, EV1, I1, IV1, D1, R1, W1, V1, v #H1

    # Warm start using bases (can also use solution)
    if LP_count - infeas_count > 1: #if v.status == GRB.OPTIMAL: 
        vbas = {i: v.getVars()[i].VBasis for i in range(v.NumVars)} # Save basis for variables
        ##psol = {i: v.getVars()[i].x for i in range(v.NumVars)}
        cbas = {i: v.getConstrs()[i].CBasis for i in range(v.NumConstrs)} # Save basis for constraints (dual)
        ##dsol = {i: v.getConstrs()[i].Pi for i in range(v.NumConstrs)}

    # Objective changes with lambda (outer loop) and tn (inner loop)
    v.setObjective((1 - nu)*D1[donor, T] + nu*D1.sum('*', T)+ l*sum(I1[a, t] for a in A_D for t in range(1, t_LP + 1))
        - (1e-7)*sum(V1[a, t] for a in A for t in range(T)), GRB.MINIMIZE)  
    
    # Some constraints change with alpha, V_cal (inner loop). Rewrite them all.
    if LP_count - infeas_count > 1:
        v.remove([constraint for constraint in v.getConstrs()]) # Remove all constraints

    v.addConstrs((V1[donor, t] <= p_k*B[t]
                 for t in range(T - T0 + 1)), "Policy_donor_limit")
    v.addConstrs((V1.sum('*', t) <= B[t] for t in range(T - T0 + 1)), "Vaccine_budget")
    v.addConstrs((V1.sum('*', t) == 0 for t in range(T - T0 + 1, T)), "No_vaccine")

    # Dynamics for start time t=1 to T-1
    ##v.addConstrs((W1[a, t+1] == W1[a, t] - alpha[a, t]/N[a]*W1[a, t] - V1[a, t] test??
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
    #v.addConstrs((H1[a, t+1] == H1[a, t] + r_d[a]*p_H*I1[a, t] + r_d[a]*p_V_H*IV1[a, t] - r_R*H1[a, t]
    #              for a in A for t in range(1, T)), "H")
    v.addConstrs((D1[a, t+1] == D1[a, t] + r_d[a]*p_D*I1[a, t] + r_d[a]*p_V_D*IV1[a, t]
                  for a in A for t in range(1, T)), "D")
    v.addConstrs((R1[a, t+1] == R1[a, t] + r_d[a]*(1 - p_D)*I1[a, t] + r_d[a]*(1 - p_V_D)*IV1[a, t]
                  for a in A for t in range(1, T)), "R")
    # Conservation constraints
    #v.addConstrs((S1[a,t]+SV1[a,t]+E1[a,t]+EV1[a,t]+I1[a,t]+IV1[a,t]+R1[a,t]+H1[a,t]+D1[a,t] == N[a]
    #              for a in A for t in range(1, T+1)), "Conserv")

    # Dynamics for start time t=0. Use global constants S0[a],SV0[a],E0[a],EV0[a],I0[a],IV0[a]  
    # constant V_cal[a, 0] but variable V1. Assume H[a, 0] = D[a, 0] = R[a, 0] = 0 (they aren't defined)
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
    #v.addConstrs((H1[a, 1] == r_d[a]*p_H*I0[a] + r_d[a]*p_V_H*IV0[a]
    #              for a in A), "H_t=0")
    v.addConstrs((D1[a, 1] == r_d[a]*p_D*I0[a] + r_d[a]*p_V_D*IV0[a]
                  for a in A), "D_t=0")
    v.addConstrs((R1[a, 1] == r_d[a]*(1 - p_D)*I0[a] + r_d[a]*(1 - p_V_D)*IV0[a]
                  for a in A), "R_t=0")

    # Regularization constraints on V_cal: constant in t
    v.addConstrs((I1[a, t] + p_e*IV1[a, t] <= V_cal[a, t] + eps #*(t/T) ##*(t/T)
                 for a in A for t in range(1, T + 1)), "V_cal_upper_bd")
    v.addConstrs((I1[a, t] + p_e*IV1[a, t] >= max(V_cal[a, t] - eps, 0) #*(t/T) ##*(t/T)
                 for a in A for t in range(1, T + 1)), "V_cal_lower_bd")

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
        #if LP_count - infeas_count == 1: 
        #    v.write("./output/lamda-" + str(l) + "solve" + str(j) + ".lp") ##test: write first LP
        return v.ObjVal
    elif v.status == GRB.INFEASIBLE:
        print("Infeasible LP at lambda =", l, "LP iteration j =", j, "exploration tol eps =", eps)
        return 100000000
    else:
        print("Optimize failed at lambda =", l, "LP iteration j =", j, "exploration tol eps =", eps)
        print("Status =", v.status)
        return 100000000

def simulate(V):
    global I, I_V
    # Define state variables, alpha, delta_E, V_star, V_plus
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
    alpha = {(a, t): 0 for a in A for t in range(T)}    #  t=0,...,T-1 b/c not computed in diff eq
    delta_E = {(a, t): 0 for a in A for t in range(T)}
    V_star = {(a, t): 0 for a in A for t in range(T)}
    V_plus = {(a, t): 0 for a in A for t in range(T)}

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
            delta_E[a, t] = min(S[a, t], alpha[a, t]*S[a, t]*V_cal[a, t]/N[a])
            if S[a, t] < 0.0000001:
                V_star[a, t] = min(W[a, t], V[a, t])
            else:
                V_star[a, t] = min(W[a, t] - W[a, t]*delta_E[a, t]/S[a, t], V[a, t])

        # Realloc: 2 areas. Recompute V_star
        if realloc_flag and n_a == 2:
            # area[0]
            if S[A[0], t] < 0.0000001:
                Wnew = W[A[0], t]
            else:
                Wnew = W[A[0], t] - W[A[0], t]*delta_E[A[0], t]/S[A[0], t]
            V_star[A[0], t] = min(Wnew, V[A[0], t] + V[A[1], t] - V_star[A[1], t]) # add avail realloc from A[1]
            # area[1]
            if S[A[1], t] < 0.0000001:
                Wnew = W[A[1], t]
            else:
                Wnew = W[A[1], t] - W[A[1], t]*delta_E[A[1], t]/S[A[1], t]
            V_star[A[1], t] = min(Wnew, V[A[1], t] + V[A[0], t] - V_star[A[0], t]) # add avail realloc from A[0]

        if realloc_flag and n_a > 2:
            print("Reaalocation for more than 2 areas not implmented")
            quit()  

        # difference equations including W
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
                (S_V[a, t]/N[a]) * \
                V_cal[a, t] - r_I*E_V[a, t]
            I[a, t + 1] = I[a, t] + r_I*E[a, t] - r_d[a]*I[a, t]
            I_V[a, t + 1] = I_V[a, t] + r_I * \
                E_V[a, t] - r_d[a]*I_V[a, t]
            #H[a, t + 1] = H[a, t] + r_d[a]*p_H*I[a, t] + \
            #    r_d[a]*p_V_H*I_V[a, t] - r_R*H[a, t]
            D[a, t + 1] = D[a, t] + r_d[a]*p_D*I[a, t] + r_d[a]*p_V_D*I_V[a, t]
            R[a, t + 1] = R[a, t] + r_d[a]*(1 - p_D)*I[a, t] + r_d[a]*(1 - p_V_D)*I_V[a, t]

        # check for the variant occuring, do not calculate if variant already emerged
        if not variant_emerge:
            I_sum = 0
            for a in A:
                if a != donor: ## donor doesn't contribute to variant
                    for t1 in range(t+1): # t1=0,...,t 
                        I_sum += I[a, t1]    # sum all infected up to current time
            # If this sum > n, variant emerges. Compute t_sim
            if I_sum > n:
                variant_emerge = True
                I_tot = 0
                for a in A:
                    if a != donor:
                        I_tot += I[a, t]
                t_sim = t + 1 - (I_sum - n)/(I_tot)
    # Compute V_cal at T (used in LP and opt output)
    for a in A:
        V_cal[a, T] = I[a, T] + p_e*I_V[a, T]

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

def import_xml(xml_path: str): # Read inputs from XML file. xml_path: path to the XML file
    root = ET.parse(xml_path).getroot()

    # read area data
    area_data = root.find("area_data")
    global A, A_D, N, rho_V, rho_I_N, delta_r, gamma, rho, donor, m, n, n_a
    A = []
    A_D = [] # areas except donor
    N = {}
    rho_V = {}
    rho_I_N = {}
    delta_r = {}
    gamma = {}
    rho = {}
    donor = area_data.find("donor").text
    m = area_data.find("m").text
    n = convert_num(area_data.find("n").text)
    for child in area_data.findall("area"):
        area = child.attrib["name"]
        A.append(area)
        if area != donor:
            A_D.append(area)
        N[area] = convert_num(child.find("N").text)
        rho_V[area] = convert_num(child.find("rho_V").text)
        rho_I_N[area] = convert_num(child.find("rho_I_N").text)
        delta_r[area] = convert_num(child.find("delta_r").text)
        gamma[area] = convert_num(child.find("gamma").text)
        rho[area] = convert_num(child.find("rho").text)
    n_a = len(A)
    # read scenario data
    scenario_data = root.find("scenario_data")
    global r_I, r_0, p_D, p_V_D, a_0, delta_a, p_e, p_r, \
        L, T_D, p, T, B_0, b_arr  #p_H, p_V_H,
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
    T = convert_num(scenario_data.find("T").text)
    B_0 = convert_num(scenario_data.find("B_0").text)
    b_arr = scenario_data.find("b").text
    if not b_arr == None:
        b_arr = b_arr.split(sep=",")
        for i in range(len(b_arr)):
            b_arr[i] = convert_num(b_arr[i])
    else:
        b_arr = []

    # read params
    params = root.find("params")
    global verbosity, simulate_only, improving, realloc_flag, nu, T0, dT, \
        t_switch, p_k, lambda_0, phi, epsilon_0, delta_I, delta, beta, \
        iter_lmt, iter_lmt_search
    verbosity = convert_num(params.find("verbosity").text)
    simulate_only = bool(convert_num(params.find("simulate_only").text))
    improving = bool(convert_num(params.find("improving").text))
    realloc_flag = bool(convert_num(params.find("realloc_flag").text))
    nu = convert_num(params.find("nu").text)
    T0 = convert_num(params.find("T0").text)
    dT = convert_num(params.find("dT").text)
    t_switch = convert_num(params.find("t_switch").text)
    p_k = convert_num(params.find("p_k").text)
    lambda_0 = convert_num(params.find("lambda_0").text)
    phi = convert_num(params.find("phi").text)
    epsilon_0 = convert_num(params.find("epsilon_0").text)
    delta_I = convert_num(params.find("delta_I").text)
    delta = convert_num(params.find("delta").text)
    beta = convert_num(params.find("beta").text)
    iter_lmt = convert_num(params.find("iter_lmt").text)
    iter_lmt_search = convert_num(params.find("iter_lmt_search").text)

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
    else:
        return int(num)

def output_state_equations(fn: TextIOWrapper,
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

# Script begins here. Prompt for input file name.
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # First positional argument (this must be present)
    parser.add_argument('input', type=str, help='Directory of input xml file')

    # Parse the command line
    args = parser.parse_args()

    input_dir = args.input
    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir,f))]
    #print(files)
    
    for f in files:
        global input_file
        global input_filename
        input_filename = f
        input_file = os.path.join(input_dir,input_filename)
        #print(input_file)
        #print(input_filename)
        main()