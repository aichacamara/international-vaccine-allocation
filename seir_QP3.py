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
    global S0, SV0, E0, EV0, I0, IV0, W0, S1, SV1, E1, EV1, I1, IV1, H1, D1, R1, W1, V1, v, iter, V_opt, V_sim, \
        t_sim0, z, donor_deaths, tot_deaths, phase, fn_base
    # read input file, compute k and B
    import_xml(xml_path=os.getcwd() + "/" + input_file)

    # initialize output directory
    try:
            os.mkdir(os.getcwd() + "/" + "output")
    except:
        pass 
    fn_base = f'./output/{input_filename.split("/")[-1][0:-4]}_T{T:03d}_nu{nu:3.1f}' # output file name uses inputs
    sys.stdout = open(fn_base + "_con_QP.out", "w") # console output redirected to file

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
        # Initialize QP Variables
        v = gp.Model("vaccine_opt")
        #v.Params.DualReductions = 0 ##test: to see if unbounded or infeasible
        v.Params.NonConvex = 2  # Suppress warning for nonconvex obj/constraints
        v.Params.LogToConsole = 1 # 1: Gurobi output turned on
        v.Params.TimeLimit = 600 # seconds for each QP
        v.Params.MIPGap = 0.02 # MIP gap for each QP
        # QP variables are S1 for state var S, etc. All are continuous and nonnegative by default.
        S1 = v.addVars(A, range(1, T + 1), name="S1")          # t = 1, ..., T
        SV1 = v.addVars(A, range(1, T + 1), name="SV1")
        E1 = v.addVars(A, range(1, T + 1), name="E1")
        EV1 = v.addVars(A, range(1, T + 1), name="EV1")
        I1 = v.addVars(A, range(1, T + 1), name="I1")
        IV1 = v.addVars(A, range(1, T + 1), name="IV1")
        H1 = v.addVars(A, range(1, T + 1), name="H1")
        D1 = v.addVars(A, range(1, T + 1), name="D1")
        R1 = v.addVars(A, range(1, T + 1), name="R1")
        W1 = v.addVars(A, range(1, T + 1), name="W1")
        V1 = v.addVars(A, range(T), name="V1")             # t = 0, ..., T-1

        # Iteration 0. Uses t_n, alpha from first sim, lambda_0. Updates t_n, alpha 
        # but z[0] is NOT used in search b/c initial V may give better z[0] than z[1], giving wrong search direction for lambda
        # Define alpha, t_n, z, l (lambda); initialize l[0], l[1], l[2]
        alpha = {(a, t): 0 for a in A for t in range(0, T)}    #  t=0,...,T-1 b/c not computed in diff eq
        t_n = {i: 0 for i in range(iter_lmt_search + 1)}
        z = {i: 0 for i in range(iter_lmt_search + 1)}  
        l = {i: 0 for i in range(iter_lmt_search + 1)}
        l[0] = lambda_0
        l[1] = lambda_0
        l[2] = lambda_0*phi
        iter = i = 0    # i is reassigned writing outputs, so store number of outer iterations in global "iter"
        phase = 1
        z[i], t_n[i], alpha, V = solve_QP(l[i], t_n[0], alpha, V)

        if iter_lmt_search > 1:
            # Iterations 1 and 2
            iter = i = 1
            z[i], t_n[i], alpha, V = solve_QP(l[i], t_n[i-1], alpha, V)
            fL = z[i]  # f at left end pt

            iter = i = 2
            z[i], t_n[i], alpha, V = solve_QP(l[i], t_n[i-1], alpha, V)
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
                z[i], t_n[i], alpha, V = solve_QP(l[i], t_n[i-1], alpha, V)
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
            z[i], t_n[i], alpha, V = solve_QP(l[i], t_n[i-1], alpha, V)

            fx = z[i]  # f(x)
            iter = i = i + 1
            l[i] = y
            z[i], t_n[i], alpha, V = solve_QP(l[i], t_n[i-1], alpha, V)
            fy = z[i]  # f(y)

            # Phase 2 loop
            while (abs(fx - fy) > delta and i < iter_lmt_search):
                iter = i = i + 1
                if fx > fy:        		    # minimum is in [a,x], so update b
                    b, fb, x, fx = (x, fx, y, fy)
                    y = b - 0.618 * (b - a)
                    l[i] = y
                    z[i], t_n[i], alpha, V = solve_QP(l[i], t_n[i-1], alpha, V)
                    fy = z[i]
                else:	                    # minimum is in [y,b], so update a
                    a, fa, y, fy = (y, fy, x, fx)
                    x = a + 0.618 * (b - a)
                    l[i] = x
                    z[i], t_n[i], alpha, V = solve_QP(l[i], t_n[i-1], alpha, V)
                    fx = z[i]
                if (fy < fx and fx > fb):
                    print("Warning: f is not unimin: y < fx and fx > fb")
                if (fy > fx and fa < fy):
                    print("Warning: f is not unimin: fy > fx and fa < fy")
        print("QP count: ", QP_count, "Infeas count: ", infeas_count)

        # Write the csv (optimize)
        with open(fn_base + "_QP_plot.csv", "w") as csv_file:
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
                     IV0[a], 0, 0, 0, S0[a], V_opt[a, 0], t_opt, L] 
                )
                for t in range(1, T + 1):
                    csv_writer.writerow(
                        [a, t, S1[a, t].x, SV1[a, t].x, E1[a, t].x, EV1[a, t].x, I1[a, t].x,
                            IV1[a, t].x, H1[a, t].x, D1[a, t].x, R1[a, t].x, W1[a, t].x, V_opt[a, t], t_opt, L] 
                        )
 
        # Write output file (optimize)
        with open(fn_base + "_QP.out", "w") as fn:
            # input echo
            fn.write("QP " + input_file + "\n\n")
            fn.write("Time Horizon: " + str(T) + "\n")
            fn.write("Donor in t_n: no\n")
            ## fn.write("Improving search: " + str(improving) + "\n")
            if n_a == 2:
                fn.write("Initial policy priority to " + str(2 - np.sign(t_switch)) + "  t_switch: " + str(t_switch) + "\n")
            fn.write("dT: " + str(dT) + "\n")
            fn.write("nu: " + str(nu) + "\n")
            fn.write("lambda: " + str(lambda_0) + "\n")
            fn.write("eps: " + str(epsilon_0) + "\n")
            fn.write("Converg. phi: " + str(phi) + "  iter_search: " + str(iter_lmt_search) + "\n\n")
            # Verbosity 0
            ##fn.write("Convergence: Min/Max change in V_cal, (sim - LP)\n\n")        
            fn.write("QP results are counted even if the time limit (600s) is reached. Suboptimality is unknown.\n\n")
            # fn.write("QP results are counted if solution has MIPgap\n")
            # fn.write("(from lower bd on opt) < 2% \n\n")
            # first sim
            fn.write( "i     lambda      deaths  tot_deaths t_n    vacc by area\n")
            fn.write("first simulation\n")
            # Compute total vacc by area (first sim)
            V_total = {a: 0 for a in A} 
            for a in A:
                for t in range(T):
                    V_total[a] += V_sim[a, t] 
            fn.write(f'0       0     {deaths_sim: 9.3f} {tot_deaths_sim: 9.3f}  {t_sim0: 6.2f}  ')
            for a in A:
                fn.write(f'{V_total[a]: 5.0f} ')                    
            fn.write("\n")

            fn.write(f'optimal (best deaths found)\n')

            # Compute total vacc by area (opt)
            V_total = {a: 0 for a in A} 
            for a in A:
                for t in range(T):
                    V_total[a] += V_opt[a, t] 
            fn.write(f'{i_opt: ^{2}} {lambda_opt: 9.4f}  ')
            fn.write(f'{z_opt: 9.3f} {tot_deaths_opt: 9.3f}  {t_opt: 6.2f}  ')
            for a in A:
                fn.write(f'{V_total[a]: 5.0f} ')                    
            fn.write("\n\n") 
     
            # Verbosity 2
            if verbosity >= 2:
                fn.write("Outer Loop\n")
                fn.write("iter   lambda    wtd deaths  subopt  t_n\n")
                
                for i in range(iter+1):
                    fn.write(f'{i: ^{2}}  {l[i]: 9.4f}  {z[i]: 9.3f} {z[i] - z_opt: 9.3f} {t_n[i]: 6.2f} ')
                    fn.write("\n")
                fn.write("\n")

            # Verbosity 1
            if verbosity >= 1:
                fn.write("\nOptimal vaccinations (V_opt) \n"
                        "  day    V_opt by area \n")
                for t in range(T - T0 + 1):
                    fn.write(f'{t: ^{7}}')
                    for a in A:
                        fn.write("  " + str(V_opt[a, t]) + "  ")
                    fn.write("\n")

    else: # Simulate only
        t_sim0, alpha, V_sim, D = simulate(V)
    ## Output not finished

def solve_QP(l, t_sim, alpha, V): 
    global deaths, donor_deaths, tot_deaths, zNLP, \
        z_opt, z_init_sim, donor_deaths_opt, tot_deaths_opt, V_opt, t_opt, lambda_opt, i_opt,  \
        deaths_sim, donor_deaths_sim, tot_deaths_sim, t_sim0, V_sim, QP_count, infeas_count

    if iter == 0:
        QP_count = infeas_count = 0
        # First simulate
        t_sim, alpha, V, D = simulate(V)
        # Compute deaths from this sim, record for output 
        deaths_sim = (1 - nu)*D[donor, T]
        donor_deaths_sim = D[donor, T]
        tot_deaths_sim = 0
        for a in A:
            deaths_sim += nu*D[a, T]
            tot_deaths_sim += D[a, T]
        t_sim0 = t_sim # t_sim0 is first sim
        V_sim = V
        # Initialize opt from this sim
        z_opt = deaths_sim
        donor_deaths_opt = donor_deaths_sim
        tot_deaths_opt = tot_deaths_sim
        t_opt = t_sim
        lambda_opt = 0 # there is no lambda for init sim
        i_opt = iter
        V_opt = {(a, t): (V[a, t] if t < T else 0) for a in A for t in range(T+1)}
                                            # t=0,...,T b/c used in output
 
    t_QP = min(math.ceil(t_sim)+dT, T) # include dT more days since t_n may increase 

    # Objective changes with lambda
    v.setObjective((1 - nu)*D1[donor, T] + nu*D1.sum('*', T)+ l*sum(I1[a, t] for a in A_D for t in range(1, t_QP + 1))
        - (1e-7)*sum(V1[a, t] for a in A for t in range(T)), GRB.MINIMIZE)   
    
    # Rewrite constraints
    if QP_count - infeas_count > 1:
        v.remove([constraint for constraint in v.getConstrs()]) # Remove all constraints

    v.addConstrs((V1[donor, t] <= p_k*B[t]
                 for t in range(T - T0 + 1)), "Policy_donor_limit")
    v.addConstrs((V1.sum('*', t) <= B[t] for t in range(T - T0 + 1)), "Vaccine_budget")
    v.addConstrs((V1.sum('*', t) == 0 for t in range(T - T0 + 1, T)), "No_vaccine")

    # Dynamics for start time t=1 to T-1
    v.addConstrs((S1[a, t+1] >= S1[a, t] - alpha[a, t]*(I1[a, t] + p_e*IV1[a, t])/N[a]*S1[a, t] - V1[a, t]
                  for a in A for t in range(1, T)), "S")
    v.addConstrs((SV1[a, t+1] >= SV1[a, t] - p_r*alpha[a, t]*(I1[a, t] + p_e*IV1[a, t])/N[a]*SV1[a, t] + V1[a, t]
                  for a in A for t in range(1, T)), "SV")
    v.addConstrs((E1[a, t+1] >= E1[a, t] + alpha[a, t]*(I1[a, t] + p_e*IV1[a, t])/N[a]*S1[a, t] - r_I*E1[a, t]
                  for a in A for t in range(1, T)), "E")
    v.addConstrs((EV1[a, t+1] >= EV1[a, t] + p_r*alpha[a, t]*(I1[a, t] + p_e*IV1[a, t])/N[a]*SV1[a, t] - r_I*EV1[a, t]
                  for a in A for t in range(1, T)), "EV")
    v.addConstrs((I1[a, t+1] >= I1[a, t] + r_I*E1[a, t] - r_d[a]*I1[a, t]
                  for a in A for t in range(1, T)), "I")
    v.addConstrs((IV1[a, t+1] >= IV1[a, t] + r_I*EV1[a, t] - r_d[a]*IV1[a, t]
                  for a in A for t in range(1, T)), "IV")
    v.addConstrs((H1[a, t+1] >= H1[a, t] + r_d[a]*p_H*I1[a, t] + r_d[a]*p_V_H*IV1[a, t] - r_R*H1[a, t]
                  for a in A for t in range(1, T)), "H")
    v.addConstrs((D1[a, t+1] >= D1[a, t] + r_R*p_D*H1[a, t]
                  for a in A for t in range(1, T)), "D")
    v.addConstrs((R1[a, t+1] == R1[a, t] + r_R*(1 - p_D)*H1[a, t] + r_d[a]*(1 - p_H)*I1[a, t] + r_d[a]*(1 - p_V_H)*IV1[a, t]
                  for a in A for t in range(1, T)), "R")
    # Conservation constraints
    #v.addConstrs((S1[a,t]+SV1[a,t]+E1[a,t]+EV1[a,t]+I1[a,t]+IV1[a,t]+R1[a,t]+H1[a,t]+D1[a,t] == N[a]
    #              for a in A for t in range(1, T+1)), "Conserv")

    # Dynamics for start time t=0. Use global constants S0[a],SV0[a],E0[a],EV0[a],I0[a],IV0[a]  
    # but variable V1. Assume H[a, 0] = D[a, 0] = R[a, 0] = 0 (they aren't defined)
    v.addConstrs((S1[a, 1] >= S0[a] - alpha[a, 0]*(I0[a] + p_e*IV0[a])/N[a]*S0[a] - V1[a, 0]
                  for a in A), "S_t=0")
    v.addConstrs((SV1[a, 1] >= SV0[a] - p_r*alpha[a, 0]*(I0[a] + p_e*IV0[a])/N[a]*SV0[a] + V1[a, 0]
                  for a in A), "SV_t=0")
    v.addConstrs((E1[a, 1] >= E0[a] + alpha[a, 0]*(I0[a] + p_e*IV0[a])/N[a]*S0[a] - r_I*E0[a]
                  for a in A), "E_t=0")
    v.addConstrs((EV1[a, 1] >= EV0[a] + p_r*alpha[a, 0]*(I0[a] + p_e*IV0[a])/N[a]*SV0[a] - r_I*EV0[a]
                  for a in A), "EV_t=0")
    v.addConstrs((I1[a, 1] >= I0[a] + r_I*E0[a] - r_d[a]*I0[a]
                  for a in A), "I_t=0")
    v.addConstrs((IV1[a, 1] >= IV0[a] + r_I*EV0[a] - r_d[a]*IV0[a]
                  for a in A), "IV_t=0")
    v.addConstrs((H1[a, 1] >= r_d[a]*p_H*I0[a] + r_d[a]*p_V_H*IV0[a]
                  for a in A), "H_t=0")
    v.addConstrs((D1[a, 1] >= 0
                  for a in A), "D_t=0")
    v.addConstrs((R1[a, 1] == r_d[a]*(1 - p_H)*I0[a] + r_d[a]*(1 - p_V_H)*IV0[a]
                  for a in A), "R_t=0")

    v.optimize()
    QP_count += 1
    if v.status == GRB.OPTIMAL or v.status == GRB.TIME_LIMIT:
#       if QP_count - infeas_count == 1: 
#           v.write("./optimization_output/lamda-" + str(l) + ".lp") ##test: write first QP

        #  sim vacc from QP
        V = {(a, t): V1[a, t].x for a in A for t in range(0, T)} # V = V1 opt values from QP
        t_sim, alpha, V, D = simulate(V)
        # Compute deaths from this sim 
        deaths = (1 - nu)*D[donor, T]
        donor_deaths = D[donor, T]
        tot_deaths = 0
        for a in A:
            deaths += nu*D[a, T]
            tot_deaths += D[a, T]

        if deaths < z_opt:  # update best deaths, t_n, lambda, vacc
            i_opt = iter
            z_opt = deaths
            donor_deaths_opt = donor_deaths
            tot_deaths_opt = tot_deaths 
            t_opt = t_sim
            lambda_opt = l
            V_opt = {(a, t): (V[a, t] if t < T else 0) for a in A for t in range(T+1)}

        return deaths, t_sim, alpha, V # return sim deaths for opt solution
    elif v.status == GRB.INFEASIBLE:
        infeas_count += 1
        print("Infeasible QP at lambda = ", l)
        return 1e7, t_sim, alpha, V # keep previous solution t_n, alpha, V
    else:
        print("Optimize failed at lambda = ", l, "  Status = ", v.status)
        return 1e7, t_sim, alpha, V # keep previous solution t_n, alpha, V

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
    H = {(a, t): 0 for a in A for t in range(T+1)}
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
            H[a, t + 1] = H[a, t] + r_d[a]*p_H*I[a, t] + \
                r_d[a]*p_V_H*I_V[a, t] - r_R*H[a, t]
            D[a, t + 1] = D[a, t] + r_R*p_D*H[a, t]
            R[a, t + 1] = R[a, t] + r_R*(1 - p_D)*H[a, t] + r_d[a]*(
                1 - p_H)*I[a, t] + r_d[a]*(1 - p_V_H)*I_V[a, t]

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
        with open("plot_" + fn_base + ".csv", "w") as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                ["area", "t", "S", "SV", "E", "EV", "I", "IV", "H", "D", "R", "W", "V", "t_n", "L"])
            csv_writer.writerow(
                [m, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_sim, L])
            for a in A:
                for t in range(T + 1):
                    if t != T:
                        csv_writer.writerow([a, t, S[a, t], S_V[a, t], E[a, t], E_V[a, t], I[a, t],
                                            I_V[a, t], H[a, t], D[a, t], R[a, t], W[a, t], V_star[a, t], t_sim, L])
                    else:
                        csv_writer.writerow([a, t, S[a, t], S_V[a, t], E[a, t], E_V[a, t], I[a, t],
                                             I_V[a, t], H[a, t], D[a, t], R[a, t], W[a, t], 0, t_sim, L])
    return t_sim, alpha, V_star, D

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
    global r_I, r_0, r_R, p_V_H, p_H, p_D, a_0, delta_a, p_e, p_r, \
        L, T_D, p, T, B_0, b_arr
    r_I = convert_num(scenario_data.find("r_I").text)
    r_0 = convert_num(scenario_data.find("r_0").text)
    r_R = convert_num(scenario_data.find("r_R").text)
    p_V_H = convert_num(scenario_data.find("p_V_H").text)
    p_H = convert_num(scenario_data.find("p_H").text)
    p_D = convert_num(scenario_data.find("p_D").text)
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

def output_state_equations(output_file: TextIOWrapper,
                           num_length: int,
                           lower_limit: int,
                           upper_limit: int,
                           state_name: str,
                           state: dict):
    """
    Helper function to generate file output listing state values over area and time in a formatted way

    Params:
        output_file: the file to have the output written to
        num_length: the length of the spacing numbers are centered inside of
        lower_limit: the lower limit of the time bounds
        upper_limit: the upper limit of the time bounds
        state_name: the name of the outputted state values
        state: the object containing the values for area and time

    Returns:
        None
    """
    output_file.write(state_name + "\n")
    output_file.write(f'{"day": ^{num_length}}')
    for a in A:
        output_file.write(f'{a: ^{num_length}}')
    output_file.write("\n")
    for t in range(lower_limit, upper_limit):
        output_file.write(f'{t: ^{num_length}}')
        for a in A:
            output_file.write(
                f'{round(state[a, t], 2): ^{num_length}}')
        output_file.write("\n")

# Script begins here. Prompt for input file name.
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # First positional argument (this must be present)
    parser.add_argument('input', type=str, help='Directory of input xml files')

    # Parse the command line
    args = parser.parse_args()

    input_dir = args.input
    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir,f))]
    print(files)

    for f in files:
        global input_file
        global input_filename
        input_filename = f
        input_file = os.path.join(input_dir,input_filename)
        print(input_file)
        print(input_filename)
        main()