from io import TextIOWrapper
import os
import math
import gurobipy as gp
from gurobipy import GRB
import xml.etree.ElementTree as ET
import argparse
import shutil
import csv

def main():
    global S0, SV0, E0, EV0, I0, IV0, W0, S1, SV1, E1, EV1, I1, IV1, H1, D1, R1, W1, V1, v, \
        I, I_V, z_opt, V_opt, V_cal_opt, lambda_opt, donor_deaths, total_deaths, t_n, t_opt, LP_count, infeas_count
    #global p_k ## test
    # read input file, compute k and B
    import_xml(xml_path=os.getcwd() + "/" + input_file)

    # initialize output files
    try:
        if simulate_only: shutil.rmtree(os.getcwd() + "/" + "simulation_output")
        else: 
            shutil.rmtree(os.getcwd() + "/" + "optimization_output")
            shutil.rmtree(os.getcwd() + "/" + "lp_output") ##test: write LP
    except:
        pass
    if simulate_only:
        try:
            os.mkdir(os.getcwd() + "/" + "simulation_output")
        except:
            pass
    else:
        try:
            os.mkdir(os.getcwd() + "/" + "optimization_output")
            os.mkdir(os.getcwd() + "/" + "lp_output")      ##test: write LP
        except:
            pass

    # Define V, z, l (lambda), t_n; initialize l[1], l[2]
    V = {(a, t): 0 for a in A for t in range(0, T)}       # t=0,...,T-1
    z = {i: 0 for i in range(0, iter_lmt_search + 1)}  
    l = {i: 0 for i in range(0, iter_lmt_search + 1)}  
    l[1] = lambda_0
    l[2] = lambda_0*phi
    t_n = {j1: 0 for j1 in range(0, iter_lmt + 1)} # t_n, iteration j of inner loop

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

    # Initialize V using S_dot = Sum of nondonor S at time 0
    S_dot = 0
    for a in A:
        if a != donor:
            S_dot += S0[a]
    #p_temp = p_k ## test
    #p_k = 0.0 ## test: initial V is "no donor"      
    for a in A:
        for t in range(0, T):
            if a != donor:
                V[a, t] = (1 - p_k) * (S0[a]/S_dot) * B[t]
            else:
                V[a, t] = p_k * B[t]
    #p_k = p_temp ## test

    # First simulate
    t_min, alpha, V_cal, V, D = simulate(V)
    # wtd deaths from this sim  
    deaths = (1 - nu)*D[donor, T] 
    for a in A:
        deaths += nu*D[a, T]
    # update best deaths, lambda, vacc, V_cal
    z_opt = deaths
    t_opt = t_min
    lambda_opt = 0 # there is no lambda for init sim
    donor_deaths = D[donor, T]
    total_deaths = 0
    for a in A:
        total_deaths += D[a, T]
    V_opt = {(a, t): (V[a, t] if t < T else 0) for a in A for t in range(0, T+1)}
    V_cal_opt = {(a, t): V_cal[a, t] for a in A for t in range(0, T+1)}  
                                            # t=0,...,T b/c used in output
    if not simulate_only:
        # Initialize outer loop (over lambda)
        phase = 1
              
        # Initialize LP Variables. LP will be updated (objective, constraints) in solve_LP
        v = gp.Model("vaccine_opt")
        v.Params.DualReductions = 0 ##test: to see if unbounded or infeasible
        v.Params.LogToConsole = 0
        v.Params.LPWarmStart = 2    ## Allows presolve with PStart, which otherwise might prevent presolve. See PStart docum. 
        # LP variables are S1 for state var S, etc. All are continuous and nonnegative by default.
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
        V1 = v.addVars(A, range(0, T), name="V1")             # t = 0, ..., T-1

        # Iteration 0. Uses t_n, alpha, V_cal from first sim, lambda_0. Updates t_n, alpha, V_cal 
        # but z[0] is NOT used in search b/c initial V may give better z[0] than z[1], giving wrong search direction for lambda
        l[0] = lambda_0
        iter = i = 0    # Init iter for use in output
        LP_count = infeas_count = 0
        z[i], t_min, alpha, V_cal = optimize_inner(l[i], t_min, alpha, V_cal)
        if iter_lmt_search > 1:
            # Iterations 1 and 2
            i = 1
            z[i], t_min, alpha, V_cal = optimize_inner(l[i], t_min, alpha, V_cal)
            fL = z[i]  # f at left end pt

            i = 2
            z[i], t_min, alpha, V_cal = optimize_inner(l[i], t_min, alpha, V_cal)
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
                z[i], t_min, alpha, V_cal = optimize_inner(l[i], t_min, alpha, V_cal)
                f3 = z[i]  # f(x3)
                if f3 > f2 or (f3 == f2 and f3 < z[1]): #tested: + delta
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
            i = i + 1
            l[i] = x
            z[i], t_min, alpha, V_cal = optimize_inner(l[i], t_min, alpha, V_cal)

            fx = z[i]  # f(x)
            i = i + 1
            l[i] = y
            z[i], t_min, alpha, V_cal = optimize_inner(l[i], t_min, alpha, V_cal)
            fy = z[i]  # f(y)

            # Phase 2 loop
            iter = i # i is reassigned writing outputs, so store number of outer iterations
            while (abs(fx - fy) > delta and i < iter_lmt_search):
                iter = i = i + 1
                if fx > fy:        		    # minimum is in [a,x], so update b
                    b, fb, x, fx = (x, fx, y, fy)
                    y = b - 0.618 * (b - a)
                    l[i] = y
                    z[i], t_min, alpha, V_cal = optimize_inner(l[i], t_min, alpha, V_cal)
                    fy = z[i]
                else:	                    # minimum is in [y,b], so update a
                    a, fa, y, fy = (y, fy, x, fx)
                    x = a + 0.618 * (b - a)
                    l[i] = x
                    z[i], t_min, alpha, V_cal = optimize_inner(l[i], t_min, alpha, V_cal)
                    fx = z[i]
                if (fy < fx and fx > fb):
                    print("Warning: f is not unimin: y < fx and fx > fb")
                if (fy > fx and fa < fy):
                    print("Warning: f is not unimin: fy > fx and fa < fy")
        print("LP count: ", LP_count, "Infeas count: ", infeas_count)

        # Write the csv (optimize)
        with open("./optimization_output/plot_" + input_file.split("/")[-1][0:-4] + ".csv", "w") as csv_file:
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
                     IV0[a], 0, 0, 0, W0[a], V_opt[a, 0], t_opt, L]
                )
                for t in range(1, T + 1):
                    csv_writer.writerow(
                        [a, t, S1[a, t].x, SV1[a, t].x, E1[a, t].x, EV1[a, t].x, I1[a, t].x,
                            IV1[a, t].x, H1[a, t].x, D1[a, t].x, R1[a, t].x, W1[a, t].x, V_opt[a, t], t_opt, L]
                        )
 
        # Write output file (optimize)
        with open("./optimization_output/output.log", "w") as output_file:
            output_file.write("Optimization\n\n")
            output_file.write("Time Horizon: " + str(T) + "\n")
            output_file.write("Optimal weighted deaths: " +  str(z_opt) + 
                "  Nondonor weight nu = " + str(nu) + "\n")
            output_file.write("Lambda opt: " + str(lambda_opt) + "\n")
            output_file.write("Donor deaths (last LP): " + str(D1[donor, T].x) + "\n")
            output_file.write("Donor deaths (sim w/ best wtd deaths): " + str(donor_deaths) + "\n")
            total_deaths_LP = 0
            for a in A:
                total_deaths_LP += D1[a, T].x
            output_file.write("Total deaths (last LP): " + str(total_deaths_LP) + "\n")
            output_file.write("Total deaths (sim w/ best wtd deaths): " + str(total_deaths) + "\n")
            if t_opt == T:
                output_file.write("Variant did not emerge\n")
            else:
                output_file.write("Day of Variant Emergence: " + str(t_opt) + "\n")
            output_file.write("Optimal vaccinations by area\n")
            total_vaccinations = 0
            for a in A:
                area_vaccinations = 0
                for t in range(0, T):
                    area_vaccinations += V_opt[a, t]
                output_file.write(a + " " + str(area_vaccinations) + "\n")
                total_vaccinations += area_vaccinations
            output_file.write("Optimal total vaccinations: " +
                              str(total_vaccinations) + "\n")
           # compute max/min change in V_cal, simulate - last LP
            dVmax = -1000
            dVmin = 1000
            dV = {(a, t): 0 for a in A for t in range(1, T+1)}     # t=1,...,T
            for t in range(1, T+1):
                for a in A:
                    # Compute difference in V_cal, last simulate - last LP
                    dV[a, t] = (I[a, t] + p_e*I_V[a, t]) - (I1[a, t].x + p_e*IV1[a, t].x)
                    dVmax = max(dVmax, dV[a, t]) # max dV so far
                    dVmin = min(dVmin, dV[a, t]) # min dV so far
            output_file.write("Max change in V_cal, (last sim - last LP): " +
                              str(dVmax) + "\n")
            output_file.write("Min change in V_cal, (last sim - last LP): " +
                              str(dVmin) + "\n")
            output_file.write("\n")
            
            # Verbosity 1
            if verbosity >= 1:
                output_file.write("Outer Loop\n")
                output_file.write(
                    "  iteration    lambda  wtd deaths  suboptimality  \n")
                for i in range(0, iter + 1):
                    output_file.write(f'{i: ^{13}}')
                    output_file.write(f'{round(l[i], 8): ^{10}}')
                    output_file.write(f'{round(z[i], 8): ^{10}}')
                    output_file.write(f'{round(z[i] - z_opt, 8): ^{17}}')
                    output_file.write("\n")
                output_file.write("\n")
                output_file.write("Inner Loop at last lambda\n")
                output_file.write("  iteration    zLP   sim zNLP    subopt from best zNLP   t_n \n")
                for i in range(1, j+1):
                    output_file.write(f'{i: ^{13}}')
                    output_file.write(f'{round(zLP[i], 8): ^{8}}' + "  ")
                    output_file.write(f'{round(zNLP[i], 8): ^{8}}' + "  ")
                    output_file.write(f'{round(zNLP[i] - zNLP_min, 8): ^{17}}' + "  ")
                    output_file.write(f'{round(t_n[i], 8): ^{8}}')                   
                    output_file.write("\n")
                output_file.write("\n")

                # Verbosity 2
                if verbosity >= 2:
                    output_file.write("\nOptimal vaccinations \n"
                        "  day    V_opt by area \n")
                    for t in range(1, T):
                        output_file.write(f'{t: ^{7}}')
                        for a in A:
                            output_file.write("  " + str(V_opt[a, t]) + "  ")
                        output_file.write("\n")

                    output_file.write("Change in V_cal in last LP\n")
                    output_file.write(
                        "  day    (V_cal last simulation) - (V_cal last LP) by area  \n")
                    for t in range(1, T + 1):
                        output_file.write(f'{t: ^{7}}')
                        for a in A:
                            output_file.write("  " + str(dV[a, t]) + "  ")
                        output_file.write("\n")

def optimize_inner(l, t_min, alpha_min, V_cal_min):
    global zNLP, zNLP_min, V_min, z_opt, V_opt, V_cal_opt, \
        lambda_opt, zLP, j, donor_deaths, total_deaths, t_n, t_opt, LP_count, infeas_count
    t_LP = min(math.ceil(t_min)+4, T) # Use fixed t_n in LP for all j, but include 2 more days since
                                    # t_n may increase  
    eps = eps_prev = epsilon_0 # eps_prev is eps at the previous best sim
    j = 0
    zLP = {j1: 0 for j1 in range(0, iter_lmt + 1)} # LP objective fcn value, iteration j
    zNLP = {j1: 0 for j1 in range(0, iter_lmt + 1)} # NLP objective from sim, iteration j
    zNLP_min = 1000000000    # initialize min ZNLP
    zNLP_prev = 2*zNLP_min   # initialize previous value of min zNLP. They must differ.
 
    while (abs(zNLP_prev - zNLP_min) > delta_I or j < 2) and j < iter_lmt:
#   while (abs(zLP[max(0,j)] - zLP[max(0,j-1)]) > delta_I or j < 2) and j < iter_lmt:
        j = j + 1
        LP_count += 1
        zLP[j] = solve_LP(l, t_LP, alpha_min, V_cal_min, eps)   # solve LP
        if v.status == GRB.OPTIMAL: # If LP infeas, update eps and iterate. All remaining j likely infeas. 
            V = {(a, t): V1[a, t].x for a in A for t in range(0, T)} # update V using optimal V1 from LP
            t_n[j], alpha, V_cal, V, D = simulate(V)   # simulate LP solution
            # compute deaths, NLP objective from sim
            deaths = (1 - nu)*D[donor, T] # wtd deaths from this sim  
            for a in A:
                deaths += nu*D[a, T]
            zNLP[j] = deaths + l*sum(I1[a, t].x for a in A for t in range(1, t_LP + 1)) \
                - .0000001*sum(V1[a, t].x for a in A for t in range(0, T))
            if zNLP[j] <= zNLP_min:    # update best deaths, obj, vacc, t_n, V_cal for this lambda
                deaths_min = deaths
                zNLP_prev = zNLP_min  # save previous value to use in stopping criterion
                zNLP_min = zNLP[j]
                V_min = {(a, t): V[a, t] for a in A for t in range(0, T)}
                V_cal_min = {(a, t): V_cal[a, t] for a in A for t in range(0, T+1)}
                t_min = t_n[j]
                eps = eps_prev      # reset eps to its value at the previous best sim objective
                eps_prev *= beta    # reduce eps b/c better sim objective achieved
            if deaths < z_opt:  # update best deaths, lambda, vacc, V_cal (for all lambda)
                z_opt = deaths
                t_opt = t_n[j]
                lambda_opt = l
                donor_deaths = D[donor, T]
                total_deaths = 0
                for a in A:
                    total_deaths += D[a, T]
                V_opt = {(a, t): (V[a, t] if t < T else 0) for a in A for t in range(0, T+1)}
                V_cal_opt = {(a, t): V_cal[a, t] for a in A for t in range(0, T+1)}
        else:
            infeas_count += 1
        eps *= beta              # reduce current eps
    #v.write("./optimization_output/lamda-" + str(l) + "solve" + str(j) + ".lp") ##test: write LP
    return deaths_min, t_min, alpha_min, V_cal_min

def solve_LP(l, t_LP, alpha, V_cal, eps):
    global S1, SV1, E1, EV1, I1, IV1, H1, D1, R1, W1, V1, v

    # Objective changes with lambda (outer loop) and tn (inner loop)
    v.setObjective((1 - nu)*D1[donor, T] + nu*D1.sum('*', T)+ l*sum(I1[a, t] for a in A for t in range(1, t_LP + 1))
        - .0000001*sum(V1[a, t] for a in A for t in range(0, T))
                       , GRB.MINIMIZE)  

    # Some constraints change with alpha, V_cal (inner loop). Rewrite them all.
    if LP_count > 1:
        vbas = {i: v.getVars()[i].VBasis for i in range(v.NumVars)} # Save basis for variables
        #psol = {i: v.getVars()[i].x for i in range(v.NumVars)}
        cbas = {i: v.getConstrs()[i].CBasis for i in range(v.NumConstrs)} # Save basis for constraints (dual)
        #dsol = {i: v.getConstrs()[i].Pi for i in range(v.NumConstrs)}
        v.remove([constraint for constraint in v.getConstrs()]) # Remove all constraints

    v.addConstrs((V1[donor, t] <= p_k*B[t] 
                 for t in range(T)), "Policy_donor_limit")
    v.addConstrs((V1.sum('*', t) <= B[t] for t in range(T)), "Vaccine_budget")

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
    v.addConstrs((H1[a, t+1] == H1[a, t] + r_d[a]*p_H*I1[a, t] + r_d[a]*p_V_H*IV1[a, t] - r_R*H1[a, t]
                  for a in A for t in range(1, T)), "H")
    v.addConstrs((D1[a, t+1] == D1[a, t] + r_R*p_D*H1[a, t]
                  for a in A for t in range(1, T)), "D")
    v.addConstrs((R1[a, t+1] == R1[a, t] + r_R*(1 - p_D)*H1[a, t] + r_d[a]*(1 - p_H)*I1[a, t] + r_d[a]*(1 - p_V_H)*IV1[a, t]
                  for a in A for t in range(1, T)), "R")

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
    v.addConstrs((H1[a, 1] == r_d[a]*p_H*I0[a] + r_d[a]*p_V_H*IV0[a]
                  for a in A), "H_t=0")
    v.addConstrs((D1[a, 1] == 0
                  for a in A), "D_t=0")
    v.addConstrs((R1[a, 1] == r_d[a]*(1 - p_H)*I0[a] + r_d[a]*(1 - p_V_H)*IV0[a]
                  for a in A), "R_t=0")

    # Regularization constraints on V_cal
    v.addConstrs((I1[a, t] + p_e*IV1[a, t] <= V_cal[a, t] + eps*t/T 
                 for a in A for t in range(1, T + 1)), "V_cal_upper_bd")
    v.addConstrs((I1[a, t] + p_e*IV1[a, t] >= max(V_cal[a, t] - eps*t/T, 0) 
                 for a in A for t in range(1, T + 1)), "V_cal_lower_bd")

    if LP_count > 1: # Warm start using VBasis/CBasis
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
    global I, I_V
    # Define state variables, alpha, delta_E, V_star, V_plus
    S = {(a, t): 0 for a in A for t in range(0, T+1)}     # t=0,...,T b/c computed in diff eq
    S_V = {(a, t): 0 for a in A for t in range(0, T+1)}
    E = {(a, t): 0 for a in A for t in range(0, T+1)}   
    E_V = {(a, t): 0 for a in A for t in range(0, T+1)}
    I = {(a, t): 0 for a in A for t in range(0, T+1)}   
    I_V = {(a, t): 0 for a in A for t in range(0, T+1)}
    R = {(a, t): 0 for a in A for t in range(0, T+1)}
    H = {(a, t): 0 for a in A for t in range(0, T+1)}
    D = {(a, t): 0 for a in A for t in range(0, T+1)}
    W = {(a, t): 0 for a in A for t in range(0, T+1)}
    V_cal = {(a, t): 0 for a in A for t in range(0, T+1)}  # t=0,...,T b/c used in output
    alpha = {(a, t): 0 for a in A for t in range(0, T)}    #  t=0,...,T-1 b/c not computed in diff eq
    delta_E = {(a, t): 0 for a in A for t in range(0, T)}
    V_star = {(a, t): 0 for a in A for t in range(0, T)}
    V_plus = {(a, t): 0 for a in A for t in range(0, T)}

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

    for t in range(0, T):   # t = 0,...,T-1
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

        # Compute V_cal, delta_E and N_dot; compute V_star and W (without realloc) 
        N_dot = 0
        for a in A:
            V_cal[a, t] = I[a, t] + p_e*I_V[a, t]
            delta_E[a, t] = min(
                S[a, t], alpha[a, t]*S[a, t]*V_cal[a, t]/N[a])
            if S[a, t] < 0.0001:
                V_star[a, t] = min(W[a, t], V[a, t])
                W[a, t + 1] = W[a, t] - V_star[a, t]
            else:
                V_star[a, t] = min(
                    W[a, t] - W[a, t]*delta_E[a, t]/S[a, t], V[a, t])
                W[a, t + 1] = W[a, t] - W[a, t] * \
                    delta_E[a, t]/S[a, t] - V_star[a, t]
            if a != donor and W[a, t + 1] > 0:
                    N_dot += N[a]

        # Realloc: compute V_plus using N_dot, then recompute V_star and W using V_plus 
        if realloc_flag:   
            for a in A:
                if a != donor and W[a, t + 1] > 0:
                    V_plus[a, t] = V[a, t] + \
                        (N[a]/N_dot) * (V[donor, t] - V_star[donor, t])
                else:
                    V_plus[a, t] = V[a, t]
                # update V_star, W with V_plus
                if S[a, t] < 0.0001:
                    V_star[a, t] = min(W[a, t], V_plus[a, t])
                    W[a, t + 1] = W[a, t] - V_star[a, t]
                else:
                    V_star[a, t] = min(
                        W[a, t] - W[a, t]*delta_E[a, t]/S[a, t], V_plus[a, t])
                    W[a, t + 1] = W[a, t] - W[a, t] * \
                        delta_E[a, t]/S[a, t] - V_star[a, t]

        # difference equations
        for a in A:
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
                for t_temp in range(0, t):
                    I_sum += I[a, t_temp] # sum all infected up to current time
            # If this sum > n, variant emerges. Compute t_sim
            if I_sum > n:
                variant_emerge = True
                I_tot = 0
                for a in A:
                    I_tot += I[a, t]
                t_sim = t - 1 + (I_sum - n)/(I_tot)
    # Compute V_cal at T (used in LP and opt output)
    for a in A:
        V_cal[a, T] = I[a, T] + p_e*I_V[a, T]

    return t_sim, alpha, V_cal, V_star, D

def import_xml(xml_path: str): # Read inputs from XML file. xml_path: path to the XML file
    root = ET.parse(xml_path).getroot()

    # read area data
    area_data = root.find("area_data")
    global A, N, rho_V, rho_I_N, delta_r, gamma, rho, donor, m, n
    A = []
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
        N[area] = convert_num(child.find("N").text)
        rho_V[area] = convert_num(child.find("rho_V").text)
        rho_I_N[area] = convert_num(child.find("rho_I_N").text)
        delta_r[area] = convert_num(child.find("delta_r").text)
        gamma[area] = convert_num(child.find("gamma").text)
        rho[area] = convert_num(child.find("rho").text)

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
        for i in range(0, len(b_arr)):
            b_arr[i] = convert_num(b_arr[i])
    else:
        b_arr = []

    # read params
    params = root.find("params")
    global verbosity, simulate_only, realloc_flag, nu, \
        p_k, lambda_0, phi, epsilon_0, delta_I, delta, beta, iter_lmt, \
        iter_lmt_search
    verbosity = convert_num(params.find("verbosity").text)
    simulate_only = bool(convert_num(params.find("simulate_only").text))
    realloc_flag = bool(convert_num(params.find("realloc_flag").text))
    nu = convert_num(params.find("nu").text)
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
    B = {(t): B_0 * 1 for t in range(0, T)}
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
    parser.add_argument('input', type=str, help='Name of input xml file')

    # Parse the command line
    args = parser.parse_args()

    global input_file
    input_file = args.input
    main()