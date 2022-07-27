# Type imports
from io import TextIOWrapper

# Packages
import numpy as np
import time
import os
import math
import gurobipy as gp
import xml.etree.ElementTree as ET
import argparse


def main(xml_path: str):
    import_xml(xml_path=xml_path)
    simulate()


# def formulate_LP():          
#     # key inputs: (I, IV, alpha) for all areas and times, tn, m, lambda[i], eps
#     # key outputs:  LP model
    
#     t_int = math.ceil(t_n)

#     # Compute alpha[a,t] here if not done in simulate()?? 
#       # Use constant rate before tn or tn + L and variable rate after

#     v = gp.Model("vaccine_opt")
#     # v.setParam("NumericFocus",1)  # Increased protection from roundoff/instability. Shouldn't need.

#     # Define LP variables, e.g., S1 for state var S. All are continuous and nonnegative by default. 

#     S1 = v.addVars(A, range(1,T + 1), name="S1")  # range(1,T+1): 1, ..., T
#     SV1 = v.addVars(A, range(1,T + 1), name="SV1")
#     E1 = v.addVars(A, range(1,T + 1), name="E1")
#     EV1 = v.addVars(A, range(1,T + 1), name="EV1")
#     I1 = v.addVars(A, range(1,T + 1), name="I1")
#     IV1 = v.addVars(A, range(1,T + 1), name="IV1")
#     H1 = v.addVars(A, range(1,T + 1), name="H1")
#     D1 = v.addVars(A, range(1,T + 1), name="D1")
#     R1 = v.addVars(A, range(1,T + 1), name="R1")
#     W1 = v.addVars(A, range(1,T + 1), name="W1")
#     V1 = v.addVars(A, range(0, T), name="V1")             # t = 0, ..., T-1


#     if non_donor_deaths: #or FALSE
#         v.setObjective(D1[donor, T] + l[i]*sum(I1[a,t] for a in areas for t in range(1, t_int + 1) ) \
#             + 1/200000*D1.sum('*', T), gp.GRB.MINIMIZE)  # include non-donor deaths w/ small weight
#                                                     # Omit constant I[a,0] from objective, so t = 1, ..., t_int
#     else:
#             v.setObjective(D1[donor, T] + lambda[i]*sum(I1[a, t] for a in areas for t in range(1, t_int + 1) ), \
#                  GRB.MINIMIZE)
                 
#     v.addConstrs((V1[donor, t] <= pk*B[t] for t in range(T)), "Policy: donor limit") # constraint must be in () if name argument used
#     v.addConstrs((V1.sum('*', t) <= B[t] for t in range(T)), "Vaccine budget")                  

#     # Dynamics for start time t=1 to T-1
#     v.addConstrs((W1[a,t+1] == W1[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*W1[a,t] - V1[a,t] \
#                         for a in areas for t in range(1,T)), "Vaccine willingness")
#     v.addConstrs((S1[a,t+1] == S1[a,t] - alpha[a,t]*V_cal[a,t]/N[a]*S1[a,t] - V1[a,t] \
#                         for a in areas for t in range(1,T)), "S")                    
#     v.addConstrs((SV1[a,t+1] == SV1[a,t] - pr*alpha[a,t]*V_cal[a,t]/N[a]*SV1[a,t] + V1[a,t] \
#                         for a in areas for t in range(1,T)), "SV")
#     v.addConstrs((E1[a,t+1] == E1[a,t] + alpha[a,t]*V_cal[a,t]/N[a]*S[a,t] - rI*E1[a,t] \
#                         for a in areas for t in range(1,T)), "E")  
#     v.addConstrs((EV1[a,t+1] == EV1[a,t] + pr*alpha[a,t]*V_cal[a,t]/N[a]*SV1[a,t] - rI*EV1[a,t] \
#                         for a in areas for t in range(1,T)), "EV")
#     v.addConstrs((I1[a,t+1] == I1[a,t] + rI*E1[a,t] - rd*I1[a,t] \
#                         for a in areas for t in range(1,T)), "I")
#     v.addConstrs((IV1[a,t+1] == IV1[a,t] + rI*EV1[a,t] - rd*IV1[a,t] \
#                         for a in areas for t in range(1,T)), "IV")
#     v.addConstrs((H1[a,t+1] == H1[a,t] + rd*pH*I1[a,t] + rd*pHV*IV1[a,t] - rR*H1[a,t] \ 
#                         for a in areas for t in range(1,T)), "H")
#     v.addConstrs((D1[a,t+1] == D1[a,t] + rR*pD*H1[a,t] \
#                         for a in areas for t in range(1,T)), "D")
#     v.addConstrs((R1[a,t+1] == R1[a,t] + rR*(1 - pD)*H1[a,t] + rd*(1 - pH)*I1[a,t] + rd*(1 - pHV)*IV1[a,t] \
#                         for a in areas for t in range(1,T)), "R")

#     # Dynamics for start time t=0. Use constant W[a,0} etc. (but variable V1). Use same constraint names??
#     v.addConstrs((W1[a,1] == W[a,0] - alpha[a,0]*V_cal[a,0]/N[a]*W[a,0] - V1[a,0] \
#                         for a in areas), "Vaccine willingness t=0")
#     v.addConstrs((S1[a,1] == S[a,0] - alpha[a,0]*V_cal[a,0]/N[a]*S[a,0] - V1[a,0] \
#                         for a in areas), "S t=0")                   
#     v.addConstrs((SV1[a,1] == SV[a,0] - pr*alpha[a,0]*V_cal[a,0]/N[a]*SV[a,0] + V1[a,0] \
#                         for a in areas), "SV t=0")
#     v.addConstrs((E1[a,1] == E[a,0] + alpha[a,t]*V_cal[a,0]/N[a]*S[a,0] - rI*E[a,0] \
#                         for a in areas), "E t=0")  
#     v.addConstrs((EV1[a,1] == EV[a,0] + pr*alpha[a,0]*V_cal[a,0]/N[a]*SV[a,0] - rI*EV[a,0] \
#                         for a in areas), "EV t=0")
#     v.addConstrs((I1[a,1] == I[a,0] + rI*E[a,0] - rd*I[a,0] \
#                         for a in areas), "I t=0")
#     v.addConstrs((IV1[a,1] == IV[a,0] + rI*EV[a,0] - rd*IV[a,0] \
#                         for a in areas), "IV t=0")
#     v.addConstrs((H1[a,1] == H[a,0] + rd*pH*I[a,0] + rd*pHV*IV[a,0] - rR*H[a,0] \ 
#                         for a in areas), "H t=0")
#     v.addConstrs((D1[a,1] == D[a,0] + rR*pD*H[a,0] \
#                         for a in areas), "D t=0")
#     v.addConstrs((R1[a,1] == R[a,0] + rR*(1 - pD)*H[a,0] + rd*(1 - pH)*I[a,0] + rd*(1 - pHV)*IV[a,0] \
#                         for a in areas), "R t=0")

#     # Regularization constraints on I, IV
#     v.addConstrs((I1[a,t] - I[a,t] <= eps for a in areas for t in range(1,T+1)), "I upper bd")
#     v.addConstrs((I1[a,t] - I[a,t] >= -eps for a in areas for t in range(1,T+1)), "I lower bd")
#     v.addConstrs((IV1[a,t] - IV[a,t] <= eps for a in areas for t in range(1,T+1)), "IV upper bd")
#     v.addConstrs((IV1[a,t] - IV[a,t] >= -eps for a in areas for t in range(1,T+1)), "IV lower bd")

def simulate():
    '''
    Preforms and SEIR simulation based on data imported from areas, scenario, 
    and run_params python scripts
    '''

    # Define state equations for each area and each time in our time horizon
    S = {(area, t): 0 for area in A for t in range(0, T)}   # Susceptible
    S_V = {(area, t): 0 for area in A for t in range(0, T)} # Susceptible Vaccinated
    E = {(area, t): 0 for area in A for t in range(0, T)}   # Exposed
    E_V = {(area, t): 0 for area in A for t in range(0, T)} # Exposed Vaccinated
    I = {(area, t): 0 for area in A for t in range(0, T)}   # Infected
    I_V = {(area, t): 0 for area in A for t in range(0, T)} # Infected Vaccinated
    R = {(area, t): 0 for area in A for t in range(0, T)}   # Recovered
    H = {(area, t): 0 for area in A for t in range(0, T)}   # Hospitalized
    D = {(area, t): 0 for area in A for t in range(0, T)}   # Dead
    W = {(area, t): 0 for area in A for t in range(0, T)}   # Willing to be vaccinated
    V_calligraphy = {(area, t): 0 for area in A for t in range(0, T)}   # Unvaccinated infectious people
    alpha = {(area, t): 0 for area in A for t in range(0, T)}   # alpha for a given t
    alpha_0 = {(area): a_0 * gamma[area] for area in A}       # Initial alpha
    delta_E = {(area, t): 0 for area in A for t in range(0, T)}
    V = {(area, t): 0 for area in A for t in range(0, T)}
    V_star = {(area, t): 0 for area in A for t in range(0, T)}
    V_plus = {(area, t): 0 for area in A for t in range(0, T)}
    N_dot = {(t): 0 for t in range(0, T)}
    r_d = {(area): r_0 + delta_r[area] for area in A}
    t_n = -1 # Initialize t_n so it can be used in output


    # Set up initial values for the state equations
    for area in A:
        E[area, 0] = ((1 - rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/r_I)
        E_V[area, 0] = ((p_r*rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/r_I)
        I[area, 0] = ((1 - rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/(r_0 + delta_r[area]))
        I_V[area, 0] = ((p_r*rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/(r_0 + delta_r[area]))
        S_V[area, 0] = rho_V[area]*N[area] - E_V[area, 0] - I_V[area, 0]
        S[area, 0] = N[area] - E[area, 0] - E_V[area, 0] - I[area, 0] - I_V[area, 0] - S_V[area, 0]
        W[area, 0] = rho[area]*N[area] - S_V[area, 0] - E_V[area, 0] - I_V[area, 0] - rho[area]*E[area, 0] - rho[area]*I[area, 0]


    if len(b) == 0:
        B = {(t): B_0 * 1 for t in range(0, T)}
    else:
        B = {(t): B_0 * b[t] for t in range(0, T)}
    
    S_dot = 0
    for area in A:
        if area != donor:
            S_dot += S[area, 0]       # Sum of nondonor S at time 0
    
    for area in A:
        for t in range(0, T):
            if area != donor:
                V[area, t] = (1 - p_k) * (S[area, 0]/S_dot) * B[t]
            else:
                V[area, t] = p_k * B[t]

    variant_emerge = False

    # loop over time
    for t in range(0, T):
        # alpha calculated (only if variant found, equations 4 and 5)
        if variant_emerge:
            a_t = a_0 + delta_a/(1 + math.exp(-k*(t - (t_n + T_D))))
            alpha[m, t] = a_t * gamma[m]
            for area in A:
                if area != m:
                    if t - L < 0:
                        alpha[area, t] = alpha_0[area]
                    else:
                        alpha[area, t] = a_0 + delta_a/(1 + math.exp(-k*(t - L - (t_n + T_D))))
        else:
            for area in A:
                alpha[area, t] = alpha_0[area]

        # must compute donor V and V_star before doing rest of areas
        if realloc_flag:
            # compute V calligraphy in order to use the donor delta e before inner for loop
            V_calligraphy[donor, t] = I[donor, t] + p_e*I_V[donor, t]
            # compute delta e (equation 7)
            delta_E[donor, t] = min(S[donor, t], alpha[donor, t]*S[donor, t]*V_calligraphy[donor, t]/N[donor])
            # v star (equation 8)
            V_star[donor, t] = min(W[donor, t] - W[donor, t]*delta_E[donor, t]/S[donor, t], V[donor, t])
            # compute equation 13
            N_dot = 0
            for a in A:
                if a != donor and W[area, t + 1] > 0:
                    N_dot += N[a]

        # loop over areas
        for area in A:
            # compute v calligraphy
            V_calligraphy[area, t] = I[area, t] + p_e*I_V[area, t]
            # compute delta e (equation 7)
            delta_E[area, t] = min(S[area, t], alpha[area, t]*S[area, t]*V_calligraphy[area, t]/N[area])
            
            if S[area, t] < 0.0001:
                # v star (equation 8)
                V_star[area, t] = min(W[area, t], V[area, t])
                # compute w (equation 9)
                W[area, t + 1] = W[area, t] - V_star[area, t]
            else:
                # v star (equation 8)
                V_star[area, t] = min(W[area, t] - W[area, t]*delta_E[area, t]/S[area, t], V[area, t])
                # compute w (equation 9
                W[area, t + 1] = W[area, t] - W[area, t]*delta_E[area, t]/S[area, t] - V_star[area, t]
            
            
            

            # if realloc flag is true
            if realloc_flag:    
                if area != donor and W[area, t + 1] > 0:
                    # compute equation 14
                    V_plus[area, t] = V[area, t] + (N[area]/N_dot) * (V[donor, t] - V_star[donor, t])
                else:
                    V_plus[area, t] = V[area, t]
                # update v star with v (equation 8)
                V_star[area, t] = min(W[area, t] - W[area, t]*delta_E[area, t]/S[area, t], V_plus[area, t])
                # update w with new v (equation 9)
                W[area, t + 1] = W[area, t] - W[area, t]*delta_E[area, t]/S[area, t] - V_star[area, t]
            
            # difference equations
            S[area, t + 1] = S[area, t] - delta_E[area, t] - V_star[area, t]
            S_V[area, t + 1] = S_V[area, t] + V_star[area, t] - p_r*alpha[area, t]*(S_V[area, t]/N[area])*V_calligraphy[area, t]
            E[area, t + 1] = E[area, t] + delta_E[area, t] - r_I*E[area, t]
            E_V[area, t + 1] = E_V[area, t] + p_r*alpha[area, t]*(S_V[area, t]/N[area])*V_calligraphy[area, t] - r_I*E_V[area, t]
            I[area, t + 1] = I[area, t] + r_I*E[area, t] - r_d[area]*I[area, t]
            I_V[area, t + 1] = I_V[area, t] + r_I*E_V[area, t] - r_d[area]*I_V[area, t]
            H[area, t + 1] = H[area, t] + r_d[area]*p_H*I[area, t] + r_d[area]*p_V_H*I_V[area, t] - r_R*H[area, t]
            D[area, t + 1] = D[area, t] + r_R*p_D*H[area, t]
            R[area, t + 1] = R[area, t] + r_R*(1 - p_D)*H[area, t] + r_d[area]*(1 - p_H)*I[area, t] + r_d[area]*(1 - p_V_H)*I_V[area, t]


        # equation 10 (check for the variant occuring, do not calculate if variant already emerged)
        if not variant_emerge:
            I_sum = 0
            # looping over area
            for area in A:
                # looping from 0 to the current time
                for t_temp in range(0, t):
                    # sum all infected up to current time
                    I_sum += I[area, t_temp]
            # compare sum from above loop to person days (parameter n)
            if I_sum > n:
                variant_emerge = True
                I_tot = 0
                for area in A:
                    I_tot += I[area, t]
                # if true calculate equation 11 (t_n)
                t_n = t - 1 + (I_sum - n)/(I_tot)

    with open("output.log", "w") as output_file:
        output_file.writelines("Simulation \n")
        output_file.writelines("Time Horizon: " + str(T) + "\n")

        output_file.writelines("Donor Deaths: " + str(D[donor, T]) + "\n")

        total_deaths = 0
        for area in A:
            total_deaths += D[area, T]
        output_file.writelines("Total Deaths: " + str(total_deaths) + "\n")

        total_vaccinations = 0
        for area in A: 
            for t in range(0, T):
                total_vaccinations += V_star[area, t]
        output_file.writelines("Total Vaccinations: " + str(total_vaccinations) + "\n")

        donor_vaccinations = 0
        for t in range(0, T):
            donor_vaccinations += V_star[donor, t]
        output_file.writelines("Donor Vaccinations: " + str(donor_vaccinations) + "\n")
        
        output_file.writelines("Variant Area: " + m + "\n")

        if t_n < 0:
            output_file.writelines("Variant did not emerge\n")
        else:
            output_file.writelines("Day of Variant Emergence: " + str(t_n) + "\n")

        output_file.writelines("\n")
        
        if verbosity >= 1:
            output_file.writelines("Vaccination Rates by Day \n")
            for t in range(0, T):
                output_file.writelines("    Day " + str(t) + "\n")
                for area in A:
                    output_file.writelines("        " + area + " " + str(round(V_star[area, t], 4)) + "\n")
                output_file.writelines("\n")
        output_file.writelines("\n\n")


        if verbosity >= 2:
            output_file.writelines("State Variables\n\n")
            lower_limit = T - 1
        
            if verbosity >= 3:
                lower_limit = 0

            population_max = 0
            for area in A:
                if N[area] > population_max:
                    population_max = N[area]
            num_length = max(len(str(population_max)) + 4, 7)

            output_state_equations(output_file, num_length, lower_limit, T, "Alpha", alpha)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Susceptible", S)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Susceptible Vaccinated", S_V)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Exposed", E)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Exposed Vaccinated", E_V)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Infected", I)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Infected Vaccinated", I_V)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Recovered", R)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Hospitalized", H)
            output_file.writelines("\n")
            output_state_equations(output_file, num_length, lower_limit + 1, T + 1, "Dead", D)
            output_file.writelines("\n")

def import_xml(xml_path: str):
    """
    Imports SEIR data from a specified XML file 

    Params:
        xml_path: the path string to where the XML file exists

    Returns:
        None
    """
    root = ET.parse(xml_path).getroot()


    area_data = root.find("area_data")
    global A, N, rho_V, rho_I_N, delta_r, gamma, rho, donor, \
        m, n
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
        area_name = child.attrib["name"]
        A.append(area_name)
        N[area_name] = convert_num(child.find("N").text)
        rho_V[area_name] = convert_num(child.find("rho_V").text)
        rho_I_N[area_name] = convert_num(child.find("rho_I_N").text)
        delta_r[area_name] = convert_num(child.find("delta_r").text)
        gamma[area_name] = convert_num(child.find("gamma").text)
        rho[area_name] = convert_num(child.find("rho").text)

    
    scenario_data = root.find("scenario_data")
    global r_I, r_0, r_R, p_V_H, p_H, p_D, a_0, delta_a, p_e, p_r, \
        L, T_D, p, T, B_0, b

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
    b = scenario_data.find("b").text
    if not b == None:
        b = b.split(sep=",")
        for i in range(0, len(b)):
            b[i] = convert_num(b[i])
    else:
        b = []

    params = root.find("params")
    global verbosity, simulate_only, realloc_flag, non_donor_deaths, \
        p_k, lambda_0, phi, epsilon_0, delta_I, delta, beta, iter_lmt, \
        iter_lmt_search

    verbosity = convert_num(params.find("verbosity").text)
    simulate_only = bool(convert_num(params.find("simulate_only").text))
    realloc_flag = bool(convert_num(params.find("realloc_flag").text))
    non_donor_deaths = bool(convert_num(params.find("non_donor_deaths").text))
    p_k = convert_num(params.find("p_k").text)
    lambda_0 = convert_num(params.find("lambda_0").text)
    phi = convert_num(params.find("phi").text)
    epsilon_0 = convert_num(params.find("epsilon_0").text)
    delta_I = convert_num(params.find("delta_I").text)
    delta = convert_num(params.find("delta").text)
    beta = convert_num(params.find("beta").text)
    iter_lmt = convert_num(params.find("iter_lmt").text)
    iter_lmt_search = convert_num(params.find("iter_lmt_search").text)

    global k
    # Calculate constant from given inputs
    k = math.log((1-p)/p)/T_D     #natural log, ln()

def convert_num(num: str):
    """
    Takes an input number string and coverts it to float or int accordingly

    Params:
        num: the number string to convert
    
    Returns:
        the number equivalent
    """
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
    output_file.writelines(state_name + "\n")
    output_file.writelines(f'{"day": ^{num_length}}')
    for area in A:
        output_file.writelines(f'{area: ^{num_length}}')
    output_file.writelines("\n")
    for t in range(lower_limit, upper_limit):
        output_file.writelines(f'{t: ^{num_length}}')
        for area in A:
            output_file.writelines(f'{round(state[area, t], 2): ^{num_length}}')
        output_file.writelines("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # First positional argument (this must be present)
    parser.add_argument('input', type=str, help='Name of input xml file')

    # Parse the command line
    args = parser.parse_args()

    main(os.getcwd() + "/" + args.input)

