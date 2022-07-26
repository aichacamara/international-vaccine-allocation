# import numpy as np
# import time
from io import TextIOWrapper
import math


from areasT1 import *
from scenarioT1_1 import *
from paramsT1 import *


# Calculate constant from given inputs
k = math.log((1-p)/p)/T_D     #natural log, ln()


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
    W = {(area, t): 0 for area in A for t in range(0, T)}
    V_calligraphy = {(area, t): 0 for area in A for t in range(0, T)}
    alpha = {(area, t): 0 for area in A for t in range(0, T)}
    alpha_0 = {(area): a_0 * gamma[area] for area in A}       # Store initial alpha
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


            # total_people = S[area, t] + S_V[area, t] + E[area, t] + E_V[area, t] + I[area, t] + I_V[area, t] + R[area, t] + H[area, t] + D[area, t]
            # outputs.writelines("            Total people on this day: " + str(total_people) + "\n")
            # outputs.writelines("\n")

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



simulate()




