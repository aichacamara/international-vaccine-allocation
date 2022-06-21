import numpy as np
import time
import math

from areas import *
from scenario import *
from run_params import *


start_time = time.time() # Start a timer to see how long the calculations take

k = math.log((1-p)/p)/T_D     #natural log, ln()

initial_pop_states = ["S", "E", "I", "SV", "IV", "EV"]

state_vars = ["S", "E", "I",
              "SV", "EV", "IV",
              "H", "D", "R"]


# Calculates the population
initial_pop_a = {(ar, compartment): rho_V[ar] if "V" in compartment else 1 - rho_V[ar] for ar in A for compartment in initial_pop_states}
initial_pop_a = {(ar, compartment): initial_pop_a[ar, compartment] * 100000 if "S" in compartment else rho_I_N[ar]*initial_pop_a[ar, compartment] for ar in A for compartment in initial_pop_states}
initial_pop_a = {(ar, compartment): initial_pop_a[ar, compartment] * 5 / 3.5 if "E" in compartment else initial_pop_a[ar, compartment] for ar in A for compartment in initial_pop_states}


# PLACEHOLDER, IT SHOULD BE DEFINED IN AREAS FILE BUT I DON'T HAVE THE DATA FOR POPULATION YET
N = {ar : sum(initial_pop_a[ar, compartment]
                for compartment in initial_pop_states)
       for ar in A}

B_t = {i: 1750 for i in range(T)}

# REMOVE THIS LINE WHEN KNOW HOW USE DELTA R
w = 0.9
t_N = T
# ------------------------------------------

areas = A[:len(A)]
global_population = sum(N.values())
r_d_t = {ar: r_0 + w * (delta_r[ar] / N[ar]) for ar in areas}


#Create a simulate object, with a new object once per iteration.
#The object holds all the values

def simulate_covid(vax_dictionary):
    global t_N, m, D
    print(f"------------------starting simulation----------------")
    print(f"Variant debug info: t_N = {t_N}, most infected area = {m}.")
    state_variables_simu = {(ar, state_var, i): 0 for ar in areas for i in range(T) for state_var in state_vars}
    weighted_I =  {(ar, i): 0 for ar in areas for i in range(T)}  # effective infectious population
    vectors = {(ar, i): 0 for ar in areas for i in range(T)}  # script V in the paper
    infected_sum = 0
    area_infected_sum = {(ar, i): 0 for ar in areas for i in range(T)}
    sus_pop = {}
    for ar in areas:
        for state_var in initial_pop_states:
            state_variables_simu[ar, state_var, 0] = initial_pop_a[ar, state_var]
        weighted_I[ar, 0] = state_variables_simu[ar, "I", 0] + (p_e) * state_variables_simu[
            ar, "IV", 0]  # initial vectors
        area_infected_sum[ar, 0] = weighted_I[ar, 0]
   # alpha = {(m, i): ALPHA_VALUE + lmbda / (1 + np.exp(-k * (i - (t_N + T_D)))) for i in range(T-1)}   # Calculates alpha(t)
    alpha = {(m, i): a_0 for i in range(T-1)}  # initialize alpha for the most infected area.
    for i in range(T-1):
        if i >= t_N:
            alpha[m, i] = a_0 + delta_a / (1 + np.exp(-k * (i - (t_N + T_D))))
    for ar in areas:
        if ar != m:
            for i in range(T-1):
                alpha[ar, i] = alpha[m, max(i - L, 0)]
    for i in range(T - 1):
        #We need to calculate the infected counts for all areas first
        #That means when we calculate the travel params, we have the info for all areas.
        for ar in areas:
           if i != 0:
                weighted_I[ar, i] = state_variables_simu[ar, "I", i] + (p_e) * state_variables_simu[ar, "IV", i]
                area_infected_sum[ar, i] = area_infected_sum[ar, i - 1] + weighted_I[ar, i]
        for ar in areas:
            vectors[ar, i] = weighted_I[ar, i] 
            deltaE = min(state_variables_simu[ar, "S", i], alpha[ar, i] * state_variables_simu[ar, "S", i] * vectors[ar, i]/N[ar])
            deltaSV = min(state_variables_simu[ar, "S", i] - deltaE, vax_dictionary[ar,i])
            deltaEV = min(state_variables_simu[ar, "SV", i] + deltaSV,\
                         p_r * alpha[ar, i] * state_variables_simu[ar, "SV", i] * vectors[ar, i] * (1 / N[ar]))

            state_variables_simu[ar, "S", i + 1] = state_variables_simu[ar, "S", i] - deltaSV - deltaE
            sus_pop[ar, i] = state_variables_simu[ar, "S", i]
            state_variables_simu[ar, "E", i + 1] = state_variables_simu[ar, "E", i] + deltaE - r_I * state_variables_simu[ar, "E", i]
            state_variables_simu[ar, "I", i + 1] = state_variables_simu[ar, "I", i] + r_I * state_variables_simu[
                ar, "E", i] - r_d_t[ar] * state_variables_simu[ar, "I", i]
            state_variables_simu[ar, "SV", i + 1] = state_variables_simu[ar, "SV", i] + deltaSV - deltaEV
            state_variables_simu[ar, "EV", i + 1] = state_variables_simu[ar, "EV", i] + deltaEV - r_I * state_variables_simu[ar, "EV", i]
            state_variables_simu[ar, "IV", i + 1] = state_variables_simu[ar, "IV", i] + r_I * state_variables_simu[
                ar, "EV", i] - r_d_t[ar] * state_variables_simu[ar, "IV", i]
            state_variables_simu[ar, "H", i + 1] = state_variables_simu[ar, "H", i] + r_d_t[ar] * p_H * state_variables_simu[
                ar, "I", i] + r_d_t[ar] * p_V_H * state_variables_simu[ar, "IV", i] - r_R * state_variables_simu[ar, "H", i]
            state_variables_simu[ar, "D", i + 1] = state_variables_simu[ar, "D", i] + r_R * p_D * state_variables_simu[ar, "H", i]
            state_variables_simu[ar, "R", i + 1] = state_variables_simu[ar, "R", i] + r_R * (1 - p_D) * state_variables_simu[ar, "H", i]\
                                                   + r_d_t[ar]*(1-p_H)*state_variables_simu[ar, "I", i]\
                                                   + r_d_t[ar]*(1-p_V_H)*state_variables_simu[ar, "IV", i]
            if infected_sum < n and infected_sum + weighted_I[ar, i] >= n:
                t_N = i
            infected_sum += weighted_I[ar, i]
        if (i % 50 == 0 and i != 0) or i == T-2:
            print(f"Simulation step for iteration {iteration_number} computed up to day {i}.")
    if min(state_variables_simu.values()) < -.001:
        print(f"Warning: In the simulation, some states became negative. The lowest state was {min(state_variables_simu.values())}.\
         You may have put in too many vaccinations.")
        if state_variables_simu[D, 'S', T-1] < -1:
            print(f"There is some evidence for that: in the donor area, {D}, the susceptible population ended with {state_variables_simu[D, 'S', T-1]}.")
        for ar in areas:
            print(f"The final sus pop in {ar} is {state_variables_simu[ar, 'S', T - 1]}.")
    min_compartment = min(state_variables_simu.values())
    print(f"min value is {min_compartment}, max value is {max(state_variables_simu.values())}")
    print(f"There were {sum(vax_dictionary.values())} (simulation) vaccinations administered.")
    printed_values = 0
    for key, value in state_variables_simu.items():
        if value < -10:
            print(key, value)
            printed_values += 1
        if printed_values > 10:
            break
    if t_N != T: #If t_N is reached, and a variant is triggered...
        for ar in areas:
            if area_infected_sum[ar, t_N] == max([area_infected_sum[l, t_N] for l in areas]):
                m = ar
        print(f"The variant is introduced on day {t_N}")
        print(f"The most infected area at day {t_N} is {m}.")
        print("\n")
    else: #Otherwise,
        print("No variant was introduced.\n")
    infected_dictionary = {(ar, "I", i): state_variables_simu[ar, "I", i] for i in range(T) for ar in areas}
    infected_V_dictionary = {(ar, "IV", i): state_variables_simu[ar, "IV", i] for i in range(T) for ar in areas}
    total_deaths = sum([state_variables_simu[ar, "D", T - 1] for ar in areas])
    print(f"The total number of deaths in this simulation was {total_deaths}.")
    print(f"The total number of deaths in the donor country was {state_variables_simu[D, 'D', T - 1]}.")
    return infected_dictionary, infected_V_dictionary, state_variables_simu[D, 'D', T - 1], alpha, t_N, m, sus_pop, min_compartment


iteration_number = 0
vax_allocation = {(iteration_number, ar, i): 0 for i in range(T-1) for ar in areas}
infections = {(iteration_number, ar, inf, i): 0 for i in range(T) for ar in areas for inf in ["I", "IV"]}
obj_function_values = {0: float("inf")}
print("\nInitialized problem.")
iteration_number += 1  # it's just 1

def intensive_simulation_preparing_vaccination(vax_dictionary):
    global T, areas, N, B_t, global_population
    infected_I, infected_IV, obj, alpha, t_N, m, sus_pop, min_compartment = simulate_covid(
        vax_dictionary=vax_dictionary)  # this run helps us find vax_dict for the next run
    for i in range(T-1):
        for ar in areas:
            if sus_pop[ar, i] < 1:
                vax_dictionary[ar, i] = 0
            else:
                vax_dictionary[ar, i] = min(B_t[i] * (N[ar] / global_population), sus_pop[ar, i])
                vax_dictionary[ar, i] = max(vax_dictionary[ar, i], 0) #removing this changes the output, somehow
    return vax_dictionary, min_compartment


# Get feasible solution set up
vax_dictionary = {(ar, i): B_t[i] * (N[ar] / global_population) for ar in areas for i in range(T-1)}
intensive_computation = 0
while intensive_computation < 3:
    print(f"Preparing to simulate for run number {intensive_computation+1} of 3.")
    vax_dictionary, min_compartment = intensive_simulation_preparing_vaccination(vax_dictionary)
    print(f"The minimum was {min_compartment}.")
    intensive_computation += 1
print(f"The minimum value is {min_compartment} (hopefully positive or near 0) so we will now do two more then optimize.")
infected_I, infected_IV, first_obj, alpha, t_N, m, sus_pop, min_compartment = simulate_covid(vax_dictionary=vax_dictionary)
for ar in areas:
    for i in range(T-1):
        vax_allocation[iteration_number, ar, i] = vax_dictionary[ar, i]
        infections[iteration_number, ar, "I", i] = infected_I[ar, "I", i]
        infections[iteration_number, ar, "IV", i] = infected_IV[ar, "IV", i]
    infections[iteration_number, ar, "I", T-1] = infected_I[ar, "I", T-1]
    infections[iteration_number, ar, "IV", T-1] = infected_IV[ar, "IV", T-1]
vaccine_allocation_totals = {iteration_number: sum(vax_dictionary.values())}
obj_function_values[iteration_number] = first_obj
print(f"\nFound feasible solution with donor deaths {first_obj}.\n")