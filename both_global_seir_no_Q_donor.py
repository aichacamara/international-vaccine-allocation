# anaconda
# cd /D C:\Users\Abraham\miniconda3\envs\snowflakes\Scripts\Global_SEIR\Real Sims
# python iterate_travel_alpha.py
four_or_not = input("Run simulation on just four areas? [y/n]")

if four_or_not == "y":
    four_binary = True
elif four_or_not == "n":
    four_binary = False
else:
    print("Try again, this time writing 'y' or 'n'.")
    exit()
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd
if four_binary:
    from user_input_params import * #T, t_N, alpha, num areas, vax days, b_t, areas
    from four_mid_params import * #only differ in their N
    from four_mid_seir_data import * #most_infected_area, state_vars,
    # m_a, vax_rates, initial pop, N_a, travel_rates, tests_per_day, T_a_l
else:
    from user_input_params_real import * #alpha
    from mid_process_df import *  #B_t, T, t_N
    from mid_params_real import * #all
    from mid_seir_data_real import * #num_areas, areas, same as four_mid_seir_data but no travel_rates
import csv

state_vars = ["S", "E", "I",
              "SV", "EV", "IV",
              "H", "D", "R"]

#Known issues:
#Too many vaccinations in the initial solution isn't feasible
#leads to negative S state. Still converges, but the objective values have to go up before they can go down.
#too high alpha leads to S state getting negative. "Infeasible" model.
#Shouldn't be an issue as pandemic shouldn't go through whole population anyway.
global_population = sum(N_a.values())
policy_params = True
if policy_params:
    f_max_a = {(ar, i): 10*B_t[i]/global_population for i in range(T) for ar in areas}
#f_max_a = {(ar, i): 1 for i in range(T) for ar in areas}
if four_binary:
    initial_pop_states = ["S", "E", "I", "SV"]
    donor = "area1"
else:
    initial_pop_states = ["S", "IV", "I", "SV"]
    donor = "USA"
    initial_pop_a = {}
    for state_var in initial_pop_states:
        for ar in areas:
            initial_pop_a[ar, state_var] = vaccine_params_df[state_var][area_lookup[ar]]

def simulate_covid(vax_dictionary):
    global t_N, most_infected_area, four_binary
    print(f"Variant debug info: t_N = {t_N}, most infected area = {most_infected_area}.")
    alpha = {}
    state_variables_simu = {}
    weighted_I = {}
    vectors = {}
    infected_sum = 0
    area_infected_sum = {}
    sus_pop = {}
    being_exposed = {}
    for ar in areas:
        for i in range(T):
            for state_var in state_vars:
                state_variables_simu[ar, state_var, i] = 0  # dictionary set to 0
            vectors[ar, i] = 0  # dictionary set to 0
            weighted_I[ar, i] = 0
            area_infected_sum[ar, i] = 0


    for state_var in initial_pop_states:
        for ar in areas:
            state_variables_simu[ar, state_var, 0] = initial_pop_a[ar, state_var]
    for ar in areas:
        weighted_I[ar, 0] = state_variables_simu[ar, "I", 0] + (1 - p_e) * state_variables_simu[
            ar, "IV", 0]  # initial vectors
        area_infected_sum[ar, 0] = weighted_I[ar, 0]
    for i in range(T):
        alpha[most_infected_area, i] = ALPHA_VALUE + lmbda / (1 + np.exp(-k * (i - (t_N + T_D))))
    for i in range(T):
        for ar in areas:
            if ar != most_infected_area:
                alpha[ar, i] = alpha[most_infected_area, max(i - L, 0)]
    for i in range(T - 1):
        for ar in areas:
            if i != 0:  # not a difference equation, just a simplifcation
                weighted_I[ar, i] = state_variables_simu[ar, "I", i] + (1 - p_e) * state_variables_simu[ar, "IV", i]
                area_infected_sum[ar, i] += area_infected_sum[ar, i - 1] + weighted_I[ar, i]

            vectors[ar, i] = weighted_I[ar, i] + sum(
                [(T_a_l[l, ar] / 365) * (weighted_I[l, i] / N_a[l]) * (1 / r_d) for l in areas]) \
                             - sum(
                [(T_a_l[ar, l] / 365) * (weighted_I[ar, i] / N_a[ar]) * (1 / r_d) for l in areas])  # travel
            state_variables_simu[ar, "S", i + 1] = state_variables_simu[ar, "S", i] - vax_dictionary[ar, i]\
                        - alpha[ar, i] *  state_variables_simu[ar, "S", i] * vectors[ar, i] / N_a[ar]

            sus_pop[ar, i] = state_variables_simu[ar, "S", i]
            being_exposed[ar, i] = alpha[ar, i] *  state_variables_simu[ar, "S", i] * vectors[ar, i] / N_a[ar]

            state_variables_simu[ar, "E", i + 1] = state_variables_simu[ar, "E", i] + alpha[ar, i] * \
                                                    state_variables_simu[
                                                       ar, "S", i] * vectors[ar, i] / N_a[ar] - r_l * \
                                                   state_variables_simu[ar, "E", i]
            state_variables_simu[ar, "I", i + 1] = state_variables_simu[ar, "I", i] + r_l * state_variables_simu[
                ar, "E", i] - r_d * state_variables_simu[ar, "I", i]
            state_variables_simu[ar, "SV", i + 1] = state_variables_simu[ar, "SV", i] + vax_dictionary[
                ar, i] - p_r * alpha[ar, i] * state_variables_simu[ar, "SV", i] * vectors[ar, i] * (
                                                                1 / N_a[ar])
                                                #1/(N_a[ar]-state_variables_simu[ar, "R", i]- state_variables_simu[ar,"D", i]))
            state_variables_simu[ar, "EV", i + 1] = state_variables_simu[ar, "EV", i] + p_r * alpha[ar, i] * \
                                                     state_variables_simu[
                                                        ar, "SV", i] * vectors[ar, i] / N_a[ar] - r_l * \
                                                    state_variables_simu[ar, "EV", i]
            state_variables_simu[ar, "IV", i + 1] = state_variables_simu[ar, "IV", i] + r_l * state_variables_simu[
                ar, "EV", i] - r_d * state_variables_simu[ar, "IV", i]
            state_variables_simu[ar, "H", i + 1] = state_variables_simu[ar, "H", i] + r_d * p_H * state_variables_simu[
                ar, "I", i] + r_d * p_V_H * state_variables_simu[ar, "IV", i] - r_R * state_variables_simu[ar, "H", i]
            state_variables_simu[ar, "D", i + 1] = state_variables_simu[ar, "D", i] + r_R * m_a[ar] * state_variables_simu[ar, "H", i]
            state_variables_simu[ar, "R", i + 1] = state_variables_simu[ar, "R", i] + r_R * (1 - m_a[ar]) * state_variables_simu[ar, "H", i]\
                                                   + r_d*(1-p_H)*state_variables_simu[ar, "I", i]\
                                                   + r_d*(1-p_V_H)*state_variables_simu[ar, "IV", i]
            if infected_sum < N and infected_sum + weighted_I[ar, i] >= N:
                t_N = i
            infected_sum += weighted_I[ar, i]
        if (i % 50 == 0 and i != 0) or i == T-2:
            print(f"Simulation step for iteration {iteration_number} computed up to day {i}.")
    for ar in areas:
        sus_pop[ar, T-1] = state_variables_simu[ar, "S", T-1]
        being_exposed[ar, T-1] = alpha[ar, T-1] * state_variables_simu[ar, "S", T-1] * vectors[ar, T-1] / N_a[ar]
    if min(state_variables_simu.values()) < -.001:
        print(f"Warning: In the simulation, some states became negative. The lowest state was {min(state_variables_simu.values())}.\
         You may have put in too many vaccinations.")
        if not four_binary:
            test_ar = "USA"
        else:
            test_ar = "area1"
        if state_variables_simu[test_ar, 'S', T-1] < -1:
            print(f"There is some evidence for that: in {test_ar}, the susceptible population ended with {state_variables_simu[test_ar, 'S', T-1]}.")
        #else:
            #print(f"sorry, here's a lot of data {state_variables_simu}")
            #exit()
    min_compartment = min(state_variables_simu.values())
    print(f"min value is {min_compartment}, max value is {max(state_variables_simu.values())}")
    if t_N != T: #If t_N is reached, and a variant is triggered...
        for ar in areas:
            if area_infected_sum[ar, t_N] == max([area_infected_sum[l, t_N] for l in areas]):
                most_infected_area = ar
        print(f"The variant is introduced on day {t_N}")
        print(f"The most infected area at day {t_N} is {most_infected_area}.")
        print("\n\n\n")
    else: #Otherwise,
        print("No variant was introduced.\n")
    infected_dictionary = {(ar, "I", i): state_variables_simu[ar, "I", i] for i in range(T) for ar in areas}
    infected_V_dictionary = {(ar, "IV", i): state_variables_simu[ar, "IV", i] for i in range(T) for ar in areas}
    total_deaths = sum([state_variables_simu[ar, "D", T - 1] for ar in areas])
    return infected_dictionary, infected_V_dictionary, total_deaths, alpha, t_N, most_infected_area, sus_pop, being_exposed, min_compartment
# Opt mod
exploration_tolerance = 1000 #1000
termination_tolerance = .5
iteration_number = 0
vax_allocation = {}
infections = {}
obj_function_values = {0: float("inf")}
for i in range(T):
    for ar in areas:
        vax_allocation[iteration_number, ar, i] = 0
        for inf in ["I", "IV"]:
            infections[iteration_number, ar, inf, i] = 0
print("\nInitialized problem.")

# Get feasible solution set up
iteration_number += 1  # it's just 1
vax_dictionary = {}

for i in range(T):
    for ar_i, ar in enumerate(areas):
        vax_dictionary[ar, i] = 0
def intensive_simulation_preparing_vaccination(vax_dictionary):
    global T, areas, sus_pop, N_a, B_t, being_exposed, global_population
    infected_I, infected_IV, obj, alpha, t_N, most_infected_area, sus_pop, being_exposed, min_compartment = simulate_covid(
        vax_dictionary=vax_dictionary)  # this run helps us find vax_dict for the next run
    for i in range(T):
        for ar in areas:
            if sus_pop[ar, i] < 1:
                vax_dictionary[ar, i] = 0
            else:
                vax_dictionary[ar, i] = min(B_t[i] * (N_a[ar] / global_population),
                                            1 * (sus_pop[ar, i] - being_exposed[ar, i]))
                vax_dictionary[ar, i] = max(vax_dictionary[ar, i], 0)
    return vax_dictionary, min_compartment

vax_dictionary, min_compartment = intensive_simulation_preparing_vaccination(vax_dictionary)
intensive_computation = 1
while intensive_computation < 5:
    print(f"Preparing to simulate for run number {intensive_computation+1} of 5.")
    vax_dictionary, min_compartment = intensive_simulation_preparing_vaccination(vax_dictionary)
    print(f"The minimum was {min_compartment}.")
    intensive_computation += 1
    #if intensive_computation > 10:
    #    print("We've done this computation too many times. Exiting.")
    #    exit()
print(f"The minimum value is positive {min_compartment} so we will now do two more then optimize.")
infected_I, infected_IV, first_obj, alpha, t_N, most_infected_area, sus_pop, being_exposed, min_compartment = simulate_covid(vax_dictionary=vax_dictionary)
for i in range(T):
    for ar in areas:
        if sus_pop[ar, i] < 1:
            vax_dictionary[ar, i] = 0
        else:
            vax_dictionary[ar, i] = min(B_t[i] * (N_a[ar] / global_population), sus_pop[ar, i] - being_exposed[ar, i])

for ar in areas:
    for i in range(T):
        vax_allocation[iteration_number, ar, i] = vax_dictionary[ar, i]
        infections[iteration_number, ar, "I", i] = infected_I[ar, "I", i]
        infections[iteration_number, ar, "IV", i] = infected_IV[ar, "IV", i]
obj_function_values[iteration_number] = first_obj
print(f"\nFound feasible solution with objective value {first_obj}.\n")

obj_difference = abs(
    obj_function_values[iteration_number] - obj_function_values[iteration_number - 1])  # right now it's infinity
infection_difference = (1 / num_areas) * sum([abs(i) for i in [
    infections[iteration_number, ar, "I", i] - infections[iteration_number - 1, ar, "I", i] for ar in areas for i in
    range(T)]])

# while max(obj_difference, infection_difference) > termination_tolerance & iteration_number < 30:
no_loop = True #Iterations can loop - in this case, terminate
if not four_binary:
    if sum(i > 0 for i in B_t.values()) == 0: #nothing to opt
        print(f"\nThe initial feasible solution had objective value {first_obj}.\n")
        print(f"Since there are no vaccinations to assign there is nothing to optimize. Stopping here.")
        exit()

while obj_difference > termination_tolerance and ((iteration_number < 30) and no_loop):
    # while iteration_number < 20:
    iteration_number += 1
    vax_dict = {(ar, i): vax_allocation[iteration_number - 1, ar, i] for ar in areas for i in range(T)}
    #infected_I, infected_IV, obj, alpha, t_N, most_infected_area, sus_pop, being_exposed, min_compartment = simulate_covid(vax_dictionary=vax_dict)
    print(f"obj is {first_obj} at iteration {iteration_number}")
    for ar in areas:
        for i in range(T):
            infections[iteration_number, ar, "I", i] = infected_I[ar, "I", i]
            infections[iteration_number, ar, "IV", i] = infected_IV[ar, "IV", i]
    # Put infections into opt Mod
    weighted_I = {}
    vectors = {}
    travel_sum = {} #Not used, just good to look at
    for i in range(T):
        for ar in areas:
            weighted_I[ar, i] = infected_I[ar, "I", i] + (1 - p_e) * infected_IV[ar, "IV", i]
        for ar in areas:
            vectors[ar, i] = weighted_I[ar, i] + sum(
                [(T_a_l[l, ar] / 365) * (weighted_I[l, i] / N_a[l]) * (1 / r_d) for l in areas]) \
                             - sum(
                [(T_a_l[ar, l] / 365) * (weighted_I[ar, i] / N_a[ar]) * (1 / r_d) for l in areas])
            travel_sum[ar, i] = vectors[ar,i] - weighted_I[ar, i] #Not used, just good to look at
    opt_Mod = gp.Model("vaccine_opt" + str(iteration_number))
    opt_Mod.setParam("NumericFocus",1)
    # opt_Mod.setParam("MIPGap", 5e-3)

    # Create variables
    #state_vars = ["S", "E", "I",
    #              "SV", "EV", "IV",
    #              "H", "D", "R"]

    # Define Variables. All are continuous and nonnegative by default
    state_variables = opt_Mod.addVars(areas, state_vars, T, name="state_variables")
    vax_variables = opt_Mod.addVars(areas, T, name="vax_variables")
    diff_infected = opt_Mod.addVars(areas, T, lb=-GRB.INFINITY,
                                    name="diff_infected")  # auxiliary variables. lower bounds because this can be negative
    abso_infected = opt_Mod.addVars(areas, T, name="abso_infected")  # auxiliary variables.
    diff_vax = opt_Mod.addVars(areas, (T - 1), lb=-GRB.INFINITY,
                               name="diff_vax")  # auxiliary variables. lower bounds because this can be negative
    abso_vax = opt_Mod.addVars(areas, T - 1, name="abso_vax")  # auxiliary variables.
    #crucial. Shows vac less than budget
    opt_Mod.addConstrs((vax_variables.sum('*', i) <= B_t[i]
                        for i in [*range(T)]), "Vax_Budget")
    #Only 95% of sus pop are willing to be vaxxed
    opt_Mod.addConstrs((vax_variables.sum(ar, '*') <= rho * (initial_pop_a[ar, "S"]+initial_pop_a[ar, "SV"]) - initial_pop_a[ar, "SV"]
                        for ar in areas), "Vax_Willingness")
    #10*budget*pop/global pop  = 10*1837*1/4 = 5000. Not an active constraint since N/global pop is so big.
    opt_Mod.addConstrs((vax_variables[ar, i] <= f_max_a[ar, i] * N_a[ar]
                        for ar in areas for i in [*range(T)]), "Vaccine_Capacity")
    opt_Mod.addConstrs((vax_variables[ar, i] <= state_variables[ar, "S", i]
                        for ar in areas for i in [*range(T)]), "Vax_Real")

    opt_Mod.addConstrs((diff_vax[ar, i] == vax_variables[ar, i + 1] - vax_variables[ar, i]
                        for ar in areas for i in [*range(T - 1)]), "diff_vax_constraint")
    opt_Mod.addConstrs((abso_vax[ar, i] == gp.abs_(diff_vax[ar, i])
                        for ar in areas for i in [*range(T - 1)]), "abso_vax_constraint")
    opt_Mod.addConstrs((abso_vax[ar, i] <= delta * f_max_a[ar, i] * N_a[ar]
                        for ar in areas for i in [*range(T - 1)]), "Vaccine_Agility")

    # Objective Function
    opt_Mod.setObjective(state_variables.sum(donor, "D", (T - 1)), GRB.MINIMIZE)  # deaths over donor
    #opt_Mod.setObjective(vax_variables["area1", 10], GRB.MINIMIZE)  # deaths over all areas

    # Add constraints
    opt_Mod.addConstrs((state_variables[ar, compartment, 0] == initial_pop_a[ar, compartment]
                        for ar in areas for compartment in initial_pop_states), "initial_states")
    opt_Mod.addConstrs((state_variables[ar, compartment, 0] == 0
                    for ar in areas for compartment in
                    [i for i in state_vars if i not in initial_pop_states]),
                   "initials")
    opt_Mod.addConstrs((state_variables[ar, compartment, T - 1] <= N_a[ar]
                        for ar in areas for compartment in state_vars), "upper_bound")
    print("\n\n\n\n\n")
    # opt_Mod.addConstrs((vectors[ar, i] == state_variables[ar, "I", i] + (1 - p_e) * state_variables[ar, "IV", i]
    #                    for ar in areas for i in [*range(T)]), "vectors_equality")
    opt_Mod.addConstrs(
        (state_variables[ar, "S", i + 1] == state_variables[ar, "S", i] - vax_variables[ar, i] - alpha[ar, i] *
         state_variables[ar, "S", i] * vectors[ar, i] / N_a[ar]
         for ar in areas for i in [*range(T - 1)]), "S_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "E", i + 1] == state_variables[ar, "E", i] + alpha[ar, i] * state_variables[
            ar, "S", i] * vectors[ar, i] / N_a[ar] - r_l * state_variables[ar, "E", i]
         for ar in areas for i in [*range(T - 1)]), "E_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "I", i + 1] == state_variables[ar, "I", i] + r_l * state_variables[ar, "E", i] - r_d *
         state_variables[ar, "I", i]
         for ar in areas for i in [*range(T - 1)]), "I_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "SV", i + 1] == state_variables[ar, "SV", i] + vax_variables[ar, i] - p_r * alpha[ar, i] *
         state_variables[ar, "SV", i] * vectors[ar, i] * (1 / N_a[ar])
         for ar in areas for i in [*range(T - 1)]), "SV_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "EV", i + 1] == state_variables[ar, "EV", i] + p_r * alpha[ar, i] *
         state_variables[
             ar, "SV", i] * vectors[ar, i] / N_a[ar] - r_l * state_variables[ar, "EV", i]
         for ar in areas for i in [*range(T - 1)]), "EV_diff")
    opt_Mod.addConstrs((state_variables[ar, "IV", i + 1] == state_variables[ar, "IV", i] + r_l * state_variables[
        ar, "EV", i] - r_d * state_variables[ar, "IV", i]
                        for ar in areas for i in [*range(T - 1)]), "IV_diff")
    opt_Mod.addConstrs((state_variables[ar, "H", i + 1] == state_variables[ar, "H", i] + r_d * p_H * state_variables[
        ar, "I", i] + r_d * p_V_H * state_variables[ar, "IV", i] - r_R * state_variables[ar, "H", i]
                        for ar in areas for i in [*range(T - 1)]), "H_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "D", i + 1] == state_variables[ar, "D", i] + r_R * m_a[ar] * state_variables[ar, "H", i]
         for ar in areas for i in [*range(T - 1)]), "D_diff")
    opt_Mod.addConstrs((state_variables[ar, "R", i + 1] == state_variables[ar, "R", i] + r_R * (1 - m_a[ar]) *
                        state_variables[ar, "H", i] \
                        + r_d * (1 - p_H) * state_variables[ar, "I", i] \
                        + r_d * (1 - p_V_H) * state_variables[ar, "IV", i]
                        for ar in areas for i in [*range(T - 1)]), "R_diff")
    opt_Mod.write("optimization_model.lp")
    opt_Mod.update()

    #opt_Mod.addConstrs((diff_infected[ar, i] == state_variables[ar, "I", i] - infections[iteration_number, ar, "I", i]
    #                    + state_variables[ar, "IV", i] - infections[iteration_number, ar, "IV", i]
    #                    for ar in areas for i in [*range(T)]), "diff_infected_constraint")
    #opt_Mod.addConstrs((abso_infected[ar, i] == gp.abs_(diff_infected[ar, i])
    #                    for ar in areas for i in [*range(T)]), "abso_infected_constraint")
    #opt_Mod.addConstrs((abso_infected[ar, i] <= exploration_tolerance
    #                    for ar in areas for i in [*range(T)]), "exploration_constraint")

    try:
        opt_Mod.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        # Solve bilinear model
        opt_Mod.Params.NonConvex = 2
        opt_Mod.optimize()
    opt_Mod.write("test.lp")
    # update Vac for sim_mod to be output of opt_mod

    for var in opt_Mod.getVars():
        var_name_list = var.VarName.split(",")
        if "vax_variables" in var_name_list[0]:
            # format is ['vax_variables[area4', '79]']
            if not four_binary:
                area_number = var_name_list[0][-3:]
            else:
                area_number = var_name_list[0][-5:]
            day = int(var_name_list[-1][:-1])
            vax_allocation[iteration_number, area_number, day] = var.x
    obj_function_values[iteration_number] = opt_Mod.getObjective().getValue()
    obj_difference = abs(obj_function_values[iteration_number] - obj_function_values[iteration_number - 1])
    infection_difference = (1 / num_areas) * sum([abs(i) for i in [
        infections[iteration_number, ar, "I", i] - infections[iteration_number - 1, ar, "I", i] for ar in areas for i in
        range(T)]])
    if iteration_number > 6: #check for loops
        #rounds the obj dictionary to be at the term tolerance
        tol = termination_tolerance/2
        rounded_obj = [round(i / tol) * tol for i in list(obj_function_values.values())[2:]]
        #rounded_obj = [i in obj_function_values.values()]
        number_of_repeated_obj_values = len(rounded_obj)-len(set(rounded_obj))
        if number_of_repeated_obj_values > 2:
            print("Warning: The algorithm is looping. Terminating due to lack of convergence.")
            no_loop = False
    print(obj_difference, infection_difference)
print(vax_dictionary)
print_solution = [[] for i in range(num_areas)]
var_names = []
var_values = []
vax_var_names = []
vax_var_values = []

for var in opt_Mod.getVars():
    var_name_list = var.VarName.split(",")
    if "vax_variables" in var_name_list[0]:
        # format is ['vax_variables[ISO', '79]']
        day = int(var_name_list[-1][:-1])
        if not four_binary:
            area_number = var_name_list[0][-3:]
            print_solution[area_lookup[area_number]].append(float(var.x))
        else:
            area_number = var_name_list[0][-5:]
            print_solution[int(area_number[-1]) - 1].append(float(var.x))
        # print(f"At {area_number} on day {day} you should give {var.x} vaccinations.")
        if var.x >= 0:
            vax_var_names.append(str(var.varName))
            vax_var_values.append(var.x)
    if str(var.varName)[:5] == "state" and var.x >= 0:
        var_names.append(str(var.varName))
        var_values.append(var.X)
vax_decision_variables = []
for i in print_solution:
    i = [0 if j < 0.0001 else j for j in i]
    print(i)
    vax_decision_variables.append(i)
print(f"Bt is: {B_t}")
print("vaccines given ", np.array(vax_decision_variables).sum(axis=0))
print(f"\nThe initial feasible solution had objective value {first_obj}.\n")
print(f"After {iteration_number} iterations, we found the solution {obj_function_values[iteration_number]}.")
print(f"The objective values were {obj_function_values}.")
    #vax_dict = {(ar, i): vax_allocation[iteration_number, ar, i] for ar in areas for i in range(T) if
    #            vax_allocation[iteration_number, ar, i] > 1}
#infect_dict = {(ar, i): infections[iteration_number, ar, i] for ar in areas for i in range(T)}
#print(f"The infections were {infected_I}")
#print(f"The final vaccinations were {vax_dict}")


if not four_binary:
    print_area = "all_countries_"
else:
    print_area = "four_countries_"
filename = f"donor_{print_area}global_simulation_{ALPHA_VALUE}_alpha.csv"
with open("simulation_data/" + filename, 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerows(zip(var_names, var_values))
    wr.writerows(zip(vax_var_names, vax_var_values))
print(f"State Variables written to {filename}.")
opt_Mod.write("solution_sol.sol")
