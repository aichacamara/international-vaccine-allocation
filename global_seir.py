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
import time
if four_binary:
    from user_input_params import * #T, t_N, alpha, num areas, vax days, b_t, areas
    from four_mid_params import * #only differ in their N
    from four_mid_seir_data import * #most_infected_area, state_vars,
    #p_D, vax_rates, initial pop, N_a, travel_rates, tests_per_day, T_a_l
else:
    from user_input_params_real import * #alpha
    from mid_process_df import *  #B_t, T, t_N
    from mid_params_real import * #all
    from mid_seir_data_real import * #num_areas, areas, same as four_mid_seir_data but no travel_rates
import csv

start_time = time.time()
initial_pop_states = ["S", "E", "I", "SV", "IV", "EV"]
many_areas = areas #Used in case areas is fewer than expected
small_run = False
if four_binary:
    donor = "area1"
else:
    if small_run:
        areas = ["USA", "GRC", "BRA", "PER", "JPN", "RUS", "POL", "IND", "GHA", "ZMB", "VEN", "PAK"]
        full_population = sum(N_a.values())
        N_a = {ar: N_a[ar] for ar in areas}
        global_population = sum(N_a.values())
        N *= global_population/full_population
        B_t = {i: B_t[i]*(global_population/full_population) for i in range(len(B_t))}
        num_areas = len(areas)
    donor = "USA"
    initial_pop_a = {}
    for state_var in initial_pop_states:
        for ar in areas:
            initial_pop_a[ar, state_var] = vaccine_params_df[state_var][many_areas.index(ar)]


global_population = sum(N_a.values())
r_d_t = {ar: r_d + w * (tests_per_day[ar] / N_a[ar]) for ar in areas}


#Create a simulate object, with a new object once per iteration.
#The object holds all the values

def simulate_covid(vax_dictionary):
    global t_N, most_infected_area, four_binary, donor
    print(f"------------------starting simulation----------------")
    print(f"Variant debug info: t_N = {t_N}, most infected area = {most_infected_area}.")
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
   # alpha = {(most_infected_area, i): ALPHA_VALUE + lmbda / (1 + np.exp(-k * (i - (t_N + T_D)))) for i in range(T-1)}   # Calculates alpha(t)
    for i in range(T-1):
        for ar in areas:
            if ar == most_infected_area:
                if i < t_N:
                    alpha[ar,i] = ALPHA_VALUE
                else:
                    alpha[ar,i] = ALPHA_VALUE + lmbda / (1 + np.exp(-k * (i - (t_N + T_D))))               
            else:
                alpha[ar, i] = alpha[most_infected_area, max(i - L, 0)]
    for i in range(T - 1):
        #We need to calculate the infected counts for all areas first
        #That means when we calculate the travel params, we have the info for all areas.
        for ar in areas:
            if i != 0:
                weighted_I[ar, i] = state_variables_simu[ar, "I", i] + (1 - p_e) * state_variables_simu[ar, "IV", i]
                area_infected_sum[ar, i] = area_infected_sum[ar, i - 1] + weighted_I[ar, i]
        for ar in areas:
            #if i != 0:  # not a difference equation, just a simplifcation
            #    weighted_I[ar, i] = state_variables_simu[ar, "I", i] + (1 - p_e) * state_variables_simu[ar, "IV", i]
            #    area_infected_sum[ar, i] = area_infected_sum[ar, i - 1] + weighted_I[ar, i]

            vectors[ar, i] = weighted_I[ar, i] + (1 / 365)*sum(
                [T_a_l[l, ar] * (weighted_I[l, i] / (N_a[l]*r_d_t[l])) for l in areas]) \
                             - (1 / 365)*sum(
                [T_a_l[ar, l] * (weighted_I[ar, i] / (N_a[ar]*r_d_t[ar])) for l in areas])  # travel
            state_variables_simu[ar, "S", i + 1] = state_variables_simu[ar, "S", i] - vax_dictionary[ar, i]\
                        - alpha[ar, i] * state_variables_simu[ar, "S", i] * vectors[ar, i] / N_a[ar]

            state_variables_simu[ar, "S", i+1] = max(state_variables_simu[ar, "S", i+1], 0)
            sus_pop[ar, i] = state_variables_simu[ar, "S", i]

            state_variables_simu[ar, "E", i + 1] = state_variables_simu[ar, "E", i] + alpha[ar, i] * \
                                                    state_variables_simu[ar, "S", i]* vectors[ar, i] / N_a[ar] - r_l * \
                                                   state_variables_simu[ar, "E", i]
            state_variables_simu[ar, "I", i + 1] = state_variables_simu[ar, "I", i] + r_l * state_variables_simu[
                ar, "E", i] - r_d_t[ar] * state_variables_simu[ar, "I", i]
            state_variables_simu[ar, "SV", i + 1] = max(state_variables_simu[ar, "SV", i] + vax_dictionary[
                ar, i] - p_r * alpha[ar, i] * state_variables_simu[ar, "SV", i] * vectors[ar, i] * (
                                                                1 / N_a[ar]),0)
                                                #1/(N_a[ar]-state_variables_simu[ar, "R", i]- state_variables_simu[ar,"D", i]))
            state_variables_simu[ar, "EV", i + 1] = state_variables_simu[ar, "EV", i] + p_r * alpha[ar, i] * \
                                                     state_variables_simu[
                                                        ar, "SV", i] * vectors[ar, i] / N_a[ar] - r_l * \
                                                    state_variables_simu[ar, "EV", i]
            state_variables_simu[ar, "IV", i + 1] = state_variables_simu[ar, "IV", i] + r_l * state_variables_simu[
                ar, "EV", i] - r_d_t[ar] * state_variables_simu[ar, "IV", i]
            state_variables_simu[ar, "H", i + 1] = state_variables_simu[ar, "H", i] + r_d_t[ar] * p_H * state_variables_simu[
                ar, "I", i] + r_d_t[ar] * p_V_H * state_variables_simu[ar, "IV", i] - r_R * state_variables_simu[ar, "H", i]
            state_variables_simu[ar, "D", i + 1] = state_variables_simu[ar, "D", i] + r_R * p_D * state_variables_simu[ar, "H", i]
            state_variables_simu[ar, "R", i + 1] = state_variables_simu[ar, "R", i] + r_R * (1 - p_D) * state_variables_simu[ar, "H", i]\
                                                   + r_d_t[ar]*(1-p_H)*state_variables_simu[ar, "I", i]\
                                                   + r_d_t[ar]*(1-p_V_H)*state_variables_simu[ar, "IV", i]
            if infected_sum < N and infected_sum + weighted_I[ar, i] >= N:
                t_N = i
            infected_sum += weighted_I[ar, i]
        if (i % 50 == 0 and i != 0) or i == T-2:
            print(f"Simulation step for iteration {iteration_number} computed up to day {i}.")
    if min(state_variables_simu.values()) < -.001:
        print(f"Warning: In the simulation, some states became negative. The lowest state was {min(state_variables_simu.values())}.\
         You may have put in too many vaccinations.")
        if state_variables_simu[donor, 'S', T-1] < -1:
            print(f"There is some evidence for that: in the donor area, {donor}, the susceptible population ended with {state_variables_simu[donor, 'S', T-1]}.")
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
                most_infected_area = ar
        print(f"The variant is introduced on day {t_N}")
        print(f"The most infected area at day {t_N} is {most_infected_area}.")
        print("\n")
    else: #Otherwise,
        print("No variant was introduced.\n")
    infected_dictionary = {(ar, "I", i): state_variables_simu[ar, "I", i] for i in range(T) for ar in areas}
    infected_V_dictionary = {(ar, "IV", i): state_variables_simu[ar, "IV", i] for i in range(T) for ar in areas}
    total_deaths = sum([state_variables_simu[ar, "D", T - 1] for ar in areas])
    print(f"The total number of deaths in this simulation was {total_deaths}.")
    print(f"The total number of deaths in the donor country was {state_variables_simu[donor, 'D', T - 1]}.")
    return infected_dictionary, infected_V_dictionary, state_variables_simu[donor, 'D', T - 1], alpha, t_N, most_infected_area, sus_pop, min_compartment


iteration_number = 0
vax_allocation = {(iteration_number, ar, i): 0 for i in range(T-1) for ar in areas}
infections = {(iteration_number, ar, inf, i): 0 for i in range(T) for ar in areas for inf in ["I", "IV"]}
obj_function_values = {0: float("inf")}
print("\nInitialized problem.")
iteration_number += 1  # it's just 1
def intensive_simulation_preparing_vaccination(vax_dictionary):
    global T, areas, N_a, B_t, global_population
    infected_I, infected_IV, obj, alpha, t_N, most_infected_area, sus_pop, min_compartment = simulate_covid(
        vax_dictionary=vax_dictionary)  # this run helps us find vax_dict for the next run
    for i in range(T-1):
        for ar in areas:
            if sus_pop[ar, i] < 1:
                vax_dictionary[ar, i] = 0
            else:
                vax_dictionary[ar, i] = min(B_t[i] * (N_a[ar] / global_population), sus_pop[ar, i])
                vax_dictionary[ar, i] = max(vax_dictionary[ar, i], 0) #removing this changes the output, somehow
    return vax_dictionary, min_compartment


# Get feasible solution set up
vax_dictionary = {(ar, i): B_t[i] * (N_a[ar] / global_population) for ar in areas for i in range(T-1)}
intensive_computation = 0
while intensive_computation < 3:
    print(f"Preparing to simulate for run number {intensive_computation+1} of 3.")
    vax_dictionary, min_compartment = intensive_simulation_preparing_vaccination(vax_dictionary)
    print(f"The minimum was {min_compartment}.")
    intensive_computation += 1
print(f"The minimum value is {min_compartment} (hopefully positive or near 0) so we will now do two more then optimize.")
infected_I, infected_IV, first_obj, alpha, t_N, most_infected_area, sus_pop, min_compartment = simulate_covid(vax_dictionary=vax_dictionary)
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


obj_difference = abs(
    obj_function_values[iteration_number] - obj_function_values[iteration_number - 1])  # right now it's infinity
infection_difference = sum([abs(i) for i in [
    infections[iteration_number, ar, "I", i] - infections[iteration_number - 1, ar, "I", i] for ar in areas for i in
    range(T)]])
opt_vax_totals = {iteration_number: 0}
no_loop = True #Iterations can loop - in this case, terminate
if not four_binary:
    if sum(i > 0 for i in B_t.values()) == 0: #nothing to opt
        print(f"\nThe initial feasible solution had objective value {first_obj}.\n")
        print(f"Since there are no vaccinations to assign there is nothing to optimize. Stopping here.")
        exit()


iteration_limit = 20
if not four_binary:
    if small_run:
        exploration_tolerance = 1000000
        termination_tolerance = 50  # .01
    else:
        exploration_tolerance = 1000000#0
        termination_tolerance = 100#.01
else:
    exploration_tolerance = 100
    termination_tolerance = .01
gamma = .9


while (obj_difference > termination_tolerance and ((iteration_number < iteration_limit) and no_loop)):
    iteration_number += 1
    opt_vax_totals[iteration_number] = 0
    vax_dict = {(ar, i): vax_allocation[iteration_number - 1, ar, i] for ar in areas for i in range(T-1)}
    infected_I, infected_IV, obj, alpha, t_N, most_infected_area, sus_pop, min_compartment = simulate_covid(vax_dictionary=vax_dict)
    print(f"Simulation obj is {obj} at iteration {iteration_number}")
    for ar in areas:
        for i in range(T):
            infections[iteration_number, ar, "I", i] = infected_I[ar, "I", i]
            infections[iteration_number, ar, "IV", i] = infected_IV[ar, "IV", i]
    weighted_I = {(ar, i): infected_I[ar, "I", i] + (1 - p_e) * infected_IV[ar, "IV", i] for i in range(T) for ar in areas}
    vectors = {(ar, i): weighted_I[ar, i] + (1 / 365) * sum(
                [T_a_l[l, ar] * (weighted_I[l, i] / (N_a[l] * r_d_t[l])) for l in areas]) \
                             - (1 / 365) * sum(
                [T_a_l[ar, l] * (weighted_I[ar, i] / (N_a[ar] * r_d_t[ar])) for l in areas]) for i in range(T-1) for ar in areas}
    travel_sum = {(ar, i): vectors[ar,i] - weighted_I[ar, i] for i in range(T-1) for ar in areas} #Not used, just good information to print
    entering_donor = {(donor, i): (1 / 365) * sum([T_a_l[l, donor] * (weighted_I[l, i] / (N_a[l] * r_d_t[l])) for l in areas]) for i in range(T-1)}

    opt_Mod = gp.Model("vaccine_opt")
    opt_Mod.setParam("NumericFocus",1)

    # Create variables
    # Define Variables. All are continuous and nonnegative by default
    state_variables = opt_Mod.addVars(areas, state_vars, T, name="state_variables")
    vax_variables = opt_Mod.addVars(areas, T-1, name="vax_variables")
    diff_infected = opt_Mod.addVars(areas, ["I", "IV"], T, lb=-GRB.INFINITY,
                                    name="diff_infected")  # auxiliary variables. lower bounds because this can be negative
    abso_infected = opt_Mod.addVars(areas, ["I", "IV"], T, name="abso_infected")  # auxiliary variables.
    diff_vax_agility = opt_Mod.addVars(areas, T - 2, lb=-GRB.INFINITY,
                               name="diff_vax_agility")  # auxiliary variables. lower bounds because this can be negative
    abso_vax_agility = opt_Mod.addVars(areas, T - 2, name="abso_vax_agility")  # auxiliary variables.

    if vaccine_incentive:
        opt_Mod.setObjective(state_variables[donor, "D", T-1] - 1/200000*vax_variables.sum('*', '*'), GRB.MINIMIZE)  # deaths over donor
    else:
        opt_Mod.setObjective(state_variables[donor, "D", T-1], GRB.MINIMIZE)  # deaths over donor

    if four_binary and strict_policy:
        if p_k == .5:
            exploration_tolerance *= 15#4.4#3.46
        elif p_k == .75:
            exploration_tolerance *= 4#1000#3.46
        elif p_k == 1:
            exploration_tolerance *= 3.46
    upper_lim = {ar: rho * N_a[ar] - initial_pop_a[ar, "SV"] for ar in areas}
    opt_Mod.update()
    if strict_policy:
        temp_variables2_W = opt_Mod.addVars(T-1, name="temporary_variables2_W")
        opt_Mod.addConstr(temp_variables2_W[0] == upper_lim[donor], "temp_2_0")
        opt_Mod.addConstrs((temp_variables2_W[i] == upper_lim[donor] - vax_variables.sum(donor, range(i))
                            for i in [*range(1,T-1)]), "temp_2_T")
        opt_Mod.addConstrs((vax_variables[donor, i] == gp.min_(temp_variables2_W[i], state_variables[donor, "S", i], constant = p_k*B_t[i])
                            for i in [*range(10, T-1)]), "Vaccine_Policy_X")
        opt_Mod.addConstrs((vax_variables[donor, i] == p_k*B_t[i]
                           for i in [*range(10)]), "Vaccine_Policy2")
    elif policy:
        opt_Mod.addConstrs((vax_variables[donor, i] <= p_k*B_t[i]
                        for i in [*range(T-1)]), "Vaccine_Policy2")

    #Full budget use, described in the appendix
    if False:
        opt_Mod.addConstrs((vax_variables.sum('*', i) <= B_t[i] for i in [*range(85,T-1)]), "Vax_Budget")
        opt_Mod.addConstrs((vax_variables.sum('*', i) == B_t[i] for i in [*range(85)]), "Vax_Budget_e")
    else:
        opt_Mod.addConstrs((vax_variables.sum('*', i) <= B_t[i] for i in [*range(T - 1)]), "Vax_Budget")


    #Only 95% of sus pop are willing to be vaxxed
    opt_Mod.addConstrs((vax_variables.sum(ar, '*') <= upper_lim[ar]
                        for ar in areas), "Vax_Willingness")
    opt_Mod.addConstrs((vax_variables[ar, i] <= state_variables[ar, "S", i]
                        for ar in areas for i in [*range(T-1)]), "Vax_Real")
    print(f"Willingness constraint ensures vaccinations are less than {upper_lim}")

    opt_Mod.addConstrs((state_variables[ar, compartment, 0] == initial_pop_a[ar, compartment]
                        for ar in areas for compartment in initial_pop_states), "initial_states")
    opt_Mod.addConstrs((state_variables[ar, compartment, 0] == 0
                    for ar in areas for compartment in
                    [i for i in state_vars if i not in initial_pop_states]),
                   "initials")
    opt_Mod.addConstrs((state_variables[ar, compartment, T - 1] <= N_a[ar]
                        for ar in areas for compartment in state_vars), "upper_bound")

    opt_Mod.addConstrs((diff_infected[ar, infection_state, i] == state_variables[ar, infection_state, i] - infections[iteration_number, ar, infection_state, i]
                        for ar in areas for infection_state in ["I", "IV"] for i in [*range(T)]), "diff_infected_constraint")
    opt_Mod.addConstrs((abso_infected[ar, infection_state, i] == gp.abs_(diff_infected[ar, infection_state, i])
                        for ar in areas for infection_state in ["I", "IV"] for i in [*range(T)]), "abso_infected_constraint")
    opt_Mod.addConstrs((abso_infected[ar, infection_state, i] <= exploration_tolerance #2.4
                        for ar in areas for infection_state in ["I", "IV"] for i in [*range(T)]), "exploration_constraint")
    opt_Mod.addConstrs((diff_vax_agility[ar, i] == vax_variables[ar, i + 1] - vax_variables[ar, i]
                        for ar in areas for i in [*range(T - 2)]), "diff_vax_agility_constraint")
    opt_Mod.addConstrs((abso_vax_agility[ar, i] == gp.abs_(diff_vax_agility[ar, i])
                        for ar in areas for i in [*range(T - 2)]), "abso_vax_agility_constraint")
    if not strict_policy:
        opt_Mod.addConstrs((abso_vax_agility[ar, i] <= B_t[i]*N_a[ar]/ global_population
                        for ar in areas for i in [*range(T - 2)]), "Vaccine_Agility")








    opt_Mod.addConstrs(
        (state_variables[ar, "S", i + 1] == state_variables[ar, "S", i] - vax_variables[ar, i] - alpha[ar, i] *
         state_variables[ar, "S", i] * (vectors[ar, i] / N_a[ar])
         for ar in areas for i in [*range(T-1)]), "S_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "E", i + 1] == state_variables[ar, "E", i] + alpha[ar, i] * state_variables[ar, "S", i] * (vectors[ar, i] / N_a[ar]) - r_l * state_variables[ar, "E", i]
         for ar in areas for i in [*range(T - 1)]), "E_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "I", i + 1] == state_variables[ar, "I", i] + r_l * state_variables[ar, "E", i] - r_d_t[ar] *
         state_variables[ar, "I", i]
         for ar in areas for i in [*range(T - 1)]), "I_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "SV", i + 1] == state_variables[ar, "SV", i] + vax_variables[ar, i] - p_r * alpha[ar, i] *
         state_variables[ar, "SV", i] * (vectors[ar, i] / N_a[ar])
         for ar in areas for i in [*range(T - 1)]), "SV_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "EV", i + 1] == state_variables[ar, "EV", i] + p_r * alpha[ar, i] *
         state_variables[ar, "SV", i] * (vectors[ar, i] / N_a[ar]) - r_l * state_variables[ar, "EV", i]
         for ar in areas for i in [*range(T - 1)]), "EV_diff")
    opt_Mod.addConstrs((state_variables[ar, "IV", i + 1] == state_variables[ar, "IV", i] + r_l * state_variables[
        ar, "EV", i] - r_d_t[ar] * state_variables[ar, "IV", i]
                        for ar in areas for i in [*range(T - 1)]), "IV_diff")
    opt_Mod.addConstrs((state_variables[ar, "H", i + 1] == state_variables[ar, "H", i] + r_d_t[ar] * p_H * state_variables[
        ar, "I", i] + r_d_t[ar] * p_V_H * state_variables[ar, "IV", i] - r_R * state_variables[ar, "H", i]
                        for ar in areas for i in [*range(T - 1)]), "H_diff")
    opt_Mod.addConstrs(
        (state_variables[ar, "D", i + 1] == state_variables[ar, "D", i] + r_R * p_D * state_variables[ar, "H", i]
         for ar in areas for i in [*range(T - 1)]), "D_diff")
    opt_Mod.addConstrs((state_variables[ar, "R", i + 1] == state_variables[ar, "R", i] + r_R * (1 - p_D) *
                        state_variables[ar, "H", i] \
                        + r_d_t[ar] * (1 - p_H) * state_variables[ar, "I", i] \
                        + r_d_t[ar] * (1 - p_V_H) * state_variables[ar, "IV", i]
                        for ar in areas for i in [*range(T - 1)]), "R_diff")
    opt_Mod.write("optimization_model.lp")
    opt_Mod.update()
    try:
        opt_Mod.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        # Solve bilinear model
        opt_Mod.Params.NonConvex = 2
        opt_Mod.optimize()
    # update Vac for sim_mod to be output of opt_mod
    toprintvar = 0
    print("Negative variables:")
    for var in opt_Mod.getVars():
        if toprintvar > 30:
            break
        if var.x < -1 and "diff_" not in var.VarName:
            print(var.Varname, var.x)
            toprintvar += 1
        var_name_list = var.VarName.split(",")
        if "vax_variables" in var_name_list[0]:
            # format is ['vax_variables[area4', '79]']
            if not four_binary:
                area_number = var_name_list[0][-3:]
            else:
                area_number = var_name_list[0][-5:]
            day = int(var_name_list[-1][:-1])
            vax_allocation[iteration_number, area_number, day] = var.x
            opt_vax_totals[iteration_number] += var.x
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

