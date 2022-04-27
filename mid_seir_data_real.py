import numpy as np
from mid_process_df import *

most_infected_area = "USA" #USA has had 600k+ deaths by June 1, this is the most.
#https://ourworldindata.org/covid-deaths
state_vars = ["S", "E", "I",
              "SV", "EV", "IV",
              "H", "D", "R"]
all_areas = vaccine_params_df["iso_code"]
areas = list(vaccine_params_df["iso_code"])
area_lookup = {ar: area_num for area_num, ar in enumerate(all_areas)}
num_areas = len(areas)
m_a = {vaccine_params_df["iso_code"][i] : vaccine_params_df["mortality_rate"][i]
       for i in range(len(vaccine_params_df))}
vax_rates_s = {vaccine_params_df["iso_code"][i] : vaccine_params_df["vax"][i]
       for i in range(len(vaccine_params_df))}
N_a = {vaccine_params_df["iso_code"][i] : vaccine_params_df["population"][i]
       for i in range(len(vaccine_params_df))}
Tr = {vaccine_params_df["iso_code"][i] : vaccine_params_df["annual_passengers"][i]
       for i in range(len(vaccine_params_df))}
Tr = vaccine_params_df["annual_passengers"]
T_a_l_list = []
for i in range(len(Tr)):
    summed_list = [(Tr[i] * j) / sum(Tr) for j in Tr]
    T_a_l_list.append(summed_list)
T_a_l = {}
for i in range(len(Tr)):
    for j in range(len(Tr)):
        T_a_l[all_areas[i], all_areas[j]] = T_a_l_list[i][j]
        if T_a_l[all_areas[i], all_areas[j]] < 2:
            T_a_l[all_areas[i], all_areas[j]] = 0
#T_1,1 is not a problem. In the equations, they cancel out! T_1,1 in leaving term = T_1,1 in arriving term
