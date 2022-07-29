import numpy as np
from mid_process_df import *

most_infected_area = "IND" #USA has had 600k+ deaths by June 1, and India had 355k.
#Since it can't be the donor, we use India. https://ourworldindata.org/covid-deaths
state_vars = ["S", "E", "I",
              "SV", "EV", "IV",
              "H", "D", "R"]
all_areas = vaccine_params_df["iso_code"]
areas = list(vaccine_params_df["iso_code"])
area_lookup = {ar: area_num for area_num, ar in enumerate(all_areas)}
num_areas = len(areas)
#m_a = {vaccine_params_df["iso_code"][i] : vaccine_params_df["mortality_rate"][i]
#       for i in range(len(vaccine_params_df))}
p_D = 0.1
vax_rates_s = {vaccine_params_df["iso_code"][i] : vaccine_params_df["vax"][i]
       for i in range(len(vaccine_params_df))}
N_a = {vaccine_params_df["iso_code"][i] : vaccine_params_df["population"][i]
       for i in range(len(vaccine_params_df))}
tests_per_day = {vaccine_params_df["iso_code"][i] : vaccine_params_df["tests_per_day"][i]
       for i in range(len(vaccine_params_df))}
Tr = {vaccine_params_df["iso_code"][i] : vaccine_params_df["annual_passengers"][i]
       for i in range(len(vaccine_params_df))}
total_passengers = sum(Tr.values())
T_a_l = {(ar, l): Tr[ar]*Tr[l]/sum(Tr.values()) if Tr[ar]*Tr[l]/total_passengers >= 2 else 0 for ar in all_areas for l in all_areas}
#T_1,1 is not a problem. In the equations, they cancel out! T_1,1 in leaving term = T_1,1 in arriving term
