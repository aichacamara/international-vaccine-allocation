import numpy as np


most_infected_area = "area3"
state_vars = ["S", "E", "I",
              "SV", "EV", "IV",
              "H", "D", "R"]
vax_rates_s = {
    "area1": 0.432,
    "area2": 0.0743,
    "area3": 0.171,
    "area4": 0.010}
initial_infected = {
    "area1": 27.1,
    "area2": 9.96,
    "area3": 32.8,
    "area4": 3.78}
tests_per_day = {
    "area1": 462.07,
    "area2": 70.02,
    "area3": 106.64,
    "area4": 21.17}
initial_pop_states =  ["S", "SV", "E", "EV", "I", "IV"]
all_areas = ["area1", "area2", "area3", "area4"]
initial_pop_a = {(ar, compartment): vax_rates_s[ar] if "V" in compartment else 1-vax_rates_s[ar] for ar in all_areas for compartment in initial_pop_states}
initial_pop_a = {(ar, compartment): initial_pop_a[ar, compartment]*100000 if "S" in compartment else initial_infected[ar]*initial_pop_a[ar, compartment] for ar in all_areas for compartment in initial_pop_states}
initial_pop_a = {(ar, compartment): initial_pop_a[ar, compartment]*5/3.5 if "E" in compartment else initial_pop_a[ar, compartment] for ar in all_areas for compartment in initial_pop_states}

N_a = {ar : sum(initial_pop_a[ar, compartment]
                for compartment in initial_pop_states)
       for ar in all_areas}
travel_rates_s = {
    "area1": 1.3274,
    "area2": 0.588,
    "area3": 0.0546,
    "area4": 0.0173}

Tr = {ar: travel_rates_s[ar]*N_a[ar] for ar in all_areas}
total_passengers = sum(Tr.values())
T_a_l = {(ar, l): Tr[ar]*Tr[l]/total_passengers if Tr[ar]*Tr[l]/total_passengers >= 2 else 0 for ar in all_areas for l in all_areas}
#T_1,1 is not a problem. In the equations, they cancel out! T_1,1 in leaving term = T_1,1 in arriving term

"""
initial_pop_a = {
    ("area1", "S"): (1-vax_rates_s["area1"])*100000,
    ("area1", "E"): 27.1*(5/3.5)*(1-vax_rates_s["area1"]),
    ("area1", "EV"): 27.1*(5/3.5)*vax_rates_s["area1"],
    ("area1", "I"): 27.1*(1-vax_rates_s["area1"]),
    ("area1", "IV"): 27.1*vax_rates_s["area1"],
    ("area1", "SV"): vax_rates_s["area1"]*100000,
    ("area2", "S"): (1-vax_rates_s["area2"])*100000,
    ("area2", "E"): 9.96*(5/3.5)*(1-vax_rates_s["area2"]),
    ("area2", "EV"): 9.96*(5/3.5)*vax_rates_s["area2"],
    ("area2", "I"): 9.96 * (1 - vax_rates_s["area2"]),
    ("area2", "IV"): 9.96 * vax_rates_s["area2"],
    ("area2", "SV"): vax_rates_s["area2"]*100000,
    ("area3", "S"): (1-vax_rates_s["area3"])*100000,
    ("area3", "E"): 32.8*(5/3.5)*(1-vax_rates_s["area3"]),
    ("area3", "EV"): 32.8 * (5 / 3.5)*vax_rates_s["area3"],
    ("area3", "I"): 32.8 * (1 - vax_rates_s["area3"]),
    ("area3", "IV"): 32.8 * vax_rates_s["area3"],
    ("area3", "SV"): vax_rates_s["area3"]*100000,
    ("area4", "S"): (1-vax_rates_s["area4"])*100000,
    ("area4", "E"): 3.78*(5/3.5)*(1-vax_rates_s["area4"]),
    ("area4", "EV"): 3.78*(5/3.5)*vax_rates_s["area4"],
    ("area4", "I"): 3.78 * (1 - vax_rates_s["area4"]),
    ("area4", "IV"): 3.78 * vax_rates_s["area4"],
    ("area4", "SV"): vax_rates_s["area4"]*100000}
"""