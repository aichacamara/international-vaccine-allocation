import numpy as np


most_infected_area = "area3"
state_vars = ["S", "E", "I",
              "SV", "EV", "IV",
              "Q", "H", "D", "R"]
m_a = {
    "area1": 0.0156,
    "area2": 0.0213,
    "area3": 0.0181,
    "area4": 0.0282}  # Also needed are m__a_H
vax_rates_s = {
    "area1": 0.4604, #.742
    "area2": 0.0887,
    "area3": 0.1859,
    "area4": 0.0117}
initial_pop_a = {
    ("area1", "S"): (1-vax_rates_s["area1"])*100000,
    ("area1", "E"): 15.1*(5/3.5),
    ("area1", "I"): 15.1,
    ("area1", "SV"): vax_rates_s["area1"]*100000,
    ("area2", "S"): (1-vax_rates_s["area2"])*100000,
    ("area2", "E"): 5.8*(5/3.5),
    ("area2", "I"): 5.8,
    ("area2", "SV"): vax_rates_s["area2"]*100000,
    ("area3", "S"): (1-vax_rates_s["area3"])*100000,
    ("area3", "E"): 18.17*(5/3.5),
    ("area3", "I"): 18.17,
    ("area3", "SV"): vax_rates_s["area3"]*100000,
    ("area4", "S"): (1-vax_rates_s["area4"])*100000,
    ("area4", "E"): 2.15*(5/3.5),
    ("area4", "I"): 2.15,
    ("area4", "SV"): vax_rates_s["area4"]*100000}
all_areas = ["area1", "area2", "area3", "area4"]
for ar in all_areas:
    initial_pop_a[ar, "SV"] = vax_rates_s[ar]*100000
    initial_pop_a[ar, "S"] = (1-vax_rates_s[ar]) * 100000
N_a = {ar : sum(initial_pop_a[ar, i]
                for i in ["S", "E", "I", "SV"])
       for ar in all_areas}
travel_rates_s = {
    "area1": 1.3268,
    "area2": 0.588,
    "area3": 0.054,
    "area4": 0.017}
tests_per_day = {
    "area1": 70090.295,
    "area2": 45169.362,
    "area3": 33743.429,
    "area4": 3938.763}

Tr = [travel_rates_s[ar]*N_a[ar] for ar in all_areas]
#Tr = [4,3,2,1]
T_a_l_list = []
for i in range(4):
    summed_list = [(Tr[i] * j) / sum(Tr) for j in Tr]
    T_a_l_list.append(summed_list)
T_a_l = {}
for i in range(4):
    for j in range(4):
        T_a_l[all_areas[i], all_areas[j]] = T_a_l_list[i][j]
#T_1,1 is not a problem. In the equations, they cancel out! T_1,1 in leaving term = T_1,1 in arriving term




"""
#T_a_l is this for Tr = [10000, 10000, 1000, 1000]
T_a_l = {
    ("area1", "area1") : 0,
    ("area1", "area2") : 8334,
    ("area1", "area3") : 833,
    ("area1", "area4") : 833,
    ("area2", "area1") : 8334,
    ("area2", "area2") : 0,
    ("area2", "area3") : 833,
    ("area2", "area4") : 833,
    ("area3", "area1") : 476,
    ("area3", "area2") : 476,
    ("area3", "area3") : 0,
    ("area3", "area4") : 48,
    ("area4", "area1") : 476,
    ("area4", "area2") : 476,
    ("area4", "area3") : 48,
    ("area4", "area4") : 0
}
T_a_l = {
    ("area1", "area1") : 0,
    ("area1", "area2") : 4240,
    ("area2", "area1") : 4240,
    ("area2", "area2") : 0,
}

#Tr = [4,3,2,1]
T_a_l_list = []
for i in range(4):
    Tr_now = Tr[:i] + Tr[i+1:]
    Tr_sum = sum(Tr_now)
    summed_list = [(Tr[i] * j) / Tr_sum for j in Tr_now]
    T_a_l_list.append(summed_list[:i] + [0] + summed_list[i:])
for i in range(4):
    Tr_sum = sum(Tr)
    summed_list = [(Tr[i] * j) / Tr_sum for j in Tr]
    #summed_list = [(Tr[i] * j) / (Tr_sum*(Tr[i]-Tr[i]*Tr[i])) for j in Tr]
    #summed_list = [(Tr[i]/(Tr[i]-(Tr[i]*Tr[i])/Tr_sum))*((Tr[i] * j) / Tr_sum) for j in Tr]
    #summed_list = [((Tr[i] * j) / (Tr_sum-Tr[i])) for j in Tr]
    #summed_list = [((Tr[i] * j) / (Tr_sum - 0.5*(Tr[i]+j))) for j in Tr]
    #summed_list[i] = 0
    T_a_l_list.append(summed_list)
summm = 0
for i in range(4):
    print(sum(T_a_l_list[i]))
    summm += sum(T_a_l_list[i])
T_a_l = {}
for i in range(4):
    for j in range(4):
        T_a_l[all_areas[i], all_areas[j]] = T_a_l_list[i][j]
"""