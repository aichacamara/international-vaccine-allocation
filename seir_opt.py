import numpy as np
import time
import math

from areas import *
from scenario import *
from run_params import *


def simulate():
    '''
    Preforms and SEIR simulation based on data imported from areas, scenario, 
    and run_params python scripts
    '''
    # Define state equations for each area and each time in our time horizon
    S = {(area, i): 0 for area in A for i in range(0, T)}   # Susceptible
    S_V = {(area, i): 0 for area in A for i in range(0, T)} # Susceptible Vaccinated
    E = {(area, i): 0 for area in A for i in range(0, T)}   # Exposed
    E_V = {(area, i): 0 for area in A for i in range(0, T)} # Exposed Vaccinated
    I = {(area, i): 0 for area in A for i in range(0, T)}   # Infected
    I_V = {(area, i): 0 for area in A for i in range(0, T)} # Infected Vaccinated
    R = {(area, i): 0 for area in A for i in range(0, T)}   # Recovered
    H = {(area, i): 0 for area in A for i in range(0, T)}   # Hospitalized
    D = {(area, i): 0 for area in A for i in range(0, T)}   # Dead

    # Set up initial values for the state equations
    for area in A:
        E[area, 0] = ((1 - rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/r_I)
        E_V[area, 0] = ((p_r*rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/r_I)
        I[area, 0] = ((1 - rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/(r_0 + delta_r[area]))
        I_V[area, 0] = ((p_r*rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/(r_0 + delta_r[area]))
        S_V[area, 0] = rho_V[area]*N[area] - E_V[area, 0] - I_V[area, 0]
        S[area, 0] = N[area] - E[area, 0] - E_V[area, 0] - I[area, 0] - I_V[area, 0] - S_V[area, 0]

    


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

r_d_t = {ar: r_0 + w * (delta_r[ar] / N[ar]) for ar in A}

simulate()




