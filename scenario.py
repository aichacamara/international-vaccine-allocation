# Scenario parameters for seir_opt script

r_I = 1 / 5     # rate out of exposed (E) state into infectious (I) state
r_0 = 1 / 3.5   # rate out of infectious state without testing
r_R =  1 / 15   # rate out of state hospitalized (H) into recovery or death
p_V_H = 0.0296  # P(Hospitalized | Infected and Vaxxed)
p_H = 0.296     # P(Hospitalized | Infected)
p_D = 0.1       # p(Dead | Hospitalized)
a_0 = 0.6       # Initial infection rate (proportion/day)
delta_a = 0.6   # Change in infection rate for a new variant (proportion/day)
p_e = 0.6       # Transmission rate from a vaccinated person (proportion of unvaccinated rate)
p_r = 0.6       # Infection rate for a vaccinated person (proportion of unvaccinated rate)
L = 20          # Lag parameter in days for the time a variant takes to reach other area
T_D = 45        # Time for variant to dominated area
p = 0.01        # proportion of people in infected state that have new variant when introduced
T = 180         # Time horizon
# b seems to be this value repeated T times, will ask if that is correct
B_0 = 1750      # Vaccine available day 1



