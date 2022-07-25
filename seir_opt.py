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
    alpha_0 = {{area}: a_0*gamma[area] for area in A}       # Store initial alpha
    delta_E = {(area, t): 0 for area in A for t in range(0, T)}
    V_star = {(area, t): 0 for area in A for t in range(0, T)}
    V_plus = {(area, t): 0 for area in A for t in range(0, T)}
    N_dot = {(t): 0 for t in range(0, T)}
    r_d = {{area}: r_0 + delta_r[area] for area in A}

    if len(b) == 0:
        B = {{t}: B_0 * 1 for t in range(0, T)}
    else:
        B = {{t}: B_0 * b[t] for t in range(0, T)}


    S_dot = S.sum(a,0) - S[donor,0]     # Sum of nondonor S at time 0
    V[a,t]  = {a,t: (1 - pk) * S[a,0]/Sdot * B[t] for a in areas for t in [0,T-1]}
    V[donor,t]  = {t: pk * B[t] for t in [0,T-1]}

    # Set up initial values for the state equations
    for area in A:
        E[area, 0] = ((1 - rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/r_I)
        E_V[area, 0] = ((p_r*rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/r_I)
        I[area, 0] = ((1 - rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/(r_0 + delta_r[area]))
        I_V[area, 0] = ((p_r*rho_V[area])/(p_r*rho_V[area] + 1 - rho_V[area]))*(rho_I_N[area]/(r_0 + delta_r[area]))
        S_V[area, 0] = rho_V[area]*N[area] - E_V[area, 0] - I_V[area, 0]
        S[area, 0] = N[area] - E[area, 0] - E_V[area, 0] - I[area, 0] - I_V[area, 0] - S_V[area, 0]
        W[area, 0] = rho[area]*N[area] - S_V[area, 0] - E_V[area, 0] - I_V[area, 0] - rho[area]*E[area, 0] - rho[area]*I[area, 0]
        V_calligraphy[area, 0] = I[area, 0] + p_e*I_V[area, 0]
        alpha[area, 0] = alpha_0[area]  # Initialize alpha in relation to time

    variant_emerge = False

    # loop over time
    for t in range(0, T):
        # loop over areas
        for area in A:
            # alpha calculated (only if variant found, equations 4 and 5)
            if variant_emerge:
                a_t = a_0 + delta_a/(1 + math.exp(-k*(t - (t_n + T_D))))
                alpha[area, t]

            # compute delta e (equation 7)
            delta_E[area, t] = min(S[area, t], alpha[area, t]*S[area, t]*V_calligraphy[area, t]/N[area])
            # v star (equation 8)
            V_star[area, t] = min(W[area, t] - W[area, t]*delta_E[area, t]/S[area, t], V)
            # compute w (equation 9)
            W[area, t + 1] = W[area, t] - W[area, t]*delta_E[area, t]/S[area, t] - V_star[area, t]

            # if realloc flag is true
            if realloc_flag:
                # compute equation 13
                for a in A[1:]: # TODO make so ignores donor, not necessarily first area
                    if W[area, t + 1] > 0:
                        N_dot += N[a]

                # compute equation 14 (update v with equation 14)
                V_plus[area, t]
                # update v star with v (equation 8)
                # update w with new v (equation 9)
            
            # difference equations
            S[area, t + 1] = S[area, t] - delta_E[area, t] - V_star[area, t]
            S_V[area, t + 1] = S_V[area, t] + V_star[area, t] - p_r*alpha[area, t]
            E[area, t + 1] = E[area, t] + delta_E[area, t] + r_I*E[area, t]
            E_V[area, t + 1] = E_V[area, t] + p_r*alpha[area, t] - r_I*E_V[area, t]
            I[area, t + 1] = I[area, t] + r_I*E[area, t] - r_d[area]*I[area, t]
            I_V[area, t + 1] = I_V[area, t] + r_I*E_V[area, t] - r_d[area]*I_V[area, t]
            H[area, t + 1] = H[area, t] + r_d[area]*p_H*I[area, t] + r_d[area]*p_V_H*I_V[area, t] - r_R*H[area, t]
            D[area, t + 1] = D[area, t] + r_R*p_D*H[area, t]
            R[area, t + 1] = R[area, t] + r_R*(1 - p_D)*H[area, t] + r_d[area]*(1 - p_H)*I[area, t] + r_d[area]*(1 - p_V_H)*I_V[area, t]


        # equation 10 (check for the variant occuring, do not calculate if variant already emerged)
        if not variant_emerge:
            cum_sum = 0
            # looping over area
            for area in A:
                # looping from 0 to the current time
                for t_star in range(0, t):
                    # sum all infected up to current time
                    cum_sum += I[area, t_star]
            # compare sum from above loop to person days (parameter n)
            if cum_sum > n:
                variant_emerge = True
                t_star_sum = 0
                for area in A:
                    t_star_sum += I[area, t]
                # if true calculate equation 11 (t_n)
                t_n = t - 1 + (cum_sum - n)/(t_star_sum)

    
# Calculate constant from given inputs
k = math.log((1-p)/p)/T_D     #natural log, ln()

simulate()




