# Area inputs for the seir_opt python script.

A = [           # Set of areas to be simulated
        "area1",
        "area2" 
    ]

donor = "area1"    # Donor area ### D -> donor

m = "area2"       # Variant area

n = 10000  	# Person infection-days until new variant

N = {       	# Population
        "area1": 100000,	# 10^5
        "area2": 100000
    }
rho_V = {       # Vaccination rate
        "area1": 0.5,
        "area2": 0.1
    }

rho_I_N = {     # Initial cases per day = rho_I * N
        "area1": 10,
        "area2": 10
    }

delta_r = {     # Testing effect on rate out of I 
        "area1": 0,
        "area2": 0
    }

gamma = {       # Behavior infection multiplier
        "area1": 1,
        "area2": 1
    }

rho = {       # Proportion willing to be vaccinated 
        "area1": 1.0,
        "area2": 1.0
    }