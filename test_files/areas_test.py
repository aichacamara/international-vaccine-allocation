# Area inputs for the seir_opt python script.

A = [           # Set of areas to be simulated
        "area1",
        "area2", 
        "area3", 
        "area4"
    ]

donor = "area1"    # Donor area ### D -> donor

m = "area3"       # Variant area

n = 50000  	# Person infection-days until new variant

N = {       	# Population
        "area1": 100000,	# 10^5
        "area2": 100000,
        "area3": 100000,
        "area4": 100000
    }
rho_V = {       # Vaccination rate
        "area1": 0.432,
        "area2": 0.0743,
        "area3": 0.171,
        "area4": 0.010
    }

rho_I_N = {     # Initial cases per day = rho_I * N
        "area1": 27.1,
        "area2": 9.96,
        "area3": 32.8,
        "area4": 3.78
    }

delta_r = {     # Testing effect on rate out of I 
        "area1": 0,
        "area2": 0,
        "area3": 0,
        "area4": 0
    }

gamma = {       # Behavior infection multiplier
        "area1": 1,
        "area2": 1,
        "area3": 1,
        "area4": 1
    }

rho = {       # Proportion willing to be vaccinated 
        "area1": 0.95,
        "area2": 0.95,
        "area3": 0.95,
        "area4": 0.95
    }