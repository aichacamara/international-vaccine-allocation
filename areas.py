# This file is defines the area inputs for the seir_opt python script.

A = [           # The set of areas to be simulated
        "area1",
        "area2", 
        "area3", 
        "area4"
    ]

D = "area1"    # The donor area

m = "IND"       # Variant area

n = 1000000000  # The number of person infection days until new variant

# N = will be propagated eventually, do not have data to propogate right now

rho_V = {       # The vaccination rates for each area
        "area1": 0.432,
        "area2": 0.0743,
        "area3": 0.171,
        "area4": 0.010
    }

rho_I_N = {     # Initial cases per day by area
        "area1": 27.1,
        "area2": 9.96,
        "area3": 32.8,
        "area4": 3.78
    }

delta_r = {     # Testing effect on rate 
        "area1": 462.07,
        "area2": 70.02,
        "area3": 106.64,
        "area4": 21.17
    }

gamma = {       # Behavior infection multiplier by area
        "area1": 0.9,
        "area2": 0.9,
        "area3": 0.9,
        "area4": 0.9
    }

rho = {       # proportion willing to be vaccinated
        "area1": 0.95,
        "area2": 0.95,
        "area3": 0.95,
        "area4": 0.95
    }





