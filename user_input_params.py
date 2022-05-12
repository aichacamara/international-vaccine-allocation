T = 600
t_N = T
alpha_string = input("Alpha? ")
#alpha_string = "0.7"
num_areas = 4
try:
    ALPHA_VALUE = int(alpha_string)
except:
    ALPHA_VALUE = float(alpha_string)
p_k_string = input("Policy restrictions? Type 0 < p^k <= 1 for yes, 'n' for no, or 'M' for Vaccine Incentives.")
vaccine_incentives = p_k_string == "M"
try:
    p_k = float(p_k_string)
    policy = True
    print(f"Running with a policy p^k = {p_k}.")
except:
    print("Running with no policy.")
    policy = False
if policy and input("Strict policy? 'y' or 'n'.") == 'y':
    strict_policy = True
else:
    strict_policy = False
lmbda = ALPHA_VALUE
#vaccine_list_initial = [1750]*vax_days#+[0]*(T-vax_days)
B_t = {i: 1750 for i in range(T)}
areas = ["area1", "area2", "area3", "area4"][:num_areas]
#if num_areas == 2:
#    areas = ["area1", "area3"]
#    print("Using areas 1,4")
#else:
#    areas = ["area1", "area2", "area3", "area4"][:num_areas]

