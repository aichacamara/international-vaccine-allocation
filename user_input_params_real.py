alpha_string = input("Alpha? ")
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