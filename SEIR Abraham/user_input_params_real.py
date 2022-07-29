alpha_string = input("Alpha? ")
try:
    ALPHA_VALUE = int(alpha_string)
except:
    ALPHA_VALUE = float(alpha_string)


p_k_string = input("Policy restrictions? Type 0 < p^k <= 1 for yes, 'n' for no, or 'M' for Vaccine Incentives.")
vaccine_incentive = p_k_string == "M"
if vaccine_incentive:
    print("Running with vaccine incentives and no policy.")
    policy = False
else:
    try:
        p_k = float(p_k_string)
        policy = True
        print(f"Running with a policy p^k = {p_k}.")
    except:
        print("Running with no policy.")
        policy = False
strict_policy = policy and input("Strict policy? 'y' or 'n'.") == 'y'