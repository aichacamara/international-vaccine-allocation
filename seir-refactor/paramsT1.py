# Run parameters for seir_opt script

verbosity = 3           # 0 least verbose output ... 3 most verbose 
simulate_only = True
realloc_flag = False    # proportional reallocation of unused vacc to non-donors in sim
non_donor_deaths = False # include non-donor deaths in objective 
p_k = 0.75       # maximum proportion of available vacc used in donor area

lambda_0 = 1    # initial Lagrange multiplier for infection
phi = 2         # exploration multiplier for lambda
epsilon_0 = 100 # exploration tolerance for I 
delta_I = 100   # termination tolerance for I
delta = 0.1     # termination tolerance for lambda
beta = .9       # convergence parameter
iter_lmt = 20   # iteration limit for I
iter_lmt_search = 10  # iteration limit for lambda

