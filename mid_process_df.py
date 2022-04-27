import pandas as pd

vax_budget_df = pd.read_csv("vax_budget.csv")

B_t = {i-342: vax_budget_df["budget"][i-343] for i in vax_budget_df["days_since_jan_2020"] if i >= 517}
B_t = {i-174: vax_budget_df["budget"][i]/10 for i in range(len(vax_budget_df["days_since_jan_2020"])) if i >= 174}
#B_t = {i: 100 for i in range(len(vax_budget_df))}
T = len(B_t)-1
#T = 300
t_N = T
vaccine_params_df = pd.read_csv("mid_vaccine_params.csv").drop('Unnamed: 0', axis = 1)


