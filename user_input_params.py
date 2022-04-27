T = int(input("T? "))
t_N = T
alpha_string = input("Alpha? ")
num_areas = int(input("How many areas do you want, 1-4? "))
vax_days = int(input("How many days do you want to give vaccinations?"))
try:
    ALPHA_VALUE = int(alpha_string)
except:
    ALPHA_VALUE = float(alpha_string)

vaccine_list_initial = [7350/4]*vax_days+[0]*(T-vax_days)
B_t = {i: vaccine_list_initial[i] for i in range(len(vaccine_list_initial))}
areas = ["area1", "area2", "area3", "area4"][:num_areas]
#if num_areas == 2:
#    areas = ["area1", "area4"]
#    print("Using areas 1,4")

