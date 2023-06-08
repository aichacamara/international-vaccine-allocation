seir_opt Documentation	M. Veatch 10/31/22
-----------------------------------------------------------------------------
updated last: A. Cha 06/08/23

Setup/Installation 
--

### LINUX
Optimization uses Gurobi's python package. To use the Gurobi model install 
Gurobi onto your local machine here: 

https://www.gurobi.com/downloads/

This will require you to create an account and if you are planning on using 
a net licence create a file called `gurobi.lic` in an accessible location
and type the following into a terminal window:

> export GRB_LICENSE_FILE=[__path__]/gurobi.lic

Install the Gurobi Python Module via terminal 
> pip install gurobipy
or
> pip3 install gurobipy

### MAC/OSX
@TODO

### WINDOWS

Assuming you have Python,  

Download Gurobi from https://gurobi.com/downloads/  You need the Gurobi Optimizer.  You will need to register. 
Set up your Gurobi license (see below). 

Install the Gurobi Python Module via terminal 
> pip install gurobipy
or
> pip3 install gurobipy

Insert the following at the beginning of your python program 

  `import gurobipy as gp` 

  `from gurobipy import GRB`

and any other modules needed. 
Now the python program can call Gurobi. 

 

Gurobi site license 
--
Create a token server client license. It is a license file (_denoted with `.lic` extension_) created with any word editor 
```
Name: gurobi.lic  
Location: a default/accesible location on your computer, e.g., C:\Users\Mike.Veatch\Documents 

Contents: 
TOKENSERVER=gurobi.gordon.edu 
  or 
TOKENSERVER=172.27.43.55 
```

(note) this IP address may not be permanent. 
You must be on the Gordon network (on campus or VPN) to use this license. 

For testing, you may want to add a line that change the timeout in case the token server is unavailable. The default value is 30 seconds. 
    SERVERTIMEOUT=10 

-----------------------------------------------------------------------------

Overview of Code
--

The code is divided into five major functions. However, like the paper, most variables are global, so they only have one meaning and are used in several functions. A few things are passed as arguments so the same code can process them. Other variables are arguments for no reason. 

 

The functions appear first. Execution begins with the code that reads the CMD line (near end of the file). 

``` 
main() # initialization, search over lambda (outer loop), compute additional outputs after optimizing, and  
  # write outputs 
  call import_xml() to read inputs 
  initialize: output filename and directory, policy, state variables, LP variables, etc. 
  if simulate_only = 0 
    # outer loop (searches over values of lambda, which is a multiplier penalizing infections) 
    initialize for outer loop 
    iteration i = 0: call optimize_inner() 
      # not used in the minimization to avoid bias from initial policy  
    iteration i = 1: call optimize_inner() 
      # same lambda as iteration 0 but updated policy 
    iteration i = 2: call optimize_inner() 
    phase 1 iterate over i (change lambda in same direction until past min) 
      call optimize_inner() 
    phase 2 iterate over i (golden ratio search for lambda in an interval) 
      call optimize_inner() 

 

    # write outputs 
    write csv file (optimize) 
  if verbosity > 0: write input echo to output file 
  if simulate_only = 0: write outputs (optimize) 

 
  else   # simulate only 
    for each area 
      set policy to give priority to this area 
      call simulate() 
      record results of simulation 
     write results of simulation 
  for each area 
    write vaccinations from simulation 
 
  opt_inner() # inner loop over j: solves LP approximation, simulates, updates approximation, iterates 
  if outer iteration i = 0 
    initialize 
    call simulate() # first sim 
    record first sim results for output 
    set t_LP (time that infections are counted in LP) 
    initialize min solution and opt solution 
  if i > 0 and phase = 1 
    initial j = 0 solution is the best solution from i = 0 
    set alpha (transmission rate), V_cal (effective number infectious), V (vaccinations) using previous.  
  # Used by LP. 
  if phase = 2 
    initial j = 0 solution is the best solution from previous i 
    set alpha, V_cal, V using min. Used by LP. 
  initialize for j loop: j = 0, set t_LP 

 

  inner loop over j (do at least 3 iterations unless iter_lmt < 3) 
    increment j, LP_count 
    call solve_LP() # improving = 1: LP uses alpha, V_cal for the best policy found so far in this inner loop 
     # improving = 0: LP uses alpha, V_cal for the policy from the last LP 

    check that optimal solution was found 
    update the policy V using the LP solution V1 
    call simulate() 
    record sim results 

    compute dVmin, dVmax (change in V_cal, used to check convergence) 
    update “min” solution (j with best zNLP for this i) 
    update “opt” solution (best deaths so far) 
  
    if i = 0: record “min” solution in “prev” # used to initializing all i in phase 1 
    solve_LP() # LP objective and constraints, call Gurobi to solve LP 
    if there is a prior solution: save bases or solution for warm start 
    define objective 
    remove all constraints 
    define constraints 

    if there is a prior solution 
      update() # updates LP in Gurobi so we can refer to new constraints  
      set bases or solution to previous values for warm start 
    optimize # calls Gurobi to solve LP 

    simulate() simulates a policy (initial policy or policy found by LP) 
    set intial states to the values computed in main(), e.g., s[a, 0] = S0[a] 
    compute initial transmission rate alpha_0[a] 
    initialize t_sim (time variant appears in sim), variant_emerge flag 
 
    for t = 0,…, T-1 # loop over day 
      if variant_emerge: compute transmission rate alpha[a, t] for this t 
      else: transmission rate is constant
  
      compute V_cal[a, t], delta_E[a, t], and V_star[a, t] (vacc w/o reallocation) for this t  
  
      reallocate unused vaccinations to other areas, up to their vacc limits Wnew, in priority order 
        (store vacc w/o realloc in V_minus[a], then update V_star[a, t]) 
      compute next states using difference equations 
      check if variant emerges at this t, using nondonor areas 

    compute V_cal at final time T 
    if simulate_only: write csv file   

```
-----------------------------------------------------------------------------


Versions
--
## Current
seir_opt.py	optimizes using iterative LPs or simulates

### Deprecated Versions
seir_opt2.py	optimizes using iterative LPs or simulates

seirQP3.py	

optimizes using QP (at each lambda) or simulates. 

Recent changes: 
```
  MIPgap=.02
  time limit=600s for each QP
  record results of QP if status == GRB.TIME_LIMIT  
  clarified output: suboptimal QP results are counted if time limit reached,
	removed unused column
  console output is redirected to output file ..._con.out using sys.stdout
  output file name changed to include values of T and nu, e.g.,
	T2_T090_nu0.0.out
	To change which parameters are in name, look for fn_base
  csv file name changed from plot_....csv to ..._plot.csv
  ```
-----------------------------------------------------------------------------
Usage
--
If the input simulate_only = 1, several policies are simulated, one with priority given to each area. 
If simulate_only = 0, LPs are solved repeatedly searching for an optimal policy. A simulation is done before each LP is solved. 

To run, give the name of the FOLDER containing input files: 

> python seir_opt input_data	 

This does multiple runs, one for each file in this folder. Old versions did a single run, given the path to the input file. To read T2.xml in the folder "input_data":  

> python seir_opt input_data/T2.xml	 

----------------------------------------------------------------------------- 

### Controlling outputs 

The subfolder "output" is created. For each input file, seir_opt writes two or three files, e.g., 

```
  C2.3_nu0.0.out		output
  C2.3_nu0.0_con.out		console output, i.e., Gurobi
  C2.3_nu0.0_plot.csv		csv file for time plots with R
```

The output file names are the input file name, followed by the value of the input "nu".  

Note that output files are overwritten if the same input and nu are used. 

To change which input is added to the file name, edit "fn_base = ..."  

To direct Gurobi output to console, remove "sys.stdout = ..." 

To turn off console (Gurobi) output, change v.Params.LogToConsole from 1 to 0 

## The simulate_only and verbosity input control what appears in the output file: 

### Optimization 
```
verbosity >= 0: results for first sim, optimal, "min" (best LP for last lambda) 
verbosity >= 1: echo inputs, vaccinations by day/area 
verbosity >= 2: outer and inner loop showing progress of algorithm 
```
### Simulate only 
```
verbosity >= 0: sim results for each priority policy 
verbosity >= 1: echo inputs, vaccinations by day/area 
```
-----------------------------------------------------------------------------
Input Format
--

All inputs and computed constants are global, defined in:import_xml (denoted  global:import_xml) 

### area_data 
```
name[a]				name of area 
N[a]				population by area  
rho_V[a]			initial proportion vacc by area 
rho_I_N[a]			initial cases per day by area 
delta_r[a]			testing effect on rate out of I by area 
gamma[a]			behavior infection multiplier by area 
rho[a]				proportion willing to be vaccinated by area 
priority			list of areas by priority for initial policy (high priority first) 
donor				donor area name 
m				variant area name “mutation” 
n				infection-days till new variant 
```
### Computed constants from area data 
```
A[ ]				set of area names		 
A_D[ ]				set of nondonor area names 
n_a				number of areas 
```
 

### scenario_data 
```
T				time horizon (days)			 
B_0				vaccine available day 0 
nu				weight for non-donor deaths in objective 
v_u				upper limit on proportion infectious due to behavior. No behavior: 0 (local) 
p_k				max prop of vacc alloc to donor area 
r_I				rate out of state E into I 
r_0				rate out of state I w/o testing (local)		 
p_D				P(death | infected, not vacc) 
p_V_D				P(death | infected, vacc) 
a_0				initial infection rate	 
delta_a				change in infection rate for variant 
p_e				proportion of transmission from a vaccinated person  
p_r				proportion of transmission to a vaccinated person  
L				lag for variant to reach other areas (days)			 
T_D				time for variant to dominate (days)		 
p				Proportion of people in state I that hav the new variant when it is introduced 
b_arr[t]			vaccine available as prop of that avail day 0.  
     # If not in file, use the default value of 1, so that B(t) = B_0 for all t. 

``` 

### Computed constants from scenario data 
```
v_l[a]				lower proportion for behavior dynamics by area			 
v_u[a]				upper proportion for behavior dynamics by area 
g[a]				offset for behavior dynamics by area 
k				decay rate of mutation model  
B[t]				vaccine available by day 
r_d[a]				rate out of I by area 
t_switch[]			switchover policy datelines (optional)
split				splits the vaccine allocation by day (optional) 
```
 

  

### params 
```
simulate_only			1 if simulate only, 0 if optimize 
lambda				Lagrange multiplier for infection (initial value) 
phi				exploration multiplier for λ 
epsilon_0			exploration tolerance for LP 
delta_I				termination tolerance for LP 
delta				termination tolerance for λ 
beta				exploration convergence parameter for LP	 
iter_lmt			iteration limit for LP 
iter_lmt_search 		iteration limit for λ 			 
dT				days after t_n[0] in Lagrangian 
verbosity 			verbosity of output					 
```
 

### Hardcoded constants (parameters) 
```
T0 = 3				days at end of scenario with no vaccinations (b/c they don’t matter) 
improving = 0			1 if LP uses alpha, V_cal for the best policy found so far in this inner loop 
				0 if LP uses alpha, V_cal for the policy from the last LP 
```
-----------------------------------------------------------------------------
Outputs
--
The subfolder "output" is created. For each input file, seir_opt2 writes two or threee files, e.g.,

  C2.3_tsw180.out		output
  C2.3_tsw180_con.out		console output, i.e., Gurobi
  C2.3_tsw180_plot.csv		csv file for time plots with R

seir_QP3 appends "_QP" to the names.
To direct Gurobi output to console, remove "sys.stdout = ..."
To turn off console (Gurobi) output, change v.Params.LogToConsole from 1 to 0

The verbosity input controls what goes in the output file:

verbosity >= 0: echo inputs, first sim, optimal, "min" (best LP for last lambda)
verbosity >= 1: vaccinations by day/area
verbosity >= 2: outer and inner loop showing progress of algorithm

For QP, the outputs are a subset of this:

verbosity >= 0: echo inputs, first sim, optimal
verbosity >= 1: vaccinations by day/area
verbosity >= 2: outer loop showing progress of algorithm

When simulate_only = 1, the initial policy is simulated and a csv file is written. Other output for this option hasn't been done. 
----------------------------------------------------------------------------
Plots
--
Open time_plots.Rmd in RStudio. It uses the packages tidyverse, tidyr, dplyr. If you haven't installed these, add the install statements 
```
  install.packages("tidyverse")
  install.packages("tidyr")
  install.packages("dplyr")
```
before the library statements.

Change the read.csv statement to give the path to your csv file, e.g.,

sim1 = read.csv("~/research/vaccine/code/output/C2.3_plot.csv")

The first chunk creates large plots, two for each area. Also run the second chunk to get smaller plots that are readable in black and white. 

### Notes 
1. You can open the data file with Import/From text (base) in the Environment window, but then you need to rename the data frame. Look in the environment to see the data frame name. Uncomment the statement, e.g.,

sim1 = C2.3_plot

2. To change which variables are plotted, edit the filter statement, e.g.,
```
filter(state == "Infectious" | state == "Deaths"   | state == "Cases" 
         | state == "Vaccinations" | state == "Susceptible") |>
```
There must be |> (the native pipe) at the end. To see the list of variables, after running the chunk look at the state column of the data frame sim_long. It includes all the state variables, plus some that were computed in R:
```
mutate(Infectious = I + IV)
mutate(Cases = E + EV + I + IV + H + R + D)
mutate(Vaccinations = cumsum(V)
```
Note that Cases, Recovered, Deaths, and Vaccinations are cumulative, while others are not. Omitting variables that take on large values, such as Cases, will zoom in on the smaller variables. To rename variables (e.g., to change "SV" to "Susceptible Vacc"), edit the recode statement:

sim_long$state = recode(sim_long$state, "H" = "Hospitalized", "R"="Recovered", 
  "D" = "Deaths", "S" = "Susceptible") 

The x-axis has tick marks every 10 days; see the comments to change this.

All areas have dashed vertical lines at
  t_n = day variant appears
  t_n + L = day variant spreads to other areas
Each area should only have one line, but I don't know how to do different lines on different plots.


  
