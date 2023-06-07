seir_opt Documentation	M. Veatch 10/31/22
-----------------------------------------------------------------------------
Setup/Installation 

LINUX

Optimization uses Gurobi's python package. To use the Gurobi model install 
Gurobi onto your local machine here: 

https://www.gurobi.com/downloads/

This will require you to create an account and if you are planning on using 
a net licence create a file called `gurobi.lic` in an accessible location
and type the following into a terminal window:

export GRB_LICENSE_FILE=[__path__]/gurobi.lic

from there:
pip install gurobipy
or
pip3 install gurobipy

MAC/OSX
@TODO

Windows
@TODO

-----------------------------------------------------------------------------
Versions

seir_opt2.py	optimizes using iterative LPs or simulates

seirQP3.py	
optimizes using QP (at each lambda) or simulates. 
Recent changes: 
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
-----------------------------------------------------------------------------
Usage

To run seir_opt2 or seir_QP3, give the name of the FOLDER containing input files:

> python seir_opt2 input_data	

This does multiple runs, one for each file in this folder. Both programs use the same input file, though not all inputs are used by QP.

To run old programs, such as seir_QP2, give the path to the input file. To read T2.xml in the folder "input_data": 

> python seir_QP2 input_data/T2.xml	
-----------------------------------------------------------------------------
Outputs

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

Open time_plots.Rmd in RStudio. It uses the packages tidyverse, tidyr, dplyr. If you haven't installed these, add the install statements 
  install.packages("tidyverse")
  install.packages("tidyr")
  install.packages("dplyr")
before the library statements.

Change the read.csv statement to give the path to your csv file, e.g.,

sim1 = read.csv("~/research/vaccine/code/output/C2.3_plot.csv")

The first chunk creates large plots, two for each area. Also run the second chunk to get smaller plots that are readable in black and white. 

Notes 
1. You can open the data file with Import/From text (base) in the Environment window, but then you need to rename the data frame. Look in the environment to see the data frame name. Uncomment the statement, e.g.,

sim1 = C2.3_plot

2. To change which variables are plotted, edit the filter statement, e.g.,

filter(state == "Infectious" | state == "Deaths"   | state == "Cases" 
         | state == "Vaccinations" | state == "Susceptible") |>

There must be |> (the native pipe) at the end. To see the list of variables, after running the chunk look at the state column of the data frame sim_long. It includes all the state variables, plus some that were computed in R:

mutate(Infectious = I + IV)
mutate(Cases = E + EV + I + IV + H + R + D)
mutate(Vaccinations = cumsum(V)

Note that Cases, Recovered, Deaths, and Vaccinations are cumulative, while others are not. Omitting variables that take on large values, such as Cases, will zoom in on the smaller variables. To rename variables (e.g., to change "SV" to "Susceptible Vacc"), edit the recode statement:

sim_long$state = recode(sim_long$state, "H" = "Hospitalized", "R"="Recovered", 
  "D" = "Deaths", "S" = "Susceptible") 

The x-axis has tick marks every 10 days; see the comments to change this.

All areas have dashed vertical lines at
  t_n = day variant appears
  t_n + L = day variant spreads to other areas
Each area should only have one line, but I don't know how to do different lines on different plots.


  
