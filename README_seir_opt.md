seir_opt Documentation	M. Veatch 6/1/23
-----------------------------------------------------------------------------
Current version:

seir_opt.py	optimizes using iterative LPs or simulates
-----------------------------------------------------------------------------
Usage

To run, give the name of the FOLDER containing input files:

> python seir_opt input_data	

This does multiple runs, one for each file in this folder. Old versions did a single run, given the path to the input file. To read T2.xml in the folder "input_data": 

> python seir_opt input_data/T2.xml	
-----------------------------------------------------------------------------
Outputs

The subfolder "output" is created. For each input file, seir_opt writes two or three files, e.g.,

  C2.3_nu0.0.out		output
  C2.3_nu0.0_con.out		console output, i.e., Gurobi
  C2.3_nu0.0_plot.csv		csv file for time plots with R

The output file names are the input file name, followed by the value of the input "nu". 
To change which input is added to the file name, edit "fn_base = ..." 
To direct Gurobi output to console, remove "sys.stdout = ..."
To turn off console (Gurobi) output, change v.Params.LogToConsole from 1 to 0

The verbosity input controls what goes in the output file:

verbosity >= 0: echo inputs, first sim, optimal, "min" (best LP for last lambda)
verbosity >= 1: vaccinations by day/area
verbosity >= 2: outer and inner loop showing progress of algorithm

When simulate_only = 1, the initial policy is simulated. 
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
