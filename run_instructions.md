# Requirements

Must have python installed along with the gurobipy package with gurobi configured on your computer.

To install gurobipy, run command "pip3 install --user gurobipy"
All other packages are built into python and should not need additional installing.
If you do run into an issue with importing a python module, run "pip3 install --user {MODULE NAME}"

# Running the code

To run the code, in the international_vaccine_allocation directory run command
"python3 seir_opt.py {INPUT FILE LOCATION}"

Input files have been stored in the "input_data" folder in the parent directory of the repository, so an example run command would be
"python3 seir_opy.py input_data/seir_dataT1.1.xml"

After the code runs, depending on whether it was optimizing or simulating, the program will create a directory called
"optimization_output" or "simulation_output" accordingly. The primary outputs in these folders are the outputs from the run as specified in the project, as well as a CSV for plotting. These folders are deleted between runs so to keep outputs from different runs you will have to move them to a different directory.