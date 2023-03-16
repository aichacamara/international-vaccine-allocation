# Warm start tests

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd
import math
import time

### Numerical indexing: size m x n, with ranges [0,m-1] x [0,n-1]. max c^Ty s.t. Ay <= b y >= 0

M = gp.Model("Numerical Indexing")
m = 2 # rows
n = 3 # columns
A = np.array([[1, 2, 3], # A is a matrix. Refer to as A[0,0], etc.
    [1, 1, 0]])
b = [4, 1] # b and c are lists
c = [1, 1, 2]
y = M.addVars(n, name="y") # y[0], ..., y[n-1]
M.setObjective(sum(c[j]*y[j] for j in range(n) ), GRB.MAXIMIZE)
M.addConstrs((sum(A[i,j]*y[j] for j in range(n) ) <= b[i] for i in range(m)), "constr") 
M.write('M.lp') # Write model in LP format
M.optimize()

# Warm start METHOD 1 PStart/DStart. After R&R all constraints, only DStart can be used. Took 2 iter.

psol = {j: M.getVars()[j].x for j in range(M.NumVars)} # Save primal sol
dsol = {i: M.getConstrs()[i].Pi for i in range(M.NumConstrs)} # Save dual sol
M.remove([constraint for constraint in M.getConstrs()]) # Remove all constraints
M.addConstrs((sum(A[i,j]*y[j] for j in range(n) ) <= b[i] for i in range(m)), "constr")
M.update() # Must update to refer to new constraints
# Must set all PStart OR all DStart for it to be used. If set both, one will be used or may "push" & use both
for j in range(M.NumVars):
    M.getVars()[j].PStart = psol[j]
for i in range(M.NumConstrs):
    M.getConstrs()[i].DStart = dsol[i]
print("Rerun after R&R constraints, PStart, DStart")
M.optimize()

# Warm start METHOD 2 VBasis/CBasis. After R&R all constraints, only VBasis can be used? DOESN'T SAY, but 0 iter.

vbas = {j: M.getVars()[j].VBasis for j in range(M.NumVars)} # Save basis for variables
cbas = {i: M.getConstrs()[i].CBasis for i in range(M.NumConstrs)} # Save basis for constraints (dual)
M.remove([constraint for constraint in M.getConstrs()]) # Remove all constraints
M.addConstrs((sum(A[i,j]*y[j] for j in range(n) ) <= b[i] for i in range(m)), "constr")
M.update() # Must update to refer to new constraints
# Must set all VBasis OR all CBasis for it to be used. 
for j in range(M.NumVars):
    M.getVars()[j].VBasis = vbas[j]
for i in range(M.NumConstrs):
    M.getConstrs()[i].CBasis = cbas[i]
print("Rerun after R&R constraints, VBasis, CBasis")
M.optimize()

""" Not used

#M.Params.UpdateMode = 0 ## test. Delays updates of vars/constraints
M.remove(M.getConstrByName("constr[0]")) # ref constraint by name
M.remove(M.getConstrs()[0]) # ref constraint by number
M.addConstr((y[0] + 2*y[1] + 3*y[2] <= 4), "constr[0]") # Placed at end of list of constraints
M.getConstrByName("constr[0]").DStart = dsol[0] # Constraints were reordered, so ref by name
M.getConstrByName("constr[1]").DStart = dsol[1]
print("#0 Dual value", M.getConstrs()[0].Pi) #check NOT WORKING: can't retrieve attrib
print("#1 Dual value", M.getConstrs()[1].Pi)
print("#0 DStart:", M.getConstrs()[0].DStart) 
print("#1 DStart:", M.getConstrs()[1].DStart)

for j in range(n):
    print(y[j].varName, '=', y[j].x)  # y is the LP variable, ".varName" gives its name, ".x" gives its optimal value
print('Obj: %g' % M.ObjVal)
"""