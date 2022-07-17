# LP indexing examples

import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd

### Numerical indexing: size m x n, with ranges [0,m-1] x [0,n-1]. max c^Ty s.t. Ay <= b y >= 0

M = gp.Model("Numerical Indexing")

m = 2 # rows
n = 3 # columns
A = np.array([[1, 2, 3], # A is a matrix, not a list of lists
    [1, 1, 0]])

b = [4, 1] # b and c are lists, not row or column vectors
c = [1, 1, 2]

# For sparse A
# A = {i,j: 0 for i in range(m) for j in range(n)} # Like a loop, sets A = 0. Not working
A[0,0] = 1
A[0,1] = 2
A[0,2] = 3
A[1,0] = 1
A[1,1] = 1

# To assign A from other data, say D
# A = {{i,j: 7*D[i,j] for i in range(m-1) for j in range(n-1)}}

y = M.addVars(range(n), name="y") # y[0], ..., y[n-1]
                                    # name= "y" is optional. If omited, name is taken from "y ="

M.setObjective(sum(c[j]*y[j] for j in range(n) ), GRB.MAXIMIZE)

M.addConstrs((sum(A[i,j]*y[j] for j in range(n) ) <= b[i] for i in range(m)), "constr") # constraint must be in () if name argument used

M.write('M.lp') # Write model in LP format

M.optimize()

for j in range(n):
    print(y[j].varName, y[j].x)  # y is the variable name, ".x" references the optimal solution of the variable
print('Obj: %g' % M.ObjVal)

# Error handling (from matrix1.py)
# except gp.GurobiError as e:
#    print('Error code ' + str(e.errno) + ": " + str(e))

# except AttributeError:
#    print('Encountered an attribute error')

### Set indexing: see diet.py, which also has some shortcut commands


### Doubly indexed variables, bounds, binary variables: Assignment problem  with set indexing. 
#   Taken from Gurobi "intro_to_modeling" example
#   max sum over i,j of C[i,j]x[i,j] s.t. sum over i of x[i,j] = 1 for all j, sum over j of x[i,j] <= 1 for all j, x binary

# The LP formulation adds the bounds 0 <= x <= 1

a = gp.Model("Assignment problem")

I = ['Carlos', 'Joe', 'Monika'] # Resources
J = ['Tester', 'JavaDeveloper', 'Architect'] # Jobs

# Matching rewards. See gp.multidict for shortcuts entering double indices or multiple values for the same indices.
C = {
  ('Carlos','Tester'): 53, # C is a list of lists.
  # Refer to it as C['Carlos','Tester']. But must refer to nested lists with numeric indices as c[0][0], which also works for arrays.
  ('Carlos','JavaDeveloper'): 27,
  ('Carlos','Architect'): 13,
  ('Joe','Tester'): 80,
  ('Joe','JavaDeveloper'): 47,
  ('Joe','Architect'): 67,
  ('Monika','Tester'): 53,
  ('Monika','JavaDeveloper'): 73,
  ('Monika','Architect'): 47}

x   = a.addVars(I, J, lb = 0,ub = 1, name='x') # 0 <= x[i,j] <= 1,  i in I, j in J. Continuous variables.
# x = a.addVars(I, J, vtype = GRB.BINARY, name='x') # x[i,j] = 0,1, i in I, j in J. Binary variables.

a.setObjective(sum(C[i,j]*x[i,j] for i in I for j in J), GRB.MAXIMIZE)

a.addConstrs((x.sum('*',j) == 1 for j in J), name='job') # constraint must be in () if name argument used
a.addConstrs((x.sum(i,'*') <= 1 for i in I), name='resource')

a.write('assignment_problem.lp')
a.optimize()

for v in a.getVars(): # Get variable names, print only positive variables
    if v.x > 1e-6:
        print(v.varName, v.x)
#print('x[Carlos,Tester]:' % x['Carlos','Tester'].X)  # x is the variable name, ".X" references the optimal solution of the variable
print('Obj: %g' % a.ObjVal)