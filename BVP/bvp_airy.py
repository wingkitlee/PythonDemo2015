import numpy as np
import scikits.bvp_solver as bvp
import matplotlib.pyplot as plt

"""
To solve u''(x) - lam*x*u(x) = 0, -1<x<1
u(-1) = u(1) = 0

u' = v
v' = lam*x*u
"""

def function( X, Y, P):
	dY = np.zeros((2))
	dY[0] = Y[1]
	dY[1] = P*X*Y[0]
	return dY

def boundary_conditions(Ya, Yb, P):
	BCa = np.zeros((2))
	BCb = np.zeros((1))

	BCa[0] = Ya[0] # u(-1)= 0
	BCa[1] = Ya[1]-1.0 #u'(-1) = 1
	BCb[0] = Yb[0] # u(1) = 0

	return BCa, BCb

problem_definition = bvp.ProblemDefinition(num_ODE = 2,
	num_parameters = 1,
	num_left_boundary_conditions = 2,
	boundary_points = (-1.0, 1.0),
	function = function,
	boundary_conditions = boundary_conditions)

sol = bvp.solve(bvp_problem = problem_definition,
	solution_guess = [1.0      ,1.0      ],
	parameter_guess = [168.0    ],
	trace=1)

print sol.parameters

fig = plt.figure()
x = sol.mesh
y = sol.solution

plt.plot(x, y[0,:], 'k-', lw=2)
plt.plot(x, y[0,:], 'bo', markersize=5)
plt.show()
