using JuMP
using HiGHS
using Random

# For reproducible random numbers
Random.seed!(1234)

# Initialize some exemplary data
n_rows, n_cols = 4, 3
A = rand(n_rows, n_cols) ./ 3 # technology matrix
B = [ # output matrix, each column is a sector
    1.0 0.0 0.0;
    0.0 1.0 0.0;
    0.0 0.0 1.0;
    1.0 0.0 0.0;
]
l = [0.4, 0.3, 0.2, 0.1] # labor input
r = 0.1 # profit rate
C = B - (1 + r) * A
d = ones(n_cols) # numeraire

lp_solver = Model(HiGHS.Optimizer) # initialize solver
@variable(lp_solver, x[1:n_rows] ≥ 0.0) # define intensities (in the implementation a column vector)
@constraint(lp_solver, x' * C .≥ d) # define the constraint
@objective(lp_solver, Min, x' * l) # define our minimization objective
optimize!(lp_solver) # solve the LP

x_trunc = value.(x) # get the intensities -> first row is not used
val = objective_value(lp_solver) # get the objective value
w_max = 1 / val # compute the maximum wage

# Truncate the technologies: We found the techniques on the envelope for profit rate r
C_trunc = C[2:end, :]
l_trunc = l[2:end]

# Sanity check: Recompute the resulting wage
1 / (d' * inv(C_trunc) * l_trunc)