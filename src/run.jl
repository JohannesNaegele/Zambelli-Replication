using MAT
using Sraffa
using JuMP
using HiGHS
using LinearAlgebra
using Plots
using Logging

include("preprocess.jl")
include("envelope.jl")

# Calculate the envelope with LP algorithm
@time results_lp = envelope(A_data, B_data, L_data, stepsize = 0.0001)
# Calculate the envelope with (mostly) piecewise switches
@time results_fine = envelope(A_data, B_data, L_data, stepsize = 0.0001, calc_all=false)

# Get info on the results
switch_info(results_lp)
switch_info(results_fine)

# Convert Zambelli's results in our dict format for plotting
b_mod = map(x -> mod(Int(x), 30), vars["SSS_1995_2011"][1:end-1, :])
b_div = map(x -> div(Int(x), 30), vars["SSS_1995_2011"][1:end-1, :])
zamb = ((b_mod .- 1) * 31 + b_div * 31 * 30) .+ (1:31)

# We copy our results and just adjust the chosen technology
results_zamb = deepcopy(results_lp)
tech = [
    [
        vars["SSS_1995_2011"][end, i] + 0.01 => zamb[:, i + 1],
        vars["SSS_1995_2011"][end, i + 1] => zamb[:, i + 2]
    ]
    for i in 1:(size(zamb, 2) - 2)
]
pushfirst!(tech, [
        0.0 => zamb[:, 1],
        vars["SSS_1995_2011"][end, 1] => zamb[:, 2]
    ]
)
results_zamb["switches"]["technology"] = tech

# Set plotting theme
theme(:bright)

# Fast visualizations
gr(size=(900, 600))
wagecurves(results_lp, switches=false, legend=false, linewidth=1, xlim=(0, 2.7))
Plots.savefig("envelope.pdf")
env(results_lp, xlim=(0, 2.7), label="LP");
env!(results_zamb, xlim=(0, 2.7), label="VFZ");
title!("Compare LP solution to VFZ")
Plots.savefig("wage_curves.pdf")

# Slow (interactive) visualizations
plotlyjs(size=(1200, 1400))
wagecurves(results_lp, xlim=(0, 2.7))
savefig("wage_curves.html")