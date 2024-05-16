using MAT
using Sraffa
using JuMP
using HiGHS
using LinearAlgebra
using Plots
using PlotlyJS

include("preprocess.jl")
include("envelope.jl")

@time results = envelope(A_data, B_data, L_data, extend=true, stepsize = 0.0001)
switch_info(results)

theme(:bright)

# Fast visualizations
gr(size=(900, 600))
wagecurves(results, switches=false, legend=false, linewidth=1, xlim=(0, 2.7))
env(results, xlim=(0, 2.7))

savefig("wage_curves.pdf")

# Slow visualizations
plotlyjs(size=(1200, 1400))
wagecurves(results, xlim=(0, 2.7))
savefig("wage_curves.html")