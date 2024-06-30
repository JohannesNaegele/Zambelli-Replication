""" Compute the envelope for all countries at once. """
function envelope(A_tensor, B_tensor, l_tensor; stepsize = 0.01, extend = true, calc_all = true)

    # Bring data in matrix shape
    A = reduce(hcat, A_tensor[:, :, i, j]' * Diagonal(safe_divide.(1.0, B_tensor[:, i, j])) for i in axes(A_tensor, 3) for j in axes(A_tensor, 4))
    n_goods = size(A, 1)
    B = reduce(hcat, I(n_goods) for i in axes(B_tensor, 2) for j in axes(B_tensor, 3))
    l = reduce(vcat, l_tensor[:, i, j] .* safe_divide.(1.0, B_tensor[:, i, j]) for i in axes(l_tensor, 2) for j in axes(l_tensor, 3))'
    d = axes(A, 1) .== 1

    # Preallocate solver variables
    lb = zeros(axes(A, 2))
    l_trunc = zeros(n_goods)'
    C_trunc = zeros(n_goods, n_goods)

    # Create LP solvers
    solver = optimizer_with_attributes(HiGHS.Optimizer, "log_to_console" => false)
    model_x = create_intensities_r_solver(solver, l, A, B, d, lb)
    model_x_trunc = create_intensities_solver(
        solver, l_trunc, C_trunc, ones(n_goods), lb[1:n_goods]
    )
    model_p = create_prices_solver(solver, l_trunc, C_trunc, ones(n_goods), lb[1:n_goods])
    R = maximum(compute_R.(real_eigvals.(A[:, axes(A, 1) .+ (i - 1) * n_goods] for i in 1:div(size(A, 2), n_goods))))
    df_q, profit_rates_to_names, profit_rates, switches = compute_envelope(
        A = A,
        B = B,
        l = l,
        d = d,
        R = R,
        stepsize = stepsize,
        model_intensities = model_x,
        model_intensities_trunc = model_x_trunc,
        model_prices = model_p,
        effects_sectors = 1:n_goods,
        extend = extend,
        calc_all = calc_all
    )
    n_switches = length(switches["technology"])
    # labeled_tech = Matrix{String}(undef, n_goods, n_switches + 1)
    for i in eachindex(switches["technology"])
        country_index = div.(switches["technology"][i][1].second .- 1, n_goods) .+ 1
        # labeled_tech[:, i] = map(j -> ids[j], country_index)
    end
    # FIXME: n_switches could be zero
    country_index = div.(switches["technology"][n_switches][2].second .- 1, n_goods) .+ 1
    # labeled_tech[:, n_switches + 1] = map(j -> ids[j], country_index)
    return Dict(
        "A" => A,
        "B" => B,
        "l" => l,
        "d" => d,
        "intensities" => df_q,
        "names" => profit_rates_to_names,
        "profit_rates" => profit_rates,
        "switches" => switches,
        # "labeled_technology" => labeled_tech,
        "max_R" => 0.0
    )
end