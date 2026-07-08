using JLD2
using Statistics
using StatsPlots
using GLM, DataFrames
using CMPSKit
using LaTeXStrings

@load "benchmark_mu=2.0_c=1.0_c12=-0.5_D=8" Ψs Es histories
sizes = length.(histories) ./ 2

@load "benchmark_mu=2.0_c=1.0_c12=-0.5_D=8_no_prec" Ψ2s E2s histories2
sizes2 = length.(histories2) ./ 2
# @load "benchmark_tens_prod_2.0_c=1.0_c12=-0.5_D=4" Ψs Es histories
# sizes2 = length.(histories) ./ 2
# @load "benchmark_penal_2.0_c=1.0_c12=-0.5_D=8" Ψs Es histories
# sizes3 = length.(histories) ./ 2

# bond_dims = [4,6,8]

# Medians
# med_with    = [median(iters_with[d])    for d in bond_dims]
# med_without = [median(iters_without[d]) for d in bond_dims]

# Linear fit function
# linfit(x, y) = begin
#     xf, yf = float.(x), float.(y)
#     xm, ym = mean(xf), mean(yf)
#     m = sum((xf .- xm) .* (yf .- ym)) / sum((xf .- xm).^2)
#     b = ym - m * xm
#     m, b
# end

# m_no, b_no = linfit(bond_dims, med_without)
# m_pc, b_pc = linfit(bond_dims, med_with)

# p = plot(xlabel="Bond dimension", ylabel="Median iterations",
#          title="Median iterations vs bond dimension", legend=:topleft)

# scatter!(p, bond_dims, med_without; marker=:circle, label="No precond.")
# plot!(p, bond_dims, m_no .* bond_dims .+ b_no;
#       label="Trend no precond. (slope=$(round(m_no, sigdigits=3)))", linestyle=:dash)

# scatter!(p, bond_dims, med_with; marker=:diamond, label="Precond.")
# plot!(p, bond_dims, m_pc .* bond_dims .+ b_pc;
#       label="Trend precond. (slope=$(round(m_pc, sigdigits=3)))", linestyle=:dot)

# display(p)
# savefig(p, "median_iters_vs_bond_dim.png")

# @load "Iterations_with_preconditioner_D10" histories

# @show sizes, sizes2, sizes3

# # println("number of iterations: ", sizes)
# println("mean: ", mean(sizes))
# println("mean: ", mean(sizes2))

# println("std: ", std(sizes))
# println("std: ", std(sizes2))

# histogram(
#            sizes;
#            bins=50,
#            xlabel="Number of iterations",
#            ylabel="Frequency",
#            title="Histogram of number of iterations with preconditioner",
#            legend=false,
#        )

# savefig("hist_w_preconditioner_D8.pdf")

# histogram(
#            sizes2;
#            bins=50,
#            xlabel="Number of iterations",
#            ylabel="Frequency",
#            title="Histogram of number of iterations without preconditioner",
#            legend=false,
#        )

# savefig("hist_wo_preconditioner_D8.pdf")

# boxplot(
#            sizes;
#            orientation=:horizontal,
#            xlabel="Number of iterations",
#            title="Boxplot of number of iterations with preconditioner",
#            legend=false,
#        )

# savefig("box_w_preconditioner_D8.pdf")

# boxplot(
#            sizes2;
#            orientation=:horizontal,
#            xlabel="Number of iterations",
#            title="Boxplot of number of iterations without preconditioner",
#            legend=false,
#        )

# savefig("box_wo_preconditioner_D8.pdf")

cap=40001

# Convert to "time" + "event" (1=converged, 0=censored)
times_pre = min.(sizes, cap)
events_pre = times_pre .< cap                       # Bool
times_tens  = min.(sizes2, cap)
events_tens = times_tens .< cap                         # Bool
# times_penal = min.(sizes3, cap)
# events_penal = times_penal .< cap

"""
Kaplan–Meier survival curve S(t)=P(not converged by t) for right-censored data.

Returns:
  tpoints::Vector{Float64}, Spoints::Vector{Float64}
Suitable for step plotting (st=:steppost).
"""
function km_curve(times::AbstractVector{<:Real}, events::AbstractVector{Bool})
    t = Float64.(times)
    e = events

    # unique event times (exclude censored)
    event_times = sort(unique(t[e]))

    S = 1.0
    tpoints = Float64[0.0]
    Spoints = Float64[1.0]

    for τ in event_times
        at_risk = count(>=(τ), t)                   # number with time >= τ
        d = count(i -> (t[i] == τ) && e[i], eachindex(t))  # events at τ
        S *= (1.0 - d / at_risk)
        push!(tpoints, τ)
        push!(Spoints, S)
    end

    return tpoints, Spoints
end

tpre, Spre = km_curve(times_pre, events_pre)
ttens, Stens = km_curve(times_tens, events_tens)
# tpen,  Spen  = km_curve(times_penal,  events_penal)

p = plot(tpre, Spre;
    st = :steppost,
    xlabel = "Iterations",
    ylabel = "S",
    # xlims = (0, cap/2),
    ylims = (0, 1),
    label = "preconditioned MDM⁻¹",
    legend = :topright
)

plot!(p, ttens, Stens; st = :steppost, label = "unpreconditioned MDM⁻¹")
# plot!(p, tpen, Spen; st = :steppost, label = "penaly term")

display(p)
savefig("Kaplan-Meier_D8_precondition.pdf")

println("Preconditioned censored: ", count(!, events_pre), " / ", length(events_pre))
# println("Penalty censored: ", count(!, events_tens), " / ", length(events_tens))
# println("Penalty censored: ", count(!, events_penal), " / ", length(events_penal))