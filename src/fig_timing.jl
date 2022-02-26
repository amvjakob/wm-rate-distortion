# fig_timing.jl
# Figure for timing effects of rate-distortion theory of working memory.

include("utils.jl")
include("pyplot.jl")

using PyPlot, LaTeXStrings

"Plot figure for set size effects."
function plot_timing_error!(ax, timing_conditions, errors; ylim = (0, 1.0))
    if isnothing(ax)
        ax = subplot(111)
    end

    for (i, condition) in enumerate(timing_conditions)
        sns.kdeplot(
            vcat(errors[i]...);
            ax = ax,
            clip = (-π, π),
            label = condition.name,
        )
    end
    ax.set_xlim(-π - 0.2, π + 0.2)
    ax.set_xticks([-π, 0, π])
    ax.set_xticklabels([L"-\pi", "0", L"\pi"],)

    ax.set_ylim(ylim...)

    ax.set_xlabel("Error")
    ax.set_ylabel("Probability density")

    return ax
end

"Plot figure for set size effects."
function plot_timing_variance!(ax, timing_conditions, errors; ylim = (0, 1.0))
    if isnothing(ax)
        ax = subplot(111)
    end

    # compute variance and std error
    x = collect(1:length(errors))
    y = map(i -> mean(cvar.(errors[i])), x)
    e = map(i -> stde(cvar.(errors[i])), x)

    for (i, condition) in enumerate(timing_conditions)
        ax.errorbar(x, y; yerr = e, fmt = "o-", color = :black, markersize = 5)
    end
    ax.set_xticks(x)
    ax.set_xticklabels(map(c -> c.name, timing_conditions))

    ax.set_ylim(ylim...)

    ax.set_xlabel("Condition")
    ax.set_ylabel("Variance")

    return ax
end