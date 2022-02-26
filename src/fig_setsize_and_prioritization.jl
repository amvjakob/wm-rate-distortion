# fig_setsize.jl
# Figure for set size and stimulus prioritization effects of 
# rate-distortion theory of working memory.

include("utils.jl")
include("pyplot.jl")

using PyPlot, LaTeXStrings

"Plot figure for set size effects."
function plot_setsize!(ax, errors)
    Ks = sort(collect(keys(errors)))

    for (i, k) in enumerate(Ks)
        sns.kdeplot(
            vcat(errors[k]...),
            ax = ax,
            clip = (-π, π),
            label = "K = $k",
            color = get_palette(i),
        )
    end
    ax.set_xlim(-π - 0.2, π + 0.2)
    ax.set_xticks([-π, 0, π])
    ax.set_xticklabels([L"-\pi", "0", L"\pi"],)

    ax.set_ylim(0, 1.0)

    ax.set_xlabel("Error")
    ax.set_ylabel("Probability density")

    return ax
end

"Plot figure for stimulus prioritization effects."
function plot_prioritization!(ax, errors_cued, errors_uncued; fun = cvar)
    # make sure that cued and uncued errors have same keys (= set sizes)
    @assert sort(collect(keys(errors_cued))) ==
            sort(collect(keys(errors_uncued)))

    Ks = sort(collect(keys(errors_cued)))

    plot_data = [
        (errors_cued, Dict(:label => "Cued")),
        (errors_uncued, Dict(:label => "Uncued")),
    ]

    for (i, el) in enumerate(plot_data)
        (data, params) = el
        y = map(k -> mean(fun.(data[k])), Ks)
        e = map(k -> stde(fun.(data[k])), Ks)
        ax.errorbar(
            Ks,
            y,
            yerr = e,
            fmt = "o-",
            color = get_palette(i + 3),
            markersize = 5;
            params...,
        )
    end
    ax.set_xlim(1.5, 12)
    ax.set_ylim(0.125, 4)

    ax.set_xscale("log")
    ax.set_xticks(Ks)
    ax.set_xticklabels(Ks)

    ax.set_yscale("log")
    ax.set_yticks([0.25, 0.5, 1, 2])
    ax.set_yticklabels([0.25, 0.5, 1, 2])

    ax.minorticks_off()

    ax.set_xlabel("Set size")

    return ax
end