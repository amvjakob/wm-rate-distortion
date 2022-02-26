# plot.jl
# General plotting tools.

include("utils.jl")
include("model.jl")
include("circ.jl")

using LaTeXStrings, Plots, StatsPlots, Colors

"Plot spike events of simulation `sim`."
function plot_spikes(sim::Simulation, seq::Sequence; save_fig = false)
    # get optimal neuron for each stimulus
    optimal_neurons = map(index -> index[1], vec(argmin(sim.D, dims = 1)))

    # draw figure
    pyplot()
    PyPlot.figure(figsize = (11, 4))

    colors = theme_palette(:auto)
    pyplotcolor(c::RGB) = (red(c), green(c), blue(c))

    # find times of change from nnz to zero
    vlines = cumsum(subseq[2] for subseq in seq)
    PyPlot.vlines(
        vlines,
        0,
        sim.nneurons,
        label = "Stimulus time",
        color = :black,
        alpha = 0.3,
    )

    for pop = 1:sim.npopulations
        params = Dict(:color => colors[(pop-1)%length(colors)+1] |> pyplotcolor)
        PyPlot.eventplot(
            sim.spikes[:, pop],
            label = "Population $(pop)";
            params...,
        )

        # find horizontal lines (= optimal neurons for stimulus)
        hlines = map(s -> s[1][pop], seq) |> unique |> nonzeros
        hlines = optimal_neurons[hlines]

        params = Dict(:alpha => 0.3, params...)
        PyPlot.hlines(
            hlines,
            0,
            maximum(maximum.(sim.spikes)),
            label = "Stimulus neuron";
            params...,
        )
    end

    PyPlot.title("Total spikes = $(sum(length.(sim.spikes)))")
    PyPlot.xlabel("Time (ms)")
    PyPlot.ylabel("Neurons")
    save_fig && PyPlot.savefig("wm_spike_train.png")

    return PyPlot.gcf()
end
plot_spikes(sim::Simulation, seq::SequenceOnePop; kwargs...) =
    plot_spikes(sim, map(s -> ([s[1]], s[2]), seq); kwargs...)

"Plot parameter `param` from `sim` vs the simulation's tunings."
function plot_param(
    sim::Simulation,
    param,
    seq::Sequence;
    paramname = "Parameter",
)
    p = plot()
    value = getfield(sim, param)

    for pop = 1:sim.npopulations
        plot!(sim.tunings, value[:, pop], label = "$paramname pop $pop")

        if !isnothing(seq)
            stim = map(s -> s[1][pop], seq) |> unique |> nonzeros
            vline!(sim.stimuli[stim], label = "Stimulus pop $pop")
        end
    end

    return p
end
plot_param(sim::Simulation, param, seq::SequenceOnePop; kwargs...) =
    plot_param(sim, param, map(s -> ([s[1]], s[2]), seq); kwargs...)

"Filter arrays."
function filter_dists_errors(dists, errors)
    indx_clean = abs.(errors) .< 10 * std(errors)
    return dists[indx_clean], errors[indx_clean]
end

"Plot time-course of R."
function plotR(Rs; kwargs...)
    plot(Rs, label = ""; kwargs...)
    hline!([sim.C], label = "Maximal capacity")
    xlabel!("Time (ms)")
    ylabel!("Information rate")
end

"Plot time-course of β."
function plotβ(βs; kwargs...)
    plot(βs, label = ""; kwargs...)
    xlabel!("Time (ms)")
    ylabel!("β")
end

"Get blank plot"
plot_blank() =
    plot(legend = false, grid = false, foreground_color_subplot = :white)

"Get rectangle"
rectangle(w, h, x, y) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])

"Add significance bar to plot `p`."
function sigbar!(p, x, y, w, h, sig; params...)
    for s in sig
        plot!(p, rectangle(w, h, x[s], y); params...)
    end

    return p
end
