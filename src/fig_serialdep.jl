# fig_serialdep.jl
# Figure for serial dependence plots of rate-distortion theory of working
# memory.

include("utils.jl")
include("circ.jl")
include("pyplot.jl")
include("serialdep.jl")

using PyPlot, LaTeXStrings, DataFrames, DataFramesMeta

"Plot serial dependence error."
function plot_serialdep_error!(
    ax,
    x,
    y;
    e = nothing,
    color = "black",
    xlabel = "Relative color of previous trial (deg)",
    ylabel = "Error on current trial (deg)",
    params...,
)
    # compute error ribbon
    e =
        isnothing(e) ? rad2deg.(vec(mapslices(cstde, y, dims = [2]))) :
        rad2deg.(vec(mapslices(cmean, e, dims = [2])))

    xplot = rad2deg.(x)
    yplot = rad2deg.(vec(mapslices(cmean, y, dims = [2])))

    # plot serial dependence
    ax.plot(xplot, yplot, color = color; params...)
    ax.fill_between(
        xplot,
        yplot .- e,
        yplot .+ e,
        alpha = 0.25,
        facecolor = color,
        edgecolor = "none",
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax
end

"Plot serial dependence error subject by subject."
function plot_serialdep_error_by_subject!(
    ax,
    x,
    y;
    xlabel = "Relative color of previous trial (deg)",
    ylabel = "Error on current trial (deg)",
    params...,
)
    # plot subject-wise serial dependence
    ax.plot(rad2deg.(x), rad2deg.(y), alpha = 0.25; params...)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax
end

"Plot DoG fit of serial dependence."
function plot_serialdep_dog!(
    ax,
    x,
    dogparams;
    xlabel = "Relative color of previous trial (deg)",
    ylabel = "Error on current trial (deg)",
    params...,
)
    # plot DoG curve
    ax.plot(rad2deg.(x), rad2deg.(dog(x, dogparams)); params...)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax
end

"Scatter plot of serial dependence error."
function scatter_serialdep_error(
    df,
    distcol,
    errcol;
    xlabel = "Relative color of previous trial (deg)",
    ylabel = "Error on current trial (deg)",
    params...,
)
    fig = figure("scatter_serialdep_error", figsize = (6, 3))
    ax = subplot(111)

    # scatter plot of error vs dist
    ax.scatter(
        rad2deg.(df[:, distcol]),
        rad2deg.(df[:, errcol]),
        color = :gray,
        opacity = 0.4,
        markersize = 1,
        markercolor = :gray,
        markeralpha = 0.1;
        params...,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # add smoothed curve
    @chain df begin
        @select(:x = $distcol, :y = $errcol)
        @by(:x, :y = mean(:y))
        @orderby(:x)
        @df ax.plot(rad2deg.(:x), rad2deg.(:y))
    end

    return fig
end

"Plot serial dependence as single number against variable `x`."
function plot_serialdep_condition!(
    ax,
    x,
    serialdeps;
    e = nothing,
    xlabel = "Parameter (s)",
    ylabel = "Serial dependence (deg)",
    ylim = (-3, 4),
    params...,
)
    ax.plot(x, serialdeps, color = "k"; params...)
    ax.axhline(y = 0, color = "k", linestyle = "--")

    if !isnothing(e)
        ax.fill_between(
            x,
            serialdeps .- e,
            serialdeps .+ e,
            facecolor = "k",
            edgecolor = "none",
            alpha = 0.1,
        )
    end

    ax.set_ylim(ylim...)

    ax.set_xticks(x)
    ax.set_xticklabels(x)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax
end

"Plot serial dependence build-up"
function plot_serialdep_buildup!(
    ax,
    x,
    y;
    color = get_palette(6),
    ylim = (-1, 3),
)
    yplot = rad2deg.(vec(mapslices(cmean, y, dims = [2])))
    eplot = rad2deg.(vec(mapslices(cstde, y, dims = [2])))

    ax.plot(x, yplot, color = color)
    ax.fill_between(
        x,
        yplot .- eplot,
        yplot .+ eplot,
        facecolor = color,
        edgecolor = "none",
        label = "none",
        alpha = 0.25,
    )
    ax.axhline(y = 0, color = "k", linestyle = "--")

    ax.set_xlabel("Trial number")
    ax.set_ylabel("Error on current trial (deg)")

    ax.set_ylim(ylim...)

    return ax
end