# serialdep.jl
# Functions related to computing serial dependence.

include("model.jl")

using LsqFit, DataFrames, DataFramesMeta, StatsPlots

"Serial dependence model with DoG: y = sqrt(2e) * x * a * w * exp(-(w * x)²)."
dog(x, p) = @. x * p[1] * p[2] * sqrt(2exp(1)) * exp(-(p[2] * x)^2)

"Fit serial dependence model."
function fit_dog(x, y; p0 = [1.0, 1.0])
    fit = curve_fit(dog, x, y, p0)

    a, w = fit.param
    serialdep = 2 * a * sign(w)

    return [a, w], serialdep
end

"""
    get_df_consecutive(cues, stimuli, reports)

Transform given `cues`, `stimuli` and `reports` into  a `DataFrame` for serial 
dependence analysis by creating rows with `pre` and `cur` values.
The vectors are assumed to be in trial order. 
"""
function get_df_consecutive(
    cues::Vector{Int64},
    stimuli::Vector{Float64},
    reports::Vector{Float64},
)
    @assert length(cues) == length(stimuli) == length(reports)

    ntrials = length(cues)

    df = DataFrame()
    df.cue_pre = cues[1:end-1]
    df.cue_cur = cues[2:end]

    df.stimulus_pre = stimuli[1:end-1]
    df.stimulus_cur = stimuli[2:end]

    df.report_pre = reports[1:end-1]
    df.report_cur = reports[2:end]

    df.trial = 1:ntrials-1
    df.triallen = fill(ntrials - 1, nrow(df))

    return df
end


"""
    get_df_consecutive(simulation, cues, rs[; naug=100])

Get a `DataFrame` for serial dependence analysis by creating rows with `pre`
and `cur values`. Extract stimuli and reports based on `cues` and `rs`. Each 
report is sampled `naug` times.
The vectors are assumed to be in trial order. 
"""
function get_df_consecutive(
    simulation::Simulation,
    cues::Vector{Int64},
    rs::Matrix{Float64};
    naug = 100,
)
    ntrials, nneurons = size(rs)

    @assert length(cues) == ntrials
    @assert simulation.nneurons == nneurons

    # sample simulation
    reports = zeros(ntrials, naug)
    for i = 1:ntrials
        reports[i, :] =
            vec(sampletunings(simulation.tunings, rs[i, :]; N = naug))
    end

    # get df
    df = DataFrame()
    df.cue_pre = repeat(cues[1:end-1], naug)
    df.cue_cur = repeat(cues[2:end], naug)

    df.stimulus_pre = simulation.stimuli[df.cue_pre]
    df.stimulus_cur = simulation.stimuli[df.cue_cur]

    df.report_pre = vec(reports[1:end-1, :])
    df.report_cur = vec(reports[2:end, :])

    df.trial = repeat(1:ntrials-1, naug)
    df.triallen = fill(ntrials - 1, nrow(df))

    return df
end

"""
get_df_multi(simulation, cues, rs[; naug=100])

Get a `DataFrame` for serial dependence analysis by creating rows with `pre`
and `cur values`. Extract stimuli and reports based on `cues` and `rs`. Each 
report is sampled `naug` times. This function is to be used when `simulation`
contains several populations.
The vectors are assumed to be in trial order. 
"""
function get_df_multi(
    simulation::Simulation,
    cues::Matrix{Int64},
    rs::Array{Float64,3};
    naug = 100,
)
    ntrials, nneurons, npopulations = size(rs)

    @assert all(size(cues) .== (ntrials, npopulations))
    @assert simulation.nneurons == nneurons
    @assert simulation.npopulations == npopulations

    # sample simulation
    reports = zeros(ntrials, npopulations * naug)
    for i = 1:ntrials
        # store reports as one augmentation after the other:
        # reports[i,1:simulation.npopulations], then 
        # reports[i,simulation.npopulations+1:simulation.npopulations]
        # and so forth
        reports[i, :] =
        # transpose important for consecutive augmentations
            vec(sampletunings(simulation.tunings, rs[i, :, :]; N = naug)')
    end

    # get df
    df = DataFrame()
    df.cue = vec(repeat(cues, 1, naug))
    df.stimulus = simulation.stimuli[df.cue]

    df.report = vec(reports)

    df.trial = vec(repeat(1:ntrials, 1, npopulations * naug))
    df.triallen = fill(ntrials, nrow(df))

    return df
end

"Get the mean reported value for each cue."
function get_mean_report(df, cuecol = :cue, reportcol = :report)
    # check that colnames are present in df
    @assert all(string(c) in names(df) for c in [cuecol, reportcol])

    # compute mean report for each cue
    return @chain df begin
        # rename cues as :x and reports as :y
        @select(:x = $cuecol, :y = $reportcol)

        @aside μ = @chain _ begin
            # compute mean y value for each x
            @by(:x, :ymean = mean(:y))
            # convert to associative Dict :x => :ymean
            @df Dict(Pair.(:x, :ymean))
        end
        # map each x value onto its mean y value
        @select(@byrow :ymean = μ[:x])
        # return mean y value
        @df vec(:ymean)
    end
end

"Add new column with mean reported value for each cue."
function add_mean_report!(df, cuecol = :cue, reportcol = :report)
    # add new col to df with mean report for each trial's cue
    df[!, "mean_"*string(reportcol)] = get_mean_report(df, cuecol, reportcol)

    return df
end

"Compute distances and errors for data given in `df`."
function compute_dist_error!(df; folded_error_distcol = :dist)
    # check that colnames are present in df
    cols = [:cue_pre, :cue_cur, :report_pre, :report_cur]
    @assert all(string(c) in names(df) for c in cols)

    # add mean report for each cue, separately for current and previous cues
    add_mean_report!(df, :cue_pre, :report_pre)
    add_mean_report!(df, :cue_cur, :report_cur)

    # newly generated columns
    cols = [:dist, :dist_report, :dist_report_to_mean, :error, :error_to_mean]

    # compute distances and errors for df
    return @chain df begin
        # add distances and errors
        @transform!(
            # distance metrics
            :dist = cdist.(:stimulus_pre, :stimulus_cur),
            :dist_report = cdist.(:report_pre, :report_cur),
            :dist_report_to_mean = cdist.(:report_pre, :mean_report_cur),
            # error metrics
            :error = cdist.(:report_cur, :stimulus_cur),
            :error_to_mean = cdist.(:report_cur, :mean_report_cur),
        )

        # add absolute variants of generated columns
        transform!(cols .=> (v -> abs.(v)) .=> "abs_" .* string.(cols))

        # add folded error
        @transform!(
            :folded_error = :error .* sign.($folded_error_distcol),
            :folded_error_to_mean =
                :error_to_mean .* sign.($folded_error_distcol),
        )
    end
end

"Compute smoothed error signal."
function smooth_serialdep_error(
    df;
    distcol = :abs_dist,
    errcol = :folded_error,
    win = π / 3,
    step = π / 30,
)
    # check that colnames are present in df
    cols = [:subject, distcol, errcol]
    @assert all(string(c) in names(df) for c in cols)

    subjects = unique(df.subject)
    x = collect(-π-win:step:π)
    y = zeros(length(x), length(subjects))
    e = zeros(length(x), length(subjects))

    # compute serial dependence for each subject
    for (i, subject) in enumerate(subjects)
        sdf = @subset(df, :subject .== subject)

        for (j, t) in enumerate(x)
            # get values in window
            inwin = @subset(sdf, t .< $distcol .< t + win)
            if nrow(inwin) > 0
                y[j, i] = cmean(inwin[:, errcol])
                e[j, i] = cstde(inwin[:, errcol])
            else
                y[j, i] = 0
                e[j, i] = 0
            end
        end
    end

    return x .+ win / 2, y, e
end

"Compute smoothed error signal for given x and y arrays."
function smooth_serialdep_error(x, y; win = π / 3, step = π / 30)
    xx = collect(-π-win:step:π)
    yy = zeros(length(xx))
    ee = zeros(length(xx))

    # compute serial dependence
    for (j, t) in enumerate(xx)
        # get values in window
        inwin = t .< x .< t + win
        if sum(inwin) > 0
            yy[j] = cmean(y[inwin])
            ee[j] = cstde(y[inwin])
        end
    end

    return xx .+ win / 2, yy, ee
end


"Compute firing rate based on spikes."
function compute_fr!(df, spikes, timing::TimingCondition)
    @assert "trial" in names(df)

    trials = unique(df.trial)
    ntrials = length(trials)

    triallen = timing.ri + timing.iti

    intervals = Dict(
        "ri" => 1:timing.ri,
        "iti" => timing.ri+1:timing.ri+timing.iti,
        "trial" => 1:timing.ri+timing.iti,
    )
    fr = Dict(name => zeros(ntrials) for (name, _) in intervals)
    frneuron =
        Dict(name => zeros(ntrials, size(spikes)...) for (name, _) in intervals)

    # compute firing rates for all trials
    for (name, range) in intervals
        for trial in trials
            start = (trial - 1) * triallen + range[1]
            stop = (trial - 1) * triallen + range[end]

            nspikes =
                map(neuronspikes -> sum(start .< neuronspikes .< stop), spikes)

            frneuron[name][trial, :, :] = 1000 .* nspikes ./ length(range)
            fr[name][trial] = mean(frneuron[name][trial, :, :])
        end
    end

    fravg = mean(fr["trial"])
    frneuronavg = dropdims(mean(frneuron["trial"], dims = 1), dims = 1)

    ir = Dict(k => log.(v) .- log.(fr["trial"]) for (k, v) in fr)
    iravg = Dict(k => log.(v) .- log(fravg) for (k, v) in fr)

    R = Dict(
        k => vec(
            mapslices(v, dims = [2, 3]) do f
                sum(xlogx.(f ./ sum(f)) .- xlogx.(frneuronavg ./ sum(frneuronavg)))
            end,
        ) for (k, v) in frneuron
    )

    return transform!(groupby(df, :trial)) do gdf
        newcols = DataFrame()
        for (name, _) in intervals
            newcols[!, "fr_"*name] = fr[name][gdf.trial]
            newcols[!, "ir_"*name] = ir[name][gdf.trial]
            newcols[!, "iravg_"*name] = iravg[name][gdf.trial]

            newcols[!, "R_"*name] = R[name][gdf.trial]
        end
        return newcols
    end
end;