# model.jl
# Main model: rate-distortion theory of working memory.

include("utils.jl")
include("circ.jl")

using Statistics, StatsFuns, StatsBase, Distributions, Random

"Container for ITI-RI timing condition."
mutable struct TimingCondition
    name::String
    iti::Integer
    ri::Integer
end

"Type of a sequence for multiple populations."
Sequence = Vector{Tuple{Vector{Int64},Int64}}

"Type of a sequence for a single population."
SequenceOnePop = Vector{Tuple{Int64,Int64}}

"""
    sampletunings(tunings, r[; N=1000])

Sample from `tunings` using weights specified by `r`.
"""
function sampletunings(tunings, r; N = 1000)
    npopulations = size(r, 2)
    samps = zeros(N, npopulations)

    for i = 1:npopulations
        sample!(tunings, Weights(r[:, i]), @view samps[:, i])
    end

    return samps
end

"Build distortion matrix from tuning preferences and stimuli."
function build_distortion(tuning, stimuli, d = abs ∘ cdist)
    D = zeros(length(tuning), length(stimuli))
    for s = 1:size(D, 2), t = 1:size(D, 1)
        D[t, s] = d(tuning[t], stimuli[s])
    end
    return D
end

"Struct containing a simulation."
mutable struct Simulation
    nneurons::Int64
    nstimuli::Int64
    npopulations::Int64

    r::Any
    spikecount::Matrix{Int64}
    t::Int64

    w::Matrix{Float64}
    w0::Float64
    wlr::Float64
    wlim::Tuple{Float64,Float64}

    β::Any
    β0::Any
    βlr::Float64
    βlim::Tuple{Float64,Float64}

    C::Any
    R::Any

    tunings::Vector{Float64}
    stimuli::Vector{Float64}
    D::Matrix{Float64}

    spikes::Matrix{Vector{Float64}}

    rng::MersenneTwister
end

"Compute distortion between `a` and `b`."
distortion(a, b) = -cos(cdist(a, b))

"Construct a simulation."
function Simulation(
    nneurons,
    nstimuli,
    npopulations;
    C = 0.9,
    w0 = 0.5,
    wlr = 3e-6,
    wlim = (-12, 0),
    β0 = 20.0,
    βlr = 1e-3,
    βlim = (0, 200),
    kwargs...,
)
    # stability condition (might not be necessary if we clip w)
    @assert wlr ≤ -wlim[1] / (w0 * exp(-wlim[1]) - 1)

    sim = Simulation(
        nneurons,
        nstimuli,
        npopulations,
        # r, spikecount, t
        zeros(0, 0),
        zeros(Int64, 0, 0),
        0,
        # w, w0, wlr, wlim
        zeros(0, 0),
        w0,
        wlr,
        wlim,
        # β, β0, βlr, βlim
        β0,
        β0,
        βlr,
        βlim,
        # C, R
        C,
        0,
        # tunings, stimuli, D,
        zeros(0),
        zeros(0),
        zeros(0, 0),
        # spikes
        map(i -> ones(0), ones(0, 0)),
        # rng
        MersenneTwister(),
    )

    return reset!(sim; kwargs...)
end

"Reset `sim`."
function reset!(
    sim::Simulation;
    tunings = nothing,
    stimuli = nothing,
    dfun = distortion,
    seed = 1234,
)
    # firing rates
    sim.r = ones(sim.nneurons, sim.npopulations)
    sim.r ./= sum(sim.r)

    # reset spikes
    sim.spikes = map(i -> ones(0), ones(sim.nneurons, sim.npopulations))
    sim.spikecount = zeros(Int64, sim.nneurons, sim.npopulations)
    sim.t = 0

    # weights
    sim.w = log.(sim.r)

    # information rate
    sim.R = 0

    # distortion matrix
    sim.tunings =
        isnothing(tunings) ? LinRange(0, 2π, sim.nneurons + 1)[1:end-1] :
        tunings
    sim.stimuli =
        isnothing(stimuli) ? LinRange(0, 2π, sim.nstimuli + 1)[1:end-1] :
        stimuli
    sim.D = build_distortion(sim.tunings, sim.stimuli, dfun)

    # β
    sim.β = sim.β0

    # rng
    sim.rng = MersenneTwister(seed)

    return sim
end

"Run `sim`."
function run!(
    sim::Simulation,
    seq::Sequence;
    fr = 20,
    Δt = 0.001,
    spike_retention = 10,
    warmup = 0,
    πprobe = uniform(sim.npopulations),
    πθ = uniform(sim.nstimuli),
    kwargs...,
)
    @assert length(πprobe) == sim.npopulations
    @assert length(πθ) == sim.nstimuli

    # add phantom spike to make sure arrays are not empty
    if any(isempty.(sim.spikes))
        push!.(sim.spikes, -spike_retention)
        sim.spikecount .+= 1
    end

    # compute gain to achieve wanted population firing rate
    γ = fr * sim.nneurons * sim.npopulations

    # spike counter
    sim.spikecount += round.(Int, warmup * γ * Δt * sim.r)

    # track variables
    len = sum(subseq[end] for subseq in seq)
    Rs = zeros(eltype(sim.C), len)
    βs = zeros(eltype(sim.β), len)

    rs = zeros(eltype(sim.C), length(seq), size(sim.r)...)

    ravg =
        sum(sim.spikecount) > 0 ? sim.spikecount / sum(sim.spikecount) :
        zeros(size(sim.r)...)

    clock = 0

    # run simulation
    for (j, subseq) in enumerate(seq)

        # reset spike counter for subsequence
        spikecount_subseq = zeros(Int, size(sim.r)...)

        stimulus, niter = subseq
        isstimulus = all(.!iszero.(stimulus))

        for iter = 1:niter
            # increase global simulation clock
            sim.t += 1
            clock += 1

            # compute new spikes according to Poisson process
            spikes = rand.(sim.rng, Poisson.(γ .* sim.r .* Δt))
            for (neuron, nspikes) in zip(sim.spikes, spikes)
                # uniformly distribute nspikes over interval 
                times = rand(sim.rng, Uniform(sim.t, sim.t + Δt), nspikes)
                append!(neuron, sort(times))
            end

            # append to spikecount
            spikecount_subseq += spikes

            # update spike-dependent plasticity
            if sim.wlr > 0
                z = spikes .> 0
                # z = last.(sim.spikes) .> sim.t - spike_retention
                Δw = sim.w0 .* exp.(-sim.w) .* z .- 1
                sim.w .+= sim.wlr .* Δw

                # clip and normalize weights
                sim.w .= min.(max.(sim.w, sim.wlim[1]), sim.wlim[2])
                sim.w .-= logsumexp(sim.w)
            end

            # update membrane potential
            if isstimulus
                u = copy(sim.w) - sim.β * sim.D[:, stimulus] .* πprobe'
                sim.r = exp.(u .- logsumexp(u))

                # compute information rate
                r =
                    sum(spikecount_subseq) > 0 ?
                    spikecount_subseq / sum(spikecount_subseq) :
                    zeros(size(sim.r)...)

                sim.R = sum(xlogx.(r) - xlogx.(ravg))
                sim.R = max(sim.R, 0)
            else
                sim.r = exp.(sim.w)
                sim.R = 0
            end

            # update beta
            if sim.βlr > 0
                sim.β += sim.βlr * (sim.C - sim.R)
                sim.β = min(max(sim.β, sim.βlim[1]), sim.βlim[2])
            end

            # store intermediate values
            Rs[clock] = sim.R
            βs[clock] = sim.β
        end

        # sample sim
        rs[j, :, :] = sim.r

        # add spikes from previous stimulus and compute new 
        # average firing rate
        sim.spikecount += spikecount_subseq
        ravg = sim.spikecount / sum(sim.spikecount)
    end

    # remove phantom spike
    # popfirst!.(sim.spikes)

    return Rs, βs, rs
end
run!(sim::Simulation, seq::SequenceOnePop; kwargs...) = run!(
    sim::Simulation,
    map(s -> (fill(s[1], sim.npopulations), s[2]), seq);
    kwargs...,
)

"""
    samplesim(sim[; N=1000])

Sample from a simulation using weights specified by the firing rates in `sim`.
"""
samplesim(sim::Simulation; N = 1000) = sampletunings(sim.tunings, sim.r; N = N)

"Get effective capacity for a simulation with given timing condition."
getcapacity(s::Simulation, tc::TimingCondition) = s.C * (tc.ri + tc.iti) / tc.ri

"Get random sequence going through `nstimuli` stimuli `repeat` times."
function getseq(nstimuli, tc::TimingCondition, cue; repeat = 1)
    @assert repeat ≥ 1

    ris = [[(i, tc.ri) for i in randperm(nstimuli)] for _ = 1:repeat]
    itis = [[(0, tc.iti) for _ = 1:nstimuli] for _ = 1:repeat]

    ris = vcat(ris...)
    itis = vcat(itis...)

    seq = collect(Iterators.flatten(zip(ris, itis)))
    return vcat(seq, (cue, tc.ri))
end
getseq(s::Simulation, tc::TimingCondition, cue; repeat = 1) =
    getseq(s.nstimuli, tc::TimingCondition, cue; repeat = repeat)