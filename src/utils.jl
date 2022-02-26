# utils.jl
# Diverse utility functions.

using Dates, DataFrames, HypothesisTests, Bootstrap, PyCall
using Pipe: @pipe

"List of seeds."
SEEDS = [
    1239,
    12494,
    12934,
    3493,
    1012,
    125,
    3023,
    45045,
    32,
    3029,
    1294,
    986,
    4983,
    585,
    876,
    34,
    5657,
    78,
    678,
    346,
    57,
    78,
    246,
    375,
    4567,
    876,
    43645,
    58455,
    4531,
    63,
];

"Log `x`."
function lg(x...)
    println("[$(now())] $(join(x, " ")...)")
    flush(stdout)
end

"Uniform discrete distribution of length `n`."
uniform(n::Int) = ones(n) ./ n

"Get nonzero elements of `a`."
nonzeros(a; kwargs...) = filter(!iszero, a; kwargs...)

"Get nonnan elements of `a`."
nonnans(a; kwargs...) = filter(!isnan, a; kwargs...)

"Compute standard error of `a`."
stde(a) = std(a) / √length(a)

"Filter df by row function."
group_and_filter_df(df, groupcol, rowfun) = begin
    return combine(groupby(df, groupcol)) do gdf
        return filter(rowfun, gdf)
    end
end

"""
    bootstrapsig(data[; fun=mean, sig=0.95, n=1000, dims=[2]])

Bootstrap `data` and compute the statistic given by applying `fun` to 
slices on the dimensions given in by `dims`.
"""
function bootstrapsig(data; fun = mean, sig = 0.95, n = 1000, dims = [2])
    # compute if each slice is signifcant
    issig = mapslices(data; dims = dims) do dataslice
        sample = bootstrap(fun, dataslice, BasicSampling(n))
        signif = BasicConfInt(sig)
        nvar = 1

        # compute confidence interval given sample
        t0, ci_low, ci_high = confint(sample, signif, nvar)

        # significant if the ci doesn't include 0
        return ci_low > 0 || ci_high < 0
    end

    # convert list of bools to the indices of the significant entries
    return findall(issig)
end

"""
    bootstrapsig(x, y[; fun=mean, sig=0.95, n=1000, dims=[2]])

Bootstrap `x` and `y`, compute the statistic given by applying `fun` to 
slices on the dimensions given in by `dims`, check if they are different.
"""
function bootstrapsig(x, y; fun = mean, sig = 0.95, n = 1000)
    # compute if diff is signifcant
    bs1 = bootstrap(fun, x, BasicSampling(n))
    bs2 = bootstrap(fun, y, BasicSampling(n))
    z = straps(bs1)[1] - straps(bs2)[1]

    issig = pvalue(OneSampleTTest(z)) ≤ (1 - sig)

    return issig * sign(mean(z))
end

"Exponential moving average of `x` using `n` samples."
function exp_moving_avg(x, n, initialval)
    k = 1 / n
    result = copy(x)

    for (i, row) in enumerate(eachrow(result))
        if i == 1
            row .= k * row + (1 - k) * initialval
        else
            row .= k * row + (1 - k) * result[i-1, :]
        end
    end

    return result
end

"Read pickle file by first converting Linux newlines to Windows newlines."
function read_pickle_windows(src)
    py"""
    import os
    import pickle

    def read_pickle_windows(src):
        try:
            return pickle.load(open(src, 'rb'), encoding='cp850')
        except:
            # fix pickle by replacing linux newline by windows newline
            dst = src + '.tmp'
            data = open(src).read().replace('\r\n', '\n')
            open(dst, 'w', newline='\n').write(data)

            # load fixed pickle file
            data = pickle.load(open(dst, 'rb'), encoding='latin1')
            os.remove(dst)

            return data
    """

    return py"read_pickle_windows"(src)
end


"Read pickle file for Linux systems."
function read_pickle_linux(src)
    py"""
    import pickle

    def read_pickle_linux(src):
        return pickle.load(open(src, 'rb'), encoding='latin1')
    """

    return py"read_pickle_linux"(src)
end
