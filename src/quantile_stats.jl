# quantile_stats.jl
# Compute statistics based on quantiles.

using DataFrames, DataFramesMeta

"Divide `x` into `n` quantiles."
function xtile(x; n = 4)
    # compute quantiles of x
    qlimit = quantile(x, LinRange(0, 1, n + 1))

    # map each value of x to its quantile index
    qindex = map(x) do v
        min(searchsortedlast(qlimit, v), n)
    end

    return qindex, qlimit
end

"Get mean and standard error of `x` based on `q` quantiles of `y`."
function quantile_stats(x, y; q = 4)
    # divide `y` into q quantiles
    qindex, qlimit = xtile(y; n = q)
    qcenter = qlimit[1:end-1] .+ diff(qlimit) ./ 2

    # add new cols to df
    df = DataFrame()
    df.x = x
    df.q = qindex
    df.qcenter = qcenter[qindex]

    # add quantile stats
    return @by(
        df,
        :q,
        :μ = mean(:x),
        :se = stde(:x),
        :qcenter = first(:qcenter),
    )
end