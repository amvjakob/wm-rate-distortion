# circ.jl
# Functions related to circular statistics.

"Circular distance between angle `x` and `y` (given in radians)."
cdist(x, y) = angle(exp(1im * x) / exp(1im * y))

"Uncentered trigonometric moment of a of order `p`."
cmoment(a, p = 1) = sum(exp.(1im .* p .* a)) / length(a)

"Circular resultant (= magnitude of first uncentered trigonometric moment)."
cresultant(a) = abs(cmoment(a, 1))

"Circular mean."
cmean(a) = angle(cmoment(a, 1))

"Circular variance."
cvar(a) = -2log(cresultant(a))
"Alternative measure of circular variance."
cvar_bounded(a) = 1 - cresultant(a)

"Circular standard deviation."
cstd(a) = sqrt(-2log(cresultant(a)))

"Circular standard error."
cstde(a) = cstd(a) / √length(a)

"Circular kurtosis."
function ckurtosis(a)
    # compute first and second uncentered trigonometric moments
    m1 = cmoment(a, 1)
    m2 = cmoment(a, 2)

    # compute circular kurtosis
    return (abs(m2) * cos(angle(m2) - 2angle(m1)) - abs(m1)^4) /
           ((1 - abs(m1))^2)
end
