# Sum of all filter coefficients
csum(coeff :: VoF) = coeff[1]  + 2sum(coeff[2:end])

# Filter evaluated at point `z`
fn(coeff :: VoF, z :: Complex) = coeff[1] + sum(c*(z^n + z^-n) for (n, c) in enumerate(coeff[2:end]))

cnonnegative(coeff :: VoF) =
    2000 * max(0, -minimum(real(fn(coeff, exp(2π * x * im))) for x in 0:0.001:1))

# The target function for coefficients of our filter incorporates the
# following:
#
# 1) We want the sum of all coefficients to be 1 so our filter is
#    low-pass.
# 2) The filter also has to satisfy conditions of Riesz-Fejér
#    factorisation lemma, i.e. map unit circle to non-negative real
#    values.
# 3) The free coefficient is fixed at 1/length because porosity is
#    everything.
#
# The coefficients found by minimizing `targetfn` are not unique.
function targetfn(coeff :: VoF)
    len = length(coeff) + 1
    coeff1 = vcat([1/len], coeff)
    return cnonnegative(coeff1) + (csum(coeff1) - 1)^2
end

"""
    filter_coeffs(length[; initial])

Get coefficients of symmetrical low-pass FIR filter of length
`2length-1` which maps autocorrelation functions to autocorrelation
functions. Argument `length` being equal to `n` is enough for
upscaling autocorrelations function `n` times.

A keyword argument `initial` can be an array of length `length - 1`
which works as an initial guess for filter coefficients.
"""
function filter_coeffs(length :: Integer; initial :: VoF = 0.5 * ones(Float64, (length - 1)))
    res = optimize(targetfn, initial, length > 2 ? NelderMead() : SimulatedAnnealing())
    return vcat([1/length], Optim.minimizer(res))
end
