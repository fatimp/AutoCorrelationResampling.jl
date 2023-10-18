"""
    filter_coeffs(length)

Get coefficients of symmetrical low-pass FIR filter of length
`2length-1` which maps autocorrelation functions to autocorrelation
functions. Argument `length` being equal to `n` is enough for
upscaling autocorrelations function `n` times.
"""
filter_coeffs(length  :: Integer) = [(length - n)/length^2 for n in 1:(length-1)]

# Filter length function
filter_length(ratio) = max(numerator(ratio), denominator(ratio))
