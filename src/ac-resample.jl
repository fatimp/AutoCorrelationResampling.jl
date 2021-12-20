"""
    ac_resample(ac, ratio[; filter])

Resample three-dimensional autocorrelation function `ac` by a ratio of
`ratio` times along the third axis, i.e. `ac` of shape `(x,y,z)`
becomes an array of shape `(x,y,ratio * z)`.
"""
function ac_resample(ac     :: AbstractArray{T, 3},
                     ratio  :: Union{Integer, Rational};
                     filter :: VoF = ratio |> filter_length |> filter_coeffs) where T <: AF
    x, y, z = let s = size(ac); s[1], s[2], s[3] end

    # Upsampling coefficient
    n = numerator(ratio)

    # Upsampled autocorrelation function
    ac_up = zeros(T, (x, y, n*z))

    ## To preserve "energy" of a signal
    filter = n * filter

    # Upsample
    for k in 1:z
        ac_up[:,:,n*(k-1)+1] = ac[:,:,k]
    end

    # Filter
    for k in 1:n:n*z
        cureven  = ac_up[:,:,mod1(k + 0, n*z)]
        nexteven = ac_up[:,:,mod1(k + n, n*z)]

        for l in 1:n-1
            ac_up[:,:,k + l] = filter[l]*cureven + filter[n-l]*nexteven
        end
    end

    # Downsample the result
    slices = (ac_up[:,:,k] for k in 1:denominator(ratio):n*z)
    return reduce((x, y) -> cat(x, y; dims = (3)), slices)
end
