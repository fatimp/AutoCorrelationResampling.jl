"""
    ac_upsample(ac, filter)

Upsample three-dimensional autocorrelation function `ac` by
`length(filter)` times along the third axis, i.e. `ac` of shape
`(x,y,z)` becomes an array of shape `(x,y,length(filter) * z)`.
"""
function ac_upsample(ac     :: AbstractArray{T, 3},
                     filter :: VoF) where T <: AF
    s = size(ac)
    x, y, z = s[1], s[2], s[3]

    # Restored autocorrelation function
    n = length(filter)
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
            ac_up[:,:,k + l] = filter[l+1]*cureven + filter[n-l+1]*nexteven
        end
    end

    # This step is actually not needed, because filter[1] is always 1.
    for k in 1:n:n*z
        ac_up[:,:,k] = filter[1] * ac_up[:,:,k]
    end

    return ac_up
end
