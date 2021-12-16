module AutoCorrelationResampling
using Optim
using Base.Iterators

AV  = AbstractVector
AF  = AbstractFloat
VoF = AV{<:AF}

include("filter-coeffs.jl")
include("ac-upsample.jl")

export filter_coeffs, ac_upsample

end # module
