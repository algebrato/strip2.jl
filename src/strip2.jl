module strip2

using Random
using FFTW
using ConjGrad

# gensky methods
export make_CMB_T_map
export make_CMB_pol_maps
export SZ_source
export PS_exp
export PS_poisson

# Beam methods
export convolve_beam

# Noise sources
export white_noise
export atmospheric_noise
export one_over_f
export filter_noise
export observe_sky
export tod_noise

# Analysis tools
export windowing!
export get_power_spectrum

include("gensky.jl")
include("beam.jl")
include("noise_map.jl")
include("filter.jl")
include("genpoints.jl")
include("mapmaker.jl")
include("analysis.jl")

end  # module strip2
