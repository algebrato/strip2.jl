module strip2

using Random
using FFTW
using ConjGrad
using Plots
using LinearAlgebra
using ConjGrad

# gensky methods
export make_CMB_T_map
export make_CMB_pol_maps
export SZ_source
export PS_exp
export PS_poisson
export Plot_Sky

# Beam methods
export convolve_beam

# Noise sources
export white_noise
export atmospheric_noise
export one_over_f
export one_over_f2
export filter_noise
export observe_sky
export get_tod
export tod_noise
export binned_map
export destriper
export conjgrad
export denoise

# Analysis tools
export windowing!
export get_power_spectrum

abstract type TOD end

struct TOD_fake_pointing <: TOD
    time_order_data::AbstractArray{<:Number}
    pointing::AbstractArray{<:Number}
    noise::AbstractArray{<:Number}
end

include("gensky.jl")
include("beam.jl")
include("noise_map.jl")
include("filter.jl")
include("genpoints.jl")
include("mapmaker.jl")
include("analysis.jl")

end  # module strip2
