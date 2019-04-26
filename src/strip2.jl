module strip2

using Random
using Healpix
using FFTW

# gensky methods
export make_CMB_T_map
export SZ_source
export PS_exp
export PS_poisson

include("gensky.jl")

end  # module strip2
