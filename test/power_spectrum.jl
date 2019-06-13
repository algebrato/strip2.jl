include("random_sky.jl")
include("gen_noise_map.jl")
using ProgressMeter
# Aailable maps:
# cmb_T_map = CMB temperature map
# sz_map  = SZ Contribute
# point_pois, point_expo = exponential sources contribute
# mapp = cmb_T_map .+ sz_map .+ point_pois .+ point_expo
# conv_map  = Cosmological signal maps convolved with instrument beam

# All noise sources maps:
# white_map
# atmos_map
# one_ove_f
# noise_map = sum(white_map, atmos_map, one_ove_f)

# The specs of every maps are in the random_sky.jl and gen_noise_map.jl

ell_max = 5000.0
delta_ell = 30.0

DlTT2 = Array{Float64, 1}(undef, Int(round(ell_max/delta_ell)))


e_true, d_true = strip2.get_power_spectrum(cmb_T_map, cmb_T_map, ell_max, delta_ell,
                                        pix_size, NN)

e_des, d_des = strip2.get_power_spectrum(map_d, map_d, ell_max, delta_ell,
                                        pix_size, NN)

a=maximum(DlTT)
b=maximum(d[3:length(d)])
ratio = a/b

plot(ell, DlTT, yaxis = :log, xlims = (0, 1000), ylims=(10, 10000))
plot!(e_true, d_true .* ratio, yaxis = :log, xlims = (0, 1000), ylims=(10,10000))
plot!(e_des, d_des .* ratio, yaxis = :log, xlims = (0, 1000), ylims=(10,10000))
