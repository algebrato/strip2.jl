include("random_sky.jl")
include("gen_noise_map.jl")

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

Map_win = strip2.windowing!(NN, cmb_T_map)
ell2, DlTT2 = strip2.get_power_spectrum(Map_win, Map_win, ell_max, delta_ell,
                                        pix_size, NN)
plot(ell, DlTT, yaxis = :log)
plot!(ell2, DlTT2, yaxis = :log)
