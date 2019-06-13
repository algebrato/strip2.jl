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
delta_ell = 20.0

DlTT2 = Array{Float64, 1}(undef, Int(round(ell_max/delta_ell)))

avg_mappa = sum(mappa ./ hit) / (NN*NN)
avg_conv  = sum(conv_map) / (NN*NN)
avg_dest  = sum(map_d) / (NN*NN)

win_map_raw = strip2.windowing!(NN, ((mappa ./ hit) .- avg_mappa) .- (conv_map .- avg_conv) )
win_map_dest = strip2.windowing!(NN, (map_d .- avg_dest) .- (conv_map .- avg_conv))

e_raw, d_raw = strip2.get_power_spectrum(win_map_raw, win_map_raw, ell_max, delta_ell,
                                        pix_size, NN)
e_des, d_des = strip2.get_power_spectrum(win_map_dest, win_map_dest, ell_max, delta_ell,
                                        pix_size, NN)

plot(e_raw, d_raw .* ratio2, yaxis = :log)
plot!(e_des, d_des .* ratio2, yaxis = :log)
