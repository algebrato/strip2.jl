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

@showprogress for i=1:16
    # dir=true
    # mappa = zeros(NN, NN)
    # hit   = zeros(NN, NN)
    # N_baselines = 10
    # dataset_baselines = Array{Main.strip2.TOD_fake_pointing, 1}(undef, N_baselines)
    #
    # @showprogress for i = 1:N_baselines
    #     dataset = strip2.observe_sky(NN, conv_map; noise=true, direction_l_r=dir)
    #     m, h = strip2.binned_map(NN, dataset.time_order_data, dataset.pointing)
    #     mappa .+= m
    #     hit   .+= h
    #     dataset_baselines[i] = dataset
    #     dir =  true ⊻ dir
    # end
    #
    # map_des = strip2.destriper(NN, dataset_baselines; n_iter=10)
    # map_d = reshape(map_des, NN, NN)
    include("gen_noise_map.jl")
    Map = conv_map .+ noise_map
    Map_win_d = strip2.windowing!(NN, Map)
    e, d = strip2.get_power_spectrum(Map_win_d, Map_win_d, ell_max, delta_ell,
                                            pix_size, NN)
    DlTT2 .+= d

end

e, d = strip2.get_power_spectrum(conv_map, conv_map, ell_max, delta_ell,
                                        pix_size, NN)

plot(ell, DlTT, yaxis = :log, xlims = (0 ,1000), ylims=(100, 1000000))
plot!(e, d, yaxis = :log, xlims = (0 ,1000), ylims=(100,1000000))
