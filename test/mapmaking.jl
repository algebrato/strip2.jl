using Plots
using BenchmarkTools

include("random_sky.jl")

dir=true

mappa = zeros(NN, NN)
hit   = zeros(NN, NN)
N_baselines = 5


dataset_baselines = Array{Main.strip2.TOD_fake_pointing, 1}(undef, N_baselines)

for i = 1:N_baselines
    dataset = strip2.observe_sky(NN, conv_map; noise=true, direction_l_r=dir)
    m, h = strip2.binned_map(NN, dataset.time_order_data, dataset.pointing)
    mappa .+= m
    hit   .+= h
    dataset_baselines[i] = dataset
    global dir =  true ⊻ dir
end

obs_map  = strip2.Plot_Sky(NN, pix_size, mappa ./ hit)
true_map = strip2.Plot_Sky(NN, pix_size, conv_map)

plot(obs_map, true_map, layout = grid(2,1),size=(590,940))

map_des = strip2.destriper(NN, dataset_baselines; n_iter=2)
map_d = reshape(map_des, NN, NN)
p2 = strip2.Plot_Sky(NN, pix_size, (mappa ./ hit) .- conv_map, cmin = -400, cmax = 400)
p1 = strip2.Plot_Sky(NN, pix_size, (map_d .- (sum(map_d)/NN*NN))  .- conv_map, cmin = -50, cmax = 50)

plot(p1, p2, layout = grid(2,1),size=(590,940))

strip2.Plot_Sky(NN, pix_size, map_d, cmin = -400, cmax = 400)

println("RMS: ",sqrt(sum((map_d .- conv_map) .^2)) / (NN*NN))
