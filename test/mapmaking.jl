using Plots
using BenchmarkTools

include("random_sky.jl")

dir=false

mappa = zeros(NN, NN)
hit   = zeros(NN, NN)
N_baselines = 10


dataset_baselines = Array{Main.strip2.TOD_fake_pointing, 1}(undef, N_baselines)

for i = 1:N_baselines
    dataset = strip2.observe_sky(NN, conv_map; noise=true, direction_l_r=dir)
    m, h = strip2.binned_map(NN, dataset.time_order_data, dataset.pointing)
    mappa .+= m
    hit   .+= h
    dataset_baselines[i] = dataset
    global dir =  true ⊻ dir
end

obs_map  = strip2.Plot_Sky(NN, pix_size, mappa ./ hit, cmin=-600, cmax=600)
true_map = strip2.Plot_Sky(NN, pix_size, conv_map)

plot(obs_map, true_map, layout = grid(2,1),size=(590,940))

map_des = strip2.destriper(NN, dataset_baselines)
map_d = reshape(map_des, NN, NN)
p2 = strip2.Plot_Sky(NN, pix_size, conv_map .- (sum(conv_map)/(NN*NN)), cmin = -400, cmax = 400, lab_title="T Map")
p3 = strip2.Plot_Sky(NN, pix_size, map_d .- (sum(map_d)/(NN*NN)) , cmin=-400, cmax=400, lab_title="Destriped Map")
p1 = strip2.Plot_Sky(NN, pix_size, map_d .- (sum(map_d)/(NN*NN)) .- (conv_map .- (sum(conv_map)/(NN*NN))) , cmin=-8, cmax=8, lab_title="Differences between destriped and true maps")

plot(p1, p2, layout = grid(2,1), size=(590,940))

#strip2.Plot_Sky(NN, pix_size, map_d, cmin = -400, cmax = 400)

residui =  map_d .- (sum(map_d)/(NN*NN)) .- (conv_map .- (sum(conv_map)/(NN*NN)))

println("RMS: ",sqrt(sum((residui) .^2) / (NN*NN)))
