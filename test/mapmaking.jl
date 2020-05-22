using Plots
using BenchmarkTools

include("random_sky.jl")

dir=false

mappa = zeros(NN, NN)
hit   = zeros(NN, NN)
N_baselines = 15


dataset_baselines = Array{Main.strip2.TOD_fake_pointing, 1}(undef, N_baselines)

for i = 1:N_baselines
    dataset = strip2.observe_sky(NN, conv_map; noise=true, direction_l_r=dir)
    m, h = strip2.binned_map(NN, dataset.time_order_data, dataset.pointing)
    mappa .+= m
    hit   .+= h
    dataset_baselines[i] = dataset
    global dir =  true ⊻ dir
end

tod = Array{Float64, 1}(undef, 0)
puntamenti = Array{Int64}(undef, 0)
for i in dataset_baselines
    append!(tod, i.time_order_data)
    append!(puntamenti, i.pointing[:,1])
end

for i in dataset_baselines
    append!(puntamenti, i.pointing[:,2])
end

puntamenti = reshape(puntamenti, NN*NN*10, 2)



#obs_map  = strip2.Plot_Sky(NN, pix_size, mappa ./ hit, cmin=-600, cmax=600)
#true_map = strip2.Plot_Sky(NN, pix_size, conv_map)

#plot(obs_map, true_map, layout = grid(2,1),size=(590,940))

map_des, a, b = strip2.destriper(NN, dataset_baselines, tolerance=1e-6)
map_d = reshape(map_des, NN, NN)
residui =  (map_d .- (sum(map_d)/(NN*NN))) .- (conv_map .- (sum(conv_map)/(NN*NN)))


#p2 = strip2.Plot_Sky(NN, pix_size, conv_map .- (sum(conv_map)/(NN*NN)), cmin = -400, cmax = 400, lab_title="T Map")
p3 = strip2.Plot_Sky(NN, pix_size, map_d .- (sum(map_d)/(NN*NN)) .-  (conv_map .- (sum(conv_map)/(NN*NN))), cmin=-20, cmax=20, lab_title="Destriped Map")
#p1 = strip2.Plot_Sky(NN, pix_size, map_d .- (sum(map_d)/(NN*NN)) .- (conv_map .- (sum(conv_map)/(NN*NN))) , cmin=-8, cmax=8, lab_title="Differences between destriped and true maps")

#plot(p1, p2, layout = grid(2,1), size=(590,940))

#strip2.Plot_Sky(NN, pix_size, map_d, cmin = -400, cmax = 400)
