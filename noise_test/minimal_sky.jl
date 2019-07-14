include("../src/strip2.jl")
NN = 1024
pix_size = 0.5

# White noise: evaluation of σ
dc = 0.35
observation = 24 # months
cmb = 1.1 # K at 95 GHz
system = 8.5 # K
atmosphere = 25 # K
bandwidth  = 31.5 # GHz
Tsys = cmb + system + atmosphere
Number_of_sensors = 1000.0
time_per_pixel = (observation * 30 * 24 * 60 * 60 * dc) / (NN^2)  # sec/pix

# Ma osservo 10 volte
time_per_pixel /=  10.0

sigma = ((Tsys / sqrt(31.5*1E9 * time_per_pixel )) / sqrt(Number_of_sensors)) * 1e6

# SZ Effect
N_SZ_Clusters  = 100
A_SZ_Clusters = 150.0
SZ_beta = 0.86
SZ_Theta_core = 1.0

sz_map = strip2.SZ_source(NN, pix_size, N_SZ_Clusters,
                   A_SZ_Clusters, SZ_beta,
                   SZ_Theta_core)

# Overall map (cosmological sources and noise) convolved with the
# telescope's beam
conv_map = sz_map

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
