
include("minimal_sky.jl")
using Plots

grad =  ColorGradient([:white, :black])
strip2.Plot_Sky(NN, pix_size, conv_map, cmin=-100, cmax=50, gradient=grad)


noise_ss = strip2.tod_noise(length(tod), 0.0083, 0.1, 3.0, 20.76)
noise_ss2 = strip2.tod_noise(length(tod), 0.0083, 1.0, 2.0, 20.76)
dataset_baselines_2 = Array{Main.strip2.TOD_fake_pointing, 1}(undef, 1)
dat = strip2.TOD_fake_pointing(tod, puntamenti, noise_ss)
dataset_baselines_2[1] = dat

# Filtro semplice - In questo caso la sorgente cosmologica viene uccisa e portata
# molto vicino al livello di wite noise
tod_denoise = strip2.denoise(tod, noise_ss)
tod_denoise2 = strip2.denoise(tod_denoise, noise_ss2)
mappa_bin, hit = strip2.binned_map(NN, tod_denoise, puntamenti)
mappa_bin2, hit2 = strip2.binned_map(NN, tod_denoise2, puntamenti)
strip2.Plot_Sky(NN, pix_size, mappa_bin ./ hit, cmin=-0.1, cmax=0.05, gradient=grad)
strip2.Plot_Sky(NN, pix_size, mappa_bin2 ./ hit2, cmin=-0.001, cmax=0.0005, gradient=grad)

# Ma se risolvo il problema del mapmaker usando la funzione strip2.denoise
# per approssimare la matrice di covarianza allora ottengo i livelli corretti
# Ci mette molto a convergere perche` non ho praticamente nessun cross-angle (worst-case)
map_des, exit_status, iter = strip2.destriper(NN, dataset_baselines_2)
map_d = reshape(map_des, NN, NN)
strip2.Plot_Sky(NN, pix_size, map_d, cmin=-100, cmax=50, gradient=grad)


# Se adotto una scanning strategy con dei cross-angles molto densa:
# Rifare i punti prima con una ss diversa.

tod_d = strip2.get_tod(map_d, puntamenti)
dataset_baselines_3 = Array{Main.strip2.TOD_fake_pointing, 1}(undef, 1)
dat2 = strip2.TOD_fake_pointing(tod_d, puntamenti, noise_ss2)
dataset_baselines_3[1] = dat2

map_des2, exit_status, iter = strip2.destriper(NN, dataset_baselines_3)
map_d2 = reshape(map_des2, NN, NN)
strip2.Plot_Sky(NN, pix_size, map_d2, cmin=-100, cmax=50, gradient=grad)
tod_d2 = strip2.get_tod(map_d2, puntamenti)

using DSP

pp = DSP.welch_pgram(tod, 5000, 120)
pp2 = DSP.welch_pgram(tod_denoise, 5000, 120)
pp3 = DSP.welch_pgram(tod_denoise2, 5000, 120)
pp4 = DSP.welch_pgram(tod_d, 5000, 120)
pp5 = DSP.welch_pgram(tod_d2, 5000, 120)

plot(pp.freq[2:end], pp.power[2:end], xscale=:log, yscale=:log, label="Raw TOD")
plot!(pp2.freq[2:end], pp2.power[2:end], xscale=:log, yscale=:log, label="Filter1")
plot!(pp3.freq[2:end], pp3.power[2:end], xscale=:log, yscale=:log, label="Filter2")
plot!(pp4.freq[2:end], pp4.power[2:end], xscale=:log, yscale=:log, label="MaxLike1")
plot!(pp5.freq[2:end], pp5.power[2:end], xscale=:log, yscale=:log, label="MaxLike2")

plot()
