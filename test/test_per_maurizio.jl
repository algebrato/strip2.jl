using Plots

grad =  ColorGradient([:white, :black])
strip2.Plot_Sky(NN, pix_size, conv_map, cmin=-100, cmax=50, gradient=grad)


noise_ss = strip2.tod_noise(length(tod), 0.0083, 0.1, 3.0, 20.76)
dataset_baselines_2 = Array{Main.strip2.TOD_fake_pointing, 1}(undef, 1)
dat = strip2.TOD_fake_pointing(tod, puntamenti, noise_ss)
dataset_baselines_2[1] = dat

# Filtro semplice - In questo caso la sorgente cosmologica viene uccisa e portata
# molto vicino al livello di wite noise
tod_denoise = strip2.denoise(tod, noise_ss)
mappa_bin, hit = strip2.binned_map(NN, tod_denoise, puntamenti)
strip2.Plot_Sky(NN, pix_size, mappa_bin ./ hit, cmin=-0.1, cmax=0.05, gradient=grad)


# Ma se risolvo il problema del mapmaker usando la funzione strip2.denoise
# per approssimare la matrice di covarianza allora ottengo i livelli corretti
# Ci mette molto a convergere perche` non ho praticamente nessun cross-angle (worst-case)
map_des, exit_status, iter = strip2.destriper(NN, dataset_baselines_2)
map_d = reshape(map_des, NN, NN)
strip2.Plot_Sky(NN, pix_size, map_d, cmin=-100, cmax=50, gradient=grad)


# Se adotto una scanning strategy con dei cross-angles molto densa:
# Rifare i punti prima con una ss diversa.
