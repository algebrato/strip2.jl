# Noise Map
using Plots
wnl = 10.0
anl = 0.15
oofnl = 0.2

white_map = strip2.white_noise(NN, pix_size, wnl)
atmos_map = strip2.atmospheric_noise(NN, pix_size, anl)
one_ove_f = strip2.one_over_f(NN, pix_size, oofnl)

noise_map = white_map .+ atmos_map .+ one_ove_f

gradient = ColorGradient([:blue, :white, :red])
heatmap( conv_map .+ noise_map, c=gradient, xlabel = "x [px]",
                         ylabel = "y [px]",
                         clims = (-400, 400),
                         size = (580, 480)
        )

N_mask = 4.0

filtered_map = strip2.filter_noise(NN, N_mask, atmos_map)

heatmap( filtered_map, c=gradient, xlabel = "x [px]",
                         ylabel = "y [px]",
                         clims=(-400, 400),
                         size = (580, 480)
        )
