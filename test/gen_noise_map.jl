# Noise Map
using Plots
wnl = 12.0
anl = 0.15
oofnl = 0.2
oofnl2 = 0.0001

white_map = strip2.white_noise(NN, pix_size, wnl)
atmos_map = strip2.atmospheric_noise(NN, pix_size, anl)
one_ove_f = strip2.one_over_f(NN, pix_size, oofnl)
one_ove_f2 = strip2.one_over_f2(NN, pix_size, oofnl2)

noise_map = white_map .+ atmos_map .+ one_ove_f

sky = conv_map .+ noise_map

N_mask = 3.0

filtered_map = strip2.filter_noise(NN, N_mask, sky)


# white noise
strip2.Plot_Sky(NN, pix_size, white_map, cmin=-30, cmax=30)

# instrumental 1/f
strip2.Plot_Sky(NN, pix_size, one_ove_f, cmin=-40, cmax=40)

# atmospheric noise
grad = ColorGradient([:blue, :white])
strip2.Plot_Sky(NN, pix_size, strip2.convolve_beam(NN, pix_size, beam_waist, atmos_map), cmin=-100, cmax=100, gradient=grad,
                lab_title="Atmospheric Correlation Map")

# Raw map
strip2.Plot_Sky(NN, pix_size, conv_map .+ atmos_map, cmin=-500, cmax=500)

# Filtered Map
strip2.Plot_Sky(NN, pix_size, filtered_map .-conv_map, cmin=-400, cmax=400)

# @gif for iter in 1:10
#   atmos_map = strip2.atmospheric_noise(NN, pix_size, anl, iter)
#   grad = ColorGradient([:blue, :white])
#   strip2.Plot_Sky(NN, pix_size, atmos_map, cmin=-100, cmax=100, gradient=grad,
#                   lab_title="Atmospheric Correlation Map")
# end every 1
#
