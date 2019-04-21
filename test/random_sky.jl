include("../src/Strip2.jl")
using DelimitedFiles
using Plots
gr()


NN = 1024
NSIDE = 64
pix = 0.5
m = readdlm("power_spectrum/CAMB_fiducial_cosmo_scalCls.dat")

ell  = m[:, 1]
DlTT = m[:, 2]

map = Strip2.make_CMB_T_map(NN, pix, ell, DlTT)
gradient = ColorGradient([:blue, :white, :red])
heatmap(map, c=gradient, xlabel = "x [px]",
                         ylabel = "y [px]",
                         clims=(-400, 400)
        )
