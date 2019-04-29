include("../src/strip2.jl")
#using strip2
using DelimitedFiles
using Plots
gr()

# General constants
NN = 1024
pix_size = 0.5
beam_waist = 10.0

# SZ Effect
N_SZ_Clusters  = 500
A_SZ_Clusters = 50.0
SZ_beta = 0.86
SZ_Theta_core = 1.0

# Point Sources - poisson distribution
N_Sources_POIS = 5000
A_Sources_POIS = 200.0

# Point Sources - exponential distrubution
N_Sources_EXP = 50
A_Sources_EXP = 1000.0

# Read the CAMB cmb power spectrum.
m = readdlm("power_spectrum/CAMB_fiducial_cosmo_scalCls.dat")

ell  = m[:, 1]
DlTT = m[:, 2]

# Create CMB temperature map
cmb_T_map = strip2.make_CMB_T_map(NN, pix_size, ell, DlTT)

# Add the SZ contribute
sz_map = strip2.SZ_source(NN, pix_size, N_SZ_Clusters,
                   A_SZ_Clusters, SZ_beta,
                   SZ_Theta_core)

# Add the point sources contribute
point_pois = strip2.PS_poisson(NN, pix_size, N_Sources_POIS, A_Sources_POIS)
point_expo = strip2.PS_exp(NN, pix_size, N_Sources_EXP, A_Sources_EXP)


mapp = cmb_T_map .+ sz_map .+ point_pois .+ point_expo
conv_map = strip2.convolve_beam(NN, pix_size, beam_waist, mapp)

# gradient = ColorGradient([:blue, :white, :red])
# heatmap( conv_map, c=gradient, xlabel = "x [px]",
#                          ylabel = "y [px]",
#                          clims=(-400, 400),
#                          size=(580,480)
#         )

tod_d1, tod_d2 = strip2.observe_sky(NN, conv_map)
plot(tod_d1, ylims=(-400, 400))
plot!(tod_d2)
