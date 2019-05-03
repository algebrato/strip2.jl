include("../src/strip2.jl")
# using strip2
using DelimitedFiles
using Plots
gr()
# Strip2: 54deg x 54deg
# General constants
NN = 1024
pix_size = 0.5
beam_waist = 9.0

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
DlEE = m[:, 3]
DlBB = m[:, 4]
DlTE = m[:, 5]

Map_T, Map_Q, Map_U, Map_E, Map_B = strip2.make_CMB_pol_maps(NN, pix_size, ell,
                                                             DlTT, DlEE, DlTE, DlBB)

gradient = ColorGradient([:blue, :white, :red])
heatmap( Map_Q, c=gradient, xlabel = "x [px]",
                   ylabel = "y [px]",
                   #clims=(-40, 4),
                   size=(580,480))
