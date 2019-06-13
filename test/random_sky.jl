include("../src/strip2.jl")
#using strip2
using DelimitedFiles
#using Plots
#gr()
# Strip2: 54deg x 54deg
# General constants
NN = 1024
pix_size = 0.5
beam_waist = 9.0

# White noise: evaluation of σ
dc = 0.35
observation = 24 # months
cmb = 1.1 # at 95 GHz
system = 18.7
atmosphere = 102
bandwidth  = 31.5 # GHz
Tsys = cmb + system + atmosphere
Number_of_sensors = 500.0
time_per_pixel = (observation * 30 * 24 * 60 * 60 * dc) / (NN^2)  # sec/pix

# Ma osservo 10 volte
time_per_pixel /=  10.0


sigma = ((Tsys / sqrt(31.5*1E9 * time_per_pixel )) / sqrt(Number_of_sensors)) * 1e6



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

# Overall map (cosmological sources and noise) convolved with the
# telescope's beam
mapp = cmb_T_map .+ sz_map .+ point_pois .+ point_expo
conv_map = strip2.convolve_beam(NN, pix_size, beam_waist, mapp)

# gradient = ColorGradient([:blue, :white, :red])
# heatmap( conv_map, c=gradient, xlabel = "x [px]",
#                          ylabel = "y [px]",
#                          clims=(-400, 400),
#                          size=(580,480)
#         )
