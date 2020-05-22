include("../src/strip2.jl")
#using strip2
using DelimitedFiles
#using Plots
#gr()
# Strip2: 54deg x 54deg
# General constants
NN = 4024
pix_size = 0.5
beam_waist = 11.0
# White noise: evaluation of σ

dc = 0.35
observation = 24 # months
cmb = 1.1 # at 95 GHz
system = 8.5
atmosphere = 25 #102
bandwidth  = 31.5 # GHz
Tsys = cmb + system + atmosphere
Number_of_sensors = 1000.0
time_per_pixel = (observation * 30 * 24 * 60 * 60 * dc) / (NN^2)  # sec/pix

# Ma osservo 10 volte
time_per_pixel /=  10.0




sigma = ((Tsys / sqrt(31.5*1E9 * time_per_pixel )) / sqrt(Number_of_sensors)) * 1e6



# SZ Effect
N_SZ_Clusters  = 100
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
m = readdlm("power_spectrum/CAMB_totalCls.dat")

ell  = m[:, 1]
DlTT = m[:, 2]
DlEE = m[:, 3]
DlBB = m[:, 4]
DlTE = m[:, 5]

Tm, Qm, Um, Em, Bm = strip2.make_CMB_pol_maps(NN, pix_size, ell,
                                              DlTT, DlEE, DlTE, DlBB)
# Add the SZ contribute
sz_map = strip2.SZ_source(NN, pix_size, N_SZ_Clusters,
                   A_SZ_Clusters, SZ_beta,
                   SZ_Theta_core)

# Add the point sources contribute
point_pois = strip2.PS_poisson(NN, pix_size, N_Sources_POIS, A_Sources_POIS)
point_expo = strip2.PS_exp(NN, pix_size, N_Sources_EXP, A_Sources_EXP)

# Overall map (cosmological sources and noise) convolved with the
# telescope's beam
# mapp = Tm .+ sz_map .+ point_pois .+ point_expo
mapp = sz_map
conv_map = strip2.convolve_beam(NN, pix_size, beam_waist, mapp)

gradient=ColorGradient([:blue, :white, :red])
p2 = strip2.Plot_Sky(NN, pix_size, conv_map, cmin = -400, cmax = 400, lab_title="T Map")
