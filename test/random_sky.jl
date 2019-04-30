include("../src/strip2.jl")
#using strip2
using DelimitedFiles
using Plots
gr()
# Strip2: 54deg x 54deg
# General constants
NN = 1024
pix_size = 0.5
beam_waist = 3.0

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


mapp = cmb_T_map #.+ sz_map .+ point_pois .+ point_expo
conv_map = strip2.convolve_beam(NN, pix_size, beam_waist, mapp)

# Noise Map
wnl = 10.0
anl = 0.15
oofnl = 0.2

white_map = strip2.white_noise(NN, pix_size, wnl)
atmos_map = strip2.atmospheric_noise(NN, pix_size, anl)
one_ove_f = strip2.one_over_f(NN, pix_size, oofnl)
noise_map = white_map .+ atmos_map .+ one_ove_f

# Combine the two maps
map_tot_cosm_noise = conv_map #.+ noise_map

# Very dummy filtering
N_mask = 4.0
filtered_map = strip2.filter_noise(NN, N_mask, map_tot_cosm_noise)

ell_max = 5000.0
delta_ell = 50.0

bin_spec_cur = zeros(Int(round(ell_max/delta_ell)))
Map_win = strip2.windowing!(NN, map_tot_cosm_noise)

# Non Funziona Cazzo!

for i = 1:1
    cmb_T_map = strip2.make_CMB_T_map(NN, pix_size, ell, DlTT)
    conv_map = strip2.convolve_beam(NN, pix_size, beam_waist, cmb_T_map)
    Map_win = strip2.windowing!(NN, conv_map)
    ell2, DlTT2 = strip2.get_power_spectrum(Map_win, Map_win, ell_max, delta_ell, pix_size, NN)
    bin_spec_cur .+= DlTT2
    println("Itration ", i," completed")
end

bin_spec_cur2 = bin_spec_cur ./ 1.00
Mult_factor = DlTT[Int.(round.(ell2)) .- 1] ./ bin_spec_cur2

#tod_d1, tod_d2 = strip2.observe_sky(NN, conv_map, true)
#plot(tod_d1)#, ylims=(-400, 400))
#plot!(tod_d2)
plot(ell, DlTT, yaxis = :log)
plot!(ell2, Mult_factor)

cmb_T_map = strip2.make_CMB_T_map(NN, pix_size, ell, DlTT)
conv_map = strip2.convolve_beam(NN, pix_size, beam_waist, cmb_T_map)
Map_win = strip2.windowing!(NN, conv_map)
ell2, DlTT2 = strip2.get_power_spectrum(Map_win, Map_win, ell_max, delta_ell, pix_size, NN)


plot!(ell2,  (DlTT2 .* Mult_factor .*316.0))




# gradient = ColorGradient([:blue, :white, :red])
# heatmap( Map_win, c=gradient, xlabel = "x [px]",
#                          ylabel = "y [px]",
#                          clims=(-400, 400),
#                          size=(580,480)
#         )
