include("../src/strip2.jl")

#using strip2
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
m = readdlm("power_spectrum/CAMB_totalCls.dat")

ell  = m[:, 1]
DlTT = m[:, 2]
DlEE = m[:, 3]
DlBB = m[:, 4]
DlTE = m[:, 5]

Map_T, Map_Q, Map_U, Map_E, Map_B = strip2.make_CMB_pol_maps(NN, pix_size, ell,
                                                             DlTT, DlEE, DlTE, DlBB)

gradient = ColorGradient([:blue, :white, :red])
heatmap( Map_T, c=gradient, xlabel = "x [px]",
                   ylabel = "y [px]",
                   #clims=(-40, 4),
                   size=(580,480))

p1 = strip2.Plot_Sky(NN, pix_size, strip2.convolve_beam(NN,pix_size, beam_waist, Map_Q), cmin = -20, cmax = 20, lab_title="Q Map")
p2 = strip2.Plot_Sky(NN, pix_size, strip2.convolve_beam(NN,pix_size, beam_waist, Map_U), cmin = -20, cmax = 20, lab_title="U Map")
p3 = strip2.Plot_Sky(NN, pix_size, Map_E, cmin = -20, cmax = 20)
p4 = strip2.Plot_Sky(NN, pix_size, Map_B, cmin = -2, cmax = 2)

plot(p1, p2, p3, p4, layout = grid(2, 2), size = (1050, 850))

Map_B_convolved = strip2.convolve_beam(NN, pix_size, beam_waist, Map_B)

strip2.Plot_Sky(NN, pix_size, Map_B_convolved, cmin=-2, cmax=2)

ell_max = 5000.0
delta_ell = 30.0

win_b = strip2.windowing!(NN, Map_B)
win_b_con = strip2.windowing!(NN, Map_B_convolved)


ell_b, DlBB_obs = strip2.get_power_spectrum(win_b, win_b,
                                            ell_max, delta_ell, pix_size, NN)

ell_conv, DlBB_conv = strip2.get_power_spectrum(win_b_con, win_b_con,
          ell_max, delta_ell, pix_size, NN)




a=plot(ell, DlTT, yaxis=:log, xaxis=:log, label="TT", ylabel="Cl(l(l+1))[muK ^2]", xlabel="l")
b=plot!(ell, DlEE, yaxis=:log, xaxis=:log, label="EE")
c=plot!(ell, DlBB, yaxis=:log, xaxis=:log, label="BB", xlims=(1,2500))

png(c,"ciao")
