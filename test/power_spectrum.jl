include("random_sky.jl")
include("gen_noise_map.jl")
using ProgressMeter

e_des = Array{Float64, 1}(undef,100)
d_des = Array{Float64, 1}(undef,100)

for i in 1:10
    global dir=false
    mappa = zeros(NN, NN)
    hit   = zeros(NN, NN)
    N_baselines = 10
    dataset_baselines = Array{Main.strip2.TOD_fake_pointing, 1}(undef, N_baselines)

    for i = 1:N_baselines
        dataset = strip2.observe_sky(NN, conv_map; noise=true, direction_l_r=dir)
        m, h = strip2.binned_map(NN, dataset.time_order_data, dataset.pointing)
        mappa .+= m
        hit   .+= h
        dataset_baselines[i] = dataset
        global dir =  true ⊻ dir
    end

    map_des = strip2.destriper(NN, dataset_baselines)
    map_d = reshape(map_des, NN, NN)

    ell_max = 5000.0
    delta_ell = 50.0

    win_map_dest = strip2.windowing!(NN, (map_d .- avg_dest))

    a, b = strip2.get_power_spectrum(win_map_dest, win_map_dest,
                                            ell_max, delta_ell, pix_size, NN)
    e_des .+= a
    d_des .+= b
end

e_des ./= 10.0
d_des ./= 10.0

Multiply = DlTT[Int.(round.(e_des)) .- 1] ./ d_des


R = 0.0001:0.0001:0.1

mag_k = ((2.0 * π) ./ (R .+ 0.01 )).^(5.0 / 3.0)
pix_to_rad = pix_size / 60. * π  / 180.  # Scale factor between pix to rad
ell_scale_factor = 2. * π / pix_to_rad
ell_atm = R .* ell_scale_factor



(1 .+ ( max.(Freq, Freq[2]) ./ fknee).^(-alpha)) .* (sigma^2)

p(f)=(1 + f/f_knee )^(-α) * σ







pp=plot(ell, DlTT./maximum(DlTT), label="CMB Temp. normalized", yaxis=:log, xlims =(0.00001, 1000), ylims=(0.1,1), xlabel="l", ylabel="Dl/max(Dl)", size=(500,500))
pp=plot!(ell, DlBB./maximum(DlBB), label="CMB BB normalized", size=(500,500))
pp=plot!(ell_atm, mag_k ./ maximum(mag_k), label="Atmospheric normalized", size=(500,500))
