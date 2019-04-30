include("fft_util.jl")


function observe_sky(N_pix::Int, sky::Array{Float64, 2})

    pixels = range(0, N_pix*N_pix, length = N_pix*N_pix)
    x_pix = Int.(round.((N_pix-1) .* abs.(sin.( pixels .*
                                                (π/(2*N_pix))) )) .+ 1)

    y_pix = Int.(round.(reshape(reshape(repeat(range(
                                   1, N_pix, length = N_pix),
                             N_pix),
                    N_pix, N_pix)',
            1, N_pix*N_pix)))
    tod_d1 = [sky[x_pix[i], y_pix[i]] for i = 1:N_pix*N_pix]
    tod_d2 = [sky[y_pix[i], x_pix[i]] for i = 1:N_pix*N_pix]

    return tod_d1, tod_d2

end


function observe_sky(N_pix::Int, sky::Array{Float64, 2},
                     noise::Bool)
    tod_d1, tod_d2  = observe_sky(N_pix, sky)

    rng  = MersenneTwister()

    rand_d1 = fft(randn(rng, length(tod_d1)))
    noise_spec_d1 = tod_noise(length(tod_d1), 0.00416, 0.1, 3.0, 40.0)
    println(size(rand_d1), size(noise_spec_d1))
    fnoise_d1 = rand_d1 .* (noise_spec_d1.^0.5)
    tod_d1 += real.(ifft(fnoise_d1))

    rand_d2 = fft(randn(rng, length(tod_d2)))
    noise_spec_d2 = tod_noise(length(tod_d2), 0.00416, 0.1, 3.0, 40.0)
    fnoise_d2 = rand_d2 .* (noise_spec_d2.^0.5)
    tod_d2 += real.(ifft(fnoise_d2))

    return tod_d1, tod_d2

end


function tod_noise(Nsamp::Int, dt::Float64, fknee::Float64, alpha::Float64,
                   sigma::Float64)
    Freq = abs.(fftfreq(Nsamp, dt))
    return (1 .+ ( max.(Freq, Freq[2]) ./ fknee).^(-alpha)) .* (sigma^2)
end
