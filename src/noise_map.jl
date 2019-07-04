function white_noise(N::Int, pix_size::Float64, wnl::Float64)

    rng = MersenneTwister()
    white_noise_map = randn(rng, N, N) .* (wnl / pix_size)

    return white_noise_map

end


function atmospheric_noise(N::Int, pix_size::Float64, anl::Float64)

    rng = MersenneTwister()
    random_realization = randn(rng, N, N)

    inds = ((0:(N-1)) .+ 0.5 .- (N/2))
    X = reshape(repeat(inds, N), N ,N)
    Y = X'
    R = ((X .^2 .+ Y .^2).^0.5) .* (pix_size / 60.0)

    mag_k = ((2.0 * π) ./ (R .+ 0.01 )).^(5.0 / 3.0)
    atm_rand = fft(random_realization)
    atm_map = ifft(atm_rand .* fftshift(mag_k))
    atm_map_norm = atm_map .* (anl / pix_size)

    return real.(atm_map_norm)

end


function one_over_f(N::Int, pix_size::Float64, oofnl::Float64)

    rng = MersenneTwister()
    random_realization = randn(rng, N, N)

    inds = ((0:(N-1)) .+ 0.5 .- (N/2))
    X = reshape(repeat(inds, N), N ,N)' .* (pix_size / 60.0)
    kx = (2.0 * π) ./ (X .+ 0.01)
    rand_ft = fft(random_realization)
    oofn_map = ifft(rand_ft .* fftshift(kx))
    oofn_map_norm = oofn_map .* (oofnl / pix_size)

    return real.(oofn_map_norm)

end

function one_over_f2(N::Int, pix_size::Float64, oofnl::Float64)

    rng = MersenneTwister()
    random_realization = randn(rng, N, N)

    inds = ((0:(N-1)) .+ 0.5 .- (N/2))
    X = reshape(repeat(inds, N), N ,N)' .* (pix_size / 60.0)
    kx = ((2.0 * π) ./ (X .+ 0.01)).^(3)
    rand_ft = fft(random_realization)
    oofn_map = ifft(rand_ft .* fftshift(kx))
    oofn_map_norm = oofn_map .* (oofnl / pix_size)

    return real.(oofn_map_norm)

end
