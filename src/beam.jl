function convolve_beam(NN::Int, pix_size::Float64, beam_waist::Float64,
                       Map::Array{Float64, 2})

    gaussian = make_2d_gaus(NN, pix_size, beam_waist)

    FT_beam = fft(fftshift(gaussian))
    FT_Map  = fft(fftshift(Map))
    conv = FT_beam .* FT_Map

    conv_map = fftshift(real.(ifft(conv)))

    return conv_map

end


function make_2d_gaus(N::Int, pix_size::Float64, beam_waist::Float64)

    inds = ((0:(N-1)) .+ 0.5 .- (N/2)) * pix_size
    X = reshape(repeat(inds, N), N ,N)
    Y = X'
    R = (X .^2 .+ Y .^2) .^0.5

    sigma = beam_waist / sqrt(8. * log(2))
    gauss = exp.(-0.5 .* ((R ./ sigma).^2.0))
    Norm_Fact = sum(gauss)
    gauss_norm = gauss ./ Norm_Fact

    return gauss_norm

end
