function filter_noise(N::Int, N_mask::Float64, Map::Array{Float64, 2})

    inds = ((0:(N-1)) .+ 0.5 .- (N/2))
    X = reshape(repeat(inds, N), N ,N)
    Y = X'
    R = ((X .^2 .+ Y .^2).^0.5)

    mask = ones(N, N)
    mask[abs.(Y) .< N_mask] .= 0

    FMap = fftshift(ifft(fftshift(Map)))
    FMap_filtered = FMap .* mask
    Map_filtered = real.(fftshift(fft(fftshift(FMap_filtered))))

    return Map_filtered

end
