function windowing!(N::Int, Map::Array{Float64, 2}; type="default")

    inds = (((0:(N-1)) .+ 0.5 .- (N/2)) ./ N ) .* π
    X = reshape(repeat(inds, N), N ,N)
    Y = X'

    window = cos.(X) .* cos.(Y)

    window .* Map

end


function get_power_spectrum(Map1::Array{Float64, 2}, Map2::Array{Float64, 2},
                            max_ell::Float64, delta_ell::Float64,
                            pix_size::Float64, N::Int)


    inds = ((0:(N-1)) .+ 0.5 .- (N/2)) ./ (N - 1.0)
    X  = reshape(repeat(range(-0.5, 0.5, length=N), N), N, N)
    Y = X'
    R = (X .^2 .+ Y .^2) .^0.5
    pix_to_rad = pix_size / 60. * π  / 180.
    lsf = 2. * π / pix_to_rad

    ell2d = R .* lsf

    N_bins = Int(round(max_ell / delta_ell))

    ell_array = Array{Float64, 1}(undef, 0)
    Cl_array  = Array{Float64, 1}(undef, 0)

    FMap1 = ifft(fftshift(Map1))
    FMap2 = ifft(fftshift(Map2))
    PSMap =  fftshift(real.(conj.(FMap1) .* FMap2))

    for i = 1:N_bins

        l = (i - 0.5) * delta_ell
        append!(ell_array, l)
        Mask = (ell2d .>= ((i-1)*delta_ell)) .* (ell2d .< (i*delta_ell))
        val_cl = PSMap[Mask]
        Cl = sum(val_cl) / length(val_cl)
        append!(Cl_array, Cl)

    end

    Dl_array = Cl_array ./ (2.0 .* π ./ (ell_array .* (ell_array .+ 1.) ))

    return ell_array, Dl_array

end
