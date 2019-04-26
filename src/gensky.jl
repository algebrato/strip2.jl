"""
N           # pixel number in linear dimension
pix_size    # size of a pixel in arcminutes
X_width     # horizontal map width in degrees
Y_width     # vertical map width in degrees
"""


function make_CMB_T_map(N::Int, pix_size::Float64,
                        ell::AbstractArray{<:Number},
                        DlTT::AbstractArray{<:Number})

    ClTT = DlTT .* 2 .* π ./ (ell .* (ell .+ 1.) )
    inds = ((0:(N-1)) .+ .5 .- ((N-1) / 2.0)) ./ (N-2.0)

    X  = reshape(repeat(range(-0.5, 0.5, length=N), N), N, N)
    Y = X'
    R = (X .^2 .+ Y .^2) .^0.5

    pix_to_rad = pix_size / 60. * π  / 180.  # Scale factor between pix to rad
    ell_scale_factor = 2. * π / pix_to_rad

    ell2d = R .* ell_scale_factor
    max_l_R = Int64(floor(maximum(ell2d))+1)

    ClTT_ex = vcat(ClTT, zeros(max_l_R - size(ClTT)[1]))
    indici = trunc.(Int64, ell2d)

    CLTT2d = ClTT_ex[indici]

    rng = MersenneTwister()
    T_rand_arr = randn(rng, N, N)
    fft_T_rand = fft(T_rand_arr)
    FT_sky = (CLTT2d .^0.5) .* fft_T_rand
    T_cmb_map = ifft(ifftshift(FT_sky))
    T_cmb_map = T_cmb_map ./ (pix_size / 60.0 * π / 180.0)

    return real.(T_cmb_map)
end


function SZ_source(N::Int, pix_size::Float64, Number_of_SZ_Clusters::Int,
                   Mean_Amplitude_of_SZ_Clusters::Float64,
                   SZ_beta::Float64, SZ_Theta_core::Float64)

    # Make a realization of a naive SZ_map
    rng = MersenneTwister()
    SZMap = zeros(N, N)
    SZCat = zeros(3, Number_of_SZ_Clusters)
    i = 0
    for i = 1:Number_of_SZ_Clusters
        pix_x = Int64(trunc((N-1)*rand(rng) + 1))
        pix_y = Int64(trunc((N-1)*rand(rng) + 1))
        pix_amplitude = (-1) * Mean_Amplitude_of_SZ_Clusters * randexp(rng)

        SZCat[1, i] += pix_x
        SZCat[2, i] += pix_y
        SZCat[3, i] += pix_amplitude

        SZMap[pix_x, pix_y] += pix_amplitude
    end

    beta = beta_function(N, pix_size, SZ_beta, SZ_Theta_core)
    FT_beta = fft(fftshift(beta))
    FT_SZMap = fft(fftshift(SZMap))
    conv = FT_beta .* FT_SZMap
    SZMap = fftshift(real.(ifft(conv)))

    return SZMap
end


function beta_function(N::Int, pix_size::Float64, SZ_beta::Float64,
                       SZ_Theta_core::Float64)

    inds = ((0:(N-1)) .+ 0.5 .- (N/2)) * pix_size
    X = reshape(repeat(inds, N), N ,N)
    Y = X'
    R = (X .^2 .+ Y .^2) .^0.5
    beta = (1 .+ (R ./ SZ_Theta_core).^2 ) .^ ((1.0 - 3.0 * SZ_beta) / 2.0)

    return beta

end


function PS_poisson(N::Int, pix_size::Float64, N_Sources::Int,
                    A_Sources::Float64)
    PS_Pois = zeros(N, N)
    rng = MersenneTwister()
    for i = 1:N_Sources
        pix_x = Int64(trunc((N-1)*rand(rng) + 1))
        pix_y = Int64(trunc((N-1)*rand(rng) + 1))
        PS_Pois[pix_x, pix_y] += A_Sources * randn(rng)
    end

    return PS_Pois

end


function PS_exp(N::Int, pix_size::Float64, N_Sources::Int,
                    A_Sources::Float64)
    PS_expo = zeros(N, N)
    rng = MersenneTwister()
    for i = 1:N_Sources
        pix_x = Int64(trunc((N-1)*rand(rng) + 1))
        pix_y = Int64(trunc((N-1)*rand(rng) + 1))
        PS_expo[pix_x, pix_y] += A_Sources * randexp(rng)
    end

    return PS_expo

end
