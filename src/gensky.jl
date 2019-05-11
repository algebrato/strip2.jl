"""
    function make_CMB_T_map(N::Int, pix_size::Float64,
                    ell::AbstractArray{<:Number},
                    DlTT::AbstractArray{<:Number})

# Brief description
Given a CAMB power spectrum return a random CMB
only-temperature map in the`Array{Float64, 2}` format.

# Development
This functions has not presented weird behaviours.
The function return a random CMB temperature anisotropies
maps, that you can draw using the `heatmap` method in
Plots julia's package.

# Args description
N           # pixel number in linear dimension
pix_size    # size of a pixel in arcminutes
ell         # Multipole value - angular scale
DlTT        # Power spectrum value Dₗ = Cₗ⋅(l⋅(l+1)) (μk)²
              Temperature  Temperature
"""
function make_CMB_T_map(N::Int, pix_size::Float64,
                        ell::AbstractArray{<:Number},
                        DlTT::AbstractArray{<:Number})

    ClTT = DlTT .* 2 .* π ./ (ell .* (ell .+ 1.) )
    ClTT[1] = 0

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


"""
function make_CMB_pol_maps(N::Int, pix_size::Float64,
                           ell::AbstractArray{<:Number},
                           DlTT::AbstractArray{<:Number},
                           DlEE::AbstractArray{<:Number},
                           DlTE::AbstractArray{<:Number},
                           DlBB::AbstractArray{<:Number})

# Brief description
Given a CAMB power spectrum, in function of the angular scale
(`ell` scale), return a complete set of maps of a hypothetical
random sky. The maps has return concern: Temperature (T_Map),
Q_Stokes_parameter (Q_Map) and U_Stokes_parametes (U_Map).

Whithin the function, the Q and U maps are turn into E and B maps.
On this last two maps you can use the power spectra analysis tools
to extract the BB and EE power_spectrum after the introduction
of the tlescope systematics errors.

# Development
Actually, the function do not work. I have decided to add it
beacouse I think that the issue that I have noticed on this
function has due to CAMB power spectrum format. More test are needed
in particular testing this function with different power-spectrum set.


# Args description
N           # pixel number in linear dimension
pix_size    # size of a pixel in arcminutes
ell         # Multipole value - angular scale
DlTT        # Power spectrum value Dₗ = Cₗ⋅(l⋅(l+1)) (μk)²
              Temperature  Temperature (μk)²
DlEE        # E-modes auto-power_spectrum  polarization  (μk)²
DlTE        # Cross-power_spectrum between Temperature and E modes (μk)²
DlBB        # B-modes auto-power_spectrum (μk)²

"""
function make_CMB_pol_maps(N::Int, pix_size::Float64,
                           ell::AbstractArray{<:Number},
                           DlTT::AbstractArray{<:Number},
                           DlEE::AbstractArray{<:Number},
                           DlTE::AbstractArray{<:Number},
                           DlBB::AbstractArray{<:Number})

    ClTT = DlTT .* 2 .* π ./ (ell .* (ell .+ 1.) )
    ClEE = DlEE .* 2 .* π ./ (ell .* (ell .+ 1.) )
    ClTE = DlTE .* 2 .* π ./ (ell .* (ell .+ 1.) )
    ClBB = DlBB .* 2 .* π ./ (ell .* (ell .+ 1.) )

    ClTT[1:2] .= 0
    ClEE[1:2] .= 0
    ClTE[1:2] .= 0
    ClBB[1:2] .= 0

    # correlation level between T and E-modes
    corr_level_E = ClTE ./ sqrt.(ClTT)
    uncorr_level_E = ClEE .- ( (ClTE .^2.0) ./ ClTT )

    corr_level_E[1:2] .= 0
    uncorr_level_E[1:2] .= 0

    X  = reshape(repeat(range(-0.5, 0.5, length=N), N), N, N)
    Y = X'
    R = (X .^2 .+ Y .^2) .^0.5
    ang = atan.(Y, X)

    pix_to_rad = pix_size / 60. * π  / 180.  # Scale factor between pix to rad
    ell_scale_factor = 2. * π / pix_to_rad

    ell2d = R .* ell_scale_factor
    max_l_R = Int64(floor(maximum(ell2d))+1)

    ClTT_ex = vcat(ClTT, zeros(max_l_R - size(ClTT)[1]))

    ClEE_uncor_ex = vcat(uncorr_level_E,
                         zeros(max_l_R - size(uncorr_level_E)[1]))

    ClEE_corr_ex = vcat(corr_level_E,
                         zeros(max_l_R - size(corr_level_E)[1]))

    ClBB_ex = vcat(ClBB, zeros(max_l_R - size(ClBB)[1]))

    indici = trunc.(Int64, ell2d)
    CLTT2d = ClTT_ex[indici]
    ClEE_uncor2d = ClEE_uncor_ex[indici]
    ClEE_cor2d = ClEE_corr_ex[indici]
    CLBB2d = ClBB_ex[indici]


    rng = MersenneTwister()
    T_rand_arr = randn(rng, N, N)
    E_rand_arr = randn(rng, N, N)
    B_rand_arr = randn(rng, N, N)

    fft_T_rand = fft(T_rand_arr)
    fft_E_rand = fft(E_rand_arr)
    fft_B_rand = fft(B_rand_arr)

    FT_sky = (CLTT2d .^0.5) .* fft_T_rand
    FE_sky = (real.((ClEE_uncor2d .+ 0im) .^0.5)) .* fft_E_rand .+ # Error!
             (real.((ClEE_uncor2d .+ 0im) .^0.5)) .* fft_T_rand    # Error!
    FB_sky = (real.((CLBB2d .+ 0im) .^0.5)) .* fft_B_rand          # Error!

    # Conver the E and B maps into Q and U maps
    FQ_sky = (FE_sky .* cos.(2.0 .* ang)) .- (FB_sky .* sin.(2.0 .* ang))
    FU_sky = (FE_sky .* sin.(2.0 .* ang)) .+ (FB_sky .* cos.(2.0 .* ang))

    T_cmb_map = ifft(ifftshift(FT_sky))
    Q_cmb_map = ifft(ifftshift(FQ_sky))
    U_cmb_map = ifft(ifftshift(FU_sky))

    E_map = ifft(ifftshift(FE_sky))
    B_map = ifft(ifftshift(FB_sky))

    T_cmb_map = T_cmb_map ./ (pix_size / 60.0 * π / 180.0)
    Q_cmb_map = Q_cmb_map ./ (pix_size / 60.0 * π / 180.0)
    U_cmb_map = U_cmb_map ./ (pix_size / 60.0 * π / 180.0)
    E_map = E_map ./ (pix_size / 60.0 * π / 180.0)
    B_map = B_map ./ (pix_size / 60.0 * π / 180.0)

    return real.(T_cmb_map), real.(Q_cmb_map), real.(U_cmb_map), real.(E_map), real.(B_map)

end


"""
    function SZ_source(N::Int, pix_size::Float64, Number_of_SZ_Clusters::Int,
                       Mean_Amplitude_of_SZ_Clusters::Float64,
                       SZ_beta::Float64, SZ_Theta_core::Float64)

# Brief description
SZ_source return a sky map Array{Float64, 2} formatted, of the only
SZ sources. The sources are placed random on the map. The pattern is
the well-know SZ-profile given by the `beta_function`. The insensity
and dimensione are normally distribuited around a given mean-amplitude.

# Development
The function has not presented weird behaviours.

# Args description
N                      # pixel number in linear dimension
pix_size               # size of a pixel in arcminutes
Number_of_SZ_Clusters  # Number of SZ sources on the map
Mean_Amplitude         # Mean intensity of those sources
SZ_beta, SZ_Theta_core # Sz spectral indexes

"""
function SZ_source(N::Int, pix_size::Float64, Number_of_SZ_Clusters::Int,
                   Mean_Amplitude_of_SZ_Clusters::Float64,
                   SZ_beta::Float64, SZ_Theta_core::Float64)

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


"""
    function PS_poisson(N::Int, pix_size::Float64, N_Sources::Int,
                        A_Sources::Float64)
# Brief description
Return a sky-map Array{Float64, 2} formatted of the point sources
distribution on the observe sky. You can chose the number of sources
that are places randomly on the map.

This specific function return N_Sources, with amplitude A_Sources.
The amplitude is assigned using a Poisson extraction.
"""
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


"""
    function PS_exp(N::Int, pix_size::Float64, N_Sources::Int,
                    A_Sources::Float64)
# Brief description
Return a sky-map Array{Float64, 2} formatted of the point sources
distribution on the observe sky. You can chose the number of sources
that are places randomly on the map.

This specific function return N_Sources, with amplitude A_Sources.
The amplitude is assigned using a exponential extraction.
"""
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
