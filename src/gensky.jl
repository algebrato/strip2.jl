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
