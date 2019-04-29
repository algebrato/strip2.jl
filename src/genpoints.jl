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
