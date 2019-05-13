function denoise(tod::Array{Float64, 1}, noise_s::Array{Float64, 1})
        ftod = fft(tod)
        ftod ./= noise_s
        tod = real.(ifft(ftod))
        return tod
end


function conjgrad(A, b; x0=0, maxiter=10)
    if x0 == 0
        x = b .* 0
        r = b
    else
        x = copy(x0)
        r = b-A(x)
    end
    n = length(b)
    z = copy(r)
    rz = dot(r,z)
    rz0 = rz
    p = z
    err = Inf
    d = 4
    i = 0
    for i = 1:maxiter
        Ap = A(p)
        alpha = rz / dot(p, Ap)
        x += alpha*p #no scalare
        r -= alpha*Ap #no scalare
        z = copy(r)
        next_rz = dot(r, z)
        err = rz/rz0
        beta = next_rz/rz
        rz = next_rz
        p = z + beta*p
    end

    return x

end


function destriper(N::Int, dataset_baselines)
    b = zeros(N*N)
    for i in dataset_baselines
            tod = strip2.denoise(i.time_order_data, i.noise)
            a, h = strip2.binned_map(N, tod, i.pointing)
            b .+= reshape(a, N*N, 1)[:, 1]
    end

    function A(x)
        res = zeros(N, N)
        for i in dataset_baselines
            tod = strip2.get_tod(reshape(x, N, N), i.pointing)
            tod = strip2.denoise(tod, i.noise)
            mat, hmat = strip2.binned_map(N, tod, i.pointing)
            res .+= mat
        end
        return reshape(res, N*N, 1)[:, 1]
    end

    return strip2.conjgrad(A, b; maxiter=10)

end
