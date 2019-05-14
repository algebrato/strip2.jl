
"""
    function destriper(N::Int, dataset_baselines; n_iter=10)
# Brief description
You have to use this function to clean a noisy-TOD (comphrensive
of white and one-over-f nosie).

Suggestion for raster scan:
In order to suppress the atmospheric noise you should use many
baselines, better if you have observed those maps with different
scanning strategies. Now you can use only Left-Right and Up-Down
scanning.

The descriper solve the noisy-tod linear system:
[Pᵀ N⁻² P] m = [Pᵀ N⁻²] d

We can reduce it at the tipical solution of
A x = b where:

A = [Pᵀ N⁻² P]
b = [Pᵀ N⁻²] d

These matris are too big, so we can you the conjugate gradient method
in order to solve the system.
"""
function destriper(N::Int, dataset_baselines; n_iter=10)
    b = zeros(N*N)
    for i in dataset_baselines
            tod = denoise(i.time_order_data, i.noise)
            a, h = binned_map(N, tod, i.pointing)
            b .+= reshape(a, N*N, 1)[:, 1]
    end

    function A(x)
        res = zeros(N, N)
        for i in dataset_baselines
            tod = get_tod(reshape(x, N, N), i.pointing)
            tod = denoise(tod, i.noise)
            mat, hmat = binned_map(N, tod, i.pointing)
            res .+= mat
        end
        return reshape(res, N*N, 1)[:, 1]
    end

    return conjgrad(A, b; maxiter=n_iter)

end

"""
    function denoise(tod::Array{Float64, 1}, noise_s::Array{Float64, 1})

# Brief Description
This function apply the N⁻² variance matrix to the TOD working in the
fourier space.
"""
function denoise(tod::Array{Float64, 1}, noise_s::Array{Float64, 1})
        ftod = fft(tod)
        ftod ./= noise_s
        tod = real.(ifft(ftod))
        return tod
end

"""
    function conjgrad(A, b; x0=0, maxiter=10)
# Brief description
"""
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
        println("Err: ", err)
    end

    return x

end
