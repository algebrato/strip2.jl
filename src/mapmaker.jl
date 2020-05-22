
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

m = [Pᵀ N⁻² P]⁻¹ [Pᵀ N⁻²] d
These matris are too big, so we can you the conjugate gradient method
in order to solve the system.
"""

function destriper(N::Int, dataset_baselines; tolerance=1e-6, n_iter=100)
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

    return ConjGrad.cg(A, b, tol=1e-10, maxIter=100)
    #return conjgrad(A, b, tol=1e-10, maxiter=100)

end

"""
    function denoise(tod::Array{Float64, 1}, noise_s::Array{Float64, 1})

# Brief Description
This function apply the N⁻² variance matrix to the TOD working in the
fourier space.
"""
function denoise(tod::Array{Float64, 1}, noise_s::Array{Float64, 1})
        FFTW.set_num_threads(12)
        ftod = fft(tod)
        ftod ./= noise_s
        tod = real.(ifft(ftod))
        return tod
end

"""
    function conjgrad(A, b; x0=0, tol=1e-6, maxiter=100)
# Brief description
"""
genblas_dot(x, y) = dot(x,y)
genblas_scal!(a, x) = x .*= a
genblas_axpy!(a, x, y) = y .+= a.*x
genblas_nrm2(x) = norm(x)
function conjgrad(A, b; tol=1e-6, maxiter=100)

    x = b .* 0
    r = b

    n = length(b)
    z = copy(r)

    p = z
    err = Inf
    i = 0

    residual_0 = genblas_nrm2(r)

    for i = 1:maxiter
        Ap = A(p)
        gamma = genblas_dot(r,z)
        alpha = gamma / genblas_dot(p, Ap)

        if alpha == Inf || alpha < 0
            return -13, iter
        end

        genblas_axpy!(alpha, p, x)
        genblas_axpy!(-alpha, Ap, r)

        residual = genblas_nrm2(r) / residual_0

        if residual <= tol
            println("Resi: ",residual)
            println("Iter: ",i)
            return x
        end

        z = copy(r)

        beta = genblas_dot(z, r) / gamma
        genblas_scal!(beta, p)
        genblas_axpy!(1.0, z, p)

        #println("Err: ", i, " ", residual)
    end

    return x

end


function make_dataset(tod::Array{Float64, 1}, pix_idx::Array{Int64, 1},
                      Num_bas::Int64
                      )
    len_tod = length(tod)
    len2    = Int(round(len_tod / Num_bas))

    mat_tod = reshape(tod, (len2, Num_bas))
    mat_point = reshape(pix_idx, (len2, Num_bas))

    dataset_baselines = Array{Tpn, 1}(undef, Num_bas)

    noise_spec = tod_noise(len2,  0.0083, 0.1, 3.0, 20.76)
    # Ora e` tutto uguale per ogni baselines, l'ideale sarebbe
    # inserire un fit, con il modello scelto, per ogni baselines

    for i in 1:Num_bas
        dat = Tpn(mat_tod[:, i], mat_point[:, i], noise_spec)
        dataset_baselines[i] = dat
    end

    return dataset_baselines

end
