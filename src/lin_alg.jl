using BenchmarkTools
using Test
using Statistics
using CUDAdrv
using CUDAnative
using CuArrays

N=2^10
x = fill(2.0f0, N)
y = fill(3.0f0, N)


function Kernel_Dot(x, y, z)
    i = blockDim().x * (blockIdx().x-1) + threadIdx().x
    tile = @cuStaticSharedMem(Float64, (blockDim().x, blockDim().y))
    z[i] = x[i] + y[i]
    return
end

a = rand(1024)
b = rand(1024)
c = zeros(1024)

d_a = CuArray(a)
d_b = CuArray(b)
d_c = CuArray(c)

@cuda blocks=(32,32) threads=(32,32) Kernel_Dot(d_a, d_b, d_c)

c = Array(d_c)
