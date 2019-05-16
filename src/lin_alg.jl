using Test

using CUDAdrv, CUDAnative
using CuArrays

dev = CuDevice(0)
ctx = CuContext(dev)


function vadd(a, b, c)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    temp = @cuStaticSharedMem(Float64, (32,1))

    temp[threadIdx().x] += a[i] * b[i]

    sync_threads()

    c[blockIdx().x] = temp[blockIdx().x]

    return
end

dims = (32, 32)
a = rand(1024)
b = rand(1024)
c = zeros(1024)

d_a = CuArray(a)
d_b = CuArray(b)
d_c = CuArray(c)
len = prod(dims)
@cuda blocks=32 threads=32 vadd(d_a, d_b, d_c)
c = Array(d_c)



# _global__ void Vector_Dot_Product ( const float *V1 , const float *V2 , float *V3   )
# {
#  __shared__ float chache[ThreadPerBlock] ;
#
#     float temp ;
#     const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x ;
#     const unsigned int chacheindex = threadIdx.x ;
#     while ( tid < N ){
#       temp += V1[tid] * V2[tid] ;
#       tid += blockDim.x * gridDim.x ;
#       }
#
#
#       chache[chacheindex] = temp ;
#       __synchthreads () ;
#       int i  = blockDim.x / 2 ;
#
#     while ( i!=0 ){
#       if ( chacheindex < i )
#         chache[chacheindex] += chache [chacheindex + i] ;
#         __synchthreads () ;
#         i/=2 ;
#       }
#
#     if ( chacheindex == 0 )
#          V3[blockIdx.x] = chache [0] ;
#
#
# }
