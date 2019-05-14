# Create the sky and the entire Obs_dataset
include("mapmaking.jl")

using FFTW
using ConjGrad
using LinearAlgebra, SparseArrays
using Test
using Plots

# GCP!!!!!!!!!!!!!!!!!!!!!!!, No λ functions !
# Costruisco FTZ
# Arrivo qui che ho il dataset delle baselines

output = Array{Float64, 2}(undef, 1024, 1024)
