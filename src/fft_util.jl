function fftfreq(samples::Int, fs::Float64)
    freq = Array{Float64, 1}(undef, 0)
    append!(freq, 0)

    f_min = 1.0 / (fs*samples)
    f_max = (1.0 / fs) * 0.5

    append!(freq, range(f_min, f_max-f_min, length=Int(floor(samples/2))-1))
    append!(freq, range(-f_max, -f_min, length = Int(floor(samples/2))))

    return freq

end
