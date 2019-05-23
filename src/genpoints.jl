# Some fft routine like fftfreq
include("fft_util.jl")


"""
    function get_points(N_pix::Int, sky::Array{Float64, 2}; direction::Bool)
# Brief Description
This function return the pointing in accord to the sky scanning-strategy.
Now return the pointing directions as pixels. You can choose between two
different scanning strategy mode: if `direction = true` the function return a
left-right scanning. If `direction = false` the return is a up-down scanning.

Future improvements: I ought code a new `get_points` function in order to
get true sky-direction in Right-ascension and Declination coordinates

# Development
The function as is, has not presented weird behaviours.

# Args description
N_pix:      Number of pixel, in linear coordinate
sky  :      The sky that you want observe, in Array{Float64, 2} format
direction:  If true = Left-Right scanning else Up-Down.
"""
function get_points(N_pix::Int, sky::Array{Float64, 2}; direction::Bool)

    points = Array{Float64, 1}(undef, 0)
    pixels = range(1, N_pix*N_pix, length = N_pix*N_pix)

    if direction == true
        x_pix = Array{Int64, 1}(undef, 0)
        r1 = Int.(round.(range(1, N_pix, length=N_pix)))
        r2 = Int.(round.(range(N_pix, 1, length=N_pix)))

        append!(x_pix, r1)
        append!(x_pix, r2)

        x_pix = repeat(x_pix, Int.(round.(N_pix/2)))

        y_pix = Int.(round.(reshape(reshape(repeat(range(
                                   1, N_pix, length = N_pix),
                             N_pix),
                    N_pix, N_pix)',
            1, N_pix*N_pix)))

        append!(points, x_pix[:,1])
        append!(points, y_pix[1,:])
    else
        x_pix = Int.(round.(reshape(reshape(repeat(range(
                                   1, N_pix, length = N_pix),
                             N_pix),
                    N_pix, N_pix)',
            1, N_pix*N_pix)))

        y_pix = Array{Int64, 1}(undef, 0)

        r1 = Int.(round.(range(1, N_pix, length=N_pix)))
        r2 = Int.(round.(range(N_pix, 1, length=N_pix)))
        append!(y_pix, r1)
        append!(y_pix, r2)

        y_pix = repeat(y_pix, Int.(round.(N_pix/2)))

        append!(points, x_pix[1,:])
        append!(points, y_pix[:,1])

    end

    points = reshape(points, Int(round(N_pix^2)), 2)

    return points

end


"""
    function get_tod(sky::Array{Float64, 2}, points::Array{Int64, 2})

# Brief description
Given a set of points, genereted either L-R or U-D scanning, return the
time order data (TOD) series detected by a sensors.

# Development
The function is very simple and has not shown weird behaviours.

# Args description
sky  :      The sky that you want observe, in Array{Float64, 2} format
points:     Time-stream of pointings of your sky scanning-strategy.

"""
function get_tod(sky::Array{Float64, 2}, points::Array{Int64, 2})

    Num_of_pixels = size(points)[1]
    tod = [sky[points[i,1], points[i,2]] for i = 1:Num_of_pixels]

    return tod

end


"""
    function observe_sky(N_pix::Int, sky::Array{Float64, 2};
                         noise::Bool, direction_l_r::Bool)
# Brief description
This function return a noisy tod. The noise is a superimposition of a
one-over-f like: include the instrumental fluctuation and atmospheric
noise, and by a second part of white noise

# Args description
N_pix:      Number of pixel, in linear coordinate
sky  :      The sky that you want observe, in Array{Float64, 2} format
direction:  If true = Left-Right scanning else Up-Down.
noise :     Is a useless args ... I ought clean this mess!

"""
function observe_sky(N_pix::Int, sky::Array{Float64, 2};
                     noise::Bool, direction_l_r::Bool)

    points  = Int.(round.(get_points(N_pix, sky, direction=direction_l_r)))
    tod = get_tod(sky, points)

    rng  = MersenneTwister()
    Frand = fft(randn(rng, length(tod)))
    noise_spec = tod_noise(length(tod), 0.0083, 0.1, 3.0, 44.196)
    fnoise = Frand .* (noise_spec .^0.5)

    tod += real.(ifft(fnoise))
    dataset = TOD_fake_pointing(tod, points, noise_spec)

    return dataset

end


function tod_noise(Nsamp::Int, dt::Float64, fknee::Float64, alpha::Float64,
                   sigma::Float64)
    Freq = abs.(fftfreq(Nsamp, dt))
    return (1 .+ ( max.(Freq, Freq[2]) ./ fknee).^(-alpha)) .* (sigma^2)
end


"""
    function binned_map(N::Int, tod::Array{Float64, 1},
                        pointing::Array{Int64, 2})
# Brief description
The function returns the baselines sum, named Sky Binned Map. If P is the
pointing matrix, returns: PᵀP

# Args description
N_pix:      Number of pixel, in linear coordinate
tod  :      Time Order Data stream, returned by get_tod function
points:     Time-stream of pointings of your sky scanning-strategy.
"""
function binned_map(N_pix::Int, tod::Array{Float64, 1},
                    points::Array{Int64, 2})
    map = zeros(N_pix, N_pix)
    hit = zeros(N_pix, N_pix)

    for i = 1:length(tod)
        map[points[i, 1], points[i, 2]] += tod[i]
        hit[points[i, 1], points[i, 2]] += 1.0
    end

    return map, hit

end
