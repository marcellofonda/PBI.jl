using LinearAlgebra
using Images
using Plots;gr()

import Base.+
import Base.*
import Base.-
import Base./

+(x::Tuple{Float64, Float64}, y::Tuple{Float64, Float64}) = x .+ y
-(x::Tuple{Float64, Float64}, y::Tuple{Float64, Float64}) = x .- y
*(x::Tuple{Float64, Float64}, y::Float64) = x .* y
*(x::Float64, y::Tuple{Float64, Float64}) = x .* y
/(x::Tuple{Float64, Float64}, y::Float64) = x ./ y


sizee::Int64 = 300
samples::Int64 = 1000
pixelsize::Float64 = 60e-6
R::Float64 = 6.
# Material parameters
δ::Float64 = 1e-6 #Adimensional
β::Float64 = 5e-9 #Adimensional

k::Float64 = 5e10

N_F = pixelsize^2 * k / (2π * R)
println("Fresnel number N_F is: $N_F")

ra::Float64 = sizee * pixelsize * .8
a::Float64 = sizee * pixelsize * .01
V0::Float64 = .002
woodsaxon(x) = V0 / (1 + exp((norm(x) - ra) / a))


sphere_radius::Float64 = pixelsize * .5
sphere(x::Tuple{Float64,Float64}) = norm(x) > sphere_radius ? 0.0 : 2sqrt(sphere_radius^2 - dot(x,x))


function grad_woodsaxon(x::Tuple{Float64,Float64})
    r = norm(x)
    common_factor = V0 * exp((r - ra) / a) / (a * r * (1 + exp((r - ra) / a))^2)

    grad_x = -x[1] * common_factor
    grad_y = -x[2] * common_factor

    return (grad_x, grad_y)
end

function grad_sphere(x::Tuple{Float64,Float64})
    if norm(x) > sphere_radius
        return (0.0, 0.0)  # Restituisci 0 fuori dal dominio di definizione
    end

    denominator = sqrt(sphere_radius^2 - dot(x,x))
    #grad_x = -2(x[1]) / denominator
    #grad_y = -2(x[2]) / denominator
    return -2. * x / denominator
    #return (grad_x, grad_y)
end


points_in_a_square(n,L) = [(L * (rand() - .5), L * (rand() - .5) ) for _ in 1:n]

scattering_centers::Vector{Tuple{Float64,Float64}} = points_in_a_square(4000, 40 * pixelsize)

function t(x::Tuple{Float64,Float64})
    th::Float64 = woodsaxon(x)
    for center in scattering_centers
        if(norm(x-center) < sphere_radius)
            th += sphere(x-center)
        end
    end
    return th
end

function grad_t(x::Tuple{Float64,Float64})
    gr::Tuple{Float64,Float64} = grad_woodsaxon(x)
    for center in scattering_centers
        if(norm(x-center) < sphere_radius)
            gr += grad_sphere(x-center)
        end
    end
    return gr
end

#image = centered(zeros(Int8, 2sizee + 1,2sizee + 1))

min_angle::Float64 = 0.

#percentage = 0
println("Starting image formation...\n")

#indexes = zeros(Threads.nthreads())
#percentages = zeros(Threads.nthreads())

function simulate_single(range)
	id = Threads.threadid()
	image = centered(zeros(UInt32, 2sizee + 1, 2sizee + 1))
	#any(x-> x<0,image) && println("Oppa")
	for i in range::UnitRange{Int64}
        # indexes[id] += 1
		# perc = div(indexes[id]*100, ((2sizee+1)^2 * samples))
	    # if (perc > percentages[id])
	    #     percentages[id]+=1
	    #     println("$(sum(percentages))%")
	    # end
		selector = rand()
	    x = ((2sizee+1) * pixelsize * (rand()-.5), (2sizee+1) * pixelsize * (rand()-.5))
	    if (exp(-2k * β * t(x)::Float64) >= selector &&
			δ * norm(grad_t(x)::Tuple{Float64,Float64}) >= min_angle )

	        index = Int.(floor.((x + R * δ * grad_t(x)::Tuple{Float64,Float64}) / pixelsize .+ .5))
	        if (abs(index[1]) <= sizee && abs(index[2]) <= sizee)
	            image[index...] += 1 #exp(-2k * β * t(x))
	        end
	    end
		#any(x-> x<0,image) && println("Oppa")
	end
	image
end

function simulate_multi(a)
	chunks = Iterators.partition(a, length(a) ÷ Threads.nthreads())
	tasks = map(chunks) do chunk
    	Threads.@spawn simulate_single(chunk)
    end
    chunk_sums = fetch.(tasks)
    return sum(chunk_sums)
end

time_simulate(a) = @time simulate_multi(a);


image = time_simulate(1:((2sizee+1)^2 * samples))

grafico = heatmap(Array(collect(image)), color=:grays, colorbar=true, aspect_ratio=1)
display(grafico)

save("montecarlo9-int.png", image)
save("montecarlo9.png", Float64.(image)/ maximum(image))

println("Done simulating image!")
# image = Float64.(image) / maximum(image)
#
# image = Array(collect(image))
