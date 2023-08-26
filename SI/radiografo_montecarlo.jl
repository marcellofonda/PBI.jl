using LinearAlgebra
using Images

import Base.+
import Base.*
import Base.-
import Base./

+(x::Tuple{Float64, Float64}, y::Tuple{Float64, Float64}) = x .+ y
-(x::Tuple{Float64, Float64}, y::Tuple{Float64, Float64}) = x .- y
*(x::Tuple{Float64, Float64}, y::Float64) = x .* y
*(x::Float64, y::Tuple{Float64, Float64}) = x .* y
/(x::Tuple{Float64, Float64}, y::Float64) = x ./ y


size = 200
samples = 100
pixelsize = 80e-6
R = 4.5
# Material parameters
δ = 1e-6 #Adimensional
β = 5e-9 #Adimensional

k = 5e10

N_F = pixelsize^2 * k / (2π * R)
println("Fresnel number N_F is: $N_F")

R = size * pixelsize * .8
a = size * pixelsize * .01
V0 = .002
woodsaxon(x) = V0 / (1 + exp((norm(x) - R) / a))


sphere_radius = pixelsize * .1
sphere(x) = norm(x) > sphere_radius ? 0.0 : sqrt(sphere_radius^2 - dot(x,x))


function grad_woodsaxon(x)
    r = norm(x)
    common_factor = V0 * exp((r - R) / a) / (a * r * (1 + exp((r - R) / a))^2)

    grad_x = -x[1] * common_factor
    grad_y = -x[2] * common_factor

    return (grad_x, grad_y)
end

function grad_sphere(x)
    if norm(x) > sphere_radius
        return (0.0, 0.0)  # Restituisci 0 fuori dal dominio di definizione
    end

    denominator = sqrt(sphere_radius^2 - dot(x,x))
    grad_x = -(x[1]) / denominator
    grad_y = -(x[2]) / denominator

    return (grad_x, grad_y)
end


points_in_a_square(n,L) = [(L * (rand() - .5), L * (rand() - .5) ) for _ in 1:n]

scattering_centers = points_in_a_square(400, 20 * pixelsize)

function t(x)
    th = woodsaxon(x)
    for center in scattering_centers
        if(norm(x-center) < sphere_radius)
            th += sphere(x-center)
        end
    end
    return th
end

function grad_t(x)
    gr = grad_woodsaxon(x)
    for center in scattering_centers
        if(norm(x-center) < sphere_radius)
            gr += grad_sphere(x-center)
        end
    end
    return gr
end

#image = centered(zeros(Int8, 2size + 1,2size + 1))

min_angle = 0

#percentage = 0
println("Starting image formation...\n")

#indexes = zeros(Threads.nthreads())
#percentages = zeros(Threads.nthreads())

function simulate_single(range)
	id = Threads.threadid()
	image = centered(zeros(UInt8, 2size + 1, 2size + 1))
	#any(x-> x<0,image) && println("Oppa")
	for i in range
		# indexes[id] += 1
		# perc = div(indexes[id]*100, ((2size+1)^2 * samples))
	    # if (perc > percentages[id])
	    #     percentages[id]+=1
	    #     println("$(sum(percentages))%")
	    # end
		selector = rand()
	    x = ((2size+1) * pixelsize * (rand()-.5), (2size+1) * pixelsize * (rand()-.5))
	    if (exp(-2k * β * t(x)) >= selector &&
			δ * norm(grad_t(x)) >= min_angle )

	        index = Int.(floor.((x + R * δ * grad_t(x)) / pixelsize .+ .5))
	        if (abs(index[1]) <= size && abs(index[2]) <= size)
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
	any(any.(x-> x<0, chunk_sums)) && println("Oppa")
    return sum(chunk_sums)
end

image = simulate_multi(1:((2size+1)^2 * samples))
heatmap(image)
println("Done simulating image!")
# image = Float64.(image) / maximum(image)
#
# image = Array(collect(image))
