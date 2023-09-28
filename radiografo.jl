using Images
using ImageMagick

###############################
#                             #
#    ABSORPTION SIMULATION    #
#                             #
###############################


function AbsoprtionRadiography(thickness, β, k, I_0)
    f(x) = exp(-(2 * k * β) * x)
	assorb = I_0 * f.(thickness)
end


###################################
#                                 #
#    TIE ANALYTICAL SIMULATION    #
#                                 #
###################################

#Given a thickness matrix, compute a Propagation-Based, phase contrast image
function LaplacianRadiography(thickness, δ, β, k, I_0, R, pixelsize)

	t = Matrix{Float64}(thickness)

	#Absorption
	f(x) = exp(-(2 * k * β) * x)
	assorb = I_0 * f.(t)

	# Define Laplacian kernel
	kernel = [0 1 0; 1 -4 1; 0 1 0]

	# Compute Laplacian using convolution
	lap = imfilter(assorb, kernel) ./ pixelsize^2
	#println(maximum(lap))

	cost = (R * δ) / (2k * β)
	img = assorb - (cost * lap)

	return img
end


#######################################
#                                     #
#    FRESNEL PROPAGATOR SIMULATION    #
#                                     #
#######################################


"""
    wave_phase(I, γ)

Retrieve the phase from the intensity in real space, given the refractive indices' ratio δ/β=γ, with n=1-δ+iβ.
"""
function wave_phase(I, γ)
	return 0.5 * γ * log(I)
end

"""
    wave_phase(I, δ, β)

Retrieve the phase from the intensity in real space, given the refractive index: n=1-δ+iβ.
"""
function wave_phase(I, δ, β)
	return 0.5 * log(I) * δ / β
end


"""
    FresnelIntegral(wave, λ, z, pixelsize)

Forward propagate the complex lightwave through free space by fresnel propagator integral. 
Uses the PhysicalOptics implementation
"""
function FresnelIntegral(wave, λ, z, pixelsize)
	return propagate(wave, size(wave)[1] * pixelsize, z; kernel=fresnel_kernel, λ=λ, n=1.)
end


"""
    FresnelPropagate(I_0, λ, γ, z, pixelsize; scale=1)

Forward propagate the image from the object plane to the image plane using Fresnel diffraction.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function FresnelPropagate(I_0, λ, γ, z, pixelsize; scale=1)

	E = abs.(ifft(I_0)) 

	if(scale != 1)
		E=imresize(E, size(E) .* scale, method=Linear())
		imshow(E);
		p = plot(E[end÷2,:])
		png(p,"plot.png")
	end
    # Extract the complex wave:
	# amplitude
    A = sqrt.(E)
	# phase
    ϕ = wave_phase.(E, γ)

	# Now actually build the complex wave
    C = A .* exp.(im .* ϕ)

	# Propagate the wave and get the intensity
	propagated_intensity = abs2.(FresnelIntegral(C, λ, z, pixelsize / scale)) # Intensity in real space
	imshow(propagated_intensity)
	#Return intensity in Fourier space
	return fft(propagated_intensity[1:scale:end,1:scale:end])
end


"""
    FresnelPropagate(I_0, λ, δ, β, z, pixelsize)

Forward propagate the image from the object plane to the image plane using Fresnel diffraction.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function FresnelPropagate(I_0, λ, δ, β, z, pixelsize)

	inverseFourier = abs.(ifft(I_0))#[ 1:size(I_0)[1]÷2, 1:size(I_0)[2]÷2 ]

	E = inverseFourier #backgroundPadding(inverseFourier)

    # Extract the complex wave:
	# amplitude
    A = sqrt.(abs.(E))
	# phase
    ϕ = wave_phase.(abs.(E), δ, β)

	# Now actually build the complex wave
    C = A .* exp.(im .* ϕ)

	# Propagate the wave and get the intensity
	propagated_intensity = abs2.(FresnelIntegral(C, λ, z, pixelsize)) # Intensity in real space

	#Return intensity in Fourier space
	return (fft(propagated_intensity))
end

function FresnelRadiography(t, δ, β, k, I_0, R, pixelsize)
    #Absorption
	attenuation(x) = exp(-(2 * k * β) * x)
	att = I_0 * attenuation.(t)

    return ifft(FresnelPropagate(fft(att), 2π/k, δ/β, R, pixelsize)) .|> abs
end



###############################
#                             #
#    MONTECARLO SIMULATION    #
#                             #
###############################



"""
    propagateRay(x::Tuple{Float64,Float64}, R::Float64, α_xy::Tuple{Float64,Float64})

Find the 2D coordinates on the image plane of a ray leaving the object plane at position `x` and propagating a distance `R` with angles `α_xy`
"""
propagateRay(x::Tuple{Float64,Float64}, R::Float64, α_xy::Tuple{Float64,Float64}) = x .+ R .* α_xy

function montecarlo_single(size_x::Integer, size_y::Integer, photons_per_pixel, t::Function, grad_t::Function, δ, β, k, R, pixelsize)

	image = centered(zeros(UInt32, size_y, size_x))

	for _ in 1:(size_x * size_y * photons_per_pixel)
		selector = rand()
	    x = (size_y * pixelsize * (rand()-.5), size_x * pixelsize * (rand()-.5))

		α_xy::Tuple{Float64,Float64} = δ .* grad_t(x)
		attenuation::Float64 = exp(-2k * β * t(x))

	    if (attenuation >= selector && norm(α_xy) >= min_angle)
	        index::Tuple{Float64,Float64} = Int.(floor.(propagateRay(x, R, α_xy) ./ pixelsize .+ .5))

	        checkbounds(Bool, index...) && image[index...] += 1
	    end
	end
	image ./ photons_per_pixel
end


function montecarlo_multi(size_x::Integer, size_y::Integer, photons_per_pixel, t::Function, grad_t::Function, δ, β, k, R, pixelsize)
	n = Threads.nthreads()
	sum([Threads.@spawn montecarlo_single(size_x, size_y, photons_per_pixel/n, t, grad_t, δ, β, k, R, pixelsize) for _ in 1:n]) / n
end

function montecarlo_radiography(size_x::Integer, size_y::Integer, photons_per_pixel, t::Function, grad_t::Function, δ, β, k, R, pixelsize; multithreading = true)
	if (multithreading)
		return montecarlo_multi(size_x, size_y, photons_per_pixel/n, t, grad_t, δ, β, k, R, pixelsize)
	else
		return montecarlo_single(size_x, size_y, photons_per_pixel/n, t, grad_t, δ, β, k, R, pixelsize)
	end
end

