using Images
using ImageMagick

#Given a matrix of numbers, linearly rescale the values so that the minimum value
#maps to 0 and the maximum to 1. Useful for saving to PNG format.
function ConvertToPng(image; min=minimum(image), max=maximum(image))
	clamper = scaleminmax(min, max)

	return Float64.(clamper.(image))
end

#Given a thickness matrix, compute a Propagation-Based, phase contrast image
function DoRadiography(thickness, delta, beta, k, I_0, z, pixelsize)

	t = Matrix{BigFloat}(thickness)

	#Absorption
	f(x) = exp(-(2 * k * beta) * x)
	assorb = I_0 * f.(t)

	# Define Laplacian kernel
	kernel = [0 1 0; 1 -4 1; 0 1 0]

	# Compute Laplacian using convolution
	lap = imfilter(assorb, kernel)
	#println(maximum(lap))

	cost = (z * delta) / (2k * beta)
	img = assorb - (cost * lap)

	return img
end
