using Images
using ImageMagick

#Given a thickness matrix, compute a Propagation-Based, phase contrast image
function DoRadiography(thickness, delta, beta, k, I_0, z, pixelsize)

	t = Matrix{BigFloat}(thickness)

	#Absorption
	f(x) = exp(-(2 * k * beta) * x)
	assorb = I_0 * f.(t)

	# Define Laplacian kernel
	kernel = [0 1 0; 1 -4 1; 0 1 0]

	# Compute Laplacian using convolution
	lap = imfilter(assorb, kernel) ./ pixelsize^2
	#println(maximum(lap))

	cost = (z * delta) / (2k * beta)
	img = assorb - (cost * lap)

	return img
end
