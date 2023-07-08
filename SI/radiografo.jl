using Images
using ImageMagick

#Given a thickness matrix, compute a Propagation-Based, phase contrast image
function DoRadiography(thickness, δ, β, k, I_0, z, pixelsize)

	t = Matrix{BigFloat}(thickness)

	#Absorption
	f(x) = exp(-(2 * k * β) * x)
	assorb = I_0 * f.(t)

	# Define Laplacian kernel
	kernel = [0 1 0; 1 -4 1; 0 1 0]

	# Compute Laplacian using convolution
	lap = imfilter(assorb, kernel) ./ pixelsize^2
	#println(maximum(lap))

	cost = (z * δ) / (2k * β)
	img = assorb - (cost * lap)

	return img
end
