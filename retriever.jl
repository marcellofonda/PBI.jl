using Images
using FFTW

filename = "sfere"


delta = 1e-6
beta = 5e-8
k=5e4
z = 20e4


function PhaseRetrieve(input_image, delta, beta, k, z; pixelsize=1)

	image = Float64.(input_image)

	I_0 = image[1,1]

	image_size = size(image)

	# Double the size of the image for FFT purposes and fill the remaining spaces
	# with zeros
	c = zeros((image_size .* 2)...)
	c .+= I_0
	c[1:(image_size[1]),1:(image_size[2])] = image
	image = c

	# Save the size of the bigger image
	transform_size = size(image)

	#Do Fourier transform
	transform = real.(fft(image ./ I_0))
	#Put origin in the center of the image
	transform = fftshift(transform)

	# Get the spatial frequency corresponding to each pixel in the new image
	freq_grid_x = fftfreq(transform_size[1], 2pi/pixelsize) |> fftshift
	freq_grid_y = fftfreq(transform_size[2], 2pi/pixelsize) |> fftshift

	cost = (z * delta) / (2k * beta)

	# Apply the phase retrieval filter in Fourier space
	for i in 1:(transform_size[1])
		for j in 1:(transform_size[2])
			transform[i,j] /= (1 + cost * (freq_grid_x[i]^2 + freq_grid_y[j]^2))
		end
	end

	# Do the antitransform and take the real part, just to be sure no accidental complex
	# numbers appear
	antitransform = real.(ifft(ifftshift(transform)))[1:image_size[1], 1:image_size[2]]

	# Apply the remaining part of the retrieval algorithm
	antitransform = -1/(2k* beta) * log.(antitransform)

	return antitransform
end
