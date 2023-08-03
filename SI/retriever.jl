using FFTW
include("utils.jl")

function PhaseRetrieve(input_image, δ, β, k, z, pixelsize)
	println("Executing phase retrieval")
	image = Float64.(input_image)

	I_0 = image[1,1]

	image_size = size(image)

	println("I_0 = $I_0")
	
	#Do Fourier transform
	transform = (fft(complex.(image ./ I_0)))
	#Put origin in the center of the image
	#transform = fftshift(transform)

	# Get the spatial frequency corresponding to each pixel in the new image
	freq_grid_x = fftfreq(image_size[1], 2pi/pixelsize) #|> fftshift
	freq_grid_y = fftfreq(image_size[2], 2pi/pixelsize) #|> fftshift

	cost = (z * δ) / (2k * β)

	# Apply the phase retrieval filter in Fourier space
	for i in 1:(image_size[1])
		for j in 1:(image_size[2])
			transform[i,j] /= (1 + cost * (freq_grid_x[i]^2 + freq_grid_y[j]^2))
		end 
	end

	# Do the antitransform and take the real part, just to be sure no accidental complex
	# numbers appear
	#antitransform = real.(ifft(ifftshift(transform)))[1:image_size[1], 1:image_size[2]]
	antitransform = real.(ifft(complex.(transform)))#[1:image_size[1], 1:image_size[2]]
	x = abs.(ifft(complex.(transform)))
	y = imag.(ifft(complex.(transform)))

	# Apply the remaining part of the retrieval algorithm
	antitransform = -1/(2k * β) * log.(abs.(antitransform))# |> topLeftQuadrant))

	println("Phase retrieval successfully executed.")

	return antitransform#, x, y
end
