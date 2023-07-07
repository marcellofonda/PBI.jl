using FFTW
using PhysicalOptics
include("utils.jl")

"""
    TIERetrieve(Input, gamma, lambda, z, pixelsize)

Backpropagate the image from image plane to object plane using the transport of intensity equation.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function TIERetrieve(Input, gamma, lambda, z, pixelsize)
    # get the size of the input image
    n, m = size(Input)
	cost = z * gamma * lambda / 4pi

    # calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)

    return Input ./ (1 .+ cost .* freq_squared)
end

"""
    TIERetrieve(Input, delta, beta, k, z, pixelsize)

Backpropagate the image from image plane to object plane using the transport of intensity equation.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function TIERetrieve(Input, delta, beta, k, z, pixelsize)
    # get the size of the input image
    n, m = size(Input)
	cost = (z * delta) / (2k * beta)

    # calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)

    return Input ./ (1 .+ cost .* freq_squared)
end


"""
    wave_phase(I, gamma)

Retrieve the phase from the intensity in real space, given the refractive indices' ratio delta/beta=gamma.
"""
function wave_phase(I, gamma)
	return 0.5 * gamma * log(I)
end

"""
    wave_phase(I, delta, beta)

Retrieve the phase from the intensity in real space, given the refractive indices' ratio delta/beta=gamma.
"""
function wave_phase(I, delta, beta)
	return 0.5 * log(I) * delta / beta
end


"""
    FresnelIntegral(wave, lambda, z, pixelsize)

Forward propagate the complex lightwave through free space by fresnel propagator integral.
"""
function FresnelIntegral(wave, lambda, z, pixelsize)
	dimensione = size(wave)[1] * pixelsize * 2
	return propagate(wave, dimensione, z; kernel=fresnel_kernel, λ=lambda, n=1.)
end


# This function takes as input the intensity image (in Fourier space) on the
# object plane and propagates it using the fresnel diffraction integral up
# to the image plane at distance z.
#  The following steps are performed:
# 1. The intensity image is calculated from the Fourier space image I_0 through
# Fourier antitransform (shifted). Only the upper left quadrant is saved, since
# we expect zero padding.
# 2. The complex wave is extracted from the intensity image:
# the amplitude is the square root of the intensity, while the phase is computed
# through the phase(I,gamma) function
# 3. The wave is propagated from object plane to image plane through
# the Fresnel integral in free space.
# 4. The intensity is calculated as the square modulus of the wave
# 6. The fourier transform of the image is returned as output
"""
    FresnelPropagate(I_0, lambda, gamma, z, pixelsize)

Forward propagate the image from the object plane to the image plane using Fresnel diffraction.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function FresnelPropagate(I_0, lambda, gamma, z, pixelsize)

	inverseFourier = real.(ifft(I_0))[ 1:div(size(I_0)[1], 2), 1:div(size(I_0)[2], 2) ]

	E = backgroundPadding(inverseFourier)

    # Extract the complex wave:
	# amplitude
    A = sqrt.(abs.(E))
	# phase
    phi = wave_phase.(abs.(E), gamma)

	# Now actually build the complex wave
    C = A .* exp.(im .* phi)

	# Propagate the wave and get the intensity
	propagated_intensity = abs2.(FresnelIntegral(C, lambda, z, pixelsize)) # Intensity in real space

	#Return intensity in Fourier space
	return real.(fft(propagated_intensity))
end

"""
    FresnelPropagate(I_0, lambda, delta, beta, z, pixelsize)

Forward propagate the image from the object plane to the image plane using Fresnel diffraction.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function FresnelPropagate(I_0, lambda, delta, beta, z, pixelsize)

	inverseFourier = real.(ifft(I_0))[ 1:div(size(I_0)[1], 2), 1:div(size(I_0)[2], 2) ]

	E = backgroundPadding(inverseFourier)

    # Extract the complex wave:
	# amplitude
    A = sqrt.(abs.(E))
	# phase
    phi = wave_phase(abs.(E), delta, beta)

	# Now actually build the complex wave
    C = A .* exp.(im .* phi)

	# Propagate the wave and get the intensity
	propagated_intensity = abs2.(FresnelIntegral(C, lambda, z, pixelsize)) # Intensity in real space

	#Return intensity in Fourier space
	return real.(fft(propagated_intensity)) .* 4 .- 3
end


"Compute the factor for Born scattering backpropagation"
bornFactor(lambda, z, gamma, nu2) =  2( cos(pi* lambda * z * nu2) + gamma * sin( pi * lambda * z * nu2 ))

"Compute the factor for Born scattering backpropagation"
bornFactor(lambda, z, delta, beta, nu2) =  2( cos(pi* lambda * z * nu2) + delta * sin( pi * lambda * z * nu2 )) / beta


# This function does the dark field phase retrieval.
function DarkFieldRetrieve(I_R_measured, delta, beta, z, lambda, pixelsize)
	println("Starting dark field retrieve...")
	image = Float64.(I_R_measured)

	gamma = delta / beta
	k = 2pi/lambda

	image_size = size(image)
	I_0 = image[1,1]

	image = backgroundPadding(image)

	# Save the size of the bigger image
	transform_size = size(image)

	#imshow(image);
	println("Performing FFT...")
	#Do Fourier transform
	I_R = real.(fft(image))  #Intensity in Fourier space


	antitransformed = real.(ifft(I_R))
	#imshow(antitransformed)# .- (antitransformed[2,1] .* correction))


	# FIRST STEP: TIE RETRIEVE

	println("Done. Performing TIE retrieve...")
	I_0_TIE = TIERetrieve(I_R, gamma, lambda, z, pixelsize)
	#I_0_TIE = TIERetrieve(I_R, delta, beta, k, z, pixelsize)


	intermedio = -real.(log.(real.(ifft(I_0_TIE))) ./ (2beta*k))

	#imshow(intermedio[1:Int(end/2), 1:Int(end/2)])
	#plot(intermedio[Int(end/2),1:end])


	#SECOND STEP: PROPAGATE FORWARD

	println("Done. Forward propagating solution...")
	I_R_TIE = FresnelPropagate(I_0_TIE, lambda, gamma, z, pixelsize)

	#TIE_DISPLAY_IMAGE = real.(hcat(ifft(I_0_TIE),ifft(I_R_TIE)))
	#TIE_IMAGE = real.(ifft(I_R_TIE))
	#TIE_IMAGE[1:10, 1:10] = zeros(10,10)
	#imshow(TIE__DISPLAY_IMAGE)


	# THIRD STEP: ITERATIVELY RECONSTRUCT DARK FIELD IMAGE

	println("Done. Entering loop.")
	n,m = transform_size
	delta_I_R_m = zeros(n,m)
	I_0_m = zeros(n,m)
	I_R_m_meno_1 = zeros(n,m)
	I_0_m_meno_1 = zeros(n,m)

	# calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)

	for m in 1:3
		println("m = ",m)
		delta_I_R_m = I_R - I_R_TIE - I_R_m_meno_1
		I_0_m = I_0_m_meno_1 .+ delta_I_R_m ./ bornFactor.(lambda, z, gamma, freq_squared)
		println("Done backpropagating. Doing it forward again...")
		I_R_m_meno_1 = FresnelPropagate(I_0_m, lambda, gamma, z, pixelsize)
		#imshow(hcat((ifft(delta_I_R_m))[1:div(end,2), 1:div(end,2)],(ifft( I_0_m))[1:div(end,2), 1:div(end,2)],(ifft(I_R_m_meno_1))[1:div(end,2), 1:div(end,2)]) .|> real, name="step no. $m")
		println("Solution propagated. Going to next iteration")
		I_0_m_meno_1 = Array(I_0_m)
	end
	println("Done.")
	return delta_I_R_m, I_0_m, I_R_TIE, intermedio
end
