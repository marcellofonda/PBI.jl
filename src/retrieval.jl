using FFTW
using ImageTransformations, Interpolations
using PhysicalOptics
include("utils.jl")
include("simulation.jl")


#############################
#                           #
#    TIE PHASE RETRIEVAL    #
#                           #
#############################


"""
    TIERetrieve(Input, γ, λ, R, pixelsize)

Backpropagate the image from image plane to object plane using the transport of intensity equation.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function TIERetrieve(Input, γ, λ, R, pixelsize)
    # get the size of the input image
    n, m = size(Input)
	cost = R * γ * λ / 4π

    # calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)

    return Input ./ (1 .+ cost .* freq_squared)
end

"""
    TIERetrieve(Input, δ, β, k, R, pixelsize)

Backpropagate the image from image plane to object plane using the transport of intensity equation.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function TIERetrieve(Input, δ, β, k, R, pixelsize)
    # get the size of the input image
    n, m = size(Input)
	cost = (R * δ) / (2k * β)

    # calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)

    return Input ./ (1 .+ cost .* freq_squared)
end



"Execute Phase Retrieval inverting the Transport of Intensity Equation on a PB image"
function PhaseRetrieve(input_image, δ, β, k, R, pixelsize)
	image = fft(input_image)

	image = ifft(TIERetrieve(image, δ, β, k, R, pixelsize))

	image = -1/(2k * β) * log.(abs.(image))
end



##############################
#                            #
#    DARK FIELD RETRIEVAL    #
#                            #
##############################

"Compute the factor for Born scattering backpropagation"
bornFactor(λ, R, γ, ν2) =  -2( cos(π* λ * R * ν2) + γ * sin( π * λ * R * ν2 ))

"Compute the factor for Born scattering backpropagation"
bornFactor(λ, R, δ, β, ν2) =  -2( cos(π * λ * R * ν2) + δ * sin( π * λ * R * ν2 )/ β) 

tychonoffFactors(λ, R, γ, ν2, ε) = cos(π*λ*R*ν2 - atan(γ))/(2sqrt(1+γ^2) * (ε + cos(π*λ*R*ν2 - atan(γ))^2))


# This function does the dark field phase retrieval.
function DarkFieldRetrieve(input_image, δ, β, R, λ, pixelsize, ε; scale = 1, save_images = false)

	images_paths = []

	println("Starting dark field retrieve...")

	γ = δ / β
	k = 2π/λ

	println("Performing FFT...")

	# Do Fourier transform
	I_R = fft(input_image) #Intensity in Fourier space



	### FIRST STEP: TIE RETRIEVE

	println("Done. Performing TIE retrieve...")
	I_0_TIE = TIERetrieve(I_R, γ, λ, R, pixelsize)
	
	save_images && printImage(abs.(ifft(I_0_TIE)), "I_0_TIE", images_paths);

	# Full phase retrieved image for comparison with original thickness
	I_TIE_PR = -(log.(abs.(ifft(I_0_TIE))) ./ (2β*k))

	save_images && printImage(I_TIE_PR, "I_TIE_PR", images_paths);

	### SECOND STEP: PROPAGATE FORWARD

	println("Done. Forward propagating solution...")
	I_R_TIE = FresnelPropagate(I_0_TIE, λ, γ, R, pixelsize; scale=scale)

	save_images && printImage(abs.(ifft(I_R_TIE)), "I_R_TIE", images_paths);

	### THIRD STEP: ITERATIVELY RECONSTRUCT DARK FIELD IMAGE

	println("Done. Entering loop.")
	n,m = size(input_image)

	ΔI_R_m = zeros(n,m)
	I_0_m = zeros(n,m)
	I_R_m_meno_1 = zeros(n,m)
	I_0_m_meno_1 = zeros(n,m)

	# calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)
	born_factors = tychonoffFactors.(λ, R, γ, freq_squared, ε)#bornFactor.(λ, R, δ, β, freq_squared)

	for m in 1:1
		println("m = ",m)
		ΔI_R_m = I_R - I_R_TIE - I_R_m_meno_1
		
		save_images && printImage(abs.(ifft(ΔI_R_m)),"Delta I_R_$m",images_paths)

		I_0_m = I_0_m_meno_1 .+ ΔI_R_m .* born_factors

		save_images && printImage(abs.(ifft(I_0_m)), "I_0_$m", images_paths)

		println("Done backpropagating. Doing it forward again...")
		I_R_m_meno_1 = FresnelPropagate(I_0_m, λ, γ, R, pixelsize; scale=scale)
		save_images && printImage(abs.(ifft(I_R_m_meno_1)), "I_R_$m", images_paths)
		
		println("Solution propagated. Going to next iteration")
		I_0_m_meno_1 = Array(I_0_m)
	end

	### Show some images

	save_images && printImage(born_factors, "BornFactor", images_paths)

	TIE_result = abs.ifft(I_0_TIE)
	A_Born = abs.ifft(I_0_m) ./ TIE_result
	
	println("Dark Field extraction completed.")
	save_images && return A_Born, I_TIE_PR, images_paths
	
	return A_Born, I_TIE_PR
end
