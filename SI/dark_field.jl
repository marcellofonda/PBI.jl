using FFTW
using PhysicalOptics
include("utils.jl")

"""
    TIERetrieve(Input, γ, λ, z, pixelsize)

Backpropagate the image from image plane to object plane using the transport of intensity equation.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function TIERetrieve(Input, γ, λ, z, pixelsize)
    # get the size of the input image
    n, m = size(Input)
	cost = z * γ * λ / 4π

    # calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)

    return Input ./ (1 .+ cost .* freq_squared)
end

"""
    TIERetrieve(Input, δ, β, k, z, pixelsize)

Backpropagate the image from image plane to object plane using the transport of intensity equation.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function TIERetrieve(Input, δ, β, k, z, pixelsize)
    # get the size of the input image
    n, m = size(Input)
	cost = (z * δ) / (2k * β)

    # calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)

    return Input ./ (1 .+ cost .* freq_squared)
end


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
    FresnelPropagate(I_0, λ, γ, z, pixelsize)

Forward propagate the image from the object plane to the image plane using Fresnel diffraction.

The input intensity is expected to be in Fourier space, and the output is also an intensity in Fourier space.
"""
function FresnelPropagate(I_0, λ, γ, z, pixelsize)

	E = abs.(ifft(I_0)) 

    # Extract the complex wave:
	# amplitude
    A = sqrt.(abs.(E))
	# phase
    ϕ = wave_phase.(abs.(E), γ)

	# Now actually build the complex wave
    C = A .* exp.(im .* ϕ)

	# Propagate the wave and get the intensity
	propagated_intensity = abs2.(FresnelIntegral(C, λ, z, pixelsize)) # Intensity in real space

	#Return intensity in Fourier space
	return fft(propagated_intensity)
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


"Compute the factor for Born scattering backpropagation"
bornFactor(λ, z, γ, ν2) =  -2( cos(π* λ * z * ν2) + γ * sin( π * λ * z * ν2 ))

"Compute the factor for Born scattering backpropagation"
bornFactor(λ, z, δ, β, ν2) =  -2( cos(π * λ * z * ν2) + δ * sin( π * λ * z * ν2 )/ β) 


# This function does the dark field phase retrieval.
function DarkFieldRetrieve(I_R_measured, δ, β, z, λ, pixelsize)
	println("Starting dark field retrieve...")

	γ = δ / β
	k = 2π/λ

	println("Performing FFT...")

	# Do Fourier transform
	I_R = fft(I_R_measured) #Intensity in Fourier space



	### FIRST STEP: TIE RETRIEVE

	println("Done. Performing TIE retrieve...")
	I_0_TIE = TIERetrieve(I_R, γ, λ, z, pixelsize)
	
	# Full phase retrieved image for comparison with original thickness
	I_TIE_PR = -abs.(log.(abs.(ifft(I_0_TIE))) ./ (2β*k))



	### SECOND STEP: PROPAGATE FORWARD

	println("Done. Forward propagating solution...")
	I_R_TIE = FresnelPropagate(I_0_TIE, λ, γ, z, pixelsize)



	### THIRD STEP: ITERATIVELY RECONSTRUCT DARK FIELD IMAGE

	println("Done. Entering loop.")
	n,m = size(I_R_measured)
	ΔI_R_m = zeros(n,m)
	I_0_m = zeros(n,m)
	I_R_m_meno_1 = zeros(n,m)
	I_0_m_meno_1 = zeros(n,m)

	# calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)
	born_factors = bornFactor.(λ, z, δ, β, freq_squared)

	imagetoshow = []

	for m in 1:3
		println("m = ",m)
		ΔI_R_m = I_R - I_R_TIE - I_R_m_meno_1
		
		I_0_m = I_0_m_meno_1 .+ ΔI_R_m ./ born_factors
		println("Done backpropagating. Doing it forward again...")
		I_R_m_meno_1 = FresnelPropagate(I_0_m, λ, γ, z, pixelsize)

		imagetoshow = Matrix{Float64}(vcat(imagetoshow, hcat(#=abs.(ifft(ΔI_R_m)),abs.(ifft( I_0_m)),=#abs.(ifft(I_R_m_meno_1)) .|> abs)))#, name="step no. $m")
		
		println("Solution propagated. Going to next iteration")
		I_0_m_meno_1 = Array(I_0_m)
	end

	### Show some images

	imshow(born_factors; name="BornFactor");
	imshow(imagetoshow);



	println("Dark Field extraction completed.")
	return ΔI_R_m, I_0_m, I_R_TIE, I_TIE_PR
end
