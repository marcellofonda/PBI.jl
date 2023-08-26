using FFTW
using ImageTransformations, Interpolations
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


"Compute the factor for Born scattering backpropagation"
bornFactor(λ, z, γ, ν2) =  -2( cos(π* λ * z * ν2) + γ * sin( π * λ * z * ν2 ))

"Compute the factor for Born scattering backpropagation"
bornFactor(λ, z, δ, β, ν2) =  -2( cos(π * λ * z * ν2) + δ * sin( π * λ * z * ν2 )/ β) 

tychonoffFactors(λ, R, γ, ν2) = cos(π*λ*R*ν2 - atan(γ))/(2sqrt(1+γ^2) * (1e-2 + cos(π*λ*R*ν2 - atan(γ))^2))


# This function does the dark field phase retrieval.
function DarkFieldRetrieve(I_R_measured, δ, β, z, λ, pixelsize)

	images_paths = []

	println("Starting dark field retrieve...")

	γ = δ / β
	k = 2π/λ

	println("Performing FFT...")

	# Do Fourier transform
	I_R = fft(I_R_measured) #Intensity in Fourier space



	### FIRST STEP: TIE RETRIEVE

	println("Done. Performing TIE retrieve...")
	I_0_TIE = TIERetrieve(I_R, γ, λ, z, pixelsize)
	
	printImage(abs.(ifft(I_0_TIE)), "I_0_TIE", images_paths);

	# Full phase retrieved image for comparison with original thickness
	I_TIE_PR = -(log.(abs.(ifft(I_0_TIE))) ./ (2β*k))

	printImage(I_TIE_PR, "I_TIE_PR", images_paths);

	### SECOND STEP: PROPAGATE FORWARD

	println("Done. Forward propagating solution...")
	I_R_TIE = FresnelPropagate(I_0_TIE, λ, γ, z, pixelsize; scale=20)

	printImage(abs.(ifft(I_R_TIE)), "I_R_TIE", images_paths);

	### THIRD STEP: ITERATIVELY RECONSTRUCT DARK FIELD IMAGE

	println("Done. Entering loop.")
	n,m = size(I_R_measured)

	ΔI_R_m = zeros(n,m)
	I_0_m = zeros(n,m)
	I_R_m_meno_1 = zeros(n,m)
	I_0_m_meno_1 = zeros(n,m)

	# calculate the frequency domain coordinates
	freq_squared = getSquaredFrequenciesGrid(n,m,pixelsize)
	born_factors = tychonoffFactors.(λ, z, γ, freq_squared)#bornFactor.(λ, z, δ, β, freq_squared)

	imagetoshow = []

	for m in 1:3
		println("m = ",m)
		ΔI_R_m = I_R - I_R_TIE - I_R_m_meno_1
		
		printImage(abs.(ifft(ΔI_R_m)),"Delta I_R_$m",images_paths)

		I_0_m = I_0_m_meno_1 .+ ΔI_R_m .* born_factors

		printImage(abs.(ifft(I_0_m)), "I_0_$m", images_paths)

		println("Done backpropagating. Doing it forward again...")
		I_R_m_meno_1 = FresnelPropagate(I_0_m, λ, γ, z, pixelsize; scale=20)
		printImage(abs.(ifft(I_R_m_meno_1)), "I_R_$m", images_paths)
		
		#imagetoshow = Matrix{Float64}(vcat(imagetoshow, hcat(#=abs.(ifft(ΔI_R_m)),abs.(ifft( I_0_m)),=#abs.(ifft(I_R_m_meno_1)) .|> abs)))#, name="step no. $m")
		
		println("Solution propagated. Going to next iteration")
		I_0_m_meno_1 = Array(I_0_m)
	end

	### Show some images

	printImage(born_factors, "BornFactor", images_paths)


	println("Dark Field extraction completed.")
	return I_0_m, images_paths
end
