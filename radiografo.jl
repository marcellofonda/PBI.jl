using Images
using ImageMagick
include("dark_field.jl")

function AbsoprtionRadiography(thickness, β, k, I_0)
    f(x) = exp(-(2 * k * β) * x)
	assorb = I_0 * f.(thickness)
end

#Given a thickness matrix, compute a Propagation-Based, phase contrast image
function LaplacianRadiography(thickness, δ, β, k, I_0, z, pixelsize)

	t = Matrix{Float64}(thickness)

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

function FresnelRadiography(t, δ, β, k, I_0, R, pixelsize)
    t = Matrix{Float64}(t)

	#Absorption
	f(x) = exp(-(2 * k * β) * x)
	assorb = I_0 * f.(t)

    return ifft(FresnelPropagate(fft(assorb), 2π/k, δ/β, R, pixelsize)) .|> abs
end


function AngleRadiography(t, δ, β, k, I_0, R, pixelsize)

    # Apply monomorphicity conditions to calculate phi
    Φ = -k * δ * t

    graf1 = plot(Φ[end÷2,1:end])
    display(graf1)

    # Calculate intensity I
    I = I_0 * exp.(-2 * k * β * t)

    # Calculate α_x and α_y using finite differences
    α_x = diff(Φ, dims=1) / (k*pixelsize)
    α_y = diff(Φ, dims=2) / (k*pixelsize)

    # Pad α_x and α_y to have the same size as phi
    α_x = vcat(α_x, zeros(1, size(α_x, 2)))
    α_y = hcat(α_y, zeros(size(α_y, 1), 1))

    
	d = R/pixelsize

    displacement_x = d * α_x
    displacement_y = d * α_y

    grafico = plot(displacement_y[end÷2, 1:end])
    display(grafico)
    # Get the dimensions of the image
    rows, cols = size(I)

    # Create a new image to store the propagated light intensities
    new_I = zeros(rows, cols)

    # Iterate over each pixel in the image
    for i in 1:rows
        for j in 1:cols
            # Calculate the new position of the pixel after light propagation
            new_i = round(Int, i + displacement_x[i, j])
            new_j = round(Int, j + displacement_y[i, j])

            
            # Check if the new position is within the image boundaries
            if new_i > 0 && new_i <= rows && new_j > 0 && new_j <= cols
                # Propagate the light intensity to the new position
                new_I[new_i, new_j] += I[i, j]
            end
        end
    end

    return new_I
end


