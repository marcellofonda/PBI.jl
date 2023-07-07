#Given a matrix of numbers, linearly rescale the values so that the minimum value
#maps to 0 and the maximum to 1. Useful for saving to PNG format.
function ConvertToPng(image; min=minimum(image), max=maximum(image))
	clamper = scaleminmax(min, max)

	return Float64.(clamper.(image))
end


function generate_calcifs(n)
	sfera = load("..\\src_img\\sfera.png")

	sferetta=sfera[1:32:end,1:32:end]

	calcificazioni = zeros(80, 400)

	posizioni_x = rand(1:72,n)
	posizioni_y = rand(1:392,n)


	for i in 1:n
		for j in 1:8
			for k in 1:8
				calcificazioni[j+posizioni_x[i],k+ posizioni_y[i]] += sferetta[j,k]
			end
		end
	end

	return calcificazioni[1:2:end,1:2:end]
end

"""
    getSquaredFrequenciesGrid(n, m, pixelsize)

Obtain the squared sum of the spacial frequencies corresponding to an image.
"""
function getSquaredFrequenciesGrid(n,m,pixelsize)
	freq_x = fftfreq(n, 2pi/pixelsize)
	freq_y = fftfreq(m, 2pi/pixelsize)
	freq_grid_x = ones(m) * freq_x'
	freq_grid_y = freq_y * ones(n)'

	return freq_squared = freq_grid_x .^ 2 .+ freq_grid_y .^ 2
end

"Double the size of an image in both directions and fill the new parts with the value of the top left pixel."
function backgroundPadding(image)
	image_size = size(image)
	# Double the size of the image for FFT purposes and fill the remaining spaces
	# with zeros
	c = zeros((image_size .* 2)...)
	c .+= image[1,1]
	c[1:(image_size[1]),1:(image_size[2])] = image
	return c
end


function small_primes_product(n::Int)
    # find the largest power of 2 less than or equal to n
    res_p2 = 0
	res_p3 = 0
	res_p5 = 0

	result = Inf

	for p2 in 0:Int(floor(log2(n)+1))
		for p3 in 0:Int(floor(log(n)/log(3) + 1))
			for p5 in 0:Int(floor(log(n)/log(5) + 1))
				temp = 2^p2 * 3^p3 * 5^p5
				if temp >= n && temp < result
					result = temp
					res_p2 = p2
					res_p3 = p3
					res_p5 = p5
				end
			end
		end
	end

    # print the factorization of the result
    println("Closest to $n is $result with factorization: 2^$(res_p2) * 3^$(res_p3) * 5^$(res_p5)")

    return result
end
