image_folder = ""

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
	freq_x = fftfreq(n, 2π/pixelsize)
	freq_y = fftfreq(m, 2π/pixelsize)
	freq_grid_x = ones(m) * freq_x'
	freq_grid_y = freq_y * ones(n)'

	return freq_grid_x .^ 2 .+ freq_grid_y .^ 2
end


"Double the size of an image in both directions and pad with zeros."
function zeroPadding(image)
	image_size = size(image)
	# Double the size of the image for FFT purposes and fill the remaining spaces
	# with zeros
	c = zeros((image_size .* 2)...)
	#c .+= image[1,1]
	c[image_size[1]÷2+1:3(image_size[1])÷2,image_size[2]÷2+1:3(image_size[2])÷2] = image
	return c
end

function removePadding(image)
	x,y = size(image)
	return image[x÷4:3x÷4,y÷4:3y÷4]
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

function printImage(image, name, images_paths)
	#image[1,1] = 0.
	#image[end,end]=1.2
	profile_plot = plot(image[:,end÷2], legend=false)
	image_plot = heatmap(image, c=:grays, ticks=nothing, border=:none, legend=false, aspect_ratio=1.)

	# Combine the image and profile
	combined_plot = plot(image_plot, profile_plot, layout=(1,2), size=(1920,1080))
	#plot!(50:67, image[50:67, end÷2], 
	#inset = (2, bbox(0.3, 0.5, 0.3, 0.3, :bottom, :right)),
    #ticks = nothing,
    #subplot = 3,
	#bg_inside=nothing,
	#legend=false
	#)
	png(combined_plot, "$image_folder\\$name.png")

	push!(images_paths, name)

	"Image printed at $image_folder\\$name.png"
end



function create_presentation(image_filenames, title::String, parameters_info)
    # Open a file for writing
    open("$title\\$title.tex", "w") do io
        # Write the preamble
        write(io, "\\documentclass{beamer}\n")
        write(io, "\\usepackage{graphicx}\n")
        write(io, "\\begin{document}\n")

        # Write the title frame
        write(io, "\\begin{frame}\n")
        write(io, "\\frametitle{\"$(title)\"}\n")
		write(io, parameters_info)
        write(io, "\\end{frame}\n")

        # Write a frame for each image filename
        for filename in image_filenames
            write(io, "\\begin{frame}\n")
            write(io, "\\frametitle{$(replace(filename, "_" => " "))}\n")
			write(io, "\\begin{figure}\n")
            write(io, "\\includegraphics[width=\\textwidth]{\"$(filename).png\"}\n")
			write(io, "\\end{figure}\n")
            write(io, "\\end{frame}\n")
        end

        # Write the end of the document
        write(io, "\\end{document}\n")
    end
	cd(title)
    run(`cmd /C pdflatex $title.tex`)
	cd("..")
end

"utils.jl successfully imported!"
