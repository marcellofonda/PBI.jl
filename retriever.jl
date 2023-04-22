using Images
using FFTW

filename = "sfere"


delta = 1e-6
beta = 5e-8
k=5e4
z = 20e4

image = load("rendering\\$filename\\sequenze\\1$filename.png")
image = Float64.(image)
I_0 = image[1,1]

transform = real.(fft(image ./ I_0))

min = minimum(transform)
max = maximum(transform)/100
println(min," ", max)
clamper = scaleminmax(min, max)

save("rendering\\$filename\\trasformate\\$(filename)_transform.png", clamper.(transform))

transform = fftshift(transform)

min = minimum(transform)
max = maximum(transform)/100
println(min," ", max)
clamper = scaleminmax(min, max)

save("rendering\\$filename\\trasformate\\$(filename)_transform_shifted.png", clamper.(transform))

freq_grid = fftfreq(size(image, 1), 1) #|> fftshift

save("rendering\\$filename\\trasformate\\$(filename)_freq.png", clamp01.(freq_grid))

freq_grid = fftfreq(size(image)[1], 1) |> fftshift

save("rendering\\$filename\\trasformate\\$(filename)_freq_shifted.png", clamp01.(freq_grid))

antitransform = real.(ifft(ifftshift(transform)))

min = minimum(antitransform)
max = maximum(antitransform)
println(min," ", max)
clamper = scaleminmax(min, max)

save("rendering\\$filename\\trasformate\\$(filename)_antitransform.png", clamper.(antitransform))



cost = z*delta /(2k* beta)


for i in 1:size(transform)[1]
	for j in 1:size(transform)[2]
		#transform[i,j] /= (1 + cost * ((i-size(transform)[1])^2 + (j-size(transform)[2])^2)/1000000)
		transform[i,j] /= (1 + cost * (freq_grid[i]^2 + freq_grid[j]^2)*100)
	end
end

min = minimum(transform)
max = maximum(transform)/100
println(min," ", max)
clamper = scaleminmax(min, max)

save("rendering\\$filename\\trasformate\\$(filename)_filtered.png", clamper.(transform))

antitransform = real.(ifft(ifftshift(transform)))

min = minimum(antitransform) -.1
max = maximum(antitransform)
println(min," ", max)
clamper = scaleminmax(min, max)

antitransform = -1/(2k* beta) * log.(clamper.(antitransform))

min = minimum(antitransform)
max = maximum(antitransform)
println(min," ", max)
clamper = scaleminmax(min, max)

save("rendering\\$filename\\trasformate\\$(filename)_retrieval.png", clamper.(antitransform))
