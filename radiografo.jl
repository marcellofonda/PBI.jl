using Images
using ImageMagick

filename = "gruppo_1024"

t = load("src_img\\$filename.png")
#Laplaciano = load("sfera_lap.png")

t = Matrix{BigFloat}(t)


delta = 1e-6
beta = 5e-8


k=5e4
I_0 = 1


f(x) = exp(-(2 * k * beta) * x)

#ca = scaleminmax(0,maximum(t))
#save("Provaprova.png", ca.(t))

assorb = I_0 * f.(t)
max = maximum(assorb)
println(max)

#clamp_abs = scaleminmax(0,max)
#save("rendering\\$filename\\assorbimento_$filename.png", Matrix{Float64}(ca.(img)./2))

# Define Laplacian kernel
#kernel = [-1 -1 -1; -1 8 -1; -1 -1 -1]
kernel = [0 1 0; 1 -4 1; 0 1 0]

# Compute Laplacian using convolution
lap = imfilter(assorb, kernel)
println(maximum(lap))

immagini = []
max = 0
min = 1e10
for z in 20e4:0.1:20e4
    global max, min, img
    cost = (z * delta) / (2k * beta)
    println(cost)
    img = assorb - (cost * lap)
    max_loc = maximum(img)
    min_loc = minimum(img)
    if (max_loc > max)
        max = max_loc
    end
    if (min_loc < min)
        min = min_loc
    end
    push!(immagini, img)
end


clamper = scaleminmax(min, max)

immagini = [clamper.(immagini[i]) for i in 1:length(immagini)]

animazione = cat(immagini..., dims=3)

save("rendering\\$filename\\assorbimento_$filename.png", Matrix{Float64}(clamper.(assorb)))

save("rendering\\$filename\\radiografie_$filename.gif", animazione, fps=10)


for i in 1:1#length(immagini)
    save("rendering\\$filename\\sequenze\\$(i)$filename.png", immagini[i])
end
