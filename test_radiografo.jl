# NOTE: export as Z up, Y forward

using ImageView#, Gtk.ShortNames
using Plots;gr()

include("integratore.jl")
include("radiografo.jl")
include("retriever.jl")

pixelsize = 1e-5 # 1 micrometro

#Material parameters
δ = 1e-6 #Adimensional
β = 5e-9 #Adimensional

#Source parameters
k = 5e10 #Rad/m
I_0 = 1 #Still have to figure this one out.

# Setup
z = .2 # m

tetta = thicknessFromFile("tetta")[1]
calcif = thicknessFromFile("calcif_generated")[1] * 4e-6


thicknes = Matrix(tetta)

thicknes[1501: 1500 + size(calcif)[1], 1501: 1500+size(calcif)[2]] += calcif

#LINEA AGGIUNTA PER CAMBIARE ODG
#thicknes .*= 10

println(typeof(thicknes), maximum(tetta), maximum(thicknes))

image = DoRadiography(thicknes, δ, β, k, I_0, z, pixelsize)
println("Done radiography!")

justatest = PhaseRetrieve(image, δ, β, k, z, pixelsize)
another = PhaseRetrieve(image2, δ, β, k, z, pixelsize)
yetanother = PhaseRetrieve(image3, δ, β, k, z, pixelsize)


#imshow(justatest; name="justatest");
#imshow(thicknes; name="thicknes");
#imshow(image; name="image");

confronto = plot(thicknes[end÷2, 1:8:end], label="thicknes")
plot!(justatest[end÷2, 1:8:end], label = "4")
plot!(another[end÷2, 1:8:end], label = "1/4")
plot!(yetanother[end÷2, 1:8:end], label = "1")

display(confronto)