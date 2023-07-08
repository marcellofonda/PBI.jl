# NOTE: export as Z up, Y forward

using ImageView#, Gtk.ShortNames
using Plots;gr()

include("integratore.jl")
include("radiografo.jl")
include("retriever.jl")
include("dark_field.jl")

pixelsize = 1e-5 # 1 micrometro

#Material parameters
δ = 1e-6 #Adimensional
β = 5e-9 #Adimensional

#Source parameters
k = 5e10 #Rad/m
I_0 = 1 #Still have to figure this one out.

# Setup
z = .2 # m

# calcif = t .|> BigFloat
#tetta = GetThickness("tetta", pixelsize, 0.001, false) .|> BigFloat
#calcif, pixelsize = GetThickness("calcif2", pixelsize, 0.001, true) .|> BigFloat

tetta = thicknessFromFile("tetta")[1]
calcif = thicknessFromFile("calcif_generated")[1] * 4e-6

thicknes = Matrix(tetta)

thicknes[1501: 1500 + size(calcif)[1], 1501: 1500+size(calcif)[2]] += calcif

#LINEA AGGIUNTA PER CAMBIARE ODG
#thicknes .*= 10

println(typeof(thicknes), maximum(tetta), maximum(thicknes))

image = DoRadiography(thicknes, δ, β, k, I_0, z, pixelsize)

println("Done radiography!")

#Gtk.showall(gui["window"])

# Define the kernel for the box filter
kernel = ones(8, 8) / 64

# Convolve the image with the kernel using the imfilter function
convolved = imfilter(image, kernel)

# Downsample the convolved image by a factor of 8 in both dimensions
downsampled = convolved[1:8:end, 1:8:end]
#downsampled += rand(size(downsampled)...) / 200

println("OK!")
justatest = PhaseRetrieve(image, δ, β, k, z, pixelsize)

I_R_m, I_0_m, I_R_TIE, intermedio =  DarkFieldRetrieve(downsampled, δ, β, z, 2π/k, pixelsize * 8)
intermedio = intermedio[1:end÷2, 1:end÷2]

imag = ifft( I_R_TIE )[1:end÷2, 1:end÷2] .|> real

phaseretrieve=PhaseRetrieve(downsampled, δ, β, k, z, 8*pixelsize)

darkfieldimage = hcat(ifft( I_R_m )[1:end÷2, 1:end÷2],ifft( I_0_m )[1:end÷2, 1:end÷2]) .|> real

imshow(justatest; name="justatest"); 
imshow(thicknes; name="thicknes");
imshow(image; name="image");
#imshow(downsampled; name="downsampled");
#imshow(imag; name="I_R_TIE_ifft");
#imshow(intermedio; name="intermedio");
imshow(phaseretrieve; name="PhaseRetrieve");
#imshow(darkfieldimage; name="darkfieldimage");

confronto = plot(thicknes[end÷2, 1:8:end], label="thicknes")
plot!(justatest[end÷2, 1:8:end], label = "upsampled PR")
plot!(phaseretrieve[end÷2, 1:end], label="downsampled PR")
display(confronto)

println("Fattotutto")
