# NOTE: export as Z up, Y forward

using ImageView, Gtk.ShortNames


include("integratore.jl")
include("radiografo.jl")
include("retriever.jl")
include("dark_field.jl")


pixelsize = 1e-6 # 1 micrometro

#Material parameters
delta = 1e-6 #Adimensional
beta = 5e-9 #Adimensional

#Source parameters
k = 5e10 #Rad/m
I_0 = 1 #Still have to figure this one out.

# Setup
z = .2 # m

#gui = imshow_gui((300, 300), (2, 2))  # 2 columns, 1 row of images (each initially 300×300)
#canvases = gui["canvas"]





# calcif = t .|> BigFloat
#tetta = GetThickness("tetta", pixelsize, 0.001, false) .|> BigFloat
#calcif, pixelsize = GetThickness("calcif2", pixelsize, 0.001, true) .|> BigFloat

tetta = thicknessFromFile("tetta")[1]
calcif = thicknessFromFile("calcif_generated")[1] * 4e-6


thicknes = Matrix(tetta)

thicknes[1501: 1500 + size(calcif)[1], 1501: 1500+size(calcif)[2]] += calcif


println(typeof(thicknes), maximum(tetta), maximum(thicknes))

image = DoRadiography(thicknes, delta, beta, k, I_0, z, pixelsize)

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

I_R_m, I_0_m, I_R_TIE, intermedio =  DarkFieldRetrieve(downsampled, delta, beta, z, 2pi/k, pixelsize * 8)

imag = ifft( I_R_TIE )[1:div(end,2), 1:div(end,2)] .|> real
phaseretrieve=PhaseRetrieve(downsampled, delta, beta, k, z, 8*pixelsize)
darkfieldimage = hcat(ifft( I_R_m )[1:div(end,2), 1:div(end,2)],ifft( I_0_m )[1:div(end,2), 1:div(end,2)]) .|> real
intermedio = intermedio[1:div(end,2), 1:div(end,2)]

imshow(thicknes; name="thicknes");
imshow(image; name="image");
imshow(downsampled; name="downsampled");
imshow(imag; name="I_R_TIE_ifft");
imshow(intermedio; name="intermedio");
imshow(phaseretrieve; name="PhaseRetrieve");
imshow(darkfieldimage; name="darkfieldimage");

println("Fattotutto")
