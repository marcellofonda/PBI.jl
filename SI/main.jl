# NOTE: export as Z up, Y forward
using Pkg;

using ImageView#, Gtk.ShortNames
using Plots;gr()

include("integratore.jl")
include("radiografo.jl")
include("retriever.jl")
include("dark_field.jl")



### Define simulation parameters

# Material parameters
δ = 1e-6 #Adimensional
β = 5e-9 #Adimensional

# Source parameters
k = 5e10 #Rad/m
I_0 = 1 #Still have to figure this one out.

# Setup parameters
R = .1 # m
pixelsize = 1e-5 # m
padding = 400 #px of border



### Compute and print the Fresnel number, to see if approximations are valid

N_F = pixelsize^2 * k / (2π * R)
println("Fresnel number N_F is: $N_F")



### Build the thickness image. From here...

# Calculate the thickness from .obj 3D meshes
#tetta = GetThickness("tetta", pixelsize, 0.001, false)
#calcif, pixelsize = GetThickness("calcif2", pixelsize, 0.001, true)

# Initialize an empty thickness matrix
thickness = zeros(2000+2padding+1,2000+2padding+1)

# Import
tetta = thicknessFromFile("tetta")[1]
calcif = thicknessFromFile("calcif_generated")[1] * 4e-6

# Sum all the contributions to thickness of sample
tetta[1501: 1500 + size(calcif)[1], 1501: 1500+size(calcif)[2]] += calcif
thickness[padding+1 : padding+2000, padding+1 : padding+2000] = tetta

println(
    "Done setting up the thickess. The tit's max thickness is: ",
    maximum(tetta),
    " m, while including calcifications gives a max of: ",
    maximum(thickness),
    " m."
    )

### ... to here, the code can be replaced with a simple image import, e.g. thickness = load("my_img.png")



### Do radiography with 2 different methods: ray-tracing and Fresnel propagator

image_laplacian = LaplacianRadiography(thickness, δ, β, k, I_0, R, pixelsize)
image_fresnel = FresnelRadiography(thickness, δ, β, k, I_0, R, pixelsize)

println("Done radiographies!")



### Convolve the image with a filter to simulate pixel effects

# Define the kernel for the PSF
kernel = ones(8, 8) / 64 #use a box PSF

# Convolve the image with the kernel using the imfilter function and 
# downsample the convolved image by a factor of 8 in both dimensions
downsampled = imfilter(image_laplacian, kernel)[1:8:end, 1:8:end]

println("Convolution and downsampling executed!")



### Try TIE-PR on original images

Laplacian_PR = PhaseRetrieve(image_laplacian, δ, β, k, R, pixelsize)
Fresnel_PR =PhaseRetrieve(image_fresnel, δ, β, k, R, pixelsize)



### Do dark-field signal extraction

# Extract images
I_R_m, I_0_m, I_R_TIE, I_TIE_PR =  DarkFieldRetrieve(downsampled, δ, β, R, 2π/k, pixelsize * 8)
# Obtain the actual dark-field info
darkfieldimage_real = ifft( I_0_m ) .|> real
darkfieldimage_abs = ifft( I_0_m ) .|> abs



### Show images

imshow(Laplacian_PR; name="Laplacian_PR"); 
imshow(thickness; name="thickness");
imshow(Fresnel_PR; name="PhaseRetrieve");
#imshow(x; name="modulo");
#imshow(y; name="parte immaginaria");
# imshow(image_laplacian; name="image_laplacian");
#imshow(downsampled; name="downsampled");
#imshow(imag; name="I_R_TIE_ifft");
#imshow(I_TIE_PR; name="I_TIE_PR");
imshow(darkfieldimage_real; name="darkfieldimage_real");
imshow(darkfieldimage_abs; name="darkfieldimage_abs");
#imshow(real.(ifft(fft(Float64.(difference)))); name="test");



### Plot the retrieved thickness

confronto = plot(thickness[end÷2, :], label="thickness")
plot!(Laplacian_PR[end÷2,:], label = "laplacian PR")
plot!(Fresnel_PR[end÷2,:], label="fresnel PR")
display(confronto)

println("Fattotutto")
