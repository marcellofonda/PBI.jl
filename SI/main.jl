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
δ = 2.3e-6 #Adimensional
β = 4.25e-9 #Adimensional

# Source parameters
k = 5e10 #Rad/m
I_0 = 1 #Still have to figure this one out.

# Setup parameters
R = 2.5 # m
pixelsize = 1e-6 # m
padding = 400 #px of border

downscale_factor = 20

### Compute and print the Fresnel number, to see if approximations are valid

N_F(pixel) = (pixel)^2 * k / (2π * R)
println("Fresnel number N_F is: $(N_F(pixelsize)) vs $(N_F(downscale_factor * pixelsize))")

### Name of the folder where to save the images (no LaTeX special characters!!)
image_folder = "montecarlo"
images_paths = []


### Thickness map parameters
radius = .002 # m
border = .00006
max_thickness = .001

calcif_thickness = 4e-6

parameters_info = """
\$\\delta\$=$δ \n\n
\$\\beta\$=$β \n\n
k = $k rad/m\n\n
\$I_0\$=$I_0 \n\n
\$R\$ = $R m\n\n
pixel size = $pixelsize m\n\n

downscale factor = $downscale_factor \n\n
\$N_F\$=$(N_F(pixelsize)) (before downscaling)\n\n
\$N_F\$=$(N_F(downscale_factor * pixelsize)) (after downscaling)\n\n
radius = $radius m\n\n
border = $border m\n\n
padding = $padding px\n\n
calcif thickess = $calcif_thickness m\n\n
breast thickness = $max_thickness m
"""



# Check if the directory exists, and create it if not
if !isdir(image_folder)
    mkpath(image_folder)
end


pixelradius = Int(ceil(radius / pixelsize))

xoffset = pixelradius-20
yoffset = pixelradius-50

pixelborder = border / pixelsize

tetta = [radius / (1+ exp((sqrt(i^2 + j^2) - pixelradius+5pixelborder)/pixelborder)) for i in -pixelradius:pixelradius, j in -pixelradius:pixelradius]
tetta = imfilter(tetta, Kernel.gaussian(20))

printImage(tetta, "tetta", images_paths)


image_size = small_primes_product(2pixelradius+2padding)
thickness = zeros(image_size+1,image_size+1)

calcif = thicknessFromFile("calcif_generated")[1] * calcif_thickness

tetta[xoffset+1-size(calcif)[1]: xoffset + size(calcif)[1], yoffset+1-size(calcif)[2]: yoffset+size(calcif)[2]] += [calcif calcif; calcif calcif]

tetta .-= tetta[1,end÷2]

thickness[padding:padding+2pixelradius, padding:padding+2pixelradius] = tetta

printImage(thickness, "thickness", images_paths)



# ### Build the thickness image. From here...

# # Calculate the thickness from .obj 3D meshes
# #tetta = GetThickness("tetta", pixelsize, 0.001, false)
# #calcif, pixelsize = GetThickness("calcif2", pixelsize, 0.001, true)

# # Initialize an empty thickness matrix
# thickness = zeros(2000+2padding+1,2000+2padding+1)

# # Import
# tetta = thicknessFromFile("tetta")[1]
# calcif = thicknessFromFile("calcif_generated")[1] * 4e-6




# # Define breast with Wood-Saxon profile
# tetta = [0.002 / (1+ exp((sqrt(i^2 + j^2) - 800)/50)) for i in -1000:1000, j in -1000:1000]
# tetta = imfilter(tetta, Kernel.gaussian(20))




# # Sum all the contributions to thickness of sample
# tetta[xoffset+1: xoffset + size(calcif)[1], yoffset+1: yoffset+size(calcif)[2]] += calcif
# thickness[padding+1 : padding+2001, padding+1 : padding+2001] = tetta



# println(
#     "Done setting up the thickess. The tit's max thickness is: ",
#     maximum(tetta),
#     " m, while including calcifications gives a max of: ",
#     maximum(thickness),
#     " m."
#     )

# ### ... to here, the code can be replaced with a simple image import, e.g. thickness = load("my_img.png")






### Do radiography with 2 different methods: ray-tracing and Fresnel propagator

image_laplacian = LaplacianRadiography(thickness, δ, β, k, I_0, R, pixelsize)
image_fresnel = FresnelRadiography(thickness, δ, β, k, I_0, R, pixelsize)

printImage(image_laplacian, "image_laplacian", images_paths)
printImage(image_fresnel, "image_fresnel", images_paths)
#Free some memory
thickness = nothing

println("Done radiographies!")



### Convolve the image with a filter to simulate pixel effects

# Define the kernel for the PSF
kernel = Kernel.gaussian(downscale_factor) #ones(8, 8) / 64 #use a box PSF

# Convolve the image with the kernel using the imfilter function and 
# downsample the convolved image by a factor of 8 in both dimensions
downsampled = imfilter(image_fresnel, kernel)[1:downscale_factor:end, 1:downscale_factor:end]

println("Convolution and downsampling executed!")




### Try TIE-PR on original images

Laplacian_PR = PhaseRetrieve(image_laplacian, δ, β, k, R, pixelsize)
Fresnel_PR =PhaseRetrieve(image_fresnel, δ, β, k, R, pixelsize)

printImage(Laplacian_PR, "Laplacian_PR", images_paths)
printImage(Fresnel_PR, "Fresnel_PR", images_paths)
printImage(downsampled, "downsampled", images_paths)


### Do dark-field signal extraction

# Extract images
I_0_m, paths =  DarkFieldRetrieve(downsampled, δ, β, R, 2π/k, pixelsize * downscale_factor)

push!(images_paths, paths...)

# Obtain the actual dark-field info
# darkfieldimage_real = ifft( I_0_m ) .|> real
# darkfieldimage_abs = ifft( I_0_m ) .|> abs


### Show images

# #imshow(Laplacian_PR; name="Laplacian_PR"); 
# #imshow(thickness; name="thickness");
# #imshow(Fresnel_PR; name="PhaseRetrieve");
# #imshow(x; name="modulo");
# #imshow(y; name="parte immaginaria");
# #imshow(image_laplacian; name="image_laplacian");
# #imshow(image_fresnel; name="image_fresnel");
# imshow(downsampled; name="downsampled");
# #imshow(imag; name="I_R_TIE_ifft");
# imshow(I_TIE_PR; name="I_TIE_PR");
# imshow(darkfieldimage_real; name="darkfieldimage_real");
# imshow(darkfieldimage_abs; name="darkfieldimage_abs");
# #imshow(real.(ifft(fft(Float64.(difference)))); name="test");
# imshow(df; name="df")


# ### Plot the retrieved thickness
# confronto = plot(thickness[end÷2, :], label="thickness")
# plot!(Laplacian_PR[end÷2,:], label = "laplacian PR")
# plot!(Fresnel_PR[end÷2,:], label="fresnel PR")
# display(confronto)

# confronto_img = plot(image_laplacian[end÷2, :], label="laplacian")
# plot!(image_fresnel[end÷2, :], label="fresnel")
# display(confronto_img)




create_presentation(images_paths, image_folder, parameters_info)

println("Fattotutto")
