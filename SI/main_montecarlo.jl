# NOTE: export as Z up, Y forward
using Pkg;

using ImageView#, Gtk.ShortNames
using Plots;gr()

include("integratore.jl")
include("radiografo.jl")
include("retriever.jl")
include("dark_field.jl")



### Define simulation parameters


# Source parameters
k = 5e10 #Rad/m
I_0 = 1 #Still have to figure this one out.

# Setup parameters
pixelsize = 80e-6
R = 4.5

# Material parameters
δ = 1e-6 #Adimensional
β = 5e-9 #Adimensional

### Compute and print the Fresnel number, to see if approximations are valid

N_F(pixel) = (pixel)^2 * k / (2π * R)
println("Fresnel number N_F is: $(N_F(pixelsize))")

### Name of the folder where to save the images (no LaTeX special characters!!)
image_folder = "montecarlo2"
images_paths = []

sizee = 200

### Thickness map parameters
radius = sizee * pixelsize * .8 # m
border = sizee * pixelsize * .01
max_thickness = .002

calcif_thickness = pixelsize * .1

parameters_info = """
\$\\delta\$=$δ \n\n
\$\\beta\$=$β \n\n
k = $k rad/m\n\n
\$I_0\$=$I_0 \n\n
\$R\$ = $R m\n\n
pixel size = $pixelsize m\n\n

\$N_F\$=$(N_F(pixelsize))\n\n
radius = $radius m\n\n
border = $border m\n\n
calcif thickess = $calcif_thickness m\n\n
breast thickness = $max_thickness m
"""



# Check if the directory exists, and create it if not
if !isdir(image_folder)
    mkpath(image_folder)
end

downsampled = load("montecarlo5.png") .|> Float64

printImage(downsampled, "downsampled", images_paths)


### Do dark-field signal extraction

# Extract images
I_0_m, paths =  DarkFieldRetrieve(downsampled, δ, β, R, 2π/k, pixelsize)

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
