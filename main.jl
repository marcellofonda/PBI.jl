include("integratore.jl")
include("radiografo.jl")
include("retriever.jl")

# NOTE: Lengths are in units of micrometers

filename = "sfere_num"

#Material parameters
delta = 1e-6 #Adimensional
beta = 5e-8 #Adimensional

#Source parameters
k = 5e4 #Rad/um
I_0 = 1 #Still have to figure this out.

pixel_size = 1. # um

# "obj" if you want to load from a .obj file.
# "image" if you want to load from a .png image.
thickness_source = "obj"

# "compute" if you want to compute the radiography from thickness info.
# "image" if you want to load from a .png image
radiography_source = "compute"
save_radiography = true








# Access thickness image
if thickness_source == "obj"
	# Compute thickness from obj file
	t = GetThickness(filename, pixel_size, 0.001)
	t = ConvertToPng(t)

	save("src_img\\$(filename)_$(size(t)[1])_$(size(t)[2]).png", ConvertToPng(t))
elseif thickness_source == "image"
	# Load thickness from a png image
	println("Loading thickness image: $filename...")
	t = load("src_img\\$filename.png")
	println("Thickness image loaded!")
end

# Access Phase Contrast image
if radiography_source == "compute"
	# Compute image from thickness information
	println("Computing thickness")
	image = DoRadiography(t, delta, beta, k, I_0, z, pixel_size)
	if save_radiography
		save("rendering\\$filename\\pc_$filename.png", ConvertToPng(image))
	end
elseif radiography_source == "image"
	# Load image from .png file
	image = load("rendering\\$filename\\pc_$filename.png")
end

# Do Phase Retrieval
retrieval = PhaseRetrieve(image, delta, beta, k, z; pixelsize=pixel_size)
save("rendering\\$filename\\trasformate\\$(filename)_retrieval.png", ConvertToPng(retrieval))
