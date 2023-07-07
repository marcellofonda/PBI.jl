# NOTE: export as Z up, Y forward

using ImageView, Gtk.ShortNames
using Plots
gr()

include("integratore.jl")
include("radiografo.jl")
include("retriever.jl")
include("dark_field.jl")


#pixelsize = 1e-6 # 1 micrometro

#Material parameters
delta = 1e-6 #Adimensional
beta = 5e-9 #Adimensional

#Source parameters
k = 5e10 #Rad/m
I_0 = 1 #Still have to figure this one out.

# Setup
z = .2 # m

downsampled = ones(300,300)
downsampled[25:274,25:274] = load("downsampled.jld", "img")


if (false)
sizee = 0
errors = []
minerror = Inf
println("OK!")
for pixelsize in 1e-7:1e-7:2e-6
	I_R_m, I_0_m, I_R_TIE, intermedio =  DarkFieldRetrieve(downsampled, delta, beta, z, 2pi/k, pixelsize * 8)
	imag = ifft( I_R_TIE )[1:div(end,2), 1:div(end,2)] .|> real
	error = sum((imag.*4 .- 3 - downsampled).^2)
	if (error < minerror)
		global sizee, minerror
		sizee = pixelsize
		minerror = error
	end
	push!(errors, error)
end
plot(errors)
println(sizee)

end

I_R_m, I_0_m, I_R_TIE, intermedio =  DarkFieldRetrieve(downsampled, delta, beta, z, 2pi/k, pixelsize * 8)

imag = (ifft( I_R_TIE )[1:div(end,2), 1:div(end,2)] .|> real) .* 4 .- 3
phaseretrieve=PhaseRetrieve(downsampled, delta, beta, k, z, 8*pixelsize)
darkfieldimage = hcat(ifft( I_R_m )[1:div(end,2), 1:div(end,2)],ifft( I_0_m )[1:div(end,2), 1:div(end,2)]) .|> real
intermedio = intermedio[1:div(end,2), 1:div(end,2)]

# imshow(thicknes; name="thicknes");
# imshow(image; name="image");
imshow(downsampled; name="downsampled");
imshow(imag ; name="I_R_TIE_ifft");
imshow(intermedio; name="intermedio");
imshow(phaseretrieve; name="PhaseRetrieve");
# imshow(darkfieldimage; name="darkfieldimage");

println("Fattotutto")
