using Images

filename = "sfera"

# Read the image from a file
img = load("src_img\\$filename.png")

# Convert to grayscale
img = Gray.(img)

# Define Laplacian kernel
#kernel = [-1 -1 -1; -1 8 -1; -1 -1 -1]
kernel = [0 1 0; 1 -4 1; 0 1 0]
# Compute Laplacian using convolution
lap = imfilter(img, kernel)

max = Float64(maximum(lap))
println(max)
clamper=scaleminmax(0,max)
lap = clamper.(lap)

# Display the original and Laplacian image side by side
save("rendering\\$(filename)_lap.png", hcat(img, lap))
