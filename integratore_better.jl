# This program should read a .obj file and project it on the yz plane.
# As a result, a matrix representing the projection is obtained,
# of desired dimensions, containing the x coordinates of the surface
# of the mesh. This should always come as an array of even dimension,
# since the mesh is expected to be closed (unless you are very unlucky
# and you bullseye a triangle containing a segment parallel to the
# x axis. This case is not currently dealt with, but should be. Have hope).

# The mesh contained in the .obj should be
# - closed
# - triangulated

using LinearAlgebra
using Images
using FileIO
using CUDA
using StaticArrays

pixelsize = 1e-6
filename = "sfere"
scale = 0.0001

function load_obj(filename::AbstractString)
    # Read file contents
    file_contents = read(filename, String)

    # Extract vertices and faces
    vertices = []
    faces = []
    for line in eachline(IOBuffer(file_contents))
        components = split(line)
        if length(components) > 0 && components[1] == "v"
            # Extract vertex coordinates
            vertex = []
            for c in components[2:4]
                push!(vertex, parse(Float64, c))
            end
			vertex = SVector(vertex...)
            push!(vertices, vertex)
        elseif length(components) > 0 && components[1] == "f"
            # Extract face vertices
            face = []
            for c in components[2:end]
                push!(face, (parse.(Int, split(c, "/"))[1]))
            end
			face = tuple(face...)
            push!(faces, face)
        end
    end

    println("Loaded mesh from $filename with $(size(vertices, 1)) vertices and $(size(faces, 1)) faces.")

    println("Example: $(vertices[1]), $(faces[1])")
    return vertices, faces
end

function projection(v)
	return (v[2], v[3])
end

function triangle_contains_point(tri1,tri2,tri3, y::Float64, z::Float64)

    # Get the 3 vectors representing, whith respect to the projection of tri1 on the
    # yz plane, the other vertices and [y,z]
    v0 = projection(tri2) .- projection(tri1)
    v1 = projection(tri3) .- projection(tri1)
    v2 = (y,z) .- projection(tri1)

    # Do the maths of using barycentric coordinates. Idk what this means, but
    # it seems like ChatGPT knows what they're doing.
    denom = v0[1]*v1[2] - v0[2]*v1[1]
    u = (v2[1]*v1[2] - v2[2]*v1[1]) / denom
    v = (v0[1]*v2[2] - v0[2]*v2[1]) / denom
    w = 1 - u - v

    # Check if P is inside the triangle
    return (0 <= u <= 1 && 0 <= v <= 1 && 0 <= w <= 1)
end


function project_face_on_pixel(v1::SVector{3,Float64},v2::SVector{3,Float64},v3::SVector{3,Float64}, min_y_tri::Float64, max_y_tri::Float64, min_z_tri::Float64, max_z_tri::Float64, is_degenerate::Bool, normal::SVector{3,Float64}, d::Float64, point::Tuple{Float64,Float64})
	y,z = point

	# Check if we are within the bounding box
	if (y < min_y_tri) || (y > max_y_tri) || (z < min_z_tri) || (z > max_z_tri) || is_degenerate
		return 0.
	end

	# Check if point is inside the face. I currently don't
	# trust this method, but it's all I have so far. Might
	# write another if it doesn't work.
	# But somehow, it works... Must study some more maths...
	if triangle_contains_point(v1,v2,v3, y, z)
		return Float64(round( (-dot(normal,SVector(0,y,z)) - d) / normal[1] ; digits=5 ) * sign(normal[1]))
	end
	return 0.
end


function project_face_on_matrix(face::Tuple{Int64, Int64, Int64}, vertices::Vector{SVector{3,Float64}}, matrix::CuArray{Tuple{Float64, Float64}})
	# Calculate plane of face
	v1 = vertices[face[1]]
	v2 = vertices[face[2]]
	v3 = vertices[face[3]]

	normal = cross(v2 - v1, v3 - v1)
	d = -dot(normal, v1)

	# If the normal vector to the triangle has no x component, then the projection
	# on the yz plane will be degenerate
	is_degenerate = (normal[1] == 0)

	# Find the bounding box of the triangle
	max_y_tri = maximum( ( v1[2], v2[2], v3[2] ) )
	min_y_tri = minimum( ( v1[2], v2[2], v3[2] ) )
	max_z_tri = maximum( ( v1[3], v2[3], v3[3] ) )
	min_z_tri = minimum( ( v1[3], v2[3], v3[3] ) )

	return map(point -> project_face_on_pixel(v1,v2,v3, min_y_tri, max_y_tri, min_z_tri, max_z_tri, is_degenerate, normal, d, point), matrix)
end
#if !(any(isapprox.(extrema, x, rtol=0.00001)))


# function reducer(x1, x2)
# 	return findall(!isnan, union(x1,x2))
# end

function project_mesh(vertices::Vector{SVector{3,Float64}}, faces::Vector{Tuple{Int64, Int64, Int64}}, pixelsize::Float64, scale::Float64)

	vertices = [v * scale for v in vertices]
    nverts = size(vertices, 1)

    # Find the bounding box of the projection of the mesh on the yz plane
    # i.e. minimum and maximum values of y and z.
    min_y = minimum([ vertices[i][2] for i in 1:nverts ])
    max_y = maximum([ vertices[i][2] for i in 1:nverts ])
    min_z = minimum([ vertices[i][3] for i in 1:nverts ])
    max_z = maximum([ vertices[i][3] for i in 1:nverts ])

	y_range = (min_y-pixelsize):pixelsize:(max_y + pixelsize)
	z_range = (min_z-pixelsize):pixelsize:(max_z + pixelsize)

	println(max_y-min_y, " ", max_z-min_z)


    # Find the normalization constant to fit the mesh into the [0,1]x[0,1] square
    # i.e. the maximum side length of the bounding box. This will be used to
    # correctly spread the pixels throughout the 3D space.
    #scaling = maximum([max_y-min_y, max_z-min_z])

	n = size(y_range)[1]
	m = size(z_range)[1]

	println(n, " ", m)

	# Initialize output matrix
    matrix = CuArray([(y,z) for z in z_range, y in y_range])

	facets = CuArray(faces)
    println("Projecting matrix ...")
    println()


	result = mapreduce(face -> project_face_on_matrix(face, vertices, matrix), +, facets)

    println("Done projecting matrix!           ")

    return Array(result)
end

#contatore = 0
# function thickness(vector)
#     thick = 0
#     leng = length(vector)
# 	if (leng == 1)
# 		return 0
# 	end
#     if (leng % 2 != 0)
#         #global contatore
#         println("Oppa, qua ne ho trovati $leng : $vector")
#         #contatore += 1
#     end
#     for i in 1:length(vector)
#         thick += (-1)^i * vector[i]
#     end
#     return thick
# end

image = project_mesh(load_obj("obj\\$filename.obj")...,pixelsize, scale)

#image = Matrix{Float64}(thickness.(projection_matrix))
#image = projection_matrix
clamper = scaleminmax(0,maximum(image))
image =  map(clamper, image)

resolution = size(image)

save("src_img\\$(filename)_$(resolution[1])_$(resolution[2]).png", image)
