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


resolution = 1024
filename = "gruppo"

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
            push!(vertices, Vector{Float64}(vertex))
        elseif length(components) > 0 && components[1] == "f"
            # Extract face vertices
            face = []
            for c in components[2:end]
                push!(face, (parse.(Int, split(c, "/"))[1]))
            end
            push!(faces, Vector{Int}(face))
        end
    end

    println("Loaded mesh from $filename with $(size(vertices, 1)) vertices and $(size(faces, 1)) faces.")

    println("Example: $(vertices[1]), $(faces[1])")
    return vertices, faces
end

function triangle_contains_point(tri1,tri2,tri3, y::Float64, z::Float64)
    # Get the x, y, and z coordinates of the triangle vertices
    x1, y1, z1 = tri1
    x2, y2, z2 = tri2
    x3, y3, z3 = tri3

    # Get the 3 vectors representing, whith respect to the projection of v1 on the
    # yz plane, the other vertices and [y,z]
    v0 = [y2 - y1, z2 - z1]
    v1 = [y3 - y1, z3 - z1]
    v2 = [y - y1, z - z1]

    # Do the maths of using barycentric coordinates. Idk what this means, but
    # it seems like ChatGPT knows what they're doing.
    denom = v0[1]*v1[2] - v0[2]*v1[1]
    u = (v2[1]*v1[2] - v2[2]*v1[1]) / denom
    v = (v0[1]*v2[2] - v0[2]*v2[1]) / denom
    w = 1 - u - v

    # Check if P is inside the triangle
    return (0 <= u <= 1 && 0 <= v <= 1 && 0 <= w <= 1)
end



function project_mesh(vertices, faces, n)
    # Initialize output matrix
    matrix = [[[] for j in 1:n] for i in 1:n]

    nverts = size(vertices, 1)

    # Find the bounding box of the projection of the mesh on the yz plane
    # i.e. minimum and maximum values of y and z.
    min_y = minimum([ vertices[i][2] for i in 1:nverts ])
    max_y = maximum([ vertices[i][2] for i in 1:nverts ])
    min_z = minimum([ vertices[i][3] for i in 1:nverts ])
    max_z = maximum([ vertices[i][3] for i in 1:nverts])

    for i in 1:nverts
        vertices[i] -= [0, min_y, min_z]
    end

    # Find the normalization constant to fit the mesh into the [0,1]x[0,1] square
    # i.e. the maximum side length of the bounding box. This will be used to
    # correctly spread the pixels throughout the 3D space.
    scaling = maximum([max_y-min_y, max_z-min_z])

    println("Projecting matrix ...")
    println()

    n_faces = length(faces)
    face_counter = 0
    for face in faces

        # Print the progress as a percentage of faces mapped
        percentage = div(face_counter*100,n_faces)
        face_counter += 1
        print("$percentage%", "   ")
        print("\r")

        # Calculate plane of face
        v1 = vertices[face[1]]
        v2 = vertices[face[2]]
        v3 = vertices[face[3]]

        normal = cross(v2 - v1, v3 - v1)
        d = -dot(normal, v1)

        # If the normal vector to the triangle has no x component, then the projection
        # on the yz plane will be degenerate
        is_degenerate = (normal[1] == 0)

        # if (is_degenerate)
        #     println("degenere")
        # end

        # Find the bounding box of the triangle
        max_y_tri = maximum( [ v1[2], v2[2], v3[2] ] )
        min_y_tri = minimum( [ v1[2], v2[2], v3[2] ] )
        max_z_tri = maximum( [ v1[3], v2[3], v3[3] ] )
        min_z_tri = minimum( [ v1[3], v2[3], v3[3] ] )

        # Iterate over rows of matrix
        Threads.@threads for i = 1:n
            # Compute the corresponding y value
            y = (i - 0.5) * scaling / n

            # Check if we are within the bounding box
            if (y < min_y_tri)
                continue
            elseif (y > max_y_tri)
                break
            end

            # Iterate over columns of matrix
            for j = 1:n
                # Compute the corresponding z value
                z = (j - 0.5) * scaling / n

                # Check if we are within the bounding box
                if (z < min_z_tri)
                    continue
                elseif (z > max_z_tri)
                    break
                end

                # Where the x coordinates of the surface of the mesh
                # will be saved.
                extrema = matrix[i][j]

                # Check if point is inside the face. I currently don't
                # trust this method, but it's all I have so far. Might
                # write another if it doesn't work.
                # But somehow, it works... Must study some more maths...
                if triangle_contains_point(v1,v2,v3, y, z)

                    # Calculate where the point with
                    # given yz coordinates is on the mesh
                    # through the equation of the plane.
                    # If the triangle is projection-degenerate,
                    # save its center of mass
                    if (is_degenerate)
                        x = sum([v1, v2, v3])[1]/3
                    else
                        x = (-dot(normal,[0,y,z]) - d) / normal[1]
                    end

                    #If it's not already there...
                    if !(any(isapprox.(extrema, x, rtol=0.00001)))

                        # If it is contained, then save it.
                        push!(extrema, x)

                        # If the projection of the triangle is degenerate, then save the projection twice
                        # This way, the thickness will be unaffected by this point.
                        if (is_degenerate)
                            push!(extrema, x)
                        end
                    end
                end
            end
            # Sort the results, which should be in even number
            # I should add a check to see if one point is tangent.
            # In that case, calculating integrals would be a pain in the ass...

        end
    end

    println("Done projecting matrix!           ")

    matrix = hcat([s for s in matrix]...)
    sort!.(matrix)

    return matrix
end

#contatore = 0
function thickness(vector)
    thick = 0
    leng = length(vector)
    if (leng % 2 != 0)
        #global contatore
        println("Oppa, qua ne ho trovati $leng : $vector")
        #contatore += 1
    end
    for i in 1:length(vector)
        thick += (-1)^i * vector[i]
    end
    return thick
end

projection_matrix = project_mesh(load_obj("obj\\$filename.obj")..., resolution)

image = Matrix{Float64}(thickness.(projection_matrix))

clamper = scaleminmax(0,max)
image =  map(clamper, image)

save("src_img\\$filename.png", image)
