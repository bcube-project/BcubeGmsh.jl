"""
Custom way to "include" a file to print infos.
"""
function custom_include(path)
    filename = split(path, "/")[end]
    print("Running test file " * filename * "...")
    include(path)
    println("done.")
end

function isapprox_arrays(a::AbstractArray, b::AbstractArray; atol::Real = eps())
    success = all(abs.(a .- b) .< atol)
    (success == false) && (@show a, b, maximum(abs.(a .- b)))
    return success
end

"""
    compare_connectivities(x::T, y::T) where T

Check if two `Bcube.Connectivity` are identical.

Note that this function works only because connectivities are (usually) composed of
"primitive" types (Int and vectors of Int), not structs.
"""
function compare_connectivities(x::Bcube.Connectivity, y::Bcube.Connectivity)
    return all(
        fname -> getproperty(x, fname) == getproperty(y, fname),
        fieldnames(typeof(x)),
    )
end

"""
    compare_meshes(mesh_a::AbstractMesh, mesh_b::AbstractMesh, tol)

Compare two meshes with the given `tol` for floating numbers. Return true if
the two meshes are identical (with respect to the tol), false otherwise.
"""
function compare_meshes(
    mesh_a::AbstractMesh,
    mesh_b::AbstractMesh,
    atol = eps();
    verbose = false,
)
    try
        identical = true

        # Compare nodes
        if !all(
            ((x, y),) -> isapprox_arrays(x, y; atol),
            zip(get_coords.(get_nodes(mesh_a)), get_coords.(get_nodes(mesh_b))),
        )
            verbose && println("Comparison of nodes failed")
            identical = false
        end
        if verbose
            d = maximum(
                ((x, y),) -> norm(x - y),
                zip(get_coords.(get_nodes(mesh_a)), get_coords.(get_nodes(mesh_b))),
            )
            println("Maximum distance between two nodes : $d")
        end

        # Compare c2n
        if !compare_connectivities(
            Bcube.connectivities_indices(mesh_a, :c2n),
            Bcube.connectivities_indices(mesh_b, :c2n),
        )
            verbose && println("Comparison of c2n failed")
            identical = false
        end

        # Compare bc tags
        bcnodes_a = Bcube.boundary_nodes(mesh_a)
        bcnodes_b = Bcube.boundary_nodes(mesh_b)
        if Set(keys(bcnodes_a)) != Set(keys(bcnodes_b))
            if verbose
                println("BC tags differ")
                @show keys(bcnodes_a)
                @show keys(bcnodes_b)
            end
            identical = false
        end

        # Compare bc nodes
        if !all(
            tag -> Bcube.boundary_nodes(mesh_a, tag) == Bcube.boundary_nodes(mesh_b, tag),
            keys(bcnodes_a),
        )
            verbose && println("BC nodes differ")
            identical = false
        end

        return identical

    catch e
        if verbose
            @show e
            println(
                "Mesh comparison failed, most likely because the two meshes are very different (ex : different number of nodes)",
            )
        end
        return false
    end
end

""" `filepath_a` should be handled with Bcube IO interface while `filepath_b` should be handled by Serialization """
function compare_meshes_helper(
    filepath_a::String,
    filepath_b::String,
    tol = eps();
    verbose = true,
)
    mesh_a = read_mesh(filepath_a)
    mesh_b = read_mesh(filepath_b)
    return compare_meshes(mesh_a, mesh_b, tol; verbose)
end

""" Point cloud comparison : avoid using this costly function on large meshes """
function compare_mesh_nodes_cloud(
    mesh_a::AbstractMesh,
    mesh_b::AbstractMesh,
    tol = eps();
    verbose = false,
)
    if nnodes(mesh_a) != nnodes(mesh_b)
        verbose && println("Different number of nodes")
        return false
    end

    coords_a = zeros(nnodes(mesh_a), Bcube.spacedim(mesh_a))
    coords_b = zeros(nnodes(mesh_b), Bcube.spacedim(mesh_b))

    for (row_a, a, row_b, b) in zip(
        eachrow(coords_a),
        get_coords.(get_nodes(mesh_a)),
        eachrow(coords_b),
        get_coords.(get_nodes(mesh_b)),
    )
        row_a .= a
        row_b .= b
    end

    verbose && println("Computing point cloud distances with $(nnodes(mesh_a)) nodes...")
    d = pairwise(Euclidean(), coords_a, coords_b; dims = 1)
    d_min = map(minimum, eachrow(d))
    max_of_min = maximum(d_min)
    success = max_of_min < tol # alternatively, `success = maximum(d_min) < tol`
    verbose && println("d_max = $(max_of_min)")
    return success
end

""" `filepath_a` should be handled with Bcube IO interface while `filepath_b` should be handled by Serialization """
function compare_mesh_nodes_cloud_helper(
    filepath_a::String,
    filepath_b::String,
    tol = eps();
    verbose = true,
)
    mesh_a = read_mesh(filepath_a)
    mesh_b = read_mesh(filepath_b)

    return compare_mesh_nodes_cloud(mesh_a, mesh_b, tol; verbose)
end
