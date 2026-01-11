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

        # Compare c2n
        if !compare_connectivities(
            Bcube.connectivities_indices(mesh_a, :c2n),
            Bcube.connectivities_indices(mesh_b, :c2n),
        )
            verbose && println("Comparison of c2n failed")
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
    @assert endswith(filepath_b, SERIALIZED_EXT)
    mesh_b = Serialization.deserialize(filepath_b)
    return compare_meshes(mesh_a, mesh_b, tol; verbose)
end