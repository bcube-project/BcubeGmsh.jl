using Test
using BcubeGmsh
using Bcube
using SHA
using DelimitedFiles
import Bcube:
    nodes,
    cells,
    Node_t,
    Bar2_t,
    Tri3_t,
    has_nodes,
    has_vertices,
    has_edges,
    has_faces,
    has_cells,
    has_entities,
    n_entities,
    nvertices,
    nedges,
    nfaces,
    boundary_nodes,
    boundary_faces,
    connectivities_indices

"""
Custom way to "include" a file to print infos.
"""
function custom_include(path)
    filename = split(path, "/")[end]
    print("Running test file " * filename * "...")
    include(path)
    println("done.")
end

# This dir will be removed at the end of the tests
tempdir = mktempdir()

# Reading sha1 checksums
filename = "checksums"
(get(ENV, "GITHUB_CI", "false") == "true") && (filename *= "-github")
f = readdlm(joinpath(@__DIR__, "$(filename).sha1"), String)
fname2sum = Dict(r[2] => r[1] for r in eachrow(f))

@testset "BcubeGmsh.jl" begin
    custom_include("./test_read.jl")
    custom_include("./test_generators.jl")
end
