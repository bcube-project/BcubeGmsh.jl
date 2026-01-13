using Test
using BcubeGmsh
using Bcube
using Serialization
using LinearAlgebra
using Distances
import Bcube:
    AbstractMesh,
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

const SERIALIZED_EXT = ".bcube-mesh.serialized"
const REF_DIR = joinpath(@__DIR__, "references")

include(joinpath(@__DIR__, "utils.jl"))

# This dir will be removed at the end of the tests
tempdir = mktempdir()

@testset "BcubeGmsh.jl" begin
    custom_include("./test_read.jl")
    custom_include("./test_generators.jl")
end
