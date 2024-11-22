struct GmshIoHandler <: Bcube.AbstractIoHandler end

Bcube._filename_to_handler(::Val{:msh}) = GmshIoHandler()

# Constants
const GMSHTYPE = Dict(
    1 => Bar2_t(),
    2 => Tri3_t(),
    3 => Quad4_t(),
    4 => Tetra4_t(),
    5 => Hexa8_t(),
    6 => Penta6_t(),
    7 => Pyra5_t(),
    9 => Tri6_t(),
    10 => Quad9_t(),
    21 => Tri9_t(),
    21 => Tri10_t(),
    36 => Quad16_t(),
)

"""
    nodes_gmsh2cgns(entity::AbstractEntityType, nodes::AbstractArray)

Reorder `nodes` of a given `entity` from the Gmsh format to CGNS format.

See https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering

# Implementation
By default, same numbering between CGNS and Gmsh is applied. Specialized the function
`nodes_gmsh2cgns(e::Type{<:T}) where {T <: Bcube.AbstractEntityType}` to secify a
different numbering
"""
function nodes_gmsh2cgns(entity::Bcube.AbstractEntityType, nodes::AbstractArray)
    map(i -> nodes[i], nodes_gmsh2cgns(entity))
end

nodes_gmsh2cgns(e) = nodes_gmsh2cgns(typeof(e))

function nodes_gmsh2cgns(e::Type{<:T}) where {T <: Bcube.AbstractEntityType}
    nodes(e)
end
function nodes_gmsh2cgns(::Type{Bcube.Hexa27_t})
    SA[
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        17,
        18,
        19,
        20,
        13,
        14,
        15,
        16,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
    ]
end

"""
Convert a cell->node connectivity with gmsh numbering convention to a cell->node connectivity
with CGNs numbering convention.
"""
function _c2n_gmsh2cgns(celltypes, c2n_gmsh)
    n = Int[]
    indices = Int[]
    for (ct, c2nᵢ) in zip(celltypes, c2n_gmsh)
        append!(n, length(c2nᵢ))
        append!(indices, nodes_gmsh2cgns(ct, c2nᵢ))
    end
    return Connectivity(n, indices)
end