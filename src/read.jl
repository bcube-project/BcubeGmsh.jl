struct GmshMetaData <: Bcube.AbstractMeshMetaData
    name2cells::Dict{String, Vector{Int}}
end

Bcube.get_zone_names(metadata::GmshMetaData, ::Bcube.AbstractMesh) = keys(metadata.names)
function Bcube.get_zone_element_indices(metadata::GmshMetaData, ::Bcube.AbstractMesh, name)
    metadata.name2cells[name]
end

function Bcube.read_mesh(
    ::GmshIoHandler,
    filepath::String;
    domains = String[],
    spacedim::Int = 0,
    verbose::Bool = false,
    kwargs...,
)
    @assert length(domains) == 0 "Reading only some domains is not supported yet (but easy to implement)"
    return read_msh(filepath, spacedim; verbose)
end

"""
    read_msh(path::String, spaceDim::Int = 0; verbose::Bool = false)

Read a .msh file designated by its `path`.

See `read_msh()` for more details.
"""
function read_msh(path::String, spaceDim::Int = 0; verbose::Bool = false)
    isfile(path) ? nothing : error("File does not exist ", path)

    # Read file using gmsh lib
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", Int(verbose))
    gmsh.open(path)

    # build mesh
    mesh = _read_msh(spaceDim, verbose)

    # free gmsh
    gmsh.finalize()

    return mesh
end

"""
    _read_msh(spaceDim::Int, verbose::Bool)

To use this function, the `gmsh` file must have been opened already (see `read_msh(path::String)`
for instance).

The number of topological dimensions is given by the highest dimension found in the file. The
number of space dimensions is deduced from the axis dimensions if `spaceDim = 0`.
If `spaceDim` is set to a positive number, this number is used as the number of space dimensions.

# Implementation
Global use of `gmsh` module. Do not try to improve this function by passing an argument
such as `gmsh` or `gmsh.model` : it leads to problems.
"""
function _read_msh(spaceDim::Int, verbose::Bool)
    # Spatial dimension of the mesh
    dim = gmsh.model.getDimension()

    # Read nodes
    ids, xyz = gmsh.model.mesh.getNodes()

    # Create a node number remapping to ensure a dense numbering
    absolute_node_indices = [convert(Int, i) for i in ids]
    _, glo2loc_node_indices = Bcube.densify(absolute_node_indices; permute_back = true)

    # Build nodes coordinates
    xyz = reshape(xyz, 3, :)
    _spaceDim = spaceDim > 0 ? spaceDim : _compute_space_dim(verbose)
    nodes = [Node(xyz[1:_spaceDim, i]) for i in axes(xyz, 2)]

    # Read cells
    elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(dim)

    # Create a cell number remapping to ensure a dense numbering
    absolute_cell_indices = Int.(reduce(vcat, elementTags))
    _, glo2loc_cell_indices = Bcube.densify(absolute_cell_indices; permute_back = true)

    # Read boundary conditions
    bc_tags = gmsh.model.getPhysicalGroups(-1)
    bc_names = [gmsh.model.getPhysicalName(_dim, _tag) for (_dim, _tag) in bc_tags]
    # keep only physical groups of dimension "dim-1" with none-empty names.
    # bc is a vector of (tag,name) for all valid boundary conditions
    bc = [
        (_tag, _name) for
        ((_dim, _tag), _name) in zip(bc_tags, bc_names) if _dim == dim - 1 && _name ≠ ""
    ]

    bc_names = Dict(convert(Int, _tag) => _name for (_tag, _name) in bc)
    bc_nodes = Dict(
        convert(Int, _tag) => Int[
            glo2loc_node_indices[i] for
            i in gmsh.model.mesh.getNodesForPhysicalGroup(dim - 1, _tag)[1]
        ] for (_tag, _name) in bc
    )

    # Fill type of each cell
    celltypes = [
        GMSHTYPE[k] for (i, k) in enumerate(elementTypes) for t in 1:length(elementTags[i])
    ]

    # Build cell->node connectivity (with Gmsh internal numbering convention)
    c2n_gmsh = Bcube.Connectivity(
        Int[nnodes(k) for k in reduce(vcat, celltypes)],
        Int[glo2loc_node_indices[k] for k in reduce(vcat, nodeTags)],
    )

    # Convert to CGNS numbering
    c2n = _c2n_gmsh2cgns(celltypes, c2n_gmsh)

    # Read volumic physical groups (build a dict tag -> name)
    el_tags = gmsh.model.getPhysicalGroups(_spaceDim)
    _el_names = [gmsh.model.getPhysicalName(_dim, _tag) for (_dim, _tag) in el_tags]
    el = [
        (_tag, _name) for ((_dim, _tag), _name) in zip(el_tags, _el_names) if
        _dim == _spaceDim && _name ≠ ""
    ]
    el_names = Dict(convert(Int, _tag) => _name for (_tag, _name) in el)
    # el_names_inv = Dict(_name => convert(Int, _tag) for (_tag, _name) in el)

    # Read cell indices associated to each volumic physical group
    el_cells = Dict{Int, Array{Int}}()
    for (_dim, _tag) in el_tags
        v = Int[]

        for iEntity in gmsh.model.getEntitiesForPhysicalGroup(_dim, _tag)
            tmpTypes, tmpTags, tmpNodeTags = gmsh.model.mesh.getElements(_dim, iEntity)

            # Notes : a PhysicalGroup "entity" can contain different types of elements.
            # So `tmpTags` is an array of the cell indices of each type in the Physical group.
            for _tmpTags in tmpTags
                v = vcat(v, Int.(_tmpTags)) # would a "push!" be a better alternative?
            end
        end
        el_cells[_tag] = v
    end

    # Create the name => cells dict
    names2cells = Dict(
        name => [glo2loc_cell_indices[i] for i in el_cells[tag]] for
        (tag, name) in el_names
    )
    metadata = GmshMetaData(names2cells)

    mesh = Bcube.Mesh(
        nodes,
        celltypes,
        c2n;
        bc_names = bc_names,
        bc_nodes = bc_nodes,
        metadata,
    )
    Bcube.add_absolute_indices!(mesh, :node, absolute_node_indices)
    Bcube.add_absolute_indices!(mesh, :cell, absolute_cell_indices)
    return mesh
end

"""
Deduce the number of space dimensions from the mesh : if one (or more) dimension of the bounding
box is way lower than the other dimensions, the number of space dimension is decreased.

Currently, having for instance (x,z) is not supported. Only (x), or (x,y), or (x,y,z).
"""
function _compute_space_dim(verbose::Bool)
    tol = 1e-15

    topodim = gmsh.model.getDimension()

    # Bounding box
    box = gmsh.model.getBoundingBox(-1, -1)
    lx = box[4] - box[1]
    ly = box[5] - box[2]
    lz = box[6] - box[3]

    return Bcube._compute_space_dim(topodim, lx, ly, lz, tol, verbose)
end

function _apply_gmsh_options(;
    split_files = false,
    create_ghosts = false,
    msh_format = 0,
    verbose = false,
)
    gmsh.option.setNumber("General.Terminal", Int(verbose))
    gmsh.option.setNumber("Mesh.PartitionSplitMeshFiles", Int(split_files))
    gmsh.option.setNumber("Mesh.PartitionCreateGhostCells", Int(create_ghosts))
    (msh_format > 0) && gmsh.option.setNumber("Mesh.MshFileVersion", msh_format)
end
