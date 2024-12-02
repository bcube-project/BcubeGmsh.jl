# BcubeGmsh.jl
Implementation of [`Bcube`](https://github.com/bcube-project/Bcube.jl) IO interface for Gmsh format. Checkout the relative `Bcube` [documentation](https://bcube-project.github.io/Bcube.jl/stable/api/io/io_interface/) for more infos.

For now, only the `read_file` interface is implemented.

## Basic usage
```julia
using Bcube
using BcubeGmsh

mesh = read_mesh("output.msh")
@show ncells(mesh)
```

## Mesh generators
This project also contains some "common" mesh generators. Here is a non-exhaustive list:
* `gen_line_mesh(
        output;
        nx = 2,
        lx = 1.0,
        xc = 0.0,
        order = 1,
        bnd_names = ("LEFT", "RIGHT"),
        n_partitions = 0,
        kwargs...
    )`
* `gen_rectangle_mesh(
        output,
        type;
        transfinite = false,
        nx = 2,
        ny = 2,
        lx = 1.0,
        ly = 1.0,
        xc = -1.0,
        yc = -1.0,
        order = 1,
        bnd_names = ("North", "South", "East", "West"),
        n_partitions = 0,
        write_geo = false,
        transfinite_lines = true,
        lc = 1e-1,
        kwargs...
    )`
* `gen_hexa_mesh(
        output,
        type;
        transfinite = false,
        nx = 2,
        ny = 2,
        nz = 2,
        lx = 1.0,
        ly = 1.0,
        lz = 1.0,
        xc = -1.0,
        yc = -1.0,
        zc = -1.0,
        order = 1,
        bnd_names = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"),
        n_partitions = 0,
        write_geo = false,
        transfinite_lines = true,
        lc = 1e-1,
        kwargs...,
    )`
* `gen_sphere_mesh(
        output;
        radius = 1.0,
        lc = 1e-1,
        order = 1,
        n_partitions = 0,
        kwargs...,
    )`
* `gen_cylinder_mesh(
    output,
    Lz,
    nz;
    radius = 1.0,
    lc = 1e-1,
    order = 1,
    n_partitions = 0,
    kwargs...,
)`
* `gen_mesh_around_disk(
        output,
        type;
        r_in = 1.0,
        r_ext = 10.0,
        nθ = 360,
        nr = 100,
        nr_prog = 1.05,
        order = 1,
        recombine = true,
        bnd_names = ("Farfield", "Wall"),
        n_partitions = 0,
        kwargs...
    )`
* `gen_disk_mesh(
        output;
        radius = 1.0,
        lc = 1e-1,
        order = 1,
        n_partitions = 0,
        kwargs...
    )`

