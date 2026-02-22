@testset "generators" begin
    basename = "gmsh_line_mesh"
    path = joinpath(tempdir, basename * ".msh")

    n_partitions = 2 # with `3`, the result is not deterministic...
    BcubeGmsh.gen_line_mesh(
        path;
        nx = 12,
        n_partitions = n_partitions,
        split_files = true,
        create_ghosts = true,
    )

    for i in 1:n_partitions
        fname = basename * "_$i.msh"
        @test compare_meshes_helper(joinpath(tempdir, fname), joinpath(ASSETS_DIR, fname))
    end

    # gen_mesh_around_disk
    fname = "gmsh_mesh_around_disk_tri.msh"
    BcubeGmsh.gen_mesh_around_disk(joinpath(tempdir, fname), :tri; nθ = 30, nr = 10)
    @test compare_meshes_helper(
        joinpath(tempdir, fname),
        joinpath(ASSETS_DIR, fname),
        1e-10,
    )
    fname = "gmsh_mesh_around_disk_quad.msh"
    BcubeGmsh.gen_mesh_around_disk(joinpath(tempdir, fname), :quad; nθ = 30, nr = 10)
    @test compare_meshes_helper(
        joinpath(tempdir, fname),
        joinpath(ASSETS_DIR, fname),
        1e-10,
    )

    # gen_mesh_around_disk
    fname = "gmsh_rectangle_mesh_with_tri_and_quad.msh"
    BcubeGmsh.gen_rectangle_mesh_with_tri_and_quad(
        joinpath(tempdir, fname);
        nx = 5,
        ny = 6,
        lx = 2,
        ly = 3,
    )
    @test compare_meshes_helper(
        joinpath(tempdir, fname),
        joinpath(ASSETS_DIR, fname),
        1e-14,
    )

    # gen_circle_mesh
    fname = "gmsh_circle_mesh.msh"
    BcubeGmsh.gen_circle_mesh(joinpath(tempdir, fname), 30)
    @test compare_meshes_helper(joinpath(tempdir, fname), joinpath(ASSETS_DIR, fname))

    # gen_disk_mesh
    fname = "gmsh_disk_mesh.msh"
    BcubeGmsh.gen_disk_mesh(joinpath(tempdir, fname))
    @test compare_meshes_helper(
        joinpath(tempdir, fname),
        joinpath(ASSETS_DIR, fname),
        1e-15,
    )

    # # gen_star_disk_mesh
    fname = "gmsh_star_disk_mesh.msh"
    BcubeGmsh.gen_star_disk_mesh(joinpath(tempdir, fname), 0.1, 7; nθ = 100)
    @test compare_meshes_helper(
        joinpath(tempdir, fname),
        joinpath(ASSETS_DIR, fname),
        1e-15,
    )

    # gen_cylinder_shell_mesh (quad)
    fname = "gmsh_cylinder_shell_mesh_quad.msh"
    BcubeGmsh.gen_cylinder_shell_mesh(joinpath(tempdir, fname); nθ = 20, nz = 8)
    @test compare_meshes_helper(
        joinpath(tempdir, fname),
        joinpath(ASSETS_DIR, fname),
        1e-15,
    )

    # gen_cylinder_shell_mesh (tri)
    # Use `compare_mesh_nodes_cloud_helper` because node ordering (and c2n) varies between machines...
    fname = "gmsh_cylinder_shell_mesh_tri.msh"
    BcubeGmsh.gen_cylinder_shell_mesh(joinpath(tempdir, fname); lc = 0.1)
    @test compare_mesh_nodes_cloud_helper(
        joinpath(tempdir, fname),
        joinpath(ASSETS_DIR, fname),
        7e-5,
    )

    # gen_cylinder_shell_mesh
    fname = "gmsh_ring_mesh.msh"
    BcubeGmsh.gen_ring_mesh(joinpath(tempdir, fname); r_int = 1.0, r_ext = 2.0, lc = 0.1)
    @test compare_meshes_helper(
        joinpath(tempdir, fname),
        joinpath(ASSETS_DIR, fname),
        1e-14,
    )
end
