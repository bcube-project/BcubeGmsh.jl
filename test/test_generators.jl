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
        @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))
    end

    # gen_mesh_around_disk
    fname = "gmsh_mesh_around_disk_tri.msh"
    BcubeGmsh.gen_mesh_around_disk(joinpath(tempdir, fname), :tri; nθ = 30, nr = 10)
    @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))
    fname = "gmsh_mesh_around_disk_quad.msh"
    BcubeGmsh.gen_mesh_around_disk(joinpath(tempdir, fname), :quad; nθ = 30, nr = 10)
    @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))

    # gen_mesh_around_disk
    fname = "gmsh_rectangle_mesh_with_tri_and_quad.msh"
    BcubeGmsh.gen_rectangle_mesh_with_tri_and_quad(
        joinpath(tempdir, fname);
        nx = 5,
        ny = 6,
        lx = 2,
        ly = 3,
    )
    @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))

    # gen_circle_mesh
    fname = "gmsh_circle_mesh.msh"
    BcubeGmsh.gen_circle_mesh(joinpath(tempdir, fname), 30)
    @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))

    # gen_disk_mesh
    fname = "gmsh_disk_mesh.msh"
    BcubeGmsh.gen_disk_mesh(joinpath(tempdir, fname))
    @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))

    # gen_star_disk_mesh
    fname = "gmsh_star_disk_mesh.msh"
    BcubeGmsh.gen_star_disk_mesh(joinpath(tempdir, fname), 0.1, 7; nθ = 100)
    @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))

    # gen_cylinder_shell_mesh
    fname = "gmsh_cylinder_shell_mesh.msh"
    BcubeGmsh.gen_cylinder_shell_mesh(joinpath(tempdir, fname), 30, 10)
    @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))

    # gen_cylinder_shell_mesh
    fname = "gmsh_ring_mesh.msh"
    BcubeGmsh.gen_ring_mesh(joinpath(tempdir, fname); r_int = 1.0, r_ext = 2.0, lc = 0.1)
    @test fname2sum[fname] == bytes2hex(open(sha1, joinpath(tempdir, fname)))
end
