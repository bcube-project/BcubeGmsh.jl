module BcubeGmsh
using Bcube
import gmsh_jll
include(gmsh_jll.gmsh_api)
import .gmsh

include("./common.jl")
include("./read.jl")
include("./generators.jl")

end
