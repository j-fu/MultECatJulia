module MultECatJulia

using ExtendableGrids
using SimplexGridFactory
using Triangulate



greet() = print("Hello World!")

include("grids.jl")
export foursquare_electrode_grid, circular_electrode_grid, circular_anisoref_electrode_grid



end # module
