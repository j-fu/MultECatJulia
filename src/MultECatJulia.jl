module MultECatJulia
using LinearAlgebra
using ExtendableGrids
using SimplexGridFactory
using Triangulate
using VoronoiFVM
using Unitful
import PhysicalConstants
using PhysicalConstants.CODATA2018


greet() = print("Hello World!")


include("equilibrium.jl")
export EquilibriumData,default_data,L_debye, update_derived!,set_molarity!
export iφ,ip,iA,iC
export create_equilibrium_system, solve_equilibrium_system,apply_voltage!
export calc_φ,calc_p, calc_cmol,calc_cnum, calc_QBL
export calc_Cdl,Cdl0, L_Debye
include("grids.jl")
export foursquare_electrode_grid, circular_electrode_grid, circular_anisoref_electrode_grid
export polycrystal_grid2d


end # module
