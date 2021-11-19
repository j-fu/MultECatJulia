module MultECatJulia

using ExtendableGrids
using SimplexGridFactory
using Triangulate
using VoronoiFVM
using Unitful
import PhysicalConstants.CODATA2018


greet() = print("Hello World!")

include("grids.jl")
export foursquare_electrode_grid, circular_electrode_grid, circular_anisoref_electrode_grid

include("equilibrium.jl")
export EquilibriumData,default_data,L_debye, update_derived!
export iφ,ip,iA,iC
export create_equilibrium_system, solve_equilibrium_system,apply_voltage!
export nm,V
export calc_φ,calc_p, calc_cmol,calc_cnum, calc_QBL, MPa, GPa
export calc_Cdl,Cdl0
end # module
