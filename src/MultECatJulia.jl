module MultECatJulia
using LinearAlgebra
using ExtendableGrids
using SimplexGridFactory
using Triangulate
using VoronoiFVM
using Unitful
using NLsolve
using Parameters
import PhysicalConstants
using PhysicalConstants.CODATA2018


greet() = print("Hello World!")


include("equilibrium.jl")
export EquilibriumData,L_debye, apply_voltage!,set_molarity!,update_derived!
export iφ,ip,iA,iC
export create_equilibrium_system, solve_equilibrium_system,create_equilibrium_pp_system
export calc_φ,calc_p, calc_cmol,calc_c0mol,calc_cnum, calc_QBL,ysum
export calc_Cdl,Cdl0, L_Debye
include("grids.jl")
export foursquare_electrode_grid, circular_electrode_grid, circular_anisoref_electrode_grid
export polycrystal_grid2d


end # module
