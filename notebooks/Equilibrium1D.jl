### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using HypertextLiteral
	using PlutoUI
	using ExtendableGrids
	using PlutoVista
	using GridVisualize
	using VoronoiFVM
	using Colors
	default_plotter!(PlutoVista)
	using MultECatJulia
end

# ╔═╡ 882dda23-63b9-4b1e-a04e-69071deff69a
md"This notebook is only relocateable together with the whole MultECatJulia project."

# ╔═╡ f36552fd-affd-44e0-83b5-3401459f0560
TableOfContents(title="")

# ╔═╡ 0f2a3eb0-9818-4bf8-9dde-ba28ea3dd2f5
htl""" Description of the equilibrium problem: 
<a href="./open?path=$(joinpath(@__DIR__,\"..\",\"src\",\"equilibrium.jl\"))" target="_blank"> here </a>"""

# ╔═╡ 5391c130-6706-4f06-8335-6474875fe210
md"""
### Geometry
"""

# ╔═╡ 7062de01-4986-45eb-8bcf-73ff23445dc2
Vmax=2*V

# ╔═╡ 6cd423b6-a6a7-485c-8531-99a51e19de14
L=90*nm

# ╔═╡ 0254a47d-33cf-4af7-bb72-c584b9bf98b6
hmin=0.01nm

# ╔═╡ 7a123ca0-d407-4dc0-9c01-f90ac53b690d
hmax=10*nm

# ╔═╡ 2fdef862-536a-4401-8b5c-bc7e80b6a224
X=ExtendableGrids.geomspace(0,L,hmin,hmax)

# ╔═╡ ac75d0a2-df48-48be-af40-dbb7b7670bf2
grid=ExtendableGrids.simplexgrid(X)

# ╔═╡ 1263e426-c510-47d4-a0c3-6023cf98c11f
md"""
### Solution for various applied voltages
"""

# ╔═╡ bb64df89-c43c-43c8-9b0f-a19cc8611414
sys=create_equilibrium_system(grid)

# ╔═╡ 1ab27251-0999-49a7-a970-29e70d8fd800
inival=unknowns(sys,inival=0)

# ╔═╡ fcfd99a5-8213-47cc-826f-95b3f3cdb4e8
vis=GridVisualizer(resolution=(600,200),legend=:rt);vis

# ╔═╡ 7899d9b8-0edc-443f-a5a9-96a01187ff74
md"""
Change applied voltage. $(@bind voltage Slider(-Vmax:0.05:Vmax,show_value=true,default=0))
"""

# ╔═╡ 9545c38c-35d4-4adf-bd80-82af84e93564
begin
	apply_voltage!(sys,voltage)
	sol,history=solve_equilibrium_system(sys,inival=inival)
	sol
end

# ╔═╡ 5e4623cc-af45-4b48-a761-666dd0d98427
history

# ╔═╡ 11e94e40-fbcf-4d03-ab86-09cc2e75c5c4
cmol=calc_cmol(sol,sys)

# ╔═╡ 10c87e6d-4485-4d17-b7f2-8ead685e17f4
q=calc_QBL(sol,sys)

# ╔═╡ 1765d933-c6a4-4302-b19e-98d27954c0b4
let
	
	scalarplot!(vis,grid,calc_φ(sol,sys),xlimits=(0,5*nm),color=:green,clear=true,label="φ/V")
	scalarplot!(vis,grid,cmol[iA,:],color=:blue,clear=false,label="c_A/(mol/L)")
	scalarplot!(vis,grid,cmol[iC,:],color=:red,clear=false,label="c_C/(mol/L)")
	scalarplot!(vis,grid,clear=false,calc_p(sol,sys)/GPa,color=:magenta,label="p/GPa")
	reveal(vis)
end

# ╔═╡ 768bc478-339f-4bb5-8736-e0377d219744
md"""
### Double layer capacitance
"""

# ╔═╡ c20f9bde-5ff5-47fe-b543-758159f8add5
molarities=[0.01,0.1,1,10]

# ╔═╡ d176f826-b8ec-4bdf-b75f-55f96e596a34
let 
	vis=GridVisualizer(resolution=(600,300),legend=:rt,clear=true)
	hmol=1/length(molarities)
	for imol=1:length(molarities)
		c=RGB(1-imol*hmol,0,imol*hmol)
		M=molarities[imol]
		V,C=calc_Cdl(sys,vmax=1,nsteps=100,molarity=M)
		cdl0=Cdl0(sys.physics.data)
	    scalarplot!(vis,V,C,color=c,clear=false,label="$(M)M",markershape=:none)
		scalarplot!(vis,[0],[cdl0],clear=false,markershape=:circle,label="")
	end
	reveal(vis)
end

# ╔═╡ Cell order:
# ╟─882dda23-63b9-4b1e-a04e-69071deff69a
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─f36552fd-affd-44e0-83b5-3401459f0560
# ╟─0f2a3eb0-9818-4bf8-9dde-ba28ea3dd2f5
# ╟─5391c130-6706-4f06-8335-6474875fe210
# ╠═7062de01-4986-45eb-8bcf-73ff23445dc2
# ╠═6cd423b6-a6a7-485c-8531-99a51e19de14
# ╠═0254a47d-33cf-4af7-bb72-c584b9bf98b6
# ╠═7a123ca0-d407-4dc0-9c01-f90ac53b690d
# ╠═2fdef862-536a-4401-8b5c-bc7e80b6a224
# ╠═ac75d0a2-df48-48be-af40-dbb7b7670bf2
# ╟─1263e426-c510-47d4-a0c3-6023cf98c11f
# ╠═bb64df89-c43c-43c8-9b0f-a19cc8611414
# ╠═1ab27251-0999-49a7-a970-29e70d8fd800
# ╠═9545c38c-35d4-4adf-bd80-82af84e93564
# ╠═5e4623cc-af45-4b48-a761-666dd0d98427
# ╠═11e94e40-fbcf-4d03-ab86-09cc2e75c5c4
# ╠═10c87e6d-4485-4d17-b7f2-8ead685e17f4
# ╟─fcfd99a5-8213-47cc-826f-95b3f3cdb4e8
# ╟─1765d933-c6a4-4302-b19e-98d27954c0b4
# ╟─7899d9b8-0edc-443f-a5a9-96a01187ff74
# ╟─768bc478-339f-4bb5-8736-e0377d219744
# ╠═c20f9bde-5ff5-47fe-b543-758159f8add5
# ╠═d176f826-b8ec-4bdf-b75f-55f96e596a34
