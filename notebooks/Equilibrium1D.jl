### A Pluto.jl notebook ###
# v0.17.2

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
	using PyPlot
	using GridVisualize
	using VoronoiFVM
	import DrWatson: plotsdir
	using Colors
	using MultECatJulia
	PyPlot.svg(true)
end

# ╔═╡ 14f8c67a-759a-4646-811c-01d03e3cf726
if isdefined(Main,:PlutoRunner)
	using PlutoVista
	default_plotter!(PlutoVista)
	
end

# ╔═╡ b43533f6-a948-418c-8539-2d54aa8e5943
begin
	using Unitful
	SI(x)=Float64(Unitful.ustrip(Unitful.upreferred(1*x)));
	const V=SI(Unitful.V)
	const eV=SI(Unitful.eV)
	const nm=SI(Unitful.nm)
	const cm=SI(Unitful.cm)
	const μF=SI(Unitful.μF)
end


# ╔═╡ 882dda23-63b9-4b1e-a04e-69071deff69a
md"This notebook is only relocateable together with the whole MultECatJulia project."

# ╔═╡ f36552fd-affd-44e0-83b5-3401459f0560
TableOfContents(title="")

# ╔═╡ 0f2a3eb0-9818-4bf8-9dde-ba28ea3dd2f5
htl""" Description of the equilibrium problem: 
<a href="./open?path=$(joinpath(@__DIR__,\"..\",\"src\",\"equilibrium.jl\"))" target="_blank"> here </a>"""

# ╔═╡ 5bce5b30-b5fd-4e13-8b20-8218edaa6c60
md"""
Results should coincide with Fuhrmann, CPC 2015
"""

# ╔═╡ 5391c130-6706-4f06-8335-6474875fe210
md"""
### Geometry
"""

# ╔═╡ 7062de01-4986-45eb-8bcf-73ff23445dc2
Vmax=2*V

# ╔═╡ 6cd423b6-a6a7-485c-8531-99a51e19de14
L=20nm

# ╔═╡ 0254a47d-33cf-4af7-bb72-c584b9bf98b6
hmin=0.001*nm

# ╔═╡ 7a123ca0-d407-4dc0-9c01-f90ac53b690d
hmax=1*nm

# ╔═╡ 2fdef862-536a-4401-8b5c-bc7e80b6a224
X=ExtendableGrids.geomspace(0,L,hmin,hmax)

# ╔═╡ ac75d0a2-df48-48be-af40-dbb7b7670bf2
grid=ExtendableGrids.simplexgrid(X)

# ╔═╡ 1263e426-c510-47d4-a0c3-6023cf98c11f
md"""
### Solution for various applied voltages
"""

# ╔═╡ 4f57312e-1b08-43a3-b038-d30bfc62e754
data=default_data();data.χ=78; data.κ.=0.0; set_molarity!(data,0.01)

# ╔═╡ d0776266-c208-4237-99c8-1364527eaeb1
Cdl0(data)/(μF/cm^2)

# ╔═╡ bb64df89-c43c-43c8-9b0f-a19cc8611414
sys=create_equilibrium_system(grid,data)

# ╔═╡ 1ab27251-0999-49a7-a970-29e70d8fd800
inival=unknowns(sys,inival=0);

# ╔═╡ fcfd99a5-8213-47cc-826f-95b3f3cdb4e8
vis=GridVisualizer(resolution=(600,200),legend=:rt,Plotter=PlutoVista,limits=(-55,55));vis

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
	scalarplot!(vis,grid,clear=false,calc_p(sol,sys)/SI(u"GPa"),color=:magenta,label="p/GPa")
	reveal(vis)
end

# ╔═╡ 768bc478-339f-4bb5-8736-e0377d219744
md"""
### Double layer capacitance
"""

# ╔═╡ c20f9bde-5ff5-47fe-b543-758159f8add5
molarities=[0.001,0.01,0.1,1]

# ╔═╡ d176f826-b8ec-4bdf-b75f-55f96e596a34
let 
	vis=GridVisualizer(resolution=(500,300),legend=:rt,clear=true,ylabel="C_dl/(μF/cm^2)",Plotter=PyPlot)
	hmol=1/length(molarities)
	for imol=1:length(molarities)
		c=RGB(1-imol*hmol,0,imol*hmol)
		volts,caps=calc_Cdl(sys,vmax=1V,nsteps=100,molarity=molarities[imol])
		cdl0=Cdl0(sys.physics.data)
	    scalarplot!(vis,volts,caps/(μF/cm^2),color=c,clear=false,label="$(molarities[imol])M",markershape=:none)
		scalarplot!(vis,[0],[cdl0]/(μF/cm^2),clear=false,markershape=:circle,label="")
	end
	save(plotsdir("1DResults.pdf"),vis)
	reveal(vis)
end

# ╔═╡ Cell order:
# ╟─882dda23-63b9-4b1e-a04e-69071deff69a
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═14f8c67a-759a-4646-811c-01d03e3cf726
# ╟─f36552fd-affd-44e0-83b5-3401459f0560
# ╟─0f2a3eb0-9818-4bf8-9dde-ba28ea3dd2f5
# ╟─5bce5b30-b5fd-4e13-8b20-8218edaa6c60
# ╠═b43533f6-a948-418c-8539-2d54aa8e5943
# ╟─5391c130-6706-4f06-8335-6474875fe210
# ╠═7062de01-4986-45eb-8bcf-73ff23445dc2
# ╠═6cd423b6-a6a7-485c-8531-99a51e19de14
# ╠═0254a47d-33cf-4af7-bb72-c584b9bf98b6
# ╠═7a123ca0-d407-4dc0-9c01-f90ac53b690d
# ╠═2fdef862-536a-4401-8b5c-bc7e80b6a224
# ╠═ac75d0a2-df48-48be-af40-dbb7b7670bf2
# ╟─1263e426-c510-47d4-a0c3-6023cf98c11f
# ╠═4f57312e-1b08-43a3-b038-d30bfc62e754
# ╠═d0776266-c208-4237-99c8-1364527eaeb1
# ╠═bb64df89-c43c-43c8-9b0f-a19cc8611414
# ╠═1ab27251-0999-49a7-a970-29e70d8fd800
# ╟─9545c38c-35d4-4adf-bd80-82af84e93564
# ╠═5e4623cc-af45-4b48-a761-666dd0d98427
# ╠═11e94e40-fbcf-4d03-ab86-09cc2e75c5c4
# ╠═10c87e6d-4485-4d17-b7f2-8ead685e17f4
# ╠═fcfd99a5-8213-47cc-826f-95b3f3cdb4e8
# ╟─1765d933-c6a4-4302-b19e-98d27954c0b4
# ╟─7899d9b8-0edc-443f-a5a9-96a01187ff74
# ╟─768bc478-339f-4bb5-8736-e0377d219744
# ╠═c20f9bde-5ff5-47fe-b543-758159f8add5
# ╠═d176f826-b8ec-4bdf-b75f-55f96e596a34
