### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using Revise
    using PlutoUI
    using ExtendableGrids
    using PlutoVista
    using DrWatson: plotsdir
    using VoronoiFVM
    using GridVisualize
    using PhysicalConstants.CODATA2018

    import CairoMakie
    using LaTeXStrings
    using Colors
    using HypertextLiteral
    default_plotter!(CairoMakie)
    using MultECatJulia
end

# ╔═╡ cd35c466-1f29-48cb-bad4-b997d63d079d
begin
    using Unitful
    SI(x) = Float64(Unitful.ustrip(Unitful.upreferred(1 * x)))
    const N_A = SI(AvogadroConstant)

    const V = SI(Unitful.V)
    const C = SI(Unitful.C)
    const eV = SI(Unitful.eV)
    const nm = SI(Unitful.nm)
    const cm = SI(Unitful.cm)
    const μF = SI(Unitful.μF)
    const L = SI(Unitful.L)

    const mol = N_A

end

# ╔═╡ 882dda23-63b9-4b1e-a04e-69071deff69a
md"This notebook is only relocateable together with the whole MultECatJulia project."

# ╔═╡ f36552fd-affd-44e0-83b5-3401459f0560
TableOfContents(title = "")

# ╔═╡ 0f2a3eb0-9818-4bf8-9dde-ba28ea3dd2f5
htl""" Description of the equilibrium problem: 
<a href="./open?path=$(joinpath(@__DIR__,\"..\",\"src\",\"equilibrium.jl\"))" target="_blank"> here </a>"""

# ╔═╡ 5391c130-6706-4f06-8335-6474875fe210
md"""
### Data
"""

# ╔═╡ 7062de01-4986-45eb-8bcf-73ff23445dc2
Vmax = 0.75 * V

# ╔═╡ 0c70dffb-d9d6-4446-b335-3899f916b06b
μ_e = 0.1 * eV

# ╔═╡ 98ebf56d-7788-426a-b681-866c91920f76
κ = 10

# ╔═╡ 12c8879e-5807-42e9-addb-20831ea1d40c
Width = 50 * nm

# ╔═╡ aabe1f5c-14d9-4b7d-9ccf-25ee565b9b7c
ngrain = 5

# ╔═╡ 3318676c-aefa-4a81-8a79-1a35edf62021
hgmin = min(0.5 * nm, Width / (15 * ngrain))

# ╔═╡ 047a8134-c309-4e00-b894-ec1ab78fb282
grid = polycrystal_grid2d(ngrain = ngrain, W = Width, hgmin = hgmin, H = 10nm, hzmax = 1nm)

# ╔═╡ 2fdef862-536a-4401-8b5c-bc7e80b6a224
gridplot(grid, resolution = (600, 300), gridscale = 1 / 1nm, zoom = 3, linewidth = 0.1)

# ╔═╡ c6067d85-41c6-43c8-be60-17cf94dff2b4
grid[CellRegions] |> extrema

# ╔═╡ 1263e426-c510-47d4-a0c3-6023cf98c11f
md"""
### Solution for various applied voltages
"""

# ╔═╡ 589ed15d-f8af-4921-8e49-82449170ef5a
begin
    data = EquilibriumData()
    data.κ .= κ
    data.z .= [-1, 1]
    data.μ_e = [ i % 2 == 0 ? μ_e : -μ_e for i in 1:ngrain]
    update_derived!(data)
    data
end

# ╔═╡ 074dcc08-0fe8-4585-a17b-9133e80b3f4f
md"""
__Pressure Poisson ?__ $(@bind ppoisson CheckBox())
"""

# ╔═╡ bb64df89-c43c-43c8-9b0f-a19cc8611414
sys = ppoisson ? create_equilibrium_pp_system(grid, data; Γ_bulk = ngrain + 1) : create_equilibrium_system(grid, data)

# ╔═╡ 1ab27251-0999-49a7-a970-29e70d8fd800
inival = unknowns(sys, inival = 0);

# ╔═╡ b3e6eabf-5e9a-41eb-9ad5-0c823edffe01
plot_fullgrid = false

# ╔═╡ 7899d9b8-0edc-443f-a5a9-96a01187ff74
md"""
Change applied voltage: $(@bind voltage Slider(-Vmax:0.05:Vmax, show_value=true,default=0.0))
"""

# ╔═╡ 9545c38c-35d4-4adf-bd80-82af84e93564
begin
    apply_voltage!(sys, voltage)
    sol = solve(sys, inival = inival, log = true, damp_initial = 0.1)
    hist = history(sol)
end

# ╔═╡ 5e4623cc-af45-4b48-a761-666dd0d98427
scalarplot(hist, resolution = (500, 200), yscale = :log)

# ╔═╡ 11e94e40-fbcf-4d03-ab86-09cc2e75c5c4
cmol = calc_cmol(sol, sys)

# ╔═╡ 10c87e6d-4485-4d17-b7f2-8ead685e17f4
q = calc_QBL(sol, sys)

# ╔═╡ ad87c830-41cf-4a77-b200-5ad08dbe1251
begin
    n = size(sol, 2)
    qdensity = zeros(n)
    for i in 1:n
        qdensity[i] = MultECatJulia.spacecharge(sol[iφ, i], sol[ip, i], data)
    end
    qdensity
end

# ╔═╡ 8b6d6ef7-282c-49ba-9ad4-c33ac72b6833
begin
    if plot_fullgrid
        plotgrid = grid
        plotϕ = calc_φ(sol, sys)
        plotqd = qdensits
    else
        plotgrid = subgrid(grid, [2])
        plotϕ = view(calc_φ(sol, sys), plotgrid)
        plotqd = view(qdensity, plotgrid)
    end
end

# ╔═╡ 1765d933-c6a4-4302-b19e-98d27954c0b4
let
    vis = GridVisualizer(resolution = (600, 400), legend = :rt, dim = 2, linewidth = 0.5, layout = (2, 1), xlabel = "x/nm", ylabel = "y/nm")
    scalarplot!(
        vis[1, 1], plotgrid, plotϕ, clear = true, label = "φ/V", gridscale = 1 / nm, zoom = 3, colormap = :bwr, levels = 9, limits = (-Vmax, Vmax), aspect = 6,
        title = L"\fontfamily{TeXGyreHeros}\mathrm{Potential}\; ϕ/V"
    )
    scalarplot!(
        vis[2, 1], plotgrid, plotqd / (C / cm^3); gridscale = 1 / nm,
        limits = (-5.0e2, 5.0e2),
        colormap = :bwr, aspect = 6,
        title = L"\fontfamily{TeXGyreHeros}\mathrm{Charge\; density\;} q/(C/cm^3)"
    )
    reveal(vis)
end

# ╔═╡ 768bc478-339f-4bb5-8736-e0377d219744
md"""
### Double layer capacitance
"""

# ╔═╡ 11a3ebb7-ef4d-4788-a042-96857dacf9c5
nsteps = 25

# ╔═╡ fd0c34b7-1b52-4430-86c6-dc4f3cb4f1b6
@bind run2d CheckBox()

# ╔═╡ c20f9bde-5ff5-47fe-b543-758159f8add5
molarity = 0.025

# ╔═╡ 34f7bee9-e858-433c-bb78-26970f5d84a0
set_molarity!(data, molarity);L_Debye(data) / nm

# ╔═╡ 2d0ca845-f24b-4970-94c0-c78c19939362
V2d, C2d = let
    if run2d
        calc_Cdl(
            sys, vmax = Vmax, nsteps = 50, molarity = molarity, verbose = false,
            δV = 1.0e-6 * V
        )
    else
        zeros(1), zeros(1)
    end
end

# ╔═╡ 944d215e-c049-4824-9694-eea8a81775c5
grid1d = simplexgrid(geomspace(0, 20nm, 0.1nm, 2nm))

# ╔═╡ 9d812aa2-ddef-486c-887a-cfc13b03dad9
begin
    molarity
    data1D = deepcopy(data)
    data1D.μ_e = [data.μ_e[1]]
    update_derived!(data1D)
end

# ╔═╡ d04af3dd-d489-4f93-bdac-32b786c1723b
sys1d = ppoisson ? create_equilibrium_pp_system(grid1d, data1D, Γ_bulk = 2) : create_equilibrium_system(grid1d, data1D)

# ╔═╡ b769149c-2a51-4460-8281-6b70069842db
μ_e1 = μ_e

# ╔═╡ 44a97c3d-057a-4bc9-9804-5f818e08a1de
Vp, Cp = let
    data1D.μ_e = [μ_e1]
    calc_Cdl(sys1d, vmax = Vmax, nsteps = 100, molarity = molarity)
end

# ╔═╡ ab490968-e69f-4935-901d-6f0ce60cc101
Vm, Cm = let
    data1D.μ_e = [-μ_e1]
    calc_Cdl(sys1d, vmax = Vmax, nsteps = nsteps = 100, molarity = molarity)
end

# ╔═╡ 8fa7773e-8d92-47d2-8480-77549e510e38
V0, C0 = let
    data1D.μ_e = [0eV]
    calc_Cdl(sys1d, vmax = Vmax, nsteps = 100, molarity = molarity)
end

# ╔═╡ d176f826-b8ec-4bdf-b75f-55f96e596a34
let
    vis = GridVisualizer(
        resolution = (500, 400), legend = :rt, Plotter = CairoMakie,
        ylabel = L"\fontfamily{TeXGyreHeros}C_{dl}/(μF/cm^2)",
        xlabel = "U/V",
        title = "Double layer capacitance"
    )
    scalarplot!(vis, Vp, Cp / (μF / cm^2), color = RGB(1.0, 0.7, 0.7), label = L"μ_E>0")
    scalarplot!(vis, Vm, Cm / (μF / cm^2), clear = false, color = RGB(0.7, 0.7, 1), label = L"μ_E<0")
    scalarplot!(vis, V0, C0 / (μF / cm^2), clear = false, color = :gray90, label = L"μ_E=0")
    scalarplot!(vis, Vm, (Cm + Cp) / 2 / (μF / cm^2), clear = false, color = RGB(1, 0.0, 1.0), label = "1d avg")
    if length(V2d) > 1
        scalarplot!(vis, V2d, (C2d / Width) / (μF / cm^2), color = :green, clear = false, label = "2d", markershape = :none)
    end
    if run2d
        GridVisualize.save(draftplotsdir("2DResults.pdf"), vis)
    end
    reveal(vis)
end

# ╔═╡ Cell order:
# ╟─882dda23-63b9-4b1e-a04e-69071deff69a
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─f36552fd-affd-44e0-83b5-3401459f0560
# ╟─0f2a3eb0-9818-4bf8-9dde-ba28ea3dd2f5
# ╟─5391c130-6706-4f06-8335-6474875fe210
# ╠═cd35c466-1f29-48cb-bad4-b997d63d079d
# ╠═7062de01-4986-45eb-8bcf-73ff23445dc2
# ╠═0c70dffb-d9d6-4446-b335-3899f916b06b
# ╠═98ebf56d-7788-426a-b681-866c91920f76
# ╠═12c8879e-5807-42e9-addb-20831ea1d40c
# ╠═aabe1f5c-14d9-4b7d-9ccf-25ee565b9b7c
# ╠═3318676c-aefa-4a81-8a79-1a35edf62021
# ╠═047a8134-c309-4e00-b894-ec1ab78fb282
# ╠═2fdef862-536a-4401-8b5c-bc7e80b6a224
# ╠═c6067d85-41c6-43c8-be60-17cf94dff2b4
# ╟─1263e426-c510-47d4-a0c3-6023cf98c11f
# ╠═589ed15d-f8af-4921-8e49-82449170ef5a
# ╟─074dcc08-0fe8-4585-a17b-9133e80b3f4f
# ╠═bb64df89-c43c-43c8-9b0f-a19cc8611414
# ╠═1ab27251-0999-49a7-a970-29e70d8fd800
# ╠═9545c38c-35d4-4adf-bd80-82af84e93564
# ╠═5e4623cc-af45-4b48-a761-666dd0d98427
# ╠═11e94e40-fbcf-4d03-ab86-09cc2e75c5c4
# ╠═10c87e6d-4485-4d17-b7f2-8ead685e17f4
# ╠═ad87c830-41cf-4a77-b200-5ad08dbe1251
# ╠═b3e6eabf-5e9a-41eb-9ad5-0c823edffe01
# ╠═8b6d6ef7-282c-49ba-9ad4-c33ac72b6833
# ╠═1765d933-c6a4-4302-b19e-98d27954c0b4
# ╟─7899d9b8-0edc-443f-a5a9-96a01187ff74
# ╟─768bc478-339f-4bb5-8736-e0377d219744
# ╠═11a3ebb7-ef4d-4788-a042-96857dacf9c5
# ╠═fd0c34b7-1b52-4430-86c6-dc4f3cb4f1b6
# ╠═c20f9bde-5ff5-47fe-b543-758159f8add5
# ╠═34f7bee9-e858-433c-bb78-26970f5d84a0
# ╠═2d0ca845-f24b-4970-94c0-c78c19939362
# ╠═d176f826-b8ec-4bdf-b75f-55f96e596a34
# ╠═944d215e-c049-4824-9694-eea8a81775c5
# ╠═9d812aa2-ddef-486c-887a-cfc13b03dad9
# ╠═d04af3dd-d489-4f93-bdac-32b786c1723b
# ╠═b769149c-2a51-4460-8281-6b70069842db
# ╠═44a97c3d-057a-4bc9-9804-5f818e08a1de
# ╠═ab490968-e69f-4935-901d-6f0ce60cc101
# ╠═8fa7773e-8d92-47d2-8480-77549e510e38
