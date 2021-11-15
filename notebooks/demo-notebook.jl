### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
	using Revise
	using PlutoUI
	using MultECatJulia
	using MySubPackage
end

# ╔═╡ 882dda23-63b9-4b1e-a04e-69071deff69a
md"This notebook is only relocateable together with the whole MultECatJulia project."

# ╔═╡ a8e37976-5db2-485f-87aa-0cf7155e8e00
MultECatJulia.greet()

# ╔═╡ 73a94305-1c3c-45e5-969f-2a245baec10d
MySubPackage.greet()

# ╔═╡ Cell order:
# ╠═882dda23-63b9-4b1e-a04e-69071deff69a
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╠═a8e37976-5db2-485f-87aa-0cf7155e8e00
# ╠═73a94305-1c3c-45e5-969f-2a245baec10d
