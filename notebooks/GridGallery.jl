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
    using GridVisualize

    using MultECatJulia
    default_plotter!(PlutoVista)
end

# ╔═╡ 58ee6727-1074-4efb-830e-d0e684b49065
md"""
## Grid with square regions at electrode
"""

# ╔═╡ 07194f25-5735-452f-9bed-cf791958d44d
@doc foursquare_electrode_grid

# ╔═╡ 9800c97e-c25b-4d08-a014-cccd746f7d71
gxyz = MultECatJulia.foursquare_electrode_grid()

# ╔═╡ 9a219111-0275-4cd9-b97e-648d3fcfcbb9
vis = GridVisualizer(dim = 3, resolution = (400, 300));vis

# ╔═╡ 85be1677-87ff-49dc-af9a-557e575bc55f
@bind z Slider(0:0.01:1, show_value = true, default = 0.5)

# ╔═╡ 33ecc78c-16ac-46ac-8d83-cbb26288a6fb
gridplot!(vis, gxyz, show = true, zplanes = [z], xplanes = [1], yplanes = [1], outlinealpha = 0.3)

# ╔═╡ 9e640256-61a4-4fb9-a449-d5f948fb2d26
md"""
## Grid with inner circle at electrode
"""

# ╔═╡ bf0004f4-b3c8-4385-b87d-df2ef1420408
md"""
### Plain triangulation
"""

# ╔═╡ 0c3f5b1c-7bcf-491e-8318-73f2c880e0b5
@doc circular_electrode_grid

# ╔═╡ 90e32b85-d1e5-49db-a100-30984c98305b
gcxyz = circular_electrode_grid()

# ╔═╡ bff0a3f1-c545-4fcc-8be8-ef460f8479bd
visc = GridVisualizer(dim = 3, resolution = (400, 300))

# ╔═╡ 66ce2939-2957-43e0-b7ea-cdd10d9fc17e
@bind zc Slider(0:0.01:1, show_value = true, default = 0.5)

# ╔═╡ a310027d-3f9e-4de5-bc9e-5457cb19eef5
gridplot!(visc, gcxyz, show = true, zplanes = [zc], xplanes = [1], yplanes = [1], outlinealpha = 0.3)

# ╔═╡ 752bdb8f-de5a-4863-b2d1-53d69aff7dcb
md"""
### Grid with anisotropic local refinement at electrode
"""

# ╔═╡ f43bd8ff-a2ed-4740-a61e-4f5d2c615c10
grxyz = circular_anisoref_electrode_grid(nref = 0)

# ╔═╡ b9ced73f-e7a0-4573-9413-c2a231d21c78
visr = GridVisualizer(dim = 3, resolution = (400, 300))

# ╔═╡ 2c285e4f-9f10-474b-9482-83a5c4cbfe09
@bind zr Slider(0:0.01:1, show_value = true, default = 0.5)

# ╔═╡ 97d8db91-16e6-4070-81a0-7bcaff6cc9f6
gridplot!(visr, grxyz, show = true, zplanes = [zr], xplanes = [1], yplanes = [1], outlinealpha = 0.3)

# ╔═╡ 78dba5b2-52c9-40cd-bc1d-d5c343271f97
html"""<hr> """

# ╔═╡ b9cc0359-7286-4c02-ba10-35303da26a50
TableOfContents()

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─58ee6727-1074-4efb-830e-d0e684b49065
# ╠═07194f25-5735-452f-9bed-cf791958d44d
# ╠═9800c97e-c25b-4d08-a014-cccd746f7d71
# ╠═9a219111-0275-4cd9-b97e-648d3fcfcbb9
# ╠═85be1677-87ff-49dc-af9a-557e575bc55f
# ╠═33ecc78c-16ac-46ac-8d83-cbb26288a6fb
# ╟─9e640256-61a4-4fb9-a449-d5f948fb2d26
# ╟─bf0004f4-b3c8-4385-b87d-df2ef1420408
# ╟─0c3f5b1c-7bcf-491e-8318-73f2c880e0b5
# ╠═90e32b85-d1e5-49db-a100-30984c98305b
# ╠═bff0a3f1-c545-4fcc-8be8-ef460f8479bd
# ╠═66ce2939-2957-43e0-b7ea-cdd10d9fc17e
# ╠═a310027d-3f9e-4de5-bc9e-5457cb19eef5
# ╟─752bdb8f-de5a-4863-b2d1-53d69aff7dcb
# ╠═f43bd8ff-a2ed-4740-a61e-4f5d2c615c10
# ╠═b9ced73f-e7a0-4573-9413-c2a231d21c78
# ╠═2c285e4f-9f10-474b-9482-83a5c4cbfe09
# ╠═97d8db91-16e6-4070-81a0-7bcaff6cc9f6
# ╟─78dba5b2-52c9-40cd-bc1d-d5c343271f97
# ╠═b9cc0359-7286-4c02-ba10-35303da26a50
