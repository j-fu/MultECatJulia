### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
	if isdefined(Main,:PlutoRunner)
	using Pkg
	Pkg.activate(joinpath(@__DIR__,".."))
    using PlutoUI
	using Unitful
    import PhysicalConstants.CODATA2018
	using VoronoiFVM
	using ExtendableGrids
	end
end

# ╔═╡ 4082c3d3-b728-4bcc-b480-cdee41d9ab99
isdefined(Main,:PlutoRunner) && TableOfContents(title="")

# ╔═╡ 6cafe7d3-a017-49d7-ba39-5b686478b18e
md"""
## Physical constants and units
"""

# ╔═╡ aaea99bb-777b-4941-8c96-501c29b34f2e
html"""
Transform unitful data to numerical values in preferred SI units.
See documentation for <a href="http://painterqubits.github.io/Unitful.jl/stable/conversion/#Unitful.upreferred" target="_blank"> upreferred </a> and 
<a href="http://painterqubits.github.io/Unitful.jl/stable/conversion/#Unitful.ustrip" target=_blank> ustrip</a>.
"""

# ╔═╡ 085fec09-8172-49ad-8773-8983e08037b8
sibase(x)=Float64(ustrip(upreferred(x)));

# ╔═╡ 1bb99234-7bbe-48ca-9fb8-be932f8fd7a3
begin
const N_A=sibase(CODATA2018.AvogadroConstant)
const R=sibase(CODATA2018.MolarGasConstant)
const e_0=sibase(CODATA2018.ElementaryCharge)
const ε_0= sibase(CODATA2018.VacuumElectricPermittivity)
const k_B= sibase(CODATA2018.BoltzmannConstant)
end;

# ╔═╡ 7fb7d7af-c14a-4ef1-b065-533d6741583b
begin
const m=sibase(1Unitful.m)
const nm=sibase(1Unitful.nm)
const cm=sibase(1Unitful.cm)
const μm=sibase(1Unitful.μm)
const K=sibase(1Unitful.K)
const Pa=sibase(1Unitful.Pa)
const MPa=sibase(1Unitful.MPa)
const GPa=sibase(1Unitful.GPa)
const μF=sibase(1Unitful.μF)
const L=sibase(1Unitful.L)
const V=sibase(1Unitful.V)
const eV=V*e_0
const mol=N_A
end;

# ╔═╡ 7d77ad32-3df6-4243-8bad-b8df4126e6ea
md"""
## Model data
"""

# ╔═╡ 86a7e69b-20f1-43d4-988f-3b9e5a603478
struct DerivedData
	v::Vector{Float64}
	y_E::Vector{Float64}
	y0_E::Float64
end

# ╔═╡ 0d825f88-cd67-4368-90b3-29f316b72e6e
mutable struct EquilibriumData
	N::Int64
	T::Float64
	p_ref::Float64
	v0::Float64
	χ::Float64
	E_ref::Float64
	μ_e::Vector{Float64}
	n_E::Vector{Float64}
	z::Vector{Float64}
	κ::Vector{Float64}
	derived::DerivedData
	EquilibriumData()=new()
end

# ╔═╡ 30c6a176-935b-423f-9447-86f78746322f
md"""
Debye length
```math
L_{Debye}=\sqrt{ \frac{(1+χ)ε_0k_BT}{e_0^2n_E}}
```
"""

# ╔═╡ 00e536dc-34aa-4a1a-93de-4eb3f5e0a348
L_Debye(data)=sqrt( (1+data.χ)*ε_0*k_B*data.T/(e_0^2*data.n_E[1]) )

# ╔═╡ a21545da-3b53-47af-b0c4-f253b37dc84f
md"""
Double layer capacitance at $φ=0$
```math
C_{dl,0}=\sqrt{\frac{2(1+χ) ε_0e_0^2 n_E}{k_BT}}
```
"""

# ╔═╡ 1d22b09e-99c1-4026-9505-07bdffc98582
Cdl0(data)=sqrt( 2*(1+data.χ)*ε_0*e_0^2*data.n_E[1]/(k_B*data.T))

# ╔═╡ 5a210961-19fc-40be-a5f6-033a80f1414d
md"""
Check with Bard/Faulkner: the value must be 22.8
"""

# ╔═╡ 9b57f6ed-02f8-48ba-afa2-0766fe8c0c4c
md"""
### Species indices
"""

# ╔═╡ 5fed71ec-35fb-4804-99ff-e1eaf18fac1b
begin
	const iφ=1
	const ip=2
	const iA=1
	const iC=2
end;

# ╔═╡ 5eca37ba-f858-45fb-a66a-3795327dfd18
md"""
## Model equations
"""

# ╔═╡ a26cf11b-0ce1-4c1d-a64d-1917178ff676
md"""
Equilibrium expression (16)
```math
y_α(φ,p)=y_α^E\exp\left(\frac{-z_αe_0}{k_BT}(φ- φ^E)-\frac{v_α}{k_BT}(p-p^E)\right)
```
"""

# ╔═╡ 188f67d8-2ae8-474c-8e58-68b8b4fde02e
function y_α(φ,p,α, data)
	η_φ=data.z[α]*e_0*(φ-data.E_ref)
	η_p=data.derived.v[α]*(p-data.p_ref)
	data.derived.y_E[α]*exp( - (η_φ+η_p )/(k_B*data.T))
end

# ╔═╡ d7531d5f-fc2d-42b2-9cf9-6a737b0f0f8d
y0(p,data)=data.derived.y0_E*exp( -data.v0*(p-data.p_ref)/(k_B*data.T))

# ╔═╡ f6f004a6-d71b-4813-a363-9f51dc37e42a
md"""
Poisson equation (32a)

```math
-∇⋅(1+χ)ε_0∇φ = q(φ,p)
```
"""

# ╔═╡ 0e2d20a1-5f26-4263-9a91-3b40b2c2996a
function poisson_flux(f,u,edge,data)
	f[iφ]=(1.0+data.χ)*ε_0*(u[iφ,1]-u[iφ,2])
end

# ╔═╡ 2e11ce81-7d0d-498f-9ddd-7d4d836ab42f
md"""
Solvated ion volumes:
```math
	v_α=(1+κ_α)v_0
```
"""

# ╔═╡ b1e062c6-f245-4edc-aa02-871e2c776998
md"""
Incompressibility condition (14)
```math
\begin{aligned}
1&=∑\limits_α v_αn_α= n ∑\limits_α v_α y_α\\
n&=\frac1{∑\limits_α v_α y_α}
\end{aligned}
```
"""

# ╔═╡ c4cc940c-74aa-45f8-a2fa-6016d7c3c145
md"""
Space charge
```math
\begin{aligned}
q(φ,p)&=e_0∑\limits_α z_αn_α = ne_0∑\limits_α z_αy_α\\
      &=e_0\frac{∑\limits_α z_αy_α(\phi,p)}{∑\limits_α v_α y_α(\phi,p)}
\end{aligned}
```
"""

# ╔═╡ a199d4dd-34f0-40a4-ac7e-2642d12dccc1
md"""
Definition of ``y_\alpha`` (32b)
```math
∑_α y_α(φ,p)=1
```
"""

# ╔═╡ 042a452a-1130-4a56-a1b9-b2674803e445
function spacecharge_and_ysum(f,u,node,data)
	φ=u[iφ]
	p=u[ip]
	y=y0(p,data)
	sumy=y
	sumyz=zero(eltype(u))
	sumyv=data.v0*y
	for α=1:data.N
		y=y_α(φ,p,α,data)
		sumy+=y
		sumyz+=data.z[α]*y
		sumyv+=data.derived.v[α]*y
	end
	f[iφ]=-e_0*sumyz/sumyv
	f[ip]=sumy-1.0
end

# ╔═╡ 13fc2859-496e-4f6e-8b22-36d9d55768b8
md"""
Calculate bulk mole fractions from incompressibiltiy:
```math
\begin{aligned}
∑\limits_αv_αn_α^E&=1\\
n_0^E&=\frac1{v_0}\left(1-∑\limits_{α>0}v_αn_α^E\right)\\
n^E&=\frac1{v_0}\left(1-∑\limits_{α>0}v_αn_α^E\right)+ ∑\limits_{α>0}n_α^E\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}(v_α-v_0)n_α^E\right)\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}((1+ κ_α)v_0-v_0)n_α^E\right)\\
   &=\frac1{v_0}\left(1-∑\limits_{α>0}κ_αv_0n_α^E\right)\\
   &=\frac1{v_0}-∑\limits_{α>0}κ_αn_α^E\\
y_α^E&=\frac{n_α^E}{n^E}
\end{aligned}
```
"""

# ╔═╡ 32db42f3-5084-4908-9b53-59291b6133c5
function update_derived!(data)
	v=(1.0.+data.κ).*data.v0
	n_E_all=1/data.v0
	N=length(data.κ)
	for α=1:N
		n_E_all-=data.κ[α]*data.n_E[α]
	end
	y_E=data.n_E/n_E_all
	y0_E=(1/data.v0)/n_E_all
	data.derived=DerivedData(v,y_E,y0_E)
	data
end

# ╔═╡ 5cd5a8a9-ad67-42d9-a998-34c038ef9688
function default_data(;n_E=0.1*mol/L,nref=55.508*mol/L)
	data=EquilibriumData()
	data.μ_e=[0.0]
	data.N=2
	data.T=298.15*K
    data.p_ref=1.0e5*Pa
	data.E_ref=0
    data.v0=1.0/nref
	data.χ=15
	data.n_E=[n_E,n_E]
	data.κ=[10,10]
	data.z=[-1,1]
	update_derived!(data)
	data
end

# ╔═╡ 0f4f3598-808d-4e01-bf4a-99ba079940a6
default_data()

# ╔═╡ 1065b3e0-60bf-497c-b7fb-c5a065737f77
L_Debye(default_data(n_E=0.01*mol/L))/nm

# ╔═╡ fe704fb4-d07c-4591-b834-d6cf2f4f7075
let
	data=default_data(n_E=0.01*mol/L)
	data.χ=78.49-1
	update_derived!(data)
    Cdl0(data)/(μF/cm^2)
end

# ╔═╡ 243d27b5-a1b8-4127-beec-d5643ad07855
md"""
Bulk boundary condition (32d)
```math
∇ φ\to 0
```
so we choose 
```math
\partial_n φ=0
```
"""


# ╔═╡ 005289e8-6979-49fe-b20f-66afd207baea
md"""
Electrode boundary condition (32c) !!! bug in paper

```math
φ|_{Σ_i}=\frac{1}{e_0} \mu_{e,i}- (E-E^{ref})
```
"""

# ╔═╡ 0c5ed337-9310-417d-a1f6-7d69dd8c377b
φ_Σ(ifacet,data,E)=data.μ_e[ifacet]/e_0-(E-data.E_ref)

# ╔═╡ 65955950-2879-4b8c-bf73-d63e07d2ad96
md"""
Poly surface charge (34)

```math
Q_s(E)= ∑_i s_i \hat Q_s(E-E^{ref} + \frac1{e_0}\mu_{e,i})
```
"""

# ╔═╡ d8f80c62-b2d6-456f-9650-e8102e968673
md"""
Single surface charge (26)

```math
\hat Q_s = \frac{∑_\alpha z_αe_0y_{s,α} +\sum_α \sum_β ν_{αβ}z_αe_0y_{s,β}}{a_V^{ref}+ ∑_{α} a_α^{ref}z_αe_0y_{s,α}} 
```

We assume that there are no surface reactions, so we assume ``\hat Q_s=0`` and ``Q_s=0``.
"""

# ╔═╡ f4b2f509-0769-4df7-956e-e8bfc9ccd89a
md"""
Boundary layer charge (35)
```math
Q^{BL}(E)=-\frac{1}{Σ} ∫_{Ω_E} q dx
```
"""

# ╔═╡ 0bbd9482-d17d-4027-8eec-450807cff792
md"""
## System setup and solution
"""

# ╔═╡ d885ac23-ddfa-495c-b93b-54032c8a5c1f
function apply_voltage!(sys,E)
	data=sys.physics.data
	nbc=num_bfaceregions(sys.grid)
	nfacets=length(data.μ_e)
	@assert nbc>nfacets
	for ifacet=1:nfacets
		boundary_dirichlet!(sys,iφ,ifacet,φ_Σ(ifacet,data,E))
	end
	sys
end

# ╔═╡ 6e3dbf34-c1c9-460b-9546-b2ee8ee99d68
function create_equilibrium_system(grid,data::EquilibriumData=default_data())
	update_derived!(data)
physics=VoronoiFVM.Physics(data=data,flux=poisson_flux,reaction=spacecharge_and_ysum)
sys=VoronoiFVM.System(grid,physics,unknown_storage=:dense)
nreg=num_cellregions(grid)	
enable_species!(sys,iφ,1:nreg)
enable_species!(sys,ip,1:nreg)
apply_voltage!(sys,0)
end

# ╔═╡ 85fa321c-dd93-4ac3-96da-1e1bff80970a
function solve_equilibrium_system(sys;inival=unknowns(sys,inival=0),damp=0.1,log=true)
	c=NewtonControl()
	c.damp_initial=damp
	VoronoiFVM.solve(inival,sys,control=c,log=log)
end

# ╔═╡ 93428d11-a3dc-4e29-ae6d-48ba37082c74
md"""
### Postprocessing
"""

# ╔═╡ f3279037-01ed-4596-8e5a-86afe4c02c5f
calc_φ(sol,sys)=sol[iφ,:]

# ╔═╡ 2afd54ca-4240-4f07-b38a-242ba0485b45
calc_p(sol,sys)=sol[ip,:]

# ╔═╡ 55bd7b9a-a191-4a0b-9c6b-13733be5023e
md"""
number concentration:
```math
	n_α=ny_α
```
"""

# ╔═╡ 3ceda3b1-bf1c-4126-b94f-2ee03e8dde99
function c_num!(c,φ,p, data)
	y=y0(p,data)
	sumyv=data.v0*y
	for α=1:data.N
		c[α]=y_α(φ,p,α,data)
		sumyv+=c[α]*data.derived.v[α]
	end
	c./=sumyv
end

# ╔═╡ 800dfed8-9f29-4138-96f8-e8bf1f2f00e6
function calc_cnum(sol,sys)
	data=sys.physics.data
	grid=sys.grid
	nnodes=num_nodes(grid)
	conc=zeros(data.N,nnodes)
	for i=1:nnodes
	  @views c_num!(conc[:,i],sol[iφ,i],sol[ip,i],data)
	end
	conc
end

# ╔═╡ 2ee34d76-7238-46c2-94d1-a40d8b017af6
calc_cmol(sol,sys)=calc_cnum(sol,sys)/(mol/L)

# ╔═╡ 49466829-9459-4dc8-85cc-c67460e290d2
calc_QBL(sol,sys)=VoronoiFVM.integrate(sys,spacecharge_and_ysum,sol)[iφ,1]

# ╔═╡ bb9d283d-cbad-4066-90fe-5b9db0f57db3
md"""
## Double layer capacitance
"""

# ╔═╡ 77f49da5-ffd2-4148-93a6-f45382ba6d91
function calc_Cdl(sys;vmax=2*V,molarity=1,nsteps=21, δV=1.0e-3*V)
	data=sys.physics.data
	data.n_E .= molarity*mol
	update_derived!(data)
	apply_voltage!(sys,0)
	inival=unknowns(sys,inival=0)
	inival=solve(inival,sys)
	vstep=vmax/(nsteps-1)
	c=VoronoiFVM.NewtonControl()
	c.damp_growth=1.2
	c.verbose=false
	c.tol_round=1.0e-10
	c.max_round=3

	function rundlcap(dir)
		V=zeros(0)
		C=zeros(0)
	sol=copy(inival)
	oldsol=copy(inival)
	volt=0.0
	for iv=1:nsteps
		apply_voltage!(sys,volt)
		
		c.damp_initial=0.1
		solve!(sol,oldsol,sys,control=c)
		oldsol.=sol
		Q=calc_QBL(sol,sys)
		apply_voltage!(sys,volt+dir*δV)
		c.damp_initial=1
	    solve!(sol,oldsol,sys,control=c)
		oldsol.=sol
		Qδ=calc_QBL(sol,sys)
	    push!(C,(Q-Qδ)/(dir*δV))
		push!(V,volt)
		volt+=dir*vstep
	end
		V,C
	end
	Vf,Cf=rundlcap(1)
	Vr,Cr=rundlcap(-1)
	vcat(reverse(Vr),Vf),vcat(reverse(Cr),Cf)
end

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─4082c3d3-b728-4bcc-b480-cdee41d9ab99
# ╟─6cafe7d3-a017-49d7-ba39-5b686478b18e
# ╟─aaea99bb-777b-4941-8c96-501c29b34f2e
# ╠═085fec09-8172-49ad-8773-8983e08037b8
# ╠═1bb99234-7bbe-48ca-9fb8-be932f8fd7a3
# ╠═7fb7d7af-c14a-4ef1-b065-533d6741583b
# ╟─7d77ad32-3df6-4243-8bad-b8df4126e6ea
# ╠═86a7e69b-20f1-43d4-988f-3b9e5a603478
# ╠═0d825f88-cd67-4368-90b3-29f316b72e6e
# ╠═5cd5a8a9-ad67-42d9-a998-34c038ef9688
# ╠═0f4f3598-808d-4e01-bf4a-99ba079940a6
# ╠═30c6a176-935b-423f-9447-86f78746322f
# ╠═00e536dc-34aa-4a1a-93de-4eb3f5e0a348
# ╠═1065b3e0-60bf-497c-b7fb-c5a065737f77
# ╟─a21545da-3b53-47af-b0c4-f253b37dc84f
# ╠═1d22b09e-99c1-4026-9505-07bdffc98582
# ╟─5a210961-19fc-40be-a5f6-033a80f1414d
# ╠═fe704fb4-d07c-4591-b834-d6cf2f4f7075
# ╟─9b57f6ed-02f8-48ba-afa2-0766fe8c0c4c
# ╠═5fed71ec-35fb-4804-99ff-e1eaf18fac1b
# ╟─5eca37ba-f858-45fb-a66a-3795327dfd18
# ╟─a26cf11b-0ce1-4c1d-a64d-1917178ff676
# ╠═188f67d8-2ae8-474c-8e58-68b8b4fde02e
# ╠═d7531d5f-fc2d-42b2-9cf9-6a737b0f0f8d
# ╟─f6f004a6-d71b-4813-a363-9f51dc37e42a
# ╠═0e2d20a1-5f26-4263-9a91-3b40b2c2996a
# ╟─2e11ce81-7d0d-498f-9ddd-7d4d836ab42f
# ╟─b1e062c6-f245-4edc-aa02-871e2c776998
# ╟─c4cc940c-74aa-45f8-a2fa-6016d7c3c145
# ╟─a199d4dd-34f0-40a4-ac7e-2642d12dccc1
# ╠═042a452a-1130-4a56-a1b9-b2674803e445
# ╟─13fc2859-496e-4f6e-8b22-36d9d55768b8
# ╠═32db42f3-5084-4908-9b53-59291b6133c5
# ╟─243d27b5-a1b8-4127-beec-d5643ad07855
# ╟─005289e8-6979-49fe-b20f-66afd207baea
# ╠═0c5ed337-9310-417d-a1f6-7d69dd8c377b
# ╟─65955950-2879-4b8c-bf73-d63e07d2ad96
# ╟─d8f80c62-b2d6-456f-9650-e8102e968673
# ╟─f4b2f509-0769-4df7-956e-e8bfc9ccd89a
# ╟─0bbd9482-d17d-4027-8eec-450807cff792
# ╠═6e3dbf34-c1c9-460b-9546-b2ee8ee99d68
# ╠═d885ac23-ddfa-495c-b93b-54032c8a5c1f
# ╠═85fa321c-dd93-4ac3-96da-1e1bff80970a
# ╟─93428d11-a3dc-4e29-ae6d-48ba37082c74
# ╠═f3279037-01ed-4596-8e5a-86afe4c02c5f
# ╠═2afd54ca-4240-4f07-b38a-242ba0485b45
# ╟─55bd7b9a-a191-4a0b-9c6b-13733be5023e
# ╠═3ceda3b1-bf1c-4126-b94f-2ee03e8dde99
# ╠═800dfed8-9f29-4138-96f8-e8bf1f2f00e6
# ╠═2ee34d76-7238-46c2-94d1-a40d8b017af6
# ╠═49466829-9459-4dc8-85cc-c67460e290d2
# ╟─bb9d283d-cbad-4066-90fe-5b9db0f57db3
# ╠═77f49da5-ffd2-4148-93a6-f45382ba6d91
