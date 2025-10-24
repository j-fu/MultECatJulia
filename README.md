MultECatJulia
=============

## Working examples

Install this project after cloning or updating via
```
$ julia --project=.
julia> using Pkg
julia> Pkg.resolvex()
```

### Pluto notebooks



- `notebooks/GridGallery.jl` : some example grids
- `notebooks/Equilibrium1D.jl`: 1D double layer capacitance calculations based on DGLM model
- `src/equilibrium.jl` Implementation of DGLM modified Poisson-Boltzmann as part of the project-package
- `notebooks/Equilibrium2D.jl`: 1D/2D double layer capacitance calculations for polycrystalline electrodes
 
Load them as follows:

```
$ julia --project=.
julia> using Pluto
julia> Pluto.run(notebook="notebooks/<Name.jl>")
```


### Scripts

- `scripts/gridgallery.jl`: some example grids

Run them as follows:
```
$ julia -project=.
julia> using Revise
julia> includet("scripts/<name.jl>
julia> ?<name>.main # (for help)
julia> <name>.main()
```


## Project code structure

For the project code structure see  https://j-fu.github.io/marginalia/julia/project-workflow/



