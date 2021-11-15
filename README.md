MultECatJulia
=============

## Working examples
Instantiate this project after cloning or updating via
```
$ julia --project=.
julia> using Pkg
julia> Pkg.instantiate()
```
### Pluto notebooks

Load them as follows:

```
$ julia --project=.
julia> using Pluto
julia> Pluto.run(notebook="notebooks/<Name.jl>")
```


- `notebooks/GridGallery.jl`

### Scripts
Run them as follows:
```
$ julia -project=.
julia> using Revise # if this is not in your startup.jl
julia> includet("scripts/<name.jl>
julia> ?<name>.main # (for help)
julia> <name>.main()
```

- `scripts/gridgallery.jl`



## Project code structure

For the project code structure see  https://j-fu.github.io/marginalia/julia/project-workflow/



