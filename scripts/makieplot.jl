using GridVisualize
using GLMakie
using MultECatJulia

function makieplot(grid)
    gridplot(grid; Plotter=GLMakie, xplanes=[0.0],zplanes=[0.3],outlinealpha=0.3)
end
