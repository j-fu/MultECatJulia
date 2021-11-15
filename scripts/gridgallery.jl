module gridgallery
using GridVisualize
using GLMakie
using MultECatJulia


"""
    main()

Show a gallery of grids
"""
function main()
    vis=GridVisualizer(Plotter=GLMakie,resolution=(1200,400), layout=(1,3))
    gridplot!(vis[1,1],foursquare_electrode_grid(), xplanes=[0.0],zplanes=[0.3],outlinealpha=0.3)
    gridplot!(vis[1,2],circular_electrode_grid(), xplanes=[0.0],zplanes=[0.3],outlinealpha=0.3)
    gridplot!(vis[1,3],circular_anisoref_electrode_grid(), xplanes=[0.0],zplanes=[0.3],outlinealpha=0.3)
    reveal(vis)
end

end
