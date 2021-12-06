"""
    foursquare_electrode_grid(;xymax=1.0,zmax=1.0,hxymin=0.1,hxymax=0.1,hzmin=0.05,hzmax=0.2)

Tensor product based simplex grid in [-xymax,xymax]x [-xymax,xymax] x [0,zmax] with electrode 
of 4 different materials in foursquare configuration at base.
"""
function foursquare_electrode_grid(;nref=0,xymax=1.0,zmax=1.0,hxymin=0.1,hxymax=0.1,hzmin=0.05,hzmax=0.2)

    qref=2.0^nref
    Xp=geomspace(0,xymax,hxymin/qref,hxymax/qref)
    Xm=-reverse(Xp)
    X=glue(Xm,Xp)
    gridxy=simplexgrid(X,X)
    fill!(gridxy[BFaceRegions],5)
    gridxy
    


    Z=geomspace(0,zmax,hzmin/qref,hzmax/qref)
    gridxyz=simplexgrid(gridxy,Z,top_offset=1)

    bfacemask!(gridxyz,[-xymax,-xymax,0],[0,0,0],1)
    bfacemask!(gridxyz,[0,0,0],[xymax,xymax,0],2)
    bfacemask!(gridxyz,[0,-xymax,0],[xymax,0,0],3)
    bfacemask!(gridxyz,[-xymax,0,0],[0,xymax,0],4)
    gridxyz
end




"""
    circular_electrode_grid(;nref=0,xymax=1.0, maxarea=0.01,rad=0.5, zmax=1.0,hzmin=0.05,hzmax=0.2)
Prism based simplexgrid grid in   [-xymax,xymax]x [-xymax,xymax] x [0,zmax] 
with circular electrode 
"""
function circular_electrode_grid(;nref=0,xymax=1.0, maxarea=0.01,rad=0.5, zmax=1.0,hzmin=0.05,hzmax=0.2)
    qref=2.0^nref
    builder=SimplexGridBuilder(Generator=Triangulate)
    
    
    regionpoint!(builder,(xymax-0.1,xymax-0.1))
    cellregion!(builder,1)
    rect2d!(builder,[-xymax,-xymax],[xymax,xymax],facetregions=1)
    
    
    facetregion!(builder,3)
    cellregion!(builder,2)
    regionpoint!(builder, (0,0))
    circle!(builder,(0,0), rad,n=20)
    gridxy=simplexgrid(builder,maxvolume=maxarea/qref^2)
    

    Z=geomspace(0,zmax,hzmin/qref,hzmax/qref)
    gridxyz=simplexgrid(gridxy,Z,bot_offset=0)
    bfacemask!(gridxyz,[-xymax,-xymax,zmax],[xymax,xymax,zmax],5,allow_new=false)
    gridxyz
    
    
end



"""
    circular_electrode_grid(;nref=0,xymax=1.0, maxarea=0.01,rad=0.5, zmax=1.0,hzmin=0.05,hzmax=0.2)
Prism based simplexgrid grid in   [-xymax,xymax]x [-xymax,xymax] x [0,zmax] 
with circular electrode 
"""
function circular_anisoref_electrode_grid(;nref=0,nang=40,nxy=15, r=0.5, dr=0.1,hrmin=0.01,hrmax=0.05,xymax=1.0, zmax=1.0,hzmin=0.05,hzmax=0.2)
    qref=2.0^nref
    
    rad1=geomspace(r-dr,r,hrmax/qref,hrmin/qref)
    rad2=geomspace(r,r+dr,hrmin/qref,hrmax/qref)
    ang=range(0,2π,length=Int(nang*qref))
    δr=2π*r/(nang*qref)
    maxvol=0.5*δr^2
    ring1=ringsector(rad1,ang)
    ring1[CellRegions].=1
    ring2=ringsector(rad2,ang)
    ring1[CellRegions].=2
    ring=glue(ring1,ring2,breg=3)
    
    binner=SimplexGridBuilder(Generator=Triangulate)
    regionpoint!(binner,(0,0))
    cellregion!(binner,2)
    
    facetregion!(binner,2)
    bregions!(binner,ring,[1])
    ginner=simplexgrid(binner,maxvolume=maxvol,nosteiner=true)
    ginner[CellRegions].=2
    
    bouter=SimplexGridBuilder(Generator=Triangulate)
    
    holepoint!(bouter,(0,0))
    
    facetregion!(bouter,3)
    bregions!(bouter,ring,[2])
    
    regionpoint!(bouter,(xymax-0.1,xymax-0.1))
    cellregion!(bouter,3)
    rect2d!(bouter,[-xymax,-xymax],[xymax,xymax],facetregions=1,nx=Int(nxy*qref),ny=Int(nxy*qref))
    
    gouter=simplexgrid(bouter;maxvolume=maxvol*2,nosteiner=true)
    
    gridxy=glue(glue(gouter,ring),ginner)


    Z=geomspace(0,zmax,hzmin/qref,hzmax/qref)
    gridxyz=simplexgrid(gridxy,Z,bot_offset=0)
    bfacemask!(gridxyz,[-xymax,-xymax,zmax],[xymax,xymax,zmax],5,allow_new=false)
    gridxyz

end



"""
	polycrystal_grid2d(;W=50*nm,
		                H=10*nm, 
                        hgmin=0.5*nm, 
                        hzmin=0.1*nm,
                        hzmax=2*nm,
                        ngrain=5)

Create 2D grid for polycrystal surface.
Keyword arguments:
- W: Domain width (electrode width) 
- H: Domain height
- hzmin: grid z resolution at electrode
- hzmax: grid z resolution in bulk
- hgmin: grid x resolution at grain boundary
- ngrain: number of grains
"""
function polycrystal_grid2d(;W=50*nm,H=10*nm, hgmin=0.5*nm, hzmin=0.1*nm,hzmax=2*nm,ngrain=5)
    if ngrain==1
	Xgrain=range(0,W,length=5)
	Wgrain=W
    else
        Wgrain=W/(ngrain-1)
	Xgrain0=geomspace(0,Wgrain/2,Wgrain/5,hgmin)
	Xgrain1=glue(Xgrain0,Wgrain.-reverse(Xgrain0))
	Xgrain=Xgrain1
	for i=1:ngrain-2
	    Xgrain=glue(Xgrain,i*Wgrain.+Xgrain1)
	end
    end
    Z=geomspace(0,H,hzmin,hzmax)
    grid=simplexgrid(Xgrain,Z)
    bfacemask!(grid,[0,0],[W,H],allow_new=false,ngrain+1)
    if ngrain==1
	bfacemask!(grid,[0,0],[W,0],allow_new=false,1)
    else
	for igrain=0:ngrain-1
	    bfacemask!(grid,[igrain*Wgrain-Wgrain/2,0],[igrain*Wgrain+Wgrain/2,0],allow_new=false,igrain+1)
	end
    end
    grid
end

