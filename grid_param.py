def grid_param(nx,ny):
    """ grid_param(nx,ny) sets the parameters, that contain all necessary grid information. """
    
    #initialize dictionary
    param = dict()
    
    # number of gridpoints in x,y
    param['nx'],param['ny'] = nx,ny
    
    #grid vectors
    param['x'] = np.linspace(0,1,nx)
    param['y'] = np.linspace(0,1,ny)
        
    #number of inner gridpoints
    param['nix'],param['niy'] = nx-2,ny-2
    
    #number of cells in x,y
    param['mx'],param['my'] = nx-1,ny-1
    
    #number of (corner) nodes
    param['nnod'] = nx*ny
    
    #number of triangles
    param['ntri'] = (nx-1)*(ny-1)*2
    
    #number of interior nodes
    param['ninod'] = (nx-2)*(ny-2)
    
    #number of boundary nodes
    param['nbnod'] = param['nnod'] - param['ninod']
    
    #domain edges
    param['x0'],param['y0'] = param['x'][0],param['y'][0]
    param['xend'],param['yend'] = param['x'][-1],param['y'][-1]
    
    #dx,dy
    param['dx'] = param['x'][1] - param['x'][0]
    param['dy'] = param['y'][1] - param['y'][0]
    
    #coordinates for pcolor plot
    param['xpc'] = np.hstack((param['x']-param['dx']/2,param['xend']+param['dx']/2))
    param['ypc'] = np.hstack((param['y']-param['dy']/2,param['yend']+param['dy']/2))
    
    return param
    

    
