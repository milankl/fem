## FEM - create param dictionary

def set_param(nx,ny):
    
    """param = set_param(nx,ny) creates a dictionary, that contains all parameters for given gridsizes nx, ny that triangulate the domain. Hence, parameters that do not need to be recomputed."""
    
    #all grid information
    param = grid_param(nx,ny)
    
    ## PARAMETERS
    
    #bubble function?
    param['bub'] = 1
    
    #number of basic functions
    param['nbf'] = 3+param['bub']
    #number of inner nodes plus bubbles nodes
    param['ninodbub'] = param['ninod'] + param['bub']*param['ntri']
    
    #uzawa parameters
    param['uz_opt'] = 1
    param['uz_alpha'] = .5 #convergence parameter
    param['uz_nmax'] = 5000 #maximum number of iterations
    param['uz_tol'] = 1e-5 #tolerance to reach convergence
    param['uz_out'] = 1 #info about uzawa algorithm
    
    #speed up with LU decomposition?
    param['lu_decomposition'] = 1
    
    #viscosity
    param['nu'] = 1e-1
    
    #nonlinear?
    param['nonlin'] = 1
    
    #newton parameters
    param['newton_tol'] = 1e-3
    param['newton_it_max'] = 15
    
    #plot result?
    param['plot'] = 1
    #output?
    param['save_plot'] = 0
    param['xwindow'] = 1 #in case of ssh and no xwindow deactivate
    
    #plotting every nth vector
    param['plot_nquiver'] = 6
    
    #contourlevels
    param['vcontour'] = np.linspace(-.5,.5,11)
    param['pcontour'] = np.linspace(0,.07,8)
    param['ncontour'] = 10
    
    return param
