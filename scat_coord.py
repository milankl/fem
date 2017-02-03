def scat_coord(param):
    """ x,y,x0,y0 = scat_coord(param) defines the coordinates used for the scatter plot. x,y are all inner nodes and the bubble nodes. x0,y0 are just inner nodes. In case of no bubble function x,y and x0,y0 are same."""
    
    #coordinates only of inner nodes
    x0 = param['nodexy'][0,:param['ninod']]
    y0 = param['nodexy'][1,:param['ninod']]
    
    #coordinates of inner nodes and bubble nodes
    x = x0.copy()
    y = y0.copy()
    
    if param['bub'] == 1:
        x = np.hstack((x,param['nodexy'][0,-param['ntri']:]))
        y = np.hstack((y,param['nodexy'][1,-param['ntri']:]))
        
        xb = param['nodexy'][0,:]
        yb = param['nodexy'][1,:]
    
    else:
        
        xb = np.hstack((x,param['nodexy'][0,:-param['ntri']]))
        yb = np.hstack((y,param['nodexy'][1,:-param['ntri']]))
    
    
    return x,y,x0,y0,xb,yb