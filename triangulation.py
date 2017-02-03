## FEM - triangulation

def triangulation(param):
    """nodexy,trinod = triangulation(param) 
    divides a defined domain into triangles and returns the following matrices:
    
    remark: nodes and triangles are numbered row first. numbering follows the partitioned case, i.e. the inner nodes are numbered from 0 to ninod, then the boundary nodes from ninod+1 to nnod, e.g.
    
    9   10  11  12  13
    24   0   1   2  14
    23   3   4   5  15
    22   6   7   8  16
    21  20  19  18  17
    
    - nodexy is a 2 by "number of nodes" matrix with x,y coordinates of all nodes
    - trinod is a 3 by "number of triangles" matrix with nodes that belong to a certain triangle, including boundary nodes
    
    this function is based on trinod_calc.py and xy2v2xy.py which are called inside.
    
    In case of the bubble function: the bubble nodes are numbered following the order of the triangles starting in the above example with 25 as all corner nodes are numbered first.
    
    """
    
    xx,yy = np.meshgrid(param['x'],param['y'])
    nodexy = np.array([xy2v(xx),xy2v(yy)])

    trinod = trinod_calc(param,v2xy(param,np.arange(param['nnod'])))
   
    if param['bub'] == 1:
        
        #bubble nodes on odd triangles
        bubox,buboy = np.meshgrid(param['x'][:-1],param['y'][:-1]) + param['dx']/3.
        
        #bubble nodes on even triangles
        bubex,bubey = np.meshgrid(param['x'][:-1],param['y'][:-1]) + 2*param['dx']/3.
        
        #interweave
        bubx = np.vstack((bubox.flatten(),bubex.flatten())).flatten('F')
        buby = np.vstack((buboy.flatten(),bubey.flatten())).flatten('F')
        nodexy = np.hstack((nodexy,np.vstack((bubx,buby))))
        
        trinod = np.vstack((trinod,np.arange(param['nnod'],param['nnod']+param['ntri'])))

    #store nodexy and trinod in param
    param['nodexy'] = nodexy
    param['trinod'] = trinod
    
    return param