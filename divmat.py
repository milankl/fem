## divmat

def divmat(param,op):
    """ B = divmat(param,op) calculates the partitioned divergence matrix. Partitioned means: B00 is stiffness for innernodes, Bbb for boundary nodes, Buu for bubble nodes, and corresponding hybrid nodes, such that
    
          ( B00     B0b )
    B =   (             )
          ( Bb0     Bbb )
        
    without bubble function, and
    
        ( B00     B0b     B0u )
    B=  (                     )
        ( Bb0     Bbb     Bbu )
            
    Please note, that these matrices actually appear twice concatenated horizontally to each other as two velocity components called 1 and 2 need to be computed. This corresponds to d/dx for 1 and d/dy for 2."""
    
    ## bubble function included? change integral matrices accordingly.
    
    if param['bub'] == 0:
    
        #predefine integral matrices for divergence
        Bij1 = (np.array([[1,-1,0],]*3) /6.).flatten()
        Bij2 = (np.array([[1,0,-1],]*3) /6.).flatten()
    
    else:
        
        #predefine integral matrices for divergence
        Bij1 = (np.array([[1,-1,0,-.05],[1,-1,0,.05],[1,-1,0,0]]) /6.).flatten()
        Bij2 = (np.array([[1,0,-1,-.05],[1,0,-1,0],[1,0,-1,.05]]) /6.).flatten()
        
    
    ## ----------------
    
    nbf = param['nbf'] #number of basic functions
    ne = 3*nbf*2 #number of entries in the integral matrices
    
    #preallocate
    bij = np.zeros((3,param['ntri']*ne))
    
    #cloning indices - used to repeat nodes indices
    c1 = np.array([range(3)]*nbf,dtype=int).flatten('F')
    c2 = np.array([range(nbf)]*3,dtype=int).flatten('C')
    
    #node offset
    nn = param['nnod'] #number of (corner)nodes
    na = param['nnod'] + param['bub']*param['ntri'] #all nodes
    
    #loop over triangles
    for i in range(param['ntri']):
        
        pinv = Pinv(param['nodexy'][:,param['trinod'][:3,i]])
        p = abs(detP(param['nodexy'][:,param['trinod'][:3,i]]))
        b1 = p*(pinv[0,0]*Bij1 + pinv[1,0]*Bij2)
        b2 = p*(pinv[0,1]*Bij1 + pinv[1,1]*Bij2)
        
        nodes = param['trinod'][:,i]    
        bij[:,i*ne:(i+1)*ne-ne/2] = np.array([b1,nodes[c1],nodes[c2]])
        bij[:,i*ne+ne/2:(i+1)*ne] = np.array([b2,nodes[c1],na+nodes[c2]])
    
    B = sparse.coo_matrix((bij[0,:],(bij[1,:],bij[2,:])),shape=(nn,2*na)).tocsc()
    
    ## cut off boundary nodes
    #cut off the Bub, Bbu, B0b and Bb0 for 1 and 2 from the example above
    
    ni = param['ninod'] #number of inner nodes
    nb = param['ntri'] #number of bubble nodes
    
    if param['bub'] == 0:
        B = sparse.hstack([B[:ni,:ni],B[:ni,nn:nn+ni]]).tocsr()
        
    else:
        B100 = B[:ni,:ni]
        B10uB200 = B[:ni,nn:nn+nb+ni]
        B20u = B[:ni,na+nn:]
        B = sparse.hstack((B100,B10uB200,B20u))
    
    #store in op
    op['B'] = B
    op['Bt'] = B.T
    
    return op