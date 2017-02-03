## FEM - stiffness matrix computation

def stiffmass(param,op):
    """ S,M = stiffmass(param,op) calculates the partitioned stiffness and mass matrix. Partitioned means: S00 is stiffness for innernodes, Sbb for boundary nodes, S0b,Sb0 hybrid nodes, such that
    
          ( S00     S0b )
    S =   (             )
          ( Sb0     Sbb )
            
    same for M. In case of bubble nodes, they are stacked such that (u for bubble)
    
        ( S00     S0b     S0u )
        (                     )
    S=  ( Sb0     Sbb     Sbu )
        (                     )
        ( Su0     Sub     Suu )
        
    Also outputs the double stiffness matrix A = (S,0; 0,S) and the mass matrices M0b and M00. All stored in dictionnary op.
    
    """
    
    ## bubble function included?
    if param['bub'] == 0: 
    
        #predefine integral matrices for mass
        M0 = ((np.eye(3)+np.ones(3))/24.).flatten()
        
        #predefine integral matrices for stiffness
        h = 1./2.
        S1 = np.array([[h,-h,0],[-h,h,0],[0,0,0]])
        S2 = np.array([[h,0,-h],[0,0,0],[-h,0,h]])
        S3 = np.array([[2*h,-h,-h],[-h,0,h],[-h,h,0]])
    
    else:
        
        #predefine integral matrix for mass
        M0 = np.ones((4,4))*1./360.
        M0[:3,:3] = (np.eye(3)+np.ones(3))/24.
        M0[3,3] = 1./5040.
        M0 = M0.flatten()
        
        #predefine integral matrices for stiffness
        h = 1./2.
        g = 1./180.
        
        S1 = np.array([[h,-h,0,0],[-h,h,0,0],[0,0,0,0],[0,0,0,g]])
        S2 = np.array([[h,0,-h,0],[0,0,0,0],[-h,0,h,0],[0,0,0,g]])
        S3 = np.array([[2*h,-h,-h,0],[-h,0,h,0],[-h,h,0,0],[0,0,0,g]])
        
    ## -------------------------     
    
    nbf2 = param['nbf']**2 #number of basic functions squared
    
    #preallocate
    sij = np.zeros((3,param['ntri']*nbf2))
    mij = np.zeros((3,param['ntri']*nbf2))
    
    #cloning indices - used to repeat nodes indices
    c = np.array([range(param['nbf'])]*param['nbf'],dtype=int)
    c1,c2 = c.flatten('F'),c.flatten('C')
    
    #loop over triangles
    for i in range(param['ntri']):
        
        B = tritrans(param['nodexy'][:,param['trinod'][:3,i]]) 
        #only the first 3 nodes define the triangle
        
        s = (B[0,0] * S1 + B[1,1] * S2 + B[1,0] * S3).flatten()
        m = abs(detP(param['nodexy'][:,param['trinod'][:3,i]]))*M0
        
        nodes = param['trinod'][:,i]    
        sij[:,i*nbf2:(i+1)*nbf2] = np.array([s,nodes[c1],nodes[c2]])
        mij[:,i*nbf2:(i+1)*nbf2] = np.array([m,nodes[c1],nodes[c2]])
    
    #size of stiff&mass matrix depended of with/without bubble
    sh = param['nnod'] + param['bub']*param['ntri']
    
    M = sparse.coo_matrix((mij[0,:],(mij[1,:],mij[2,:])),shape=(sh,sh)).tocsc()
    S = sparse.coo_matrix((sij[0,:],(sij[1,:],sij[2,:])),shape=(sh,sh)).tocsc()
    
    ## cut off the boundary nodes
    
    ni = param['ninod'] #number of inner nodes
    nb = param['ntri'] #number of bubble functions
    
    M00 = M[:ni,:ni]
    
    if param['bub'] == 1:
        #stiffness
        S00S0u = sparse.hstack((S[:ni,:ni],S[:ni,-nb:]))
        Su0Suu = sparse.hstack((S[-nb:,:ni],S[-nb:,-nb:]))
        S = sparse.vstack((S00S0u,Su0Suu)).tocsc()
        
        #mass incl. boundary nodes (for forcing f)
        M0b = sparse.vstack((M[:ni,:],M[-nb:,:])).tocsc()
        
        #mass
        M00M0u = sparse.hstack((M[:ni,:ni],M[:ni,-nb:]))
        Mu0Muu = sparse.hstack((M[-nb:,:ni],M[-nb:,-nb:]))
        M = sparse.vstack((M00M0u,Mu0Muu)).tocsc()
        
    else:
        #inner nodes
        S = S[:ni,:ni]
        M0b = M[:ni,:]
        M = M[:ni,:ni]
    
    #double stiffness matrix to make it applicable on vectors
    #multioply by viscosity parameter
    Z = sparse.csr_matrix((param['ninodbub'],param['ninodbub']))    
    A = sparse.vstack([sparse.hstack([S,Z]),sparse.hstack([Z,S])]).tocsc()
    
    #store matrices in operator dictionnary op
    op['Z'] = Z
    op['A'] = param['nu']*A
    op['S'] = param['nu']*S
    op['M'] = M
    op['M0b'] = M0b
    op['M00'] = M00
    
    return op
    