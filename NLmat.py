def NLmat(param,op,vars):
    """ computes the nonlinear matrix with given velocity v. This is similar to the set up of the other operator matrices but involves another sum (to transform the 3rd order tensor to a matrix), done here with the einsum function for readability. """
    
    nbf = param['nbf']
    ne = nbf**2
    nn = param['nnod'] + param['ntri']
    nb = param['ntri']
    ni = param['ninod']
    
    #split v in components v1, v2
    v1 = vars['v'][:param['ninodbub']]
    v2 = vars['v'][-param['ninodbub']:]
    
    #reinsert boundary velocities for correct numbering
    v1 = np.hstack((v1[:param['ninod']],np.zeros(param['nbnod']),v1[param['ninod']:]))
    v2 = np.hstack((v2[:param['ninod']],np.zeros(param['nbnod']),v2[param['ninod']:]))
    
    #cloning indices
    cs = np.array([range(nbf)]*nbf,dtype=int).flatten('C')
    ct = np.array([range(nbf)]*nbf,dtype=int).flatten('F')
    
    #preallocate
    Nij = np.zeros((3,param['ntri']*nbf**2))
    
    for i in range(param['ntri']):
        
        pinv = Pinv(param['nodexy'][:,param['trinod'][:3,i]])
        p = abs(detP(param['nodexy'][:,param['trinod'][:3,i]]))
        N1 = (pinv[0,0]*op['Nijk1'] + pinv[1,0]*op['Nijk2'])
        N2 = (pinv[0,1]*op['Nijk1'] + pinv[1,1]*op['Nijk2'])
        
        nodes = param['trinod'][:,i]
        N1 = p*np.einsum('kij,k',N1,v1[nodes]).flatten('F')
        N2 = p*np.einsum('kij,k',N2,v2[nodes]).flatten('F')
    
        Nij[:,i*ne:(i+1)*ne] = np.array([N1+N2,nodes[ct],nodes[cs]])
    
    NL = sparse.coo_matrix((Nij[0,:],(Nij[1,:],Nij[2,:])),shape=(nn,nn)).tocsc()
    
    if param['bub'] == 1:
        N00N0u = sparse.hstack((NL[:ni,:ni],NL[:ni,-nb:]))
        Nu0Nuu = sparse.hstack((NL[-nb:,:ni],NL[-nb:,-nb:]))
        NL = sparse.vstack((N00N0u,Nu0Nuu)).tocsc()
        
    else:
        NL = NL[:ni,:ni]  #inner nodes
    
    #store
    op['NL'] = NL
    
    return op