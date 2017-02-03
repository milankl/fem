def xy2v(A):
    """v = xy2v(A) calculates vector representation of nodes from matrix representation given in A following the partitioned numbering."""
    
    v1 = A[1:-1,1:-1].flatten() #inner
    v2 = A[0,:].flatten() #down
    v3 = A[1:-1,-1].flatten() #right
    v4 = A[-1,-1::-1].flatten() #up
    v5 = A[-2:0:-1,0].flatten()
    
    return np.hstack((v1,v2,v3,v4,v5))
    
def v2xy(param,v):
    """A = v2xy(param,v) calculates matrix representation of nodes from vector representation given in v following the partitioned numbering."""
    
    if len(v) == param['ninod']:
        nix,niy = param['nx']-2,param['ny']-2
        nixy = nix*niy
        
        Ai = v[:nixy].reshape((niy,nix)) #inner
        #repad with zeros
        A = np.pad(Ai,1,'constant',constant_values=0)
    
    elif len(v) == param['nnod']:
        nx,ny = param['nx'],param['ny']
        nix,niy = nx-2,ny-2
        nixy = nix*niy
        
        Ai = v[:nixy].reshape((niy,nix)) #inner
        Ad = v[nixy:nixy+nx] #down
        Ar = v[nixy+nx:nixy+nx+niy].reshape(niy,1) #right
        Au = v[nixy+nx+niy:nixy+2*nx+niy] #up
        Al = v[-niy:].reshape(niy,1) #left
        
        #concatenate
        Alir = np.hstack((Al[-1::-1],Ai,Ar))
        A = np.vstack((Ad,Alir,Au[-1::-1]))
        
    else:
        raise ValueError('Input v is inappropriate in size.')
        
    return A
    
    