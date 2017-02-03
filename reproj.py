def reproj(param,vars):
    """ vars = reproj(param,op,vars) computes the reprojection matrix R and transforms v,f back into normal space."""
    
    
    tic = time.time()
    if param['bub'] == 1:
        
        R = reprojmat(param)
        #v = (v1,v2)
        vars['v1'],vars['v2'] = R.dot(vars['v'][:param['ninodbub']]),R.dot(vars['v'][-param['ninodbub']:])
        vars['f1'],vars['f2'] = R.dot(vars['f1']),R.dot(vars['f2'])
    
    else:
        vars['v1'],vars['v2'] = vars['v'][:param['ninodbub']],vars['v'][-param['ninodbub']:]

    print('reprojection in '+str(time.time() - tic)[:5]+'s.')
    return vars


def reprojmat(param):
    """ R = reprojmat(param,trinod) produces a matrix that reprojects the discrete variables that are partly defined on bubble nodes into the normal space. This is done via averageing over the neighbouring corner nodes for each bubble node, respectively.
    
    Hence, computing:
    u = sum_i u_i phi_i
    
    Example, for x=1/3,y=1/3 in the standard triangle:
    
    u = u_1*1/3 + u_2*1/3 + u_3*1/3 + u_b*1/27
    
    Hence transforming the coefficients (u_1...u_b) from finite element space back into normal space (u). This is only necessary for the bubble nodes, as the interior nodes are not zero where the bubble function peaks.
    """
    nb = param['ntri']
    ni = param['ninod']
    
    #indices
    c1 = np.ones(3*(nb))/3.
    c2 = ni+np.array([range(nb)]*3,dtype=int).flatten('F')

    #interior nodes onto bubble nodes
    R1 = sparse.coo_matrix((c1,(c2,param['trinod'][:3,:].flatten('F'))),shape=(nb+ni,nb+ni)).tocsc()
    R1 = sparse.hstack((R1[:,:ni],sparse.csc_matrix((ni+nb,nb))))
    
    #scale bubble nodes down as they peak with 1/27
    R2 = 1/27.*sparse.eye(nb+ni,nb+ni).tocsc()
    
    #identity matrix for inner to inner
    R3 = 26/27.*sparse.coo_matrix((np.ones(ni),(np.arange(ni),np.arange(ni))),shape=(nb+ni,nb+ni)).tocsc()
    
    return R1 + R2 + R3