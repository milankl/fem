## SET UP INITIAL CONDITIONS AND FORCING

def set_vars(param,op):
    
    vars = dict()
    
    #guess for p
    vars['p'] = np.zeros(param['ninod'])
    
    #guess for v
    vars['v'] = np.zeros(param['ninodbub']*2)
    
    #get forcing and analytic solutions for v1,v2
    vars = forcing(param,op,vars)
    
    return vars
    
def forcing(param,op,vars):
    """forcing(param,nodexy,M0b) sets up the forcing f and the analytic solution for v."""
    x = param['nodexy'][0,:param['nnod']]
    y = param['nodexy'][1,:param['nnod']]
    
    if param['bub'] == 1:
        x = np.hstack((x,param['nodexy'][0,-param['ntri']:]))
        y = np.hstack((y,param['nodexy'][1,-param['ntri']:]))
    
    pi = np.pi
    
    f1 = (1-2*x)*y*(1-y) - param['nu']*\
    2*pi**2*(np.cos(pi*x)**2 - 3*np.sin(pi*x)**2)*np.sin(pi*y)*np.cos(pi*y)
    
    f2 = (1-2*y)*x*(1-x) + param['nu']*\
    2*pi**2*(np.cos(pi*y)**2 - 3*np.sin(pi*y)**2)*np.sin(pi*x)*np.cos(pi*x)
    
    if param['nonlin']:
        f1 = f1 + pi*np.sin(pi*x)**3*np.sin(pi*y)**2*np.cos(pi*x)
        f2 = f2 + pi*np.sin(pi*y)**3*np.sin(pi*x)**2*np.cos(pi*y)
    
    vars['fn1'] = f1
    vars['fn2'] = f2
    
    ## analytic solution for v
    
    x = param['nodexy'][0,:param['ninod']]
    y = param['nodexy'][1,:param['ninod']]
    
    if param['bub'] == 1:
        x = np.hstack((x,param['nodexy'][0,-param['ntri']:]))
        y = np.hstack((y,param['nodexy'][1,-param['ntri']:]))
    
    v1ana = np.sin(pi*x)**2 * np.sin(pi*y)*np.cos(pi*y)
    v2ana = -np.sin(pi*y)**2 * np.sin(pi*x)*np.cos(pi*x)
    
    vars['v1ana'] = v1ana
    vars['v2ana'] = v2ana
    
    ## discretize f
    #by projecting the function values of f at positions x,y into FE space with Q and multiplying with the mass matrix 
    
    op = projmat(param,op)
    
    f1 = op['M0b'].dot(op['Q'].dot(f1))
    f2 = op['M0b'].dot(op['Q'].dot(f2))
    
    vars['f1'] = f1
    vars['f2'] = f2    
    vars['f'] = np.hstack((f1,f2))
    
    return vars
    
def projmat(param,op):
    """Q = projmat(param,op) defines the projection matrix, that projects a function f into the finite element space. For all nodes except bubble functions, this is simply the identity matrix. For bubble nodes at position m (= 1/3,1/3 in the standard triangle) we seek the coefficient f_b of the bubble function phi_4 via
    
    f(m) = f_1 * phi_1(m) + f_2 * phi_2(m) + f_3 * phi_3(m) + f_b * phi_4(m)

as phi_1-3(m) = 1/3 and phi_4(m) = 1/27 we get

    f_b = (f(m) - (f_1 + f_2 + f_3)/3)*27.

which is now written as f_i = Qf. f being the function values at all nodes, Q the projection matrix and f_i the coefficients for each basis function."""

    nt = param['ntri']
    ni = param['ninod']
    nb = param['nbnod']
    
    #indices
    c1 = np.ones(3*(nt))*(-9)
    c2 = ni+nb+np.array([range(nt)]*3,dtype=int).flatten('F')

    #interior&boundary nodes onto bubble nodes
    Q1 = sparse.coo_matrix((c1,(c2,param['trinod'][:3,:].flatten('F'))),shape=(nt+ni+nb,nt+ni+nb)).tocsc()
    
    #identity matrix for inner&boundary nodes
    Q2 = sparse.eye(nb+ni+nt,nb+ni+nt).tocsc()
    
    #f(m)*27, 26 from Q3 and 1 from Q2 is 27
    Q3 = 26.*sparse.coo_matrix((np.ones(nt),(np.arange(nt)+ni+nb,np.arange(nt)+ni+nb)),shape=(nb+ni+nt,nb+ni+nt)).tocsc()
    
    op['Q'] = Q1 + Q2 + Q3
    return op
    
    
    
    
    
    