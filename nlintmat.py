def nlintmat(param):
    """ Nijk1,Nijk2 = nlintmat() sets up the integral matrices for the nonlinear terms."""
    
    ## Nijk1 - partial x
    
    cl1 = np.array([-1/12.,-1/24.,-1/24.,-1/360.])
    cl2 = (-1)*cl1.copy()
    cl3 = np.zeros(4)
    cl4 = np.array([1/180.,0,1/360.,1/10080.])
    
    T1 = np.vstack((cl1,cl2,cl3,cl4)).T
    
    cl1 = np.array([-1/24.,-1/12.,-1/24.,-1/360.])
    cl2 = (-1)*cl1.copy()
    cl4 = np.array([0,-1/180.,-1/360.,0])
    
    T2 = np.vstack((cl1,cl2,cl3,cl4)).T
    
    cl1 = np.array([-1/24.,-1/24.,-1/12.,-1/360.])
    cl2 = (-1)*cl1.copy()
    cl4 = np.array([1/360.,-1/360.,0,0])
    
    T3 = np.vstack((cl1,cl2,cl3,cl4)).T
    
    cl1 = np.array([-1/360.,-1/360.,-1/360.,-1/5040.])
    cl2 = (-1)*cl1.copy()
    cl4 = np.array([1/10080.,0,0,0])
    
    T4 = np.vstack((cl1,cl2,cl3,cl4)).T
    
    Nijk1 = np.zeros((4,4,4))
    
    Nijk1[:,:,0] = T1
    Nijk1[:,:,1] = T2
    Nijk1[:,:,2] = T3
    Nijk1[:,:,3] = T4
    
    
    ## Nijk2 - partial y
    
    cl1 = np.array([-1/12.,-1/24.,-1/24.,-1/360.])
    cl2 = np.zeros(4)
    cl3 = (-1)*cl1.copy()
    cl4 = np.array([1/180.,1/360.,0,1/10080.])
    
    T1 = np.vstack((cl1,cl2,cl3,cl4)).T
    
    cl1 = np.array([-1/24.,-1/12.,-1/24.,-1/360.])
    cl3 = (-1)*cl1.copy()
    cl4 = np.array([1/360.,0,-1/360.,-1/10080.])
    
    T2 = np.vstack((cl1,cl2,cl3,cl4)).T
    
    cl1 = np.array([-1/24.,-1/24.,-1/12.,-1/360.])
    cl3 = (-1)*cl1.copy()
    cl4 = np.array([0,-1/360.,-1/180.,-1/10080.])
    
    T3 = np.vstack((cl1,cl2,cl3,cl4)).T
    
    cl1 = np.array([-1/360.,-1/360.,-1/360.,-1/5040.])
    cl3 = (-1)*cl1.copy()
    cl4 = np.array([1/10080.,-1/10080.,-1/10080.,0])
    
    T4 = np.vstack((cl1,cl2,cl3,cl4)).T
    
    Nijk2 = np.zeros((4,4,4))
    
    Nijk2[:,:,0] = T1
    Nijk2[:,:,1] = T2
    Nijk2[:,:,2] = T3
    Nijk2[:,:,3] = T4
    
    if param['bub'] == 0:
        Nijk1 = Nijk1[:3,:3,:3]
        Nijk2 = Nijk2[:3,:3,:3]
        
    #store Nijk1,Nijk2 in op dictionnary
    op = dict()
    op['Nijk1'] = Nijk1
    op['Nijk2'] = Nijk2
    
    return op
    