def detP(pqr):
    """ detP = detP(pqr) is the determinant of the transformation matrix P between the standard triangle and an arbitrary triangle defined by coordinates pqr.
    
    pqr is a 2x3 matrix with {x,y}x{p,q,r}."""
    
    P = np.array([pqr[0,1:]-pqr[0,0],pqr[1,1:]-pqr[1,0]])
    return P[0,0]*P[1,1]-P[0,1]*P[1,0]