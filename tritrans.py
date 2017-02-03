def tritrans(pqr):
    """B = tritrans(pqr)
computes the transformation onto the triangle defined by corners p,q,r as B = det(P) * Pinv * Pinv.T with P being the affin transformation and Pinv its inverse."""
    
    P = np.array([pqr[0,1:]-pqr[0,0],pqr[1,1:]-pqr[1,0]])
    detP = P[0,0]*P[1,1]-P[0,1]*P[1,0]
    #Pinv is without the 1/det -factor
    Pinv = np.array([[P[1,1],-P[0,1]],[-P[1,0],P[0,0]]])
    
    return np.dot(Pinv,Pinv.T) / abs(detP)