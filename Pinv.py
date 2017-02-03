def Pinv(pqr):
    """Inverse of the triangle transformation matrix P."""
    
    P = np.array([pqr[0,1:]-pqr[0,0],pqr[1,1:]-pqr[1,0]])
    detP = P[0,0]*P[1,1]-P[0,1]*P[1,0]
    return np.array([[P[1,1],-P[0,1]],[-P[1,0],P[0,0]]])/detP # =Pinv