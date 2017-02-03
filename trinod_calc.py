def trinod_calc(param,allnodes):
        """this function calculates the associated nodes for all triangles. allnodes as input is a matrix of all nodes (in the shape of the domain) with arbitrary numbering. Main idea is to cut off two of the four framing rows/columns which remainings are the associated nodes. Then rearranging to follow the numbering of the triangles."""
        
        #cut off 2 sides to obtain lower left (ll), lower right (lr), upper left (ul), and upper right (ur) matrix
        ll = allnodes[:-1,:-1].flatten()
        lr = allnodes[:-1,1:].flatten()
        ul = allnodes[1:,:-1].flatten()
        ur = allnodes[1:,1:].flatten()

        #rearrange for even and odd triangles
        oddtris = np.array([ll,lr,ul])
        eventris = np.array([ur,ul,lr])
        
        #vstack in order to interweave odd and even triangles
        return np.reshape(np.vstack((oddtris,eventris)),(3,param['ntri']),'F')