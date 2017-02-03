def call_uzawa_alg(param,op,vars,eval):
    """ Outer script to call the algorithm. This script distinguishes between linear and non-linear equations and calculates defects."""
    
    #save per newton iteration: n iterations of uzawa and time
    eval['it'] = np.zeros((param['newton_it_max'],2))
    
    #precompute LU decomposition of mass matrix
    if param['lu_decomposition']:
        op['luM'] = sparse.linalg.factorized(op['M00'])
    
    if param['nonlin'] == 1:    
        for j in range(param['newton_it_max']):
            eval['j'] = j
            
            tic = time.time()
            op = NLmat(param,op,vars)
            eval['it'][j,1] = time.time() - tic
            print('\n'+str(j+1)+'. NL in '+str(time.time() - tic)[:5]+'s.')
            
            op['A'] = sparse.vstack([sparse.hstack([op['S']+op['NL'],op['Z']]),\
            sparse.hstack([op['Z'],op['S']+op['NL']])]).tocsc()
            
            if param['lu_decomposition']:
                # do LU decomposition of stiffness matrix once per Newton step
                tic = time.time()
                op['luA'] = sparse.linalg.factorized(op['S']+op['NL'])
                print('LU decomposition in '+str(time.time() - tic)[:5]+'s.')
                
            f_new = op['A'].dot(vars['v'])+op['Bt'].dot(vars['p'])-vars['f']
            d,e,eval = uzawa_alg(param,op,vars,f_new,eval)
            
            vars['v'] = vars['v'] - d
            vars['p'] = vars['p'] - e
            
            print('defects, e='+format(np.linalg.norm(e),'.2e')+\
            ', d='+format(np.linalg.norm(d),'.2e'))
            
            if (np.linalg.norm(e) < param['newton_tol']) \
            and (np.linalg.norm(d) < param['newton_tol']):
                print('Newton convergence reached.')
                break
            elif j == (param['newton_it_max']-1):
                print('Newton reached maximum iterations.')
        
    else:
        eval['j'] = 0
        v,p,eval = uzawa_alg(param,op,vars,vars['f'],eval)
        vars['v'] = v
        vars['p'] = p
    
    return vars,eval


def uzawa_alg(param,op,vars,f,eval):
    """ Solves the stokes problem via uzawa (optimal) algorithm with
        stiffness matrix A (may include the nonlinear matrix N), 
        which consists of [A 0; 0 A]
        mass matrix M00
        divergence matrix B
        initial guess p for the pressure
        forcing f.
        
        the parameters for convergence are defined in param."""
    
    if param['uz_opt'] == 1:
        
        tic = time.time()
        p = vars['p']
            
        v = sparse.linalg.spsolve(op['A'],f-op['Bt'].dot(p))
        
        for k in range(param['uz_nmax']):
            d = op['B'].dot(v)
            
            if param['lu_decomposition']:
                e = op['luM'](d)
            else:
                e = sparse.linalg.spsolve(op['M00'],d)
                
            b = op['Bt'].dot(e)
            
            if param['lu_decomposition']:
                w1 = op['luA'](-b[:param['ninodbub']])
                w2 = op['luA'](-b[param['ninodbub']:])
                w = np.hstack((w1,w2))
            else:
                w = sparse.linalg.spsolve(op['A'],-b)
            
            alpha = np.dot(d,e) / np.dot(-w,b)
            
            p = p + alpha * e
            v = v + alpha * w
            
            norm1 = alpha*np.linalg.norm(e)
            norm2 = norm1
            
            if k % 100 == 0 and param['uz_out']: #output
                print('it='+str(k)+', alpha='+\
                str(np.round(alpha,decimals=4))[:4]+', pnorm='+format(norm2,'.2e'))
                
            if (norm2 < param['uz_tol']):
                break
        
        eval['it'][eval['j'],0] = k+1
        eval['it'][eval['j'],1] += time.time() - tic
        print('uz opt: conv in '+str(time.time() - tic)[:5]+'s and '+str(k)+' it')
            
    else:
        tic = time.time()
        
        for k in range(param['uz_nmax']):
            v = sparse.linalg.spsolve(A,f-Bt.dot(p))
            d = B.dot(v)
            e = sparse.linalg.spsolve(M00,d)
            p = p + param['uz_alpha'] * e
            norm1 = param['uz_alpha'] * np.linalg.norm(e)
            norm2 = norm1
            
            if k % 100 == 0 and param['uz_out']: #info about convergence
                print([k,norm2])
            if norm2 < param['uz_tol']:
                break
        
        eval['it'][eval['j'],0] = k+1
        eval['it'][eval['j'],1] += time.time() - tic
        print('uz: conv in '+str(time.time() - tic)[:5]+'s and '+str(k)+' it')

    return v,p,eval