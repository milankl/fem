## ERROR CALCULATION

def calc_error(param,vars,eval):
    """ Calculates relative and absolute error and stores them in eval."""

    xx0,yy0 = np.meshgrid(param['x'],param['y'])
    
    p = vars['p']
    v1 = vars['v1']
    v2 = vars['v2']
    v1ana = vars['v1ana']
    v2ana = vars['v2ana']
    pana = xx0*yy0*(1-xx0)*(1-yy0)
    
    #relative error, cut off boundary for p, as p_analytic=0 there.
    er1 = ((v1-v1ana)/v1ana)
    er2 = ((v2-v2ana)/v2ana)
    er3 = ((v2xy(param,p)[1:-1,1:-1]-pana[1:-1,1:-1])/pana[1:-1,1:-1])
    
    #absolute error
    er1abs = (v1-v1ana)
    er2abs = (v2-v2ana)
    er3abs = v2xy(param,p)-pana
    
    #root mean square error
    eval['er_v'] = np.sqrt(((er1abs**2) + (er2abs**2)).mean()) / np.sqrt((v1ana**2 + v2ana**2).mean())
    eval['er_p'] = np.sqrt((er3abs**2).mean()) / np.sqrt((pana**2).mean())
    
    eval['er_v_abs'] = np.sqrt(((er1abs**2) + (er2abs**2)).mean())
    eval['er_p_abs'] = np.sqrt((er3abs**2).mean())
    
    eval['er_v_norm'] = np.linalg.norm(np.hstack((er1abs,er2abs))) / np.linalg.norm(np.hstack((v1ana,v2ana)))
    
    eval['er_p_norm'] = np.linalg.norm(er3abs) / np.linalg.norm(pana)
    
    # separate for bubble and interior nodes
    
    if param['bub']:
        eval['er_v_abs_inod'] = np.sqrt(((er1abs[:param['ninod']]**2) +\
        (er2abs[:param['ninod']]**2)).mean())
        
        eval['er_v_abs_bub'] = np.sqrt(((er1abs[-param['ntri']:]**2) +\
        (er2abs[-param['ntri']:]**2)).mean())
        
        eval['er_v_bub'] = eval['er_v_abs_bub'] / \
        np.sqrt((v1ana[-param['ntri']:]**2 + v2ana[-param['ntri']:]**2).mean())
        
        eval['er_v_inod'] = eval['er_v_abs_inod'] / \
        np.sqrt((v1ana[:param['ninod']]**2 + v2ana[:param['ninod']]**2).mean())
    
    print('RMSE relative in v = %2.4f' % eval['er_v'])
    print('RMSE abs in v = %2.4f' % eval['er_v_abs'])
    print('RMSE relative in p = %2.4f' % eval['er_p'])
    print('RMSE abs in p = %2.4f' % eval['er_p_abs'])

    if param['bub']:
        print('RMSE relative in v_bub = %2.4f' % eval['er_v_bub'])
        print('RMSE relative in v_inod = %2.4f' % eval['er_v_inod'])
    
    return eval