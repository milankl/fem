    ## PLOTTING
    
def plot_uzawa(param,vars,eval):
    """ Plots the results and compares them to the analytic solution."""
    
    if param['xwindow']:
    
        x,y,x0,y0,xb,yb = scat_coord(param)
        xx0,yy0 = np.meshgrid(param['x'],param['y'])
        
        n = param['plot_nquiver']
        p = vars['p']
        v1 = vars['v1']
        v2 = vars['v2']
        v1ana = vars['v1ana']
        v2ana = vars['v2ana']
        fx = vars['fn1'] # plot f from normal space
        fy = vars['fn2']
        ni = param['ninod']
        
        import matplotlib.tri as tri
        triang = tri.Triangulation(x, y)
        
        # string for output    
        outstr = str(param['nx'])+'x'+str(param['ny'])+'_%1.0e.pdf' % param['nu']
        
        #--------------------------------
        # compare to analytic solution
        pana = xx0*yy0*(1-xx0)*(1-yy0)
        
        #redefine x,y to account for pcolor shift
        xx0p,yy0p = np.meshgrid(param['xpc'],param['ypc'])
        
        #for quiverkey
        fquimax = float('%1.e' % (fx).max())
        
        fig2, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True,sharey=True)
        Q1 = ax1.quiver(x[:ni:n],y[:ni:n],fx[:ni:n],fy[:ni:n],scale=35,width=0.008)
        ax1.set_title(r'(a) $f$')
        ax1.set_ylim(0,1)
        ax1.set_xlim(0,1)
        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.add_patch(plt.Rectangle([0.75,0.82],0.24,0.17,color='white'))
        ax1.add_patch(plt.Rectangle([0.75,0.82],0.24,0.17,fill=False,color='k'))
        ax1.quiverkey(Q1, 0.87,0.85, fquimax, '%1.e' % fquimax, zorder=5)
        ax1.set_ylabel('y')
        
        impana = ax2.contourf(xx0,yy0,pana,param['pcontour'])
        ax2.set_title(r'(b) $p$ analytic')
        plt.colorbar(impana,ax=(ax1,ax2))
        impana.set_clim(0,.5**4)
        
        imv1a = ax3.tricontourf(triang,v1ana,param['vcontour'])
        ax3.set_title(r'(c) $v_1$ analytic')
        imv1a.set_clim(-.5,.5)
        ax3.set_xlabel('x')
        ax3.set_ylabel('y')
        
        imv2a = ax4.tricontourf(triang,v2ana,param['vcontour'])
        ax4.set_title(r'(d) $v_2$ analytic')
        imv2a.set_clim(-.5,.5)
        plt.colorbar(imv2a,ax=(ax3,ax4))
        ax4.set_xlabel('x')
        
        if param['save_plot']:
            plt.savefig(path+'graphs/analytic_'+outstr)
        
        #-----------------
        #for quiverkey
        vquimax = float('%1.e' % (v1).max())
        ni = param['ninod']
        
        fig1, ((ax11,ax12),(ax13,ax14)) = plt.subplots(2,2,sharex=True,sharey=True)
        Q2 = ax11.quiver(x[:ni:n],y[:ni:n],v1[:ni:n],v2[:ni:n],scale=5,width=0.008)
        ax11.set_title(r'(a) $v$ computed')
        ax11.set_ylim(0,1)
        ax11.set_xlim(0,1)
        ax11.set_xticklabels([])
        ax11.set_yticklabels([])
        ax11.add_patch(plt.Rectangle([0.75,0.82],0.24,0.17,color='white'))
        ax11.add_patch(plt.Rectangle([0.75,0.82],0.24,0.17,fill=False,color='k'))
        ax11.quiverkey(Q2, 0.87,0.85, vquimax, '%1.e' % vquimax, zorder=5)
        ax11.set_ylabel('y')
        
        impcomp = ax12.contourf(xx0,yy0,v2xy(param,p),param['pcontour'],extend='min')
        ax12.set_title(r'(b) $p$ computed')
        c1 = plt.colorbar(impcomp,ax=(ax11,ax12))
        impcomp.set_clim(0,.5**4)
        
        imv1comp = ax13.tricontourf(triang,v1,param['vcontour'],extend='both')
        ax13.set_title(r'(c) $v_1$ computed')
        imv1comp.set_clim(-.5,.5)
        ax13.set_xlabel('x')
        ax13.set_ylabel('y')
        
        imv2comp = ax14.tricontourf(triang,v2,param['vcontour'],extend='both')
        ax14.set_title(r'(d) $v_2$ computed')
        imv2comp.set_clim(-.5,.5)
        c2 = plt.colorbar(imv2comp,ax=(ax13,ax14))
        ax14.set_xlabel('x')
        
        
        if param['save_plot']:
            plt.savefig(path+'graphs/computed_'+outstr)
        
        # ----------------------- 
        er1 = ((v1-v1ana)/v1ana)
        er2 = ((v2-v2ana)/v2ana)
        er3 = ((v2xy(param,p)[1:-1,1:-1]-pana[1:-1,1:-1])/pana[1:-1,1:-1])
        
        er1 = er1.clip(-1,1)
        er2 = er2.clip(-1,1)
        
        fig2, ((ax24,ax23),(ax21,ax22)) = plt.subplots(2,2,sharex=True,sharey=True)
        imer1 = ax21.tricontourf(triang,er1,param['ncontour'],cmap='RdBu_r',extend='both')
        ax21.set_title(r'(c) relative error $v_1$')
        ax21.set_xticklabels([])
        ax21.set_yticklabels([])
        ax21.set_ylim(0,1)
        ax21.set_xlim(0,1)
        imer1.set_clim(-1,1)
        ax21.set_xlabel('x')
        ax21.set_ylabel('y')
        
        imer2 = ax22.tricontourf(triang,er2,param['ncontour'],cmap='RdBu_r',extend='both')
        ax22.set_title(r'(d) relative error $v_2$')
        imer2.set_clim(-1,1)
        ax22.set_xlabel('x')
        
        imer3 = ax23.pcolor(xx0p[1:-1,1:-1],yy0p[1:-1,1:-1],er3,cmap='RdBu_r')
        ax23.set_title(r'(b) relative error $p$')
        plt.colorbar(imer2,ax=(ax24,ax23,ax21,ax22))
        imer3.set_clim(-1,1)
        
        equimax = float('%1.e' % (vars['v1']-vars['v1ana']).max())
        
        #plot inner node velocities
        Q3 = ax24.quiver(x[:ni:n],y[:ni:n],(v1-v1ana)[:ni:n],(v2-v2ana)[:ni:n],color='k',scale=3e-2,width=0.008)            
        ax24.set_title(r'(a) absolute error $v$')
        ax24.add_patch(plt.Rectangle([0.75,0.82],0.24,0.17,color='white'))
        ax24.add_patch(plt.Rectangle([0.75,0.82],0.24,0.17,fill=False,color='k'))
        ax24.quiverkey(Q3, 0.87,0.85, equimax, '%1.e' % equimax, zorder=5)
        ax24.set_ylabel('y')
        
        if param['save_plot']:
            plt.savefig(path+'graphs/error_'+outstr)
            print('All plots printed to file.')
        
        if param['plot']:
            plt.show()