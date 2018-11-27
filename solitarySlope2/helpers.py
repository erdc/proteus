def CreateFig():
    from tables import  openFile
    archive = openFile('nonlinear_waves.h5','r')
    import nonlinear_waves as tank
    import nonlinear_waves_so as tank_so
    import matplotlib.tri as mtri
    from matplotlib import pyplot as  plt
    import numpy as np
    domain = tank.domain
    nodes = archive.getNode("/nodesSpatial_Domain0")
    x=nodes[:,0]
    y=nodes[:,1]
    elements = archive.getNode("/elementsSpatial_Domain0")
    triang = mtri.Triangulation(x, y, elements)
    domain.L=tank.tank_dim
    domain.x=[0.,0.,0.]
    xg = np.linspace(0, domain.L[0], 20)
    yg = np.linspace(0, domain.L[1], 20)
    xi, yi = np.meshgrid(xg,yg)
    plt.figure()
    for it,t in enumerate(tank_so.tnList[:]):
        phi = archive.getNode("/phi_t"+`it`)
        vof = archive.getNode("/vof_t"+`it`)
        wvof = np.ones(vof.shape,'d')
        wvof -= vof
        u = archive.getNode("/u_t"+`it`)
        v = archive.getNode("/v_t"+`it`)
        plt.clf()
        plt.xlabel(r'z[m]')
        plt.ylabel(r'x[m]')
        colors = ['b','g','r','c','m','y','k','w']
        #plt.xlim(domain.x[0]-0.1*domain.L[0],domain.x[0]+domain.L[0]+0.1*domain.L[0])    
        for si,s in enumerate(domain.segments):
            plt.plot([domain.vertices[s[0]][0],
                     domain.vertices[s[1]][0]],
                    [domain.vertices[s[0]][1],
                     domain.vertices[s[1]][1]],
                    color=colors[domain.segmentFlags[si]-1],
                    linewidth=2,
                    marker='o')
        plt.tricontourf(x,y,elements,wvof*np.sqrt(u[:]**2 + v[:]**2))
        plt.tricontour(x,y,elements,phi,[0], linewidth=4)
        u_interp_lin = mtri.LinearTriInterpolator(triang, u[:])
        v_interp_lin = mtri.LinearTriInterpolator(triang, v[:])
        u_lin = u_interp_lin(xi, yi)
        v_lin = v_interp_lin(xi, yi)
        plt.streamplot(xg, yg, u_lin, v_lin,color='k')
        plt.title('T=%2.2f' % (t,))
        plt.xlim((0,domain.L[0]))
        plt.ylim(0,domain.L[1])
        plt.savefig('phi%4.4d.png' % (it,))