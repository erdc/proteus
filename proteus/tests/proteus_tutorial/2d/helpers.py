import matplotlib

def CreateFig(dt_output=0.1):
    from tables import  open_file
    archive = open_file('dambreak.h5','r')
    import dambreak
    dambreak.outputStepping.dt_output=dt_output
    dambreak.outputStepping.nDTout=None
    dambreak.outputStepping.setOutputStepping()
    dambreak.myTpFlowProblem.initializeAll()
    import matplotlib.tri as mtri
    from matplotlib import pyplot as  plt
    import numpy as np
    domain = dambreak.domain
    domain.L = dambreak.tank_dim
    domain.x = (0.,0.,0.)
    nodes = archive.get_node("/nodesSpatial_Domain0")
    x=nodes[:,0]
    y=nodes[:,1]
    elements = archive.get_node("/elementsSpatial_Domain0")
    triang = mtri.Triangulation(x, y, elements)
    xg = np.linspace(0, domain.L[0], 20)
    yg = np.linspace(0, domain.L[1], 20)
    xi, yi = np.meshgrid(xg,yg)
    plt.figure()
    for it,t in enumerate(dambreak.myTpFlowProblem.so.tnList[:]):
        phi = archive.get_node("/phi_t{0:d}".format(it))
        vof = archive.get_node("/vof_t{0:d}".format(it))
        wvof = np.ones(vof.shape,'d')
        wvof -= vof
        u = archive.get_node("/u_t{0:d}".format(it))
        v = archive.get_node("/v_t{0:d}".format(it))
        plt.clf()
        plt.xlabel(r'z[m]')
        plt.ylabel(r'x[m]')
        colors = ['w','b', 'g','r','c','m','y','k']*(max(domain.segmentFlags)//8 + 1)
        plt.xlim(domain.x[0]-0.1*domain.L[0],domain.x[0]+domain.L[0]+0.1*domain.L[0])    
        for si,s in enumerate(domain.segments):
            plt.plot([domain.vertices[s[0]][0],
                     domain.vertices[s[1]][0]],
                    [domain.vertices[s[0]][1],
                     domain.vertices[s[1]][1]],
                    color=colors[domain.segmentFlags[si]-1],
                    linewidth=2,
                    marker='o')
        plt.tricontourf(x,y,elements,wvof*np.sqrt(u[:]**2 + v[:]**2))
        plt.tricontour(x,y,elements,phi,[0], linewidths=4)
        u_interp_lin = mtri.LinearTriInterpolator(triang, u[:])
        v_interp_lin = mtri.LinearTriInterpolator(triang, v[:])
        u_lin = u_interp_lin(xi, yi)
        v_lin = v_interp_lin(xi, yi)
        plt.streamplot(xg, yg, u_lin, v_lin,color='k')
        plt.title('T={0:2.2f}'.format(t))
        plt.axis('equal')
        plt.xlim((0,domain.L[0]))
        plt.savefig('phi{0:04d}.png'.format(it))
