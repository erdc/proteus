from matplotlib  import pyplot as  plt

def plot_domain(domain):
    plt.figure()
    plt.rcParams['figure.figsize'] = (10.0, 8.0)
    plt.title('Bathmetry')
    plt.xlabel(r'z[m]')
    plt.ylabel(r'x[m]')
    colors = ['b','g','r','c','m','y','k','w']
    plt.xlim(domain.x[0]-0.1*domain.L[0],domain.x[0]+domain.L[0]+0.1*domain.L[0])
    for si,s in enumerate(domain.segments):
        plt.plot([domain.vertices[s[0]][0],
                  domain.vertices[s[1]][0]],
                 [domain.vertices[s[0]][1],
                  domain.vertices[s[1]][1]],
                 color=colors[domain.segmentFlags[si]-1],
                 linewidth=2,
                 marker='o')
    return plt
