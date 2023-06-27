import numpy as np
import sys, os
import matplotlib.pyplot as plt


# vertical line at 0.4, split kh and rrm, /% for 0.4


confcut = 0.127
confcut2 = 0.2

latstats = np.genfromtxt(sys.argv[1], dtype = str)
protnames = np.genfromtxt(sys.argv[2], dtype = str, delimiter = ' ', comments = None)
speciclade = np.genfromtxt(sys.argv[3], dtype = str)

for p, prot in enumerate(protnames):
    protnames[p, -1] = prot[-1].split('-')[0]

domains = ['RRM', 'KH']


reccolors = [(.9,0.6,0.3,.9), (.3,0.7,0.3,.9),(.6,0.2,0.2,.9), (.4,0.1,0.4,.9)]

xspeci = [[],[]]
kings = np.unique(speciclade[:,-1])
kings = ['Metazoa', 'Plants','Fungi', 'Protists']
for k in kings:
    for d, dom in enumerate(domains):
        speckings = speciclade[np.isin(speciclade[:,-1],k),0]
        protking = protnames[np.isin(protnames[:,-2],speckings),0]
        protdom = protnames[protnames[:,-1] == dom,0]
        protking = np.intersect1d(protking, protdom)
        xspeci[d].append(latstats[np.isin(latstats[:,0],protking),-1].astype(float))


x = [np.sort(np.concatenate(xspeci[0])), np.sort(np.concatenate(xspeci[1]))]
cumdist = [np.arange(1,len(x[0])+1)/float(len(x[0])), np.arange(1,len(x[1])+1)/float(len(x[1]))]

confy = [cumdist[0][x[0]>=confcut][0], cumdist[1][x[1]>=confcut][0]]
confy2 = [cumdist[0][x[0]>=confcut2][0], cumdist[1][x[1]>=confcut2][0]]


fig2 = plt.figure(figsize = (8,4))
conys = []
for d, xid in enumerate(xspeci):
    ax2 = fig2.add_subplot(1,2,int(d+1))
    for s, xi in enumerate(xid):
        xi = np.sort(xi)
        ax2.set_title(domains[d])
        cumdisti = np.arange(1,len(xi)+1)/float(len(xi))
        confi = cumdisti[xi>=confcut][0]*100.
        confi2 = cumdisti[xi>=confcut2][0]*100.
        ax2.plot(xi, cumdisti, c = reccolors[s], label = kings[s]+' '+str(int(confi))+'/'+str(int(confi2))+'%')
        ax2.fill_between(xi, 0, cumdisti, color = reccolors[s] , alpha = 0.4)
        conys.append(confi)


    ax2.plot(x[d], cumdist[d], c = 'k', label = 'Total '+str(int(confy[d]*100.))+'/'+str(int(confy2[d]*100.))+'%', ls = '--')

    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xlabel('JPLE e-dist to measured RBP')
    ax2.set_ylabel('Percentage Eukaryotic RBPs')
    ax2.set_yticks(np.linspace(0,1,11))
    ax2.set_yticklabels([0,10,20,30,40,50,60,70,80,90,100])
    ax2.set_xticks(np.linspace(0,1,11))
    ax2.set_xticks([confcut], minor = True)
    ax2.grid()
    ax2.set_xlim([0,np.amax(x[d])])
    ax2.set_ylim([0,1])
    ax2.legend()

    ax2.plot([confcut, confcut],[0, 1], c= 'r', ls = '--' )
    ax2.plot([confcut2, confcut2],[0, 1], c= 'r', ls = ':' )
    #ax2.text(confcut, confy+0.01, str(int(confy*100.))+'%', va = 'bottom', ha = 'right', color = 'red')

if '--savefig' in sys.argv:
    print 'Saved as', os.path.splitext(sys.argv[1])[0]+'-cummulativeclade.svg'
    fig2.savefig(os.path.splitext(sys.argv[1])[0]+'-cummulativeclade.svg', dpi = 300, bbox_inches = 'tight')
else:
    fig2.tight_layout()
    plt.show()



