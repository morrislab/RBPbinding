import numpy as np
import sys, os
import matplotlib.pyplot as plt



zscores = np.genfromtxt(sys.argv[1], delimiter = ' ', dtype = str)
rows = np.array(sys.argv[2].split(','), dtype = int)
pnames = np.genfromtxt(sys.argv[3], dtype = str)

znames = zscores[:,0]
x = zscores[:,rows[0]].astype(float)
y = zscores[:,rows[1]].astype(float)
significant = zscores[:,rows[2]] == 'True'

color = np.zeros(len(x), dtype = int)
color[significant * (x > 0) * (y > 0)] = 1
color[significant * (x < 0) * (y < 0)] = 2
colors = np.array(['grey', 'coral', 'cornflowerblue'])
col = colors[color]

#for z, zname in enumerate(znames):
    #print pnames[list(pnames[:,0]).index(zname),1], x[z], y[z], significant[z]



fig = plt.figure(figsize = (4.5,4.5), dpi= 200)
ax = fig.add_subplot(111)
ax.plot([-np.amax(np.absolute(x))-.6, np.amax(np.absolute(x))+.6], [0,0], ls = '--', c = 'grey', alpha = 0.5)
ax.plot([0,0], [-np.amax(np.absolute(y))*1.05, 1.05*np.amax(np.absolute(y))],  ls = '--', c = 'grey', alpha = 0.5)
ax.set_xlim([-np.amax(np.absolute(x))*1.05, np.amax(np.absolute(x))*1.05])
ax.set_ylim([-np.amax(np.absolute(y))*1.05, np.amax(np.absolute(y))*1.05])

ax.scatter(x, y, c = col, s = 100)
ax.set_xlabel('Z-score correlation')
ax.set_ylabel('Z-score of Z-score')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.savefig(os.path.splitext(sys.argv[1])[0]+'.jpg', bbox_inches = 'tight', dpi = 300)
addname = ''
if '--legend' in sys.argv:
    addname = '-legend'
    for s, sig in enumerate(significant):
        if x[s] > 0:
            va = 'bottom'
            ha = 'left'
        else:
            va = 'top'
            ha = 'right'
        if sig:
            if '--' in znames[s]:
                ax.text(x[s], y[s], pnames[list(pnames[:,0]).index(znames[s].split('--')[0]), 3]+'-'+pnames[list(pnames[:,0]).index(znames[s].split('--')[1]), 3], va = va, ha = ha)
            else:
                ax.text(x[s], y[s], pnames[list(pnames[:,0]).index(znames[s]), 3], va = va, ha = ha)

if '--outdir' in sys.argv:
    outdir = sys.argv[sys.argv.index('--outdir')+1]
else:
    outdir = os.path.split(sys.argv[1])[0]

if '--savefig' in sys.argv:
    fig.savefig(outdir+os.path.splitext(os.path.split(sys.argv[1])[1])[0]+'row'+str(rows[-1])+'_'+addname+'-zversusz.jpg', dpi = 300, bbox_inches = 'tight')
else:
    plt.show()




'''
significantlist = np.array(significantlist)
sortsig =np.argsort(-significantlist[:, 1].astype(float))
significantlist = significantlist[sortsig]
uniq, unid = np.unique(significantlist[:,0], return_index = True)
significantlist = significantlist[unid]
sortsig = np.argsort(-significantlist[:,1].astype(float))
obj = open(os.path.splitext(sys.argv[1])[0]+'-significant.txt', 'w')
for u in sortsig:
    print ' '.join(significantlist[u])
    obj.write(' '.join(significantlist[u])+'\n')
'''
