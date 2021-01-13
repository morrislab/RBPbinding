import numpy as np
import sys, os
import seaborn as sns
import matplotlib.pyplot as plt

numsets = int(sys.argv[1])

x = []
y = []
pairname = []
add = ''    
for n in range(numsets):
    pearson = np.genfromtxt(sys.argv[2+2*n])
    pnames = open(sys.argv[2+2*n], 'r').readline().strip().split()[1:]

    if os.path.splitext(sys.argv[3+2*n]) == '.npz':
        ident = np.load(sys.argv[3+2*n])['identmat']
        inames = np.load(sys.argv[3+2*n])['names']
    else:
        ident = np.genfromtxt(sys.argv[3+2*n])
        inames = open(sys.argv[3+2*n], 'r').readline().strip().split()[1:]

    for i, iname in enumerate(inames):
        inames[i] = iname.split('||')[1]


    if '--protlist' in sys.argv:
        add += '-set'+os.path.splitext(os.path.split(sys.argv[sys.argv.index('--protlist')+1+n])[1])[0]
        plist = np.genfromtxt(sys.argv[sys.argv.index('--protlist')+1+n], dtype = str)
        if len(np.shape(plist)) == 2:
            plist = plist[:, 1]
    else:
        plist = np.intersect1d(inames, pnames)
        
    for p, prot1 in enumerate(plist):
        for q in range(p + 1, len(plist)):
            xp = inames.index(prot1)
            xq = inames.index(plist[q])
            yp = pnames.index(prot1)
            yq = pnames.index(plist[q])
            x.append(ident[xp, xq])
            y.append(pearson[yp,yq])
            pairname.append([prot1, plist[q]])

x = np.array(x)
y = np.array(y)




if '--printpair' in sys.argv:
    pearscut = float(sys.argv[sys.argv.index('--printpair')+1])
    idcut = float(sys.argv[sys.argv.index('--printpair')+2])
    pairids = np.where((x >= idcut) & (y <= pearscut))[0]
    for pairid in pairids:
        print np.sort(pairname[pairid]), x[pairid], y[pairid]


fig = plt.figure(figsize = (5.5,4.5), dpi = 300)
ax = fig.add_subplot(121)
if '--selfcorr' in sys.argv:
    ax.set_position([0.15,0.2,0.7,0.7])
    add +='-selfAB'
    scorrfile = sys.argv[sys.argv.index('--selfcorr')+1]
    scrow = int(sys.argv[sys.argv.index('--selfcorr')+2])
    selfcorrs = np.genfromtxt(scorrfile)[:, scrow]
    axself = fig.add_subplot(122)
    axself.set_position([0.875, 0.2, 0.075, 0.7])
    
    axself.spines['top'].set_visible(False)
    axself.spines['right'].set_visible(False)
    axself.spines['left'].set_visible(False)
    axself.tick_params(left = False, labelleft = False)
    bplot = axself.boxplot(selfcorrs, patch_artist = True)
    bplot['boxes'][0].set_facecolor('grey')
    axself.set_xticks([1])
    axself.set_xticklabels(['Control'], rotation = 60)
    axself.set_xlim([0.85, 1.15])
    
    patch_artist=True
else:
    ax.set_position([0.15,0.1,0.8,0.8])

boxcuts = [0,10,20,30,40,50,60,70,80,90,100]

if '--boxplot' in sys.argv:
    boxdata = []
    for i in range(len(boxcuts)-1):
        boxdata.append(y[(x>=boxcuts[i])*(x<boxcuts[i+1])])
    ax.boxplot(boxdata)
    add += '-boxplot'
if '--scatter' in sys.argv:
    ax.scatter(x/10.+0.5, y, alpha = 0.2, color = 'grey')
    add+= '-scatter'
if '--ylabel' in sys.argv:
    ax.set_ylabel(sys.argv[sys.argv.index('--ylabel')+1])



ax.set_xticks(np.arange(11)+0.5)
ax.set_xlim([.6,10.5])
ax.set_xticklabels(np.array(boxcuts))
yticks = np.arange(int(np.amin(y)/0.2)*0.2, 1.2, 0.2)
ax.set_yticks(yticks) #[0.,0.2,0.4,0.6,0.8,1.])
ax.set_xlabel('Sequence identity %')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)    

if '--setylim' in sys.argv:
    miny = float(sys.argv[sys.argv.index('--setylim')+1])
    maxy = float(sys.argv[sys.argv.index('--setylim')+2])
    ax.set_ylim([miny, maxy])

ylims = ax.get_ylim()
if '--selfcorr' in sys.argv:
    axself.set_ylim(ylims)

fig.savefig(os.path.splitext(os.path.split(sys.argv[3])[1])[0]+'-vs-'+os.path.splitext(os.path.split(sys.argv[2])[1])[0]+add+'.jpg', bbox_inches = 'tight', transparent=True, dpi = 300)
print os.path.splitext(os.path.split(sys.argv[3])[1])[0]+'-vs-'+os.path.splitext(os.path.split(sys.argv[2])[1])[0]+add+'.jpg'
#fig.tight_layout()
plt.show()



