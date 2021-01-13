import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


masterfile = open(sys.argv[1], 'r').readlines()
if len(sys.argv) >=3:
    outdir = sys.argv[2]
else:
    outdir = ''

header = masterfile[0].strip().split('\t')
stats = []
for l, line in enumerate(masterfile):
    if l > 0:
        stats.append(line.strip().split('\t'))
stats = np.array(stats)

clades = stats[:, list(header).index('Kingdom')]
project = stats[:, list(header).index('Project')]

cladetotal, total = np.unique(clades, return_counts = True)
cladeold, old = np.unique(clades[project == 'Nature Paper 2013'], return_counts = True)
cladenew, new = np.unique(clades[project != 'Nature Paper 2013'], return_counts = True)

sort = np.argsort(-total)
cladetotal = cladetotal[sort]
old = old[sort]
new = new[sort]
total = total[sort]
oldpie = np.append(old, [np.sum(total)-np.sum(old)])
newpie = np.append([np.sum(total)-np.sum(new)], new)
oldnew = [np.sum(total)-np.sum(new),np.sum(total)-np.sum(old)]

print oldpie, newpie, oldnew, old, new, total, cladetotal

colorsold = [(.9,0.6,0.3,1), (.4,0.1,0.4,1),(.3,0.7,0.3,1),(.6,0.2,0.2,1),(1.,1.,1.,0)]
colorsnew = [(1,1,1,0),(.9,0.6,0.3,1), (.4,0.1,0.4,1),(.3,0.7,0.3,1),(.6,0.2,0.2,1)]

fig = plt.figure(figsize = (1.5,1.5), dpi = 400)
ax = fig.add_subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(which = 'both', left = False, bottom = False, labelbottom = False, labelleft = False)

ax.pie(oldnew, radius = 1.3, center = (0.1,0.02), labels = None, startangle=90, colors = ['white', 'dimgrey'])
ax.pie(newpie, center = (0.1,0.02),labels = None, startangle=90, colors = colorsnew)
ax.pie(oldpie, labels = None, startangle=90, colors = colorsold)
ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

fig.savefig(outdir+'Piekingdom.jpg', dpi = 300, bbox_inches = 'tight')

plt.show()

reccolors = [(.9,0.6,0.3,.9), (.4,0.1,0.4,.9), (.3,0.7,0.3,.9), (.6,0.2,0.2,.9)]
recognizedtypes = [cladetotal[i]+' '+str(old[i])+'+'+str(new[i])+'='+str(total[i]) for i in range(len(total))]
figlegend = plt.figure()
axf = figlegend.add_subplot(111)
axf.spines['top'].set_visible(False)
axf.spines['bottom'].set_visible(False)
axf.spines['left'].set_visible(False)
axf.spines['right'].set_visible(False)
axf.tick_params(which = 'both', left = False, bottom = False, labelbottom = False, labelleft = False)
for r, rdt in enumerate(recognizedtypes):
    axf.text(0.007, r*-0.1+0.04, rdt, fontsize=32, horizontalalignment='left', verticalalignment='center')
    axf.add_patch(Rectangle(xy=(0, r*-0.1), width=0.006,height=0.08, facecolor=reccolors[r]))
    axf.set_xlim([-0.01,0.06])
    axf.set_ylim([-0.5, 0.1])

figlegend.savefig(outdir+'Kingdomlegend.jpg', dpi = 300, bbox_inches = 'tight')




plt.show()




