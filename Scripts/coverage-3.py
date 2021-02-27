## coverage.py
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import glob
import pandas as pd


statsfile = sys.argv[1]
stats = np.genfromtxt(statsfile, dtype = str, delimiter = ',', skip_header = 1)
scol = open(statsfile, 'r').readline().strip().split(', ')[1:]
print scol
statspecies = stats[:, 0]
stats = stats[:,1:].astype(float)

if '--savefig' in sys.argv:
    plotname = sys.argv[sys.argv.index('--savefig')+1]
    plotfmt = plotname.split('.')[-1]


if '--specieslist' in sys.argv:
    lspecies = sys.argv[sys.argv.index('--specieslist')+1]
    if ',' in lspecies:
        mspec = lspecies.split(',')
    else:
        mspec = np.genfromtxt(lspecies, dtype =str)

else:
    print 'Provide list of species'
    sys.exit()
    
    
king = False
average = ''
if '--kingdoms' in sys.argv:
    king = True
    kingdoms = np.genfromtxt(sys.argv[sys.argv.index('--kingdoms')+1], dtype = str)
    average = sys.argv[sys.argv.index('--kingdoms')+2] # fraction, total
    if ~np.array_equal(np.sort(kingdoms[:,0]), mspec):
        resort = []
        for m, sp in enumerate(mspec):
            if sp in kingdoms[:,0]:
                resort.append([m, list(kingdoms[:,0]).index(sp)])
            else:
                print sp, 'not in kingdomfile'
        resort = np.array(resort).T
        mspec = mspec[resort[0]]
        kingdoms = kingdoms[resort[1]]

elif '--kingcolor' in sys.argv:
    kingcol = True
    kingdoms = np.genfromtxt(sys.argv[sys.argv.index('--kingcolor')+1], dtype = str)
    kingrow = int(sys.argv[sys.argv.index('--kingcolor')+2])
    if ~np.array_equal(np.sort(kingdoms[:,0]), mspec):
        resort = []
        for m, sp in enumerate(mspec):
            if sp in kingdoms[:,0]:
                resort.append([m, list(kingdoms[:,0]).index(sp)])
            else:
                print sp, 'not in kingdomfile'
        resort = np.array(resort).T
        mspec = mspec[resort[0]]
        kingdoms = kingdoms[resort[1]]   
    unking = np.unique(kingdoms[:, kingrow])
    kingcolors = np.zeros(len(kingdoms), dtype = int)
    for u, uk in enumerate(unking):
        kingcolors[kingdoms[:, kingrow] == uk] = u
        
    


statsort = []
for m, sp in enumerate(mspec):
    statsort.append(list(statspecies).index(sp))

print np.shape(stats)
stats = stats[statsort]
statspecies = statspecies[statsort]


totspeci = np.zeros((len(stats), 4))
totspeci[:,0] = stats[:,scol.index('Nrbpprot')]
totspeci[:,2] = stats[:,scol.index('NprotRRM')]
totspeci[:,3] = stats[:,scol.index('NprotKH')]
totspeci[:,1] = stats[:,scol.index('NprotRRM')] + stats[:,scol.index('NprotKH')]

inputs = np.array([scol.index('Nmeas'), scol.index('Nrectot'), scol.index('Nrecjple'), scol.index('Nrecid'),scol.index('Nrecinter'), scol.index('NmeasRRM'), scol.index('NrectotRRM'), scol.index('NrecjpleRRM'), scol.index('NrecidRRM'),scol.index('NrecinterRRM'), scol.index('NmeasKH'), scol.index('NrectotKH'), scol.index('NrecjpleKH'), scol.index('NrecidKH'),scol.index('NrecinterKH')])
#inputs = np.array([7, 10, 11, 12, 13, 8, 18, 19, 20, 21, 9, 26, 27, 28, 29])-1
stats = stats[:, inputs]





if king:
    kings = np.unique(kingdoms[:,1])
    print kings
    
    emeans = np.zeros(4)
    estats = np.zeros(15)
    
    kmeans = np.zeros((len(kings), 4))
    kstats = np.zeros((len(kings), 15))
    
    ktotdist = []
    
    emeans = np.mean(totspeci, axis = 0)
    
    for k, king in enumerate(kings): 
        kmeans[k] = np.mean(totspeci[kingdoms[:,1]==king], axis = 0)
        #print np.shape(totspeci[kingdoms[:,1]==king].T)
        ktotdist.append(totspeci[kingdoms[:,1]==king].T[[0,2,3]])
    print kmeans
    print emeans
    
    if average == 'fraction':
        for l in range(15):
            #print int(l/5)+1
            efraction = stats[:,l]/totspeci[:,int(l/5) + 1]
            efraction[np.isnan(efraction)] = 0.
            efraction[np.isinf(efraction)] = 0.
            estats[l] = np.mean(efraction)
            for g in range(len(kings)):
                kfraction = stats[kingdoms[:,1] == kings[g], l].T/totspeci[kingdoms[:,1] == kings[g],int(l/5) + 1]
                kfraction[np.isnan(kfraction)] = 0.
                kfraction[np.isinf(kfraction)] = 0.
                kstats[g, l] = np.mean(kfraction)
            
    elif average == 'total':
        for l in range(15):
            estats[l] = np.sum(stats[:,l])/np.sum(totspeci[:,int(l/5) + 1])
            for g in range(len(kings)):
                #print l, g, np.sum(stats[kingdoms[:,1] == kings[g], l]), np.sum(totspeci[kingdoms[:,1] == kings[g],int(l/5) + 1]), np.sum(stats[kingdoms[:,1] == kings[g], l])/np.sum(totspeci[kingdoms[:,1] == kings[g],int(l/5) + 1])
                kstats[g, l] = np.sum(stats[kingdoms[:,1] == kings[g], l])/np.sum(totspeci[kingdoms[:,1] == kings[g],int(l/5) + 1])
    else:
        print 'average not understood'
        sys.exit()
    print kstats
    stats = kstats
    mspec = kings
    totspeci = kmeans
    if '--sortkingdom' in sys.argv:
        sortk = np.array(sys.argv[sys.argv.index('--sortkingdom')+1].split(','), dtype = int)
        stats = stats[sortk]
        mspec = mspec[sortk]
        totspeci = totspeci[sortk]
    else:
        sortk = np.arange(len(kings))
    ktotbox = []
    for i in range(3):
        ktbox = []
        for k in sortk:
            ktbox.append(ktotdist[k][i])
        ktotbox.append(ktbox)

else:
    for l in range(15):
        stats[:, l] = stats[:,l].T/totspeci[:,int(l/5) + 1]

#print stats[:5, [12,13]]
#sys.exit()

if '--figsize' in sys.argv:
    xs, ys = sys.argv[sys.argv.index('--figsize')+1].split(',')
    ys = float(ys)
    xs = float(xs)
else:
    ys = 0.2*len(mspec)
    xs = 2.5

for m, ms in enumerate(mspec):
    mspec[m] = ms.replace('_', ' ')



if '--scattercomparison' in sys.argv:
    figscat = plt.figure(figsize = (xs,ys))
    axs = figscat.add_subplot(111)
    axs.set_xlabel('% Covered ID')
    axs.set_ylabel('% Covered JPLE')
    #axs.patch.set_facecolor('whitesmoke')
    if king:
        markers = ['s', 'o']
        colors = ['firebrick', 'darkgoldenrod', 'seagreen', 'darkslateblue']
        for m, ms in enumerate(mspec):
            axs.scatter([stats[m, 8]], [stats[m, 7]], alpha=.8, color=colors[m], s = 150, marker = markers[0], label = ms)
            axs.scatter([stats[m, 13]], [stats[m, 12]], alpha=.8, color=colors[m], s = 150, marker = markers[1])
    elif kingcol:
        markers = ['s', 'o']
        colors = np.array(['firebrick', 'darkgoldenrod', 'seagreen', 'darkslateblue'])
        #axs.scatter(stats[:, 8], stats[:, 7], alpha=.8, color=colors[kingcolors], s = 100, marker = markers[0])
        #axs.scatter(stats[:, 13], stats[:, 12], alpha=.8, color=colors[kingcolors], s = 100, marker = markers[1])
        for c, col in enumerate(unking):
            if '--rrm' in sys.argv:
                axs.scatter(stats[kingcolors == c, 8], stats[kingcolors == c, 7], alpha=.8, color=colors[c], s = 100, marker = markers[0], label = col)
            elif '--kh' in sys.argv:
                axs.scatter(stats[kingcolors == c, 13], stats[kingcolors == c, 12], alpha=.8, color=colors[c], s = 100, marker = markers[1], label = col)
            else:
                axs.scatter(stats[kingcolors == c, 3], stats[kingcolors == c, 2], alpha=.8, color=colors[c], s = 100, marker = markers[1], label = col)
    else:
        axs.scatter(stats[:, 8], stats[:, 7], alpha=.8, color='darkslateblue', s = 100, label = 'RRM')
        axs.scatter(stats[:, 13], stats[:, 12], alpha=.8, color='firebrick', s = 100, label = 'KH')
    axs.plot([0,1], [0,1], color = 'grey', alpha = 0.5)
    major_xticks = np.linspace(0.,1.,11)
    axs.set_xticks(major_xticks)
    axs.set_yticks(major_xticks)
    axs.set_xticklabels(np.linspace(0,100,11).astype(int).astype(str), rotation = 90)
    axs.set_yticklabels(np.linspace(0,100,11).astype(int).astype(str), rotation = 0)
    axs.set_xlim([-0.05,.8])
    axs.set_ylim([-0.05,.8])
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    axs.legend()
    if '--savefig' in sys.argv:
        figscat.savefig(os.path.splitext(plotname)[0]+'-scatteridvsjple.'+plotfmt, format=plotfmt, dpi = 300, bbox_inches = 'tight')
    else:
        plt.show()
else:

        
    fig = plt.figure(figsize = (xs,ys))
    ax = fig.add_subplot(111)
    ax.set_xlabel('% Covered')
    ax.patch.set_facecolor('whitesmoke')
    ax.barh(np.arange(1, len(mspec)+1, 1), stats[:,1], alpha=1., color='limegreen', label = 'JPLE only')
    ax.barh(np.arange(1, len(mspec)+1, 1), stats[:,3], alpha=1., color='dimgrey', label = 'ID only')
    ax.barh(np.arange(1, len(mspec)+1, 1), stats[:,4], alpha=1., color='darkgrey', label = 'JPLE '+r'$\cap$'+' ID')
    ax.barh(np.arange(1, len(mspec)+1, 1), stats[:,0], alpha=1., color='black', label = 'Measured')
    if king:
        # eukaryotes
        ax.barh([len(mspec)+1.5], [estats[1]], alpha=1., color='limegreen')
        ax.barh([len(mspec)+1.5], [estats[3]], alpha=1., color='dimgrey')
        ax.barh([len(mspec)+1.5], [estats[4]], alpha=1., color='darkgrey')
        ax.barh([len(mspec)+1.5], [estats[0]], alpha=1., color='black')

        ax.set_yticks(np.append(np.arange(1., len(mspec)+1., 1), [len(mspec)+1.5]))
        ax.set_yticklabels(np.append(mspec, ['Eukaryotes']))

    else:
        ax.set_yticks(np.arange(1., len(mspec)+1., 1))
        ax.set_yticklabels(mspec)
    major_xticks = np.linspace(0.,1.,6)
    minor_xticks = np.linspace(0.1,1.1,6)
    ax.set_xticks(major_xticks)
    ax.set_xticks(minor_xticks, minor = True)
    ax.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
    ax.grid(b=True, which='minor', axis='x', alpha = 0.7, linestyle = '--')
    ax.set_xticklabels(np.linspace(0,100,6).astype(int).astype(str), rotation = 60)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_title('Reconstruction RBP ('+average+')')
    #ax.spines['bottom'].set_linewidth(2.)
    if '--legend' in sys.argv:
        ax.legend()




    figb = plt.figure(figsize = (xs,ys))
    axb = figb.add_subplot(111)
    axb.set_xlabel('% Covered')
    axb.patch.set_facecolor('whitesmoke')
    axb.barh(np.arange(1, len(mspec)+1, 1), stats[:,6], alpha=1., color='limegreen', label = 'JPLE only')
    axb.barh(np.arange(1, len(mspec)+1, 1), stats[:,8], alpha=1., color='dimgrey', label = 'ID only')
    axb.barh(np.arange(1, len(mspec)+1, 1), stats[:,9], alpha=1., color='darkgrey', label = 'JPLE+ID')
    axb.barh(np.arange(1, len(mspec)+1, 1), stats[:,5], alpha=1., color='black', label = 'Measured')
    if king:
        # eukaryotes
        axb.barh([len(mspec)+1.5], [estats[6]], alpha=1., color='limegreen')
        axb.barh([len(mspec)+1.5], [estats[8]], alpha=1., color='dimgrey')
        axb.barh([len(mspec)+1.5], [estats[9]], alpha=1., color='darkgrey')
        axb.barh([len(mspec)+1.5], [estats[5]], alpha=1., color='black')

        axb.set_yticks(np.append(np.arange(1., len(mspec)+1., 1), [len(mspec)+1.5]))
        axb.set_yticklabels(np.append(mspec, ['Eukaryotes']))

    else:
        axb.set_yticks(np.arange(1., len(mspec)+1., 1))
        axb.set_yticklabels(mspec)
    major_xticks = np.linspace(0.,1.,11)
    minor_xticks = np.linspace(0.05,1.05,11)
    axb.set_xticks(major_xticks)
    #axb.set_xticks(minor_xticks, minor = True)
    axb.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
    #axb.grid(b=True, which='minor', axis='x', alpha = 0.5, linestyle = ':')
    axb.set_xticklabels(np.linspace(0,100,11).astype(int).astype(str), rotation = 60)
    axb.spines['right'].set_visible(False)
    axb.spines['top'].set_visible(False)
    axb.spines['left'].set_visible(False)
    axb.set_title('Reconstruction RRM ('+average+')')
    #axb.spines['bottom'].set_linewidth(2.)
    if '--legend' in sys.argv:
        axb.legend()
        
        
        
    figc = plt.figure(figsize = (xs,ys))
    axc = figc.add_subplot(111)
    axc.set_xlabel('% Covered')
    axc.patch.set_facecolor('whitesmoke')
    axc.barh(np.arange(1, len(mspec)+1, 1), stats[:,11], alpha=1., color='limegreen', label = 'JPLE only')
    axc.barh(np.arange(1, len(mspec)+1, 1), stats[:,13], alpha=1., color='dimgrey', label = 'ID only')
    axc.barh(np.arange(1, len(mspec)+1, 1), stats[:,14], alpha=1., color='darkgrey', label = 'JPLE+ID')
    axc.barh(np.arange(1, len(mspec)+1, 1), stats[:,10], alpha=1., color='black', label = 'Measured')
    if king:
        # eukaryotes
        axc.barh([len(mspec)+1.5], [estats[11]], alpha=1., color='limegreen')
        axc.barh([len(mspec)+1.5], [estats[13]], alpha=1., color='dimgrey')
        axc.barh([len(mspec)+1.5], [estats[14]], alpha=1., color='darkgrey')
        axc.barh([len(mspec)+1.5], [estats[10]], alpha=1., color='black')

        axc.set_yticks(np.append(np.arange(1., len(mspec)+1., 1), [len(mspec)+1.5]))
        axc.set_yticklabels(np.append(mspec, ['Eukaryotes']))

    else:
        axc.set_yticks(np.arange(1., len(mspec)+1., 1))
        axc.set_yticklabels(mspec)
    major_xticks = np.linspace(0.,1.,11)
    minor_xticks = np.linspace(0.05,1.05,11)
    axc.set_xticks(major_xticks)
    #ax.set_xticks(minor_xticks, minor = True)
    axc.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
    #ax.grid(b=True, which='minor', axis='x', alpha = 0.5, linestyle = ':')
    axc.set_xticklabels(np.linspace(0,100,11).astype(int).astype(str), rotation = 60)
    axc.spines['right'].set_visible(False)
    axc.spines['top'].set_visible(False)
    axc.spines['left'].set_visible(False)
    axc.set_title('Reconstruction KH ('+average+')')
    #ax.spines['bottom'].set_linewidth(2.)
    if '--legend' in sys.argv:
        axc.legend()





    fig2 = plt.figure(figsize = (xs,ys))
    ax2 = fig2.add_subplot(111)
    #ax2.set_title('Total numbers')
    ax2.set_xlabel('Number RBPs')
    ax2.patch.set_facecolor('whitesmoke')
    ax2.barh(np.arange(1, len(mspec)+1, 1), totspeci[:,0], alpha=1., color='k', edgecolor = 'k', linewidth = 1., label = 'cRBPs')
    ax2.barh(np.arange(1, len(mspec)+1, 1), totspeci[:,1], alpha=1., color='dimgrey', edgecolor = 'dimgrey', linewidth = 1., label = 'RRM' )
    ax2.barh(np.arange(1, len(mspec)+1, 1), totspeci[:,3], alpha=1., color='darkgrey', edgecolor = 'darkgrey', linewidth = 1., label = 'KH')

    if king:
        # eukaryotes
        ax2.barh([len(mspec)+1.5], emeans[0], alpha=1., color='k')
        ax2.barh([len(mspec)+1.5], emeans[1], alpha=1., color='dimgrey')
        ax2.barh([len(mspec)+1.5], emeans[3], alpha=1., color='darkgrey')
        
        ax2.set_yticks(np.append(np.arange(1., len(mspec)+1., 1), [len(mspec)+1.5]))
        ax2.set_yticklabels(np.append(mspec, ['Eukaryotes']))

    else:
        ax2.set_yticks(np.arange(1., len(mspec)+1., 1))
        ax2.set_yticklabels(mspec)
    ax2.set_xticks([0,100,200,300,400,500])
    ax2.set_xticklabels([0,100,200,300,400,500], rotation = 60)
    ax2.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
    #ax2.grid(b=True, which='minor', axis='x', alpha = 0.5, linestyle = ':')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    if '--legend' in sys.argv:
        ax2.legend() #bbox_to_anchor = (0.6, 0.6))

    if king:
        fig3 = plt.figure(figsize = (xs*3.,ys))
        ax3a= fig3.add_subplot(131)
        ax3b = fig3.add_subplot(132)
        ax3c = fig3.add_subplot(133)
        #ax3.set_title('Total numbers')
        ax3a.set_xlabel('Number RBPs')
        ax3b.set_xlabel('Number RBPs')
        ax3c.set_xlabel('Number RBPs')
        ax3a.set_title('canonical RBPs')
        ax3b.set_title('RRMs')
        ax3c.set_title('KHs')
        ax3a.patch.set_facecolor('whitesmoke')
        ax3b.patch.set_facecolor('whitesmoke')
        ax3c.patch.set_facecolor('whitesmoke')
        
        ax3a.boxplot(ktotbox[0], vert=False, patch_artist=True, labels=mspec, boxprops=dict(facecolor='k', alpha = 0.7))#, boxprops=dict(alpha=0.7, facecolor='k', color = k))
        ax3b.boxplot(ktotbox[1], vert=False, patch_artist=True, labels=mspec, boxprops=dict(facecolor='dimgrey', alpha = 0.7))#,  boxprops=dict(alpha=0.7, facecolor='dimgrey', color = k))
        ax3c.boxplot(ktotbox[2], vert=False, patch_artist=True, labels=mspec, boxprops=dict(facecolor='darkgrey', alpha = 0.7))#,  boxprops=dict(alpha=0.7, facecolor='darkgrey', color = k))

        
        ax3a.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
        ax3b.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
        ax3c.grid(b=True, which='major', axis='x', alpha = 0.7, linestyle = '-')
        #ax3.grid(b=True, which='minor', axis='x', alpha = 0.5, linestyle = ':')
        ax3a.spines['right'].set_visible(False)
        ax3a.spines['top'].set_visible(False)
        ax3a.spines['left'].set_visible(False)
        ax3b.spines['right'].set_visible(False)
        ax3b.spines['top'].set_visible(False)
        ax3b.spines['left'].set_visible(False)    
        ax3c.spines['right'].set_visible(False)
        ax3c.spines['top'].set_visible(False)
        ax3c.spines['left'].set_visible(False)

        fig3.tight_layout()
    
    if '--savefig' in sys.argv:
        fig.savefig(plotname, format=plotfmt, dpi = 300, bbox_inches = 'tight')
        figb.savefig(os.path.splitext(plotname)[0]+'-RRM.'+plotfmt, format=plotfmt, dpi = 300, bbox_inches = 'tight')
        figc.savefig(os.path.splitext(plotname)[0]+'-KH.'+plotfmt, format=plotfmt, dpi = 300, bbox_inches = 'tight')
        fig2.savefig(os.path.splitext(plotname)[0]+'-total.'+plotfmt, format=plotfmt, dpi = 300, bbox_inches = 'tight')
        if king:
            fig3.savefig(os.path.splitext(plotname)[0]+'-totalbox.'+plotfmt, format=plotfmt, dpi = 300, bbox_inches = 'tight')
    else:
        fig.tight_layout()
        figb.tight_layout()
        figc.tight_layout()
        if king:
            fig3.tight_layout()
        fig2.tight_layout()
        plt.show() 




