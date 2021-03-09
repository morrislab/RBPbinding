# precision_recall_single_prediction.py
import numpy as np
import sys, os 
from functools import reduce
import matplotlib.pyplot as plt
import matplotlib.cm as cm 

colors =['darkgrey', 'limegreen','blueviolet','steelblue', 'goldenrod', 'black', 'sienna', 'firebrick', 'turquoise']

if '--definecolors' in sys.argv:
    colors = np.array(colors)[np.array(sys.argv[sys.argv.index('--definecolors')+1].split(','), dtype = int)]


## determine how many sets you want to compare
numsets = int(sys.argv[1])

# determine the scores you want to sort by
# - before the name reverses scores
# ,int(column) after the name defines the column in the text file with the score
allnames = []
if '--scores' in sys.argv:
    scorefiles = []
    scorenames = []
    scores = []
    for n in range(numsets):
        scorename = sys.argv[sys.argv.index('--scores')+1+n]
        if scorename[0] == '-':
            scorename = scorename[1:]
            mfac = -1.
        else:
            mfac = 1.
        if ',' in scorename:
            scorename = scorename.rsplit(',',1)
            scorerow = int(scorename[-1])
            scorename = scorename[0]
        else:
            scorerow = 1
        scorefiles.append(os.path.split(scorename)[0])
        nset = np.genfromtxt(scorename, dtype = str)
        scorenames.append(nset[:, 0])
        scores.append(mfac*nset[:, scorerow].astype(float))
        allnames.append(nset[:, 0])

# provide the reconstruction performance for every rbp
# ,int(column) after the name defines the column in the text file with the reconstruction
if '--realval' in sys.argv:
    realnames = []
    realsim = []
    for n in range(numsets):
        simname = sys.argv[sys.argv.index('--realval')+1+n]
        if ',' in simname:
            simname = simname.rsplit(',',1)
            simrow = int(simname[-1])
            simname = simname[0]
        else:
            simrow = -1
        nset = np.genfromtxt(simname, dtype = str)
        #print nset
        realnames.append(nset[:, 0])
        realsim.append(nset[:, simrow].astype(float))
        allnames.append(nset[:, 0])

# use to give names to the data sets 
if '--setnames' in sys.argv:
    scorefiles = sys.argv[sys.argv.index('--setnames')+1].split(',')

if '--savefig' in sys.argv:
    figname = sys.argv[sys.argv.index('--savefig')+1]


# align scores and realvals to each other on commonnames
commonnames = np.sort(reduce(np.intersect1d, allnames))
print len(commonnames)
scoresort = []
for s, scorename in enumerate(scorenames):
    scsort = []
    for commonname in commonnames:
        scsort.append(scores[s][list(scorename).index(commonname)])
    scoresort.append(scsort)
    
realsort = []
for s, realname in enumerate(realnames):
    resort = []
    for commonname in commonnames:
        resort.append(realsim[s][list(realname).index(commonname)])
    realsort.append(resort)


scores = np.nan_to_num(np.array(scoresort))
realsim = np.nan_to_num(np.array(realsort))


sortscore = np.argsort(-scores, axis = 1)
#print 'SET0'
#print scores[0, sortscore[0,:20]]
#print realsim[0, sortscore[0,:20]]
#print 'SET1'
#print scores[1, sortscore[1,:20]]
#print realsim[1, sortscore[1,:20]]


# can make the reconstructions binary if certain threshold defines correct or incorrect
if '--cutoff' in sys.argv:
    simcut = float(sys.argv[sys.argv.index('--cutoff')+1])
    realsim[realsim >= simcut] = 1
    realsim[realsim < simcut] = 0
    
precision = []
precvar = []
avprec = []
recall = []
xlabels = []
xticklocs = []
sortscores = []
sortsims = []
splitnum = 8 # number of steps on x-axis 
for n in range(numsets):
    # sort reconstruction measures to scores
    sortvalues = np.argsort(-scores[n])
    sortscores.append(scores[n][sortvalues])
    sortsims.append(realsim[n][sortvalues])
    print scorefiles[n]
    print ' '.join(commonnames[sortvalues[:20]].astype(str))
    print ' '.join(scores[n][sortvalues[:20]].astype(str))
    print ' '.join(realsim[n][sortvalues[:20]].astype(str))
    # addjust scale of scores shown in the figure
    scale = np.mean(np.diff(scores[n][sortvalues]))
    if np.absolute(np.log10(scale)) > 3:
        factor10 = 10.**round(scale)
    else:
        factor10 = 1.
    print 'Factor', factor10
    # make labels for x-axis 
    xlabel = [round(factor10 * scores[n][sortvalues][0],3)]
    xtickl = [0]
    steps = int(float(len(sortvalues))/float(splitnum))+int(float(len(sortvalues))%float(splitnum)!=0)
    for i in range(1,splitnum):
        #print steps, i
        xlabel.append(round(factor10 * scores[n][sortvalues][steps*i],3))
        xtickl.append(steps*i)
        
    #print 'Xlabel', xlabel
    xlabels.append(xlabel)
    xticklocs.append(xtickl)
    
    # cummulative mean performance
    meannorm = np.arange(2, len(sortvalues)+1)
    cumprecision = np.cumsum(realsim[n][sortvalues])[1:]/meannorm
    cumprecision = np.append(cumprecision[[0]], cumprecision)
    precision.append(cumprecision)
    # add standard errors
    pvar = []
    for rs in range(len(sortvalues)):
        if rs > 1:
            xche = realsim[n][sortvalues][:rs]
            xcstd = np.std(xche)/np.sqrt(len(xche)-1.)
            q75 = precision[-1][rs] + xcstd
            q25 = precision[-1][rs] - xcstd
            pvar.append([q75, q25])
    adv = precision[-1][:2]
    precvar.append(np.append([[adv[0], adv[0]], [adv[1], adv[1]]], pvar, axis = 0))
    
    # Average precision in bins of 'dis'
    aprec = []
    avar = []
    dis = 25
    offset = 5
    for rs, rr in enumerate(sortvalues):
        if rs < dis:
            xche = realsim[n][sortvalues][:max(offset,rs)]
        else:
            xche = realsim[n][sortvalues][rs -dis: rs+1]
        
        xcmean = np.mean(xche)
        xcstd = np.std(xche)
        q75 = xcmean + xcstd/np.sqrt(len(xche)-1.)
        q25 = xcmean - xcstd/np.sqrt(len(xche)-1.)
        avar.append([q75, q25])
        aprec.append(np.mean(xche))
    avprec.append([np.array(aprec), np.array(avar)])
    recall.append(np.cumsum(realsim[n][sortvalues])/float(len(sortvalues)+1))


# add line at a certain average precision
if '--ycutsum' in sys.argv:
    
    ycutsum = float(sys.argv[sys.argv.index('--ycutsum')+1])
    print '\nYcutsum', ycutsum
    xloc = []
    xscore = []
    for n in range(numsets):
        for m in range(len(precision[0])):
            if precision[n][-1-m] >= ycutsum:
                xloc.append(len(precision[n])-1-m)
                xscore.append(sortscores[n][-1-m])
                break
            if m == len(precision[0]) -1:
                xloc.append(0)
                xscore.append(sortscores[n][0])
        print xloc[-1], '\t', round(xloc[-1]*100./float(len(precision[n])),1),'\t',xscore[-1],'\t', scorefiles[n]
            

## can just print some statistics in terminal instead of making plots
# if not given plot will occur
if '--print_stats' not in sys.argv:
    if np.amax(xlabels[1]) > 90:
        xlabels[1] = np.array(xlabels[1], dtype = int)
    elif np.amax(xlabels[1]) <= 0:
        xlabels[1] = -np.around(xlabels[1], 3)
    elif np.amax(xlabels[1]) < 1.1:
        xlabels[1] = np.around(xlabels[1], 3)
    
    xlabels[1] = list(xlabels[1])
    if np.amax(xlabels[0]) > 90:
        xlabels[0] = np.array(xlabels[0], dtype = int)
    elif np.amax(xlabels[0]) <= 0:
        xlabels[0] = -np.around(xlabels[0], 3)
    elif np.amax(xlabels[0]) < 1.1:
        xlabels[0] = np.around(xlabels[0], 3)

    xlabels[0] = list(xlabels[0])
    
    fig = plt.figure(figsize = (4.5,4.5), dpi= 100) #figsize=(5,5), dpi = 100)
    ax = fig.add_subplot(111)
    for n in range(numsets):
        xax = np.linspace(0.,1.,len(precision[n]))
        print len(precision[n]), len(scorefiles[n]), colors[n]
        ax.plot(xax, precision[n], label = scorefiles[n], c = colors[n])
        ax.fill_between(np.linspace(0.,1.,len(precision[n])), precvar[n][:,0], precvar[n][:,1], alpha = 0.3, facecolor = colors[n]) 
    
    ylim = [np.amin(np.array(precision)[:, -1])-0.05, 1.]
    ax.set_ylim(ylim)
    if '--ycutsum' in sys.argv:
        ax.plot([0., 1.], [ycutsum,ycutsum], c = 'crimson', alpha = 0.6, linestyle = '--')
        for n in range(2):
            ax.plot([xax[xloc[n]], xax[xloc[n]]], [ylim[n], ycutsum], c = colors[n], alpha = 0.6, linestyle = '--')
            ax.annotate(str(np.around(xscore[n],3)), (xax[xloc[n]]+0.02, np.absolute(ylim[n]+ycutsum)/2.), color = colors[n])
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    #ax.spines['left'].set_color('dimgray')
    #ax.spines['bottom'].set_color('dimgray')
    #ax.tick_params(color = 'dimgray', width = 1.)
    
    
    ax2 = ax.twiny()
    ax2.spines['right'].set_visible(False)
    ax2.set_xlim(ax.get_xlim())
    ax2.spines["top"].set_position(("axes", 1.05))
    ax2.set_xticks(xax[xticklocs[0]])
    ax2.tick_params(color = colors[0], width = 2.)
    ax2.spines['top'].set_color(colors[0])
    ax2.spines['top'].set_linewidth(2.)
    ax2.set_xticklabels(xlabels[0], rotation = 60, color =colors[0], horizontalalignment = 'left')
    
    ax2b = ax.twiny()
    ax2b.spines['right'].set_visible(False)
    ax2b.set_xlim(ax.get_xlim())
    ax2b.spines["top"].set_position(("axes", 1.2))
    ax2b.set_xticks(xax[xticklocs[1]])
    ax2b.tick_params(color = colors[1], width = 2.)
    ax2b.spines['top'].set_color(colors[1])
    ax2b.spines['top'].set_linewidth(2.)
    ax2b.set_xticklabels(xlabels[1], rotation = 60, color =colors[1], horizontalalignment = 'left')
    
    ax.set_xticks(np.linspace(0.,1., 6))
    ax.set_xticklabels([0,20,40,60,80,100])
    ax.grid(which = 'both', linestyle = '--')
    
    if '--legend' in sys.argv:
        ax.legend()








if '--ycut' in sys.argv:
    
    ycut = float(sys.argv[sys.argv.index('--ycut')+1])
    print '\nycut', ycut
    xloc2 = []
    xscore2 = []
    xloc2f = []
    xscoref = []
    xcross = []
    for n in range(numsets):
        cross = 0
        for m in range(len(avprec[0][0])):
            if avprec[n][0][-1-m] >= ycut:
                xloc2.append(len(avprec[n][0])-m-1)
                xscore2.append(np.around(sortscores[n][-1-m],4))
                break
        for m in range(len(avprec[0][0])):
            if avprec[n][0][m] <= ycut:
                xloc2f.append(m)
                xscoref.append(np.around(sortscores[n][m],4))
                break
        for mn in range(xloc2f[-1], xloc2[-1]):
            if np.sum(avprec[n][0][mn:mn+2] >= ycut) == 1:
                cross += 1
        print 'Firstcross, recall, scorefirst | scoremid | Lastcross, recall, scorelast | Ncross, Method'
        print xloc2f[-1], '\t', round(xloc2f[-1]*100./float(len(precision[n])),1), '\t', xscoref[-1], ' | ', np.around(sortscores[n][int((xloc2f[-1]+xloc2[-1])/2)],4), ' | ', xloc2[-1], '\t', round(xloc2[-1]*100./float(len(precision[n])),1), '\t', xscore2[-1],' | ', cross, '\t', scorefiles[n]


if '--print_stats' not in sys.argv:
    fig2 = plt.figure(figsize = (4.5,4.5), dpi = 100)
    ax3 = fig2.add_subplot(111)
    for n in range(numsets):
        xax = np.linspace(0.,1.,len(avprec[n][0]))
        ax3.plot(xax, avprec[n][0], label = scorefiles[n], c = colors[n]) #, marker = 'o')
        #print avprec[n][1][:,0]
        ax3.fill_between(np.linspace(0.,1.,len(precision[n])), avprec[n][1][:,0], avprec[n][1][:,1], alpha = 0.3, facecolor = colors[n])
    ylim = ax3.get_ylim()
    if '--ycut' in sys.argv:
        ax3.plot([0., 1.], [ycut,ycut], c = 'crimson', alpha = 0.6, linestyle = '--')
        for n in range(2):
            ax3.plot([xax[xloc2f[n]], xax[xloc2f[n]]], [ylim[n], ycut], c = colors[n], alpha = 0.6, linestyle = '--')
            ax3.annotate(str(np.absolute(np.around(xscoref[n],3))), (xax[xloc2f[n]]+0.02, np.absolute(n+ycut)/2.), color = colors[n])

    ax3.set_ylim([0,1])
    ax3.set_yticks([0.2,0.4,0.6,0.8,1.])
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['left'].set_linewidth(1.5)
    ax3.spines['bottom'].set_linewidth(1.5)
    #ax3.spines['left'].set_color('dimgray')
    #ax3.spines['bottom'].set_color('dimgray')
    #ax3.tick_params(color = 'dimgray', width = 1.)
    
    
    
    ax4 = ax3.twiny()
    ax4.spines['right'].set_visible(False)
    ax4.set_xlim(ax3.get_xlim())
    ax4.spines["top"].set_position(("axes", 1.05))
    ax4.set_xticks(xax[xticklocs[0]])
    ax4.tick_params(color = colors[0], width = 2.)
    ax4.spines['top'].set_color(colors[0])
    ax4.spines['top'].set_linewidth(2.)
    ax4.set_xticklabels(xlabels[0], rotation = 60, color =colors[0], horizontalalignment = 'left')
    
    ax4b = ax3.twiny()
    ax4b.spines['right'].set_visible(False)
    ax4b.set_xlim(ax3.get_xlim())
    ax4b.spines["top"].set_position(("axes", 1.2))
    ax4b.set_xticks(xax[xticklocs[1]])
    ax4b.tick_params(color = colors[1], width = 2.)
    ax4b.spines['top'].set_color(colors[1])
    ax4b.spines['top'].set_linewidth(2.)
    ax4b.set_xticklabels(xlabels[1], rotation = 60, color =colors[1], horizontalalignment = 'left')
    
    ax3.set_xticks(np.linspace(0.,1., 6))
    ax3.set_xticklabels([0,20,40,60,80,100])
    ax3.grid(which = 'both', linestyle = '--')
    if '--legend' in sys.argv:
        ax3.legend()




    fig3 = plt.figure(figsize = (4.5,4.5), dpi = 100)
    ax5 = fig3.add_subplot(111)
    fig3.subplots_adjust(top=0.85)
    
    xlabels = []
    xticks = []
    for n in range(numsets):
        lscores = len(sortsims[n])
        width = 0.07/numsets
        loc = 1./numsets
        xlab = [sortscores[n][lscores*(1+j*2)/20] for j in range(10)]
        xlabels.append(xlab)
        xticks.append(np.array([(i+n*10.*width)/10.+width/2. for i in range(10)]))
        
        for i in range(10):
            ax5.boxplot(sortsims[n][i*lscores/10:(i+1)*lscores/10], widths = width, positions = [(i+n*10.*width)/10.+width/2.], patch_artist = True, boxprops=dict(facecolor = colors[n]))
    
    if np.amax(xlabels[1]) > 90:
        xlabels[1] = np.array(xlabels[1], dtype = int)
    elif np.amax(xlabels[1]) <= 0:
        xlabels[1] = -np.around(xlabels[1], 3)
    elif np.amax(xlabels[1]) < 1.1:
        xlabels[1] = np.around(xlabels[1], 3)
    
    xlabels[1] = list(xlabels[1])
    if np.amax(xlabels[0]) > 90:
        xlabels[0] = np.array(xlabels[0], dtype = int)
    elif np.amax(xlabels[0]) <= 0:
        xlabels[0] = -np.around(xlabels[0], 3)
    elif np.amax(xlabels[0]) < 1.1:
        xlabels[0] = np.around(xlabels[0], 3)
    
    
    #ax5.set_xticks(xticks[0])
    xlim =[-0.05,1.0]
    ax5.set_xlim(xlim)
    ylim = ax5.get_ylim()
    #if '--ycut' in sys.argv:
        #ax3.plot([0., 1.], [ycut,ycut], c = 'crimson', alpha = 0.6, linestyle = '--')
        #for n in range(2):
            #ax3.plot([xax[xloc2[n]], xax[xloc2[n]]], [ylim[n], ycut], c = colors[n], alpha = 0.6, linestyle = '--')
            #ax3.annotate(str(np.around(xscore2[n],3)), (xax[xloc2[n]]+0.02, np.absolute(n+ycut)/2.), color = colors[n])
    ax5.spines['right'].set_visible(False)
    ax5.spines['top'].set_visible(False)
    ax5.spines['left'].set_linewidth(1.5)
    ax5.spines['bottom'].set_linewidth(1.5)
    #ax5.spines['left'].set_color('dimgray')
    #ax5.spines['bottom'].set_color('dimgray')
    #ax5.tick_params(width = 1.5)
    
    
    
    ax6 = ax5.twiny()
    ax6.spines['right'].set_visible(False)
    #ax6.spines['left'].set_visible(False)
    ax6.set_xlim(ax5.get_xlim())
    ax6.spines["top"].set_position(("axes", 1.05))
    ax6.set_xticks(xticks[0])
    ax6.tick_params(color = colors[0], width = 2.)
    ax6.spines['top'].set_color(colors[0])
    ax6.spines['top'].set_linewidth(2.)
    ax6.set_xticklabels(xlabels[0], rotation = 90, color =colors[0])
    
    ax7 = ax5.twiny()
    ax7.spines['right'].set_visible(False)
    #ax7.spines['left'].set_visible(False)
    ax7.set_xlim(ax5.get_xlim())
    ax7.spines["top"].set_position(("axes", 1.25))
    ax7.set_xticks(xticks[1])
    ax7.tick_params(color = colors[1], width = 2.)
    ax7.spines['top'].set_color(colors[1])
    ax7.spines['top'].set_linewidth(2.)
    ax7.set_xticklabels(xlabels[1], rotation = 90, color =colors[1])
    
    
    scorecor = np.array([int()])
    ax5.set_xticks(np.mean(xticks, axis = 0))
    ax5.set_xticklabels(['0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '90-100', '100'], rotation = 60, ha='right', va = 'top')
    ax5.set_xlabel('Recall')
    
    
    
    
    
    
    
    
    
    
    fig.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    if '--savefig' in sys.argv:
        fig.savefig(os.path.splitext(figname)[0]+'-cummeanprecision'+os.path.splitext(figname)[-1], dpi = 300)
        fig2.savefig(os.path.splitext(figname)[0]+'-Meanprecision'+os.path.splitext(figname)[-1], dpi = 300)
        fig3.savefig(os.path.splitext(figname)[0]+'-Meanprecisionboxplot'+os.path.splitext(figname)[-1], dpi = 300)
    else:
        plt.show()





