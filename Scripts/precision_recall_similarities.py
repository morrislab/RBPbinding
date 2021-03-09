# precision_recall_single_prediction.py
import numpy as np
import sys, os 
from functools import reduce
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
import seaborn as sns

    
def readinsimfile(simfile):
    sobj = open(simfile, 'r').readlines()
    t = 0
    outnames = []
    outtops = []
    trainnames = []
    for sline in sobj:
        if sline[0] == '#':
            if t > 0:
                outnames.append(tnames)
                tctop = np.array(tc, dtype = float).T
                for f in range(len(tnames)):
                    if tnames[f] in tcnames:
                        tradd = np.array(tcnames)
                        traddmask = tradd != tnames[f]
                        trainnames.append(tradd[traddmask])
                        outtops.append(np.array(tctop[f][traddmask], dtype = float))
                    else:
                        trainnames.append(tcnames)
                        outtops.append(np.array(tctop[f], dtype = float))
            tnames = sline.strip().split()[1:]
            tc = []
            tcnames = []
            t += 1
        else:
            tc.append(sline.strip().split()[1:])
            tcnames.append(sline.strip().split()[0])
    outnames.append(tnames)
    tctop = np.array(tc, dtype = float).T
    for f in range(len(tnames)):
        if tnames[f] in tcnames:
            tradd = np.array(tcnames)
            traddmask = tradd != tnames[f]
            trainnames.append(tradd[traddmask])
            outtops.append(np.array(tctop[f][traddmask], dtype = float))
        else:
            trainnames.append(tcnames)
            outtops.append(np.array(tctop[f], dtype = float))
            
    #outnames = np.concatenate(outnames)
    return outnames, trainnames, outtops

def readinclustermat(clustermat):
    clustnames = open(clustermat, 'r').readline().strip().split()[1:]
    if '||' in clustnames[0]:
        for c, cname in enumerate(clustnames):
            cname = cname.split('||')
            clustnames[c] = cname[-1]
    clusteringmat = np.genfromtxt(clustermat)
    return np.array(clustnames), clusteringmat


def opennpz(nfile, dtype, row):
    objm = np.load(nfile)
    testp = objm['testprots'].astype(str)
    if dtype == 'both':
        trpn = np.append(objm['testprots'].astype(str), objm['trainprots'].astype(str))
        if row == 'euc':
            tscore = np.append(objm['testtesteuc'], objm['testeuc'], axis = 1)
        elif row == 'cosine':
            tscore = np.append(objm['testtestcosine'], objm['testcosine'], axis = 1)
        elif row == 'pearson':
            tscore = np.append(objm['testtestpearson'], objm['testpearson'], axis = 1)
        elif row == 'dot':
            tscore = np.append(objm['testtestdot'], objm['testdot'], axis = 1)
    elif dtype == 'training':
        trpn = objm['trainprots'].astype(str)
        if row == 'euc':
            tscore = objm['testeuc']
        elif row == 'cosine':
            tscore = objm['testcosine']
        elif row == 'pearson':
            tscore = objm['testpearson']
        elif row == 'dot':
            tscore = objm['testdot']
        elif row == 'predeuc':
            tscore = objm['Yeucpred']
        elif row == 'predcosine':
            tscore = objm['Ycosinepred']
        elif row == 'predpearson':
            tscore = objm['Ypearsonpred']
        elif row == 'predsim':
            tscore = objm['Ysimpred']
    elif dtype == 'test':
        trpn = objm['testprots'].astype(str)
        if row == 'euc':
            tscore = objm['testtesteuc']
        elif row == 'cosine':
            tscore = objm['testtestcosine']
        elif row == 'pearson':
            tscore = objm['testtestpearson']
        elif row == 'dot':
            tscore = objm['testtestdot']
        elif row == 'predeuc':
            tscore = objm['Yeucpredtest']
        elif row == 'predcosine':
            tscore = objm['Ycosinepredtest']
        elif row == 'predpearson':
            tscore = objm['Ypearsonpredtest']
        elif row == 'predsim':
            tscore = objm['Ysimpredtest']
    return testp, trpn, tscore


def readinnpz(fstart, fend, numcv, row, mfac, dtype):
    outnames = []
    trainnames = []
    outtops = []
    if numcv == 0:
        testp, trpn, tscore = opennpz(fstart+fend, dtype, row)
        for t, tetp in enumerate(testp):
            trainmask = trpn != tetp
            trainnames.append(trpn[trainmask])
            outnames.append([tetp])
            outtops.append(np.array([mfac*tscore[t, trainmask]]))
    else:
        for nr in range(numcv):
            testp, trpn, tscore = opennpz(fstart+str(nr)+fend, dtype, row)
            for t, tetp in enumerate(testp):
                trainmask = trpn != tetp
                trainnames.append(trpn[trainmask])
                outnames.append([tetp])
                outtops.append(np.array([mfac*tscore[t, trainmask]]))
    outnames = np.array(outnames)
    #print outtops
    return outnames, np.array(trainnames), outtops
    






numsets = int(sys.argv[1])
colors =['limegreen', 'darkgrey', 'blueviolet','steelblue', 'goldenrod',  'firebrick', 'turquoise']
if numsets > len(colors) + 1 or '--rainbowc' in sys.argv:
    colors = cm.rainbow(np.arange(numsets+1)/float(numsets+1))
if '--definecolors' in sys.argv:
    colors = np.array(colors)[np.array(sys.argv[sys.argv.index('--definecolors')+1].split(','), dtype = int)]


allnames = []
scorefiles = []
scorenames1 = []
scorenames2 = []
scores = []
if '--scores' in sys.argv:
    for n in range(numsets):
        scorename = sys.argv[sys.argv.index('--scores')+1+n]
        if scorename[0] == '-':
            scorename = scorename[1:]
            mfac = -1.
        else:
            mfac = 1.
        scorefiles.append(os.path.splitext(os.path.split(scorename)[1])[0])
        sname1, sname2, sscore = readinsimfile(scorename)
        scorenames1.append(sname1)
        scorenames2.append(sname2)
        scores.append(sscore)
        allnames.append(sname1)


elif '--scorenpz' in sys.argv:
    numcv = int(sys.argv[sys.argv.index('--scorenpz')+1])
    scotype = sys.argv[sys.argv.index('--scorenpz')+2]
    for n in range(numsets):
        scorename, scorefname, scorerow = sys.argv[sys.argv.index('--scorenpz')+3+n].split(':')
        # to change the order of scores to -scores
        if scorename[0] == '-':
            scorename = scorename[1:]
            mfac = -1.
        else:
            mfac = 1.
        
        scorefiles.append(os.path.split(scorename)[0])
        sname1, sname2, sscore = readinnpz(scorename, scorefname, numcv, scorerow, mfac, scotype)
        scorenames1.append(sname1) ## test sets
        scorenames2.append(sname2) ## train sets
        scores.append(sscore)
        allnames.append(np.concatenate(sname1))



commonnames = np.sort(reduce(np.intersect1d, allnames))
print 'Commonnames', len(commonnames)


if '--Identities' in sys.argv:
    print 'read in Identitties'
    numsets += 1
    identityfile= sys.argv[sys.argv.index('--Identities')+1]
    
    scorefiles.append(os.path.splitext(identityfile)[0])
    idnames, idmat = readinclustermat(identityfile)
    
    if '>' in idnames[0]:
        for i, idname in enumerate(idnames):
            idnames[i] = idname[1:]
    
    keep = []
    for i, comname in enumerate(commonnames):
        keep.append(list(idnames).index(comname))
    idnames = idnames[keep]
    idmat = idmat[keep]
    idmat = idmat[:, keep]
    
    itest = []
    itrain = []
    iscore = []
    for i, idn in enumerate(idnames):
        itest.append([idn])
        itrain.append(idnames[idnames != idn])
        iscore.append([idmat[i, idnames!= idn]])
    
    scorenames1.append(itest)
    scorenames2.append(itrain)
    scores.append(iscore)
    allnames.append(idnames)

if '--realval' in sys.argv:
    simfile = sys.argv[sys.argv.index('--realval')+1]
    realnames, realsim = readinclustermat(simfile)


##### sort out pair, so that identity works on the same set of pairs!!! 

# Build 1d elements that can be tested! 
commonnames = np.sort(reduce(np.intersect1d, allnames))
print 'Commonnames', len(commonnames)


scoresorta = []
simisorta = []
for c, commonname in enumerate(commonnames):
    commontrain =[]
    commontrainscore = []
    for s, scorename in enumerate(scorenames1):
        for t, sconam in enumerate(scorename):
            if commonname in sconam:
                commontrain.append(list(scorenames2[s][t])) 
                commontrainscore.append(scores[s][t][0])
                break
    commontraincomb = reduce(np.intersect1d, commontrain)
    scsort = []
    sisort = []
    rj = list(realnames).index(commonname)
    for s, cscores in enumerate(commontrainscore):
        sc = []
        si = []
        for d, ctr in enumerate(commontraincomb):
            if ctr in commontrain[s] and ctr in realnames:
                sc.append(cscores[commontrain[s].index(ctr)])
                si.append(realsim[rj, list(realnames).index(ctr)])
        scsort.append(sc)
        sisort.append(si)
    scoresorta.append(scsort)
    simisorta.append(sisort)

scoresort = []
simisort = []    
simsims = []
scoresims = []
for s in range(len(scorenames1)):
    scsort = []
    sisort = []
    for t in range(len(scoresorta)):
        scsort.append(scoresorta[t][s])
        sisort.append(simisorta[t][s])
    scoresort.append(scsort)
    simisort.append(sisort)
    simsims.append(np.concatenate(sisort))
    scoresims.append(np.concatenate(scsort))

'''
scoresort = []
simisort = []
for s, scorename in enumerate(scorenames1):
    scsort = []
    sisort = []
    for su, subset in enumerate(scorename):
        trainsort = []
        for d, trsname in enumerate(scorenames2[s][su]): 
            rind2 = list(realnames).index(trsname)
            trainsort.append(rind2)
        for tp, tprot in enumerate(subset):
            if tprot in commonnames:
                csort = scores[s][su][tp]
                isort = realsim[list(realnames).index(tprot), trainsort]
                scsort.append(csort)
                sisort.append(isort)
    print len(scsort), len(scsort[0])
    print len(sisort), len(sisort[0])
    scoresort.append(scsort)
    simisort.append(sisort)

'''


if '--cutoff' in sys.argv:
    simcut = float(sys.argv[sys.argv.index('--cutoff')+1])
    cutout = '-'+str(simcut)
    for s, sisort in enumerate(simisort):
        for c, isort in enumerate(sisort):
            isort = np.array(isort)
            isort[isort >= simcut] = 1
            isort[isort < simcut] = 0
            simisort[s][c] = isort
 
else:
    cutout = ''
    for s, sisort in enumerate(simisort):
        for c, isort in enumerate(sisort):
            isort = np.array(isort)
            isort[isort < 0] = 0
            simisort[s][c] = isort
### Single PR, combined mean precision recall and individual   

# combined
precision = []
precvar = []
recall = []
xlabels = []
xticklocs = []
sortscores = []
sortsims = []
splitnum = 8
for n in range(numsets):
    scoresn = np.concatenate(scoresort[n])
    simn = np.concatenate(simisort[n])
    sortvalues = np.argsort(-scoresn)
    sortscores.append(scoresn[sortvalues])
    sortsims.append(simn[sortvalues])
    reca = np.cumsum(simn[sortvalues])/np.sum(simn)
    # all for xticks
    if scoresn[sortvalues][0] < 0.01 and scoresn[sortvalues][0] > 0.:
        factor10 = round(-np.log10(scoresn[sortvalues][-1]))
        factor10 = 10.**factor10
    else:
        factor10 = 1.
    xlabel = [round(factor10 * scoresn[sortvalues][0],3)]
    xtickl = [0]
    steps = 1./float(splitnum)
    for i in range(1,splitnum):
        pi = np.amin(np.where(reca >= i*steps)[0])
        print i*steps
        print reca[pi], scoresn[sortvalues][pi]
        
        xlabel.append(round(factor10 * scoresn[sortvalues][pi],3))
        xtickl.append(pi)
        
    #xlabel.append(round(factor10 * scores[n][sortvalues][-1],2))
    #xtickl.append(-1)
    print xlabel
    xlabels.append(xlabel)
    xticklocs.append(xtickl)
    precision.append(np.cumsum(simn[sortvalues])/np.arange(1, len(sortvalues)+1))
    if '--plotvar' in sys.argv:
        pvar = []
        sosim = simn[sortvalues]
        for i in range(len(sortvalues)):
            if i > 1:
                pchx = sosim[:i]
                pme = np.mean(pchx)
                pst = np.std(pchx)/np.sqrt(len(pchx))
                pvar.append([pme-pst, pme+pst])
        pvar = np.append([[precision[-1][0],precision[-1][0]],[precision[-1][1],precision[-1][1]]], pvar, axis = 0)
        precvar.append(pvar)
    recall.append(reca)
    

avsims = []
avscores = []
xs=[]
for n in range(numsets):
    ov = 300
    stsize = 20
    avsim = [np.mean(sortsims[n][:20])]
    varsim = []
    x = [recall[n][0]]
    avscore = [sortscores[n][0]]
    for s in range(1,ov/stsize):
        x.append(recall[n][s*stsize])
        avsim.append(np.mean(sortsims[n][:s*stsize]))
        avscore.append(sortscores[n][s*stsize])
    for s in range(len(recall[n])/stsize):
        if s > (len(recall[n])-ov)/stsize:
            avsim.append(np.mean(sortsims[n][-ov:]))
            x.append(recall[n][stsize*s])
            avscore.append(sortscores[n][stsize*s])
        elif s > ov/stsize:
            avsim.append(np.mean(sortsims[n][s*stsize-ov/2:s*stsize+ov/2]))
            x.append(recall[n][stsize*s])
            avscore.append(sortscores[n][stsize*s])
    xs.append(x)
    avsims.append(avsim)
    avscores.append(avscore)

if '--savestats' in sys.argv:
    for s, scorefile in enumerate(scorefiles):
        statname = scorefile+cutout+'-avprecicion.stats.txt'
        np.savetxt(statname, np.array([avscores[s], avsims[s]]).T)

if '--setnames' in sys.argv:
    scorefiles = sys.argv[sys.argv.index('--setnames')+1].split(',')

if '--savefig' in sys.argv:
    figname = sys.argv[sys.argv.index('--savefig')+1]
    




                  


reccuts = [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6]
print '\n# Recall 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6 Method'
for s, scf  in enumerate(scorefiles):
    recvals = []
    recscores = []
    for r, rec in enumerate(reccuts):
        wprec = np.where(precision[s]>rec)[0]
        if len(wprec) > 0:
            wprec = wprec[-1]
        else:
            wprec = 0
        recvals.append(recall[s][wprec])
        recscores.append(np.around(sortscores[s][wprec],4))
    print '\t'.join(np.around(np.array(recvals)*100., 1).astype(str)), '\t'+scf
    #print '%g\t,%g\t,%g\t,%g\t,%g\t,%g\t,%g' % recscores, '\t'+scf
    print '\t'.join(np.array(recscores).astype(str)), '\t'+scf

reccuts = [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5]
print '\n# Average_precision 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5 Method'
for s, scf  in enumerate(scorefiles):
    recvals = []
    recscores = []
    for r, rec in enumerate(reccuts):
        wprec = np.where(np.array(avsims[s])>rec)[0]
        if len(wprec) > 0:
            wprec = wprec[-1]
        else:
            wprec = 0
        recvals.append(xs[s][wprec])
        recscores.append(np.around(avscores[s][wprec],4))
    print '\t'.join(np.around(np.array(recvals)*100., 1).astype(str)), '\t'+scf
    #print '%g\t,%g\t,%g\t,%g\t,%g\t,%g\t,%g' % recscores, '\t'+scf
    print '\t'.join(np.array(recscores).astype(str)), '\t'+scf

    
if '--supplot' not in sys.argv:
    fig = plt.figure(figsize = (6,6), dpi = 100)
    ax = fig.add_subplot(111)
    for n in range(numsets):
        ax.plot(recall[n], precision[n], label = scorefiles[n], c = colors[n])
        if '--plotvar' in sys.argv:
            ax.fill_between(recall[n], precvar[n][:,0], precvar[n][:,1], alpha = 0.3, facecolor = colors[n])
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    #ax.spines['left'].set_color('dimgray')
    #ax.spines['bottom'].set_color('dimgray')
    #ax.tick_params(color = 'dimgray', width = 1.)
    ax.set_ylim([0,1.])
    
    ax2 = ax.twiny()
    ax2.spines['right'].set_visible(False)
    ax2.set_xlim(ax.get_xlim())
    ax2.spines["top"].set_position(("axes", 1.2))
    ax2.set_xticks(recall[0][xticklocs[0]])
    ax2.tick_params(color = colors[0], width = 2.)
    ax2.spines['top'].set_color(colors[0])
    ax2.spines['top'].set_linewidth(2.)
    ax2.set_xticklabels(np.absolute(np.around(xlabels[0], 3)), rotation = 60, color =colors[0], horizontalalignment = 'left')
    
    ax2b = ax.twiny()
    ax2b.spines['right'].set_visible(False)
    ax2b.set_xlim(ax.get_xlim())
    ax2b.spines["top"].set_position(("axes", 1.05))
    ax2b.set_xticks(recall[-1][xticklocs[1]])
    ax2b.tick_params(color = colors[-1], width = 2.)
    ax2b.spines['top'].set_color(colors[-1])
    ax2b.spines['top'].set_linewidth(2.)
    ax2b.set_xticklabels(np.array(xlabels[-1],dtype= int), rotation = 60, color =colors[-1], horizontalalignment = 'left')
    
    ax.set_xticks(np.linspace(0.,1., 6))
    ax.set_xticklabels([0,20,40,60,80,100])
    ax.grid(which = 'both', linestyle = '--')
    ax.set_ylabel('Fraction correct')
    ax.set_xlabel('Recall')
    ax.legend()




    fig4 = plt.figure(figsize = (6,6), dpi = 100)
    ax4 = fig4.add_subplot(111)
    for n in range(numsets):
        ax4.plot(xs[n], avsims[n], label = scorefiles[n], c = colors[n], marker = '.',)
    
    ax4.set_ylim([0,1.])
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.spines['left'].set_linewidth(1.5)
    ax4.spines['bottom'].set_linewidth(1.5)
    
    ax5 = ax4.twiny()
    ax5.spines['right'].set_visible(False)
    ax5.set_xlim(ax4.get_xlim())
    ax5.spines["top"].set_position(("axes", 1.2))
    ax5.set_xticks(recall[0][xticklocs[0]])
    ax5.tick_params(color = colors[0], width = 2.)
    ax5.spines['top'].set_color(colors[0])
    ax5.spines['top'].set_linewidth(2.)
    ax5.set_xticklabels(np.absolute(np.around(xlabels[0], 3)), rotation = 60, color =colors[0], horizontalalignment = 'left')
    
    ax5b = ax4.twiny()
    ax5b.spines['right'].set_visible(False)
    ax5b.set_xlim(ax4.get_xlim())
    ax5b.spines["top"].set_position(("axes", 1.05))
    ax5b.set_xticks(recall[-1][xticklocs[1]])
    ax5b.tick_params(color = colors[-1], width = 2.)
    ax5b.spines['top'].set_color(colors[-1])
    ax5b.spines['top'].set_linewidth(2.)
    ax5b.set_xticklabels(np.array(xlabels[-1],dtype= int), rotation = 60, color =colors[-1], horizontalalignment = 'left')
    
    ax4.set_xticks(np.linspace(0.,1., 6))
    ax4.set_xticklabels([0,20,40,60,80,100])
    ax4.grid(which = 'both', linestyle = '--')
    ax4.set_ylabel('Fraction correct')
    ax4.set_xlabel('Recall')
    
    ax4.legend()
    
    
    #ax4.set_xticks(np.linspace(0.,1.,9))
    #ax5 = ax4.twiny()
    #ax5.set_xlim(ax4.get_xlim())
    #ylim = ax4.get_ylim()
    #ax5.set_xticks(np.linspace(0.,1.,9))
    #ax4.set_xticklabels(np.array(xlabels[0]).astype(str), rotation = 90, color =colors[0])
    #ax5.set_xticklabels(np.array(xlabels[-1]).astype(str), rotation = 90, color =colors[numsets-1])
    #ax4.grid(which = 'both', linestyle = '--')
    #percentage = ['0', '25', '50', '75', '100']
    #for p, perc in enumerate(percentage):
        #ax4.annotate(perc, (float(p)/float(len(percentage)-1), ylim[0]+0.01))
    #ax4.legend()
    
    
    
    
    
    
    fig.tight_layout()
    fig4.tight_layout()
    if '--savefig' in sys.argv:
        fig4.savefig(os.path.splitext(figname)[0]+'-Meanprecision'+os.path.splitext(figname)[-1], dpi = 300)
        fig.savefig(os.path.splitext(figname)[0]+'-PR'+os.path.splitext(figname)[-1], dpi = 300)
    else:
        plt.show()





    if '--densities' in sys.argv:
        simcut = float(sys.argv[sys.argv.index('--densities')+1])
        simcut2 = float(sys.argv[sys.argv.index('--densities')+2])
        sortclasses = []
        for s, sisort in enumerate(simsims):
            sisort = np.array(sisort)
            clsort = np.zeros(len(sisort))
            clsort[sisort >= simcut] = 1
            #clsort[(sisort < simcut) & (sisort >= simcut2)] = 0
            clsort[sisort < simcut2] = -1
            sortclasses.append(clsort)

        bins = 41
        rows = int(np.sqrt(numsets))+int((numsets%int(np.sqrt(numsets)))>0)
        colums = int(numsets/rows)
        print colums, rows
        
        fd = 3.5
        
        figd = plt.figure('Density', figsize = (fd*colums/rows, fd))
        figd2 = plt.figure('Densitynonorm', figsize = (fd*colums/rows, fd))
        for n in range(numsets):
            axd = figd.add_subplot(rows, colums, n+1)
            axd2 = figd2.add_subplot(rows, colums, n+1)
            axd.spines['top'].set_visible(False)
            axd.spines['left'].set_visible(False)
            axd.spines['right'].set_visible(False)
            axd2.spines['top'].set_visible(False)
            axd2.spines['left'].set_visible(False)
            axd2.spines['right'].set_visible(False)
            #axd.spines['bottom'].set_linewidth(2.)
            #axd2.spines['bottom'].set_linewidth(2.)
            axd.tick_params(left = False, labelleft = False)#, width = 2.)
            axd2.tick_params(left = False, labelleft = False)#, width = 2.)
            als = [0.4,0.6, 0.7]
            linecol = ['grey', 'dimgrey', 'k'] #, colors[n]]
            iis = [0, -1, 1]
            for j in range(3):
                i = iis[j]
                print n, i, len(scoresims[n][sortclasses[n] == i]), np.mean(scoresims[n][sortclasses[n] == i]), np.std(scoresims[n][sortclasses[n] == i])
                if len(scoresims[n][sortclasses[n] == i]) > 0:
                    sns.distplot(scoresims[n][sortclasses[n] == i], kde = False, bins = bins, color = colors[n], hist_kws=dict(alpha=als[j], linestyle = '-', linewidth = 1), norm_hist = True, ax = axd)#, label = scorefiles[n]+' class('+str(i)+')')
                    sns.distplot(scoresims[n][sortclasses[n] == i], kde = False, bins = bins, color = linecol[j], hist_kws=dict(alpha=1., linestyle = '-',  linewidth = 2, histtype = 'step'),  norm_hist = True, ax = axd)#, label = scorefiles[n]+' class('+str(i)+')')
                    sns.distplot(scoresims[n][sortclasses[n] == i], kde = False, bins = bins, color = colors[n], hist_kws=dict(alpha=als[j], linestyle = '-', linewidth = 1), norm_hist = False, ax = axd2)#, label = scorefiles[n]+' class('+str(i)+')')
                    sns.distplot(scoresims[n][sortclasses[n] == i], kde = False, bins = bins, color = linecol[j], hist_kws=dict(alpha=1., linestyle = '-',  linewidth = 2, histtype = 'step'),  norm_hist = False, ax = axd2)#, label = scorefiles[n]+' class('+str(i)+')')
            #axd.legend()
            #axd2.legend()
            

        figav = plt.figure('Average', figsize = (fd*colums/rows, fd))
        fig2av = plt.figure('Averagenorm', figsize = (fd*colums/rows, fd))
        figdav = plt.figure('Average+density', figsize = (fd*colums/rows, fd))
        for n in range(numsets):
            
            axav = figav.add_subplot(rows, colums, n+1)
            ax2av = fig2av.add_subplot(rows, colums, n+1)
            
            axav.spines['top'].set_visible(False)
            axav.spines['right'].set_visible(False)
            ax2av.spines['top'].set_visible(False)
            ax2av.spines['right'].set_visible(False)
            
            #axav.spines['bottom'].set_linewidth(2.)
            #ax2av.spines['bottom'].set_linewidth(2.)
            #axav.spines['left'].set_linewidth(2.)
            #ax2av.spines['left'].set_linewidth(2.)
            
            #axav.tick_params(width = 2.)
            #ax2av.tick_params(width = 2.)
            axav.set_ylim([0.,1.02])
            ax2av.set_ylim([0.,1.02])
            
            axdav2 = figdav.add_subplot(rows, colums, n+1)
            axdav = axdav2.twinx()
            
            axdav.spines['top'].set_visible(False)
            axdav.spines['right'].set_visible(False)
            
            #axdav.spines['bottom'].set_linewidth(2.)
            #axdav.spines['left'].set_linewidth(2.)
            axdav.tick_params(labelleft = True, labelright = False, right = False, left = True)
            axdav.set_ylim([0.,1.04])
            
            
            axdav2.spines['top'].set_visible(False)
            axdav2.spines['right'].set_visible(False)
            axdav2.spines['left'].set_visible(False)
            axdav2.spines['bottom'].set_visible(False)
            
            axdav2.tick_params(left = False, right = False, labelleft = False,labelright = False)
            
            
            sns.distplot(scoresims[n][sortclasses[n] == -1], kde = False, bins = bins, color = colors[n], hist_kws=dict(alpha=als[1]*0.7, linestyle = '-', linewidth = 1), norm_hist = True, ax = axdav2)#, label = scorefiles[n]+' class('+str(i)+')')
            sns.distplot(scoresims[n][sortclasses[n] == -1], kde = False, bins = bins, color = linecol[1], hist_kws=dict(alpha=.7, linestyle = '-',  linewidth = 2, histtype = 'step'),  norm_hist = True, ax = axdav2)#, label = scorefiles[n]+' class('+str(i)+')')
                    
            
            als = [0.4,0.6, 0.7]
            linecol = ['white', 'grey', 'k'] #, colors[n]]
            iis = [0, 1]
            for j in range(2):
                i = iis[j]
                if len(scoresims[n][sortclasses[n] == i]) > 0:
                    xbarr = np.linspace(np.amin(scoresims[n]), np.amax(scoresims[n]), bins)
                    xplot = []
                    ybarr = []
                    ybarrnorm = []
                    for xi in range(len(xbarr)-1):
                        xplot.append((xbarr[xi]+xbarr[xi+1])/2.)
                        
                        ybarr.append(np.sum((sortclasses[n] >= i)[(scoresims[n] < xbarr[xi+1]) & (scoresims[n] > xbarr[xi])])/float(len(sortclasses[n][(scoresims[n] < xbarr[xi+1]) & (scoresims[n] > xbarr[xi])])))
                        
                        ybarrnorm.append((float(np.sum((sortclasses[n] >= i)[(scoresims[n] < xbarr[xi+1]) & (scoresims[n] > xbarr[xi])]))/float(np.sum(sortclasses[n] >= i)))/((np.sum((sortclasses[n] >= i)[(scoresims[n] < xbarr[xi+1]) & (scoresims[n] > xbarr[xi])])/float(np.sum(sortclasses[n] >= i)))+(np.sum((sortclasses[n] < i)[(scoresims[n] < xbarr[xi+1]) & (scoresims[n] > xbarr[xi])])/float(np.sum(sortclasses[n] < i)))))
                    
                    axav.plot(xplot, ybarr, linestyle = None, marker = 'o', color = colors[n], markeredgecolor = linecol[j+1],label = scorefiles[n])
                    ax2av.plot(xplot, ybarrnorm, linestyle = None, marker = 'o', color = colors[n], markeredgecolor = linecol[j+1], label = scorefiles[n])
                    sns.distplot(scoresims[n][sortclasses[n] == i], kde = False, bins = bins, color = colors[n], hist_kws=dict(alpha=als[j+1]*0.7, linestyle = '-', linewidth = 1), norm_hist = True, ax = axdav2)#, label = scorefiles[n]+' class('+str(i)+')')
                    sns.distplot(scoresims[n][sortclasses[n] == i], kde = False, bins = bins, color = linecol[j+1], hist_kws=dict(alpha=.7, linestyle = '-',  linewidth = 2, histtype = 'step'),  norm_hist = True, ax = axdav2)#, label = scorefiles[n]+' class('+str(i)+')')
                    axdav.plot(xplot, ybarr, linestyle = '-', color = 'k', linewidth = 2.)
                    axdav.plot(xplot, ybarr, linestyle = '', marker = 'o', markersize = 10., alpha = 1., color = colors[n], markeredgecolor = 'k',label = scorefiles[n])
                    
                    
            #ax2av.legend()
            if n == 0:
                axav.set_ylabel('Fraction correct')
                ax2av.set_ylabel('Fraction of fractions correct')
                axdav.yaxis.set_label_position("left")
                axdav.set_ylabel('Fraction correct')
            
        if '--savefig' in sys.argv:
            figav.savefig(os.path.splitext(figname)[0]+'-Scoreprobability'+os.path.splitext(figname)[-1], dpi = 300)
            figdav.savefig(os.path.splitext(figname)[0]+'-Scoredistribution_probability'+os.path.splitext(figname)[-1], dpi = 300)
            fig2av.savefig(os.path.splitext(figname)[0]+'-Scoreprobabilitynormedclass'+os.path.splitext(figname)[-1], dpi = 300)
            figd.savefig(os.path.splitext(figname)[0]+'-ScoreNormedDistribution'+os.path.splitext(figname)[-1], dpi = 300)
            figd2.savefig(os.path.splitext(figname)[0]+'-ScoreDistribution'+os.path.splitext(figname)[-1], dpi = 300)
        else:
            plt.show()




'''
# single
precision = []
recall = []
precisionvar = []

xlabels = []
xticks = []
sortscores = []
sortsims = []
for n in range(numsets):
    scoresn = scoresort[n]
    simn = simisort[n]
    precisiona = []
    recalla = []
    avnum = len(simn)
    for m in range(len(scoresn)):
        scoresnm = scoresn[m]
        simnm = simn[m]
        sortvalues = np.argsort(-scoresnm)
        sortscores.append(scoresnm[sortvalues])
        sortsims.append(simnm[sortvalues])
        # Compute precision recall
        precisiona.append(np.cumsum(simnm[sortvalues])/np.arange(1, len(sortvalues)+1))
        recalla.append(np.cumsum(simnm[sortvalues])/np.sum(simnm))
    recalla = np.concatenate(recalla)
    precisiona = np.concatenate(precisiona)
    recallsort = np.argsort(recalla)
    recalla = recalla[recallsort]
    precisiona = precisiona[recallsort]
    precision_mean = []
    recall_mean = []
    precision_var = []
    for i in range(len(precisiona)/10 - avnum/10):
        precision_mean.append(np.mean(precisiona[(i*10):i*10+avnum]))
        precision_var.append(np.std(precisiona[(i*10):i*10+avnum]))          
        recall_mean.append(np.mean(recalla[(i*10):i*10+avnum]))
    precision.append(precision_mean)
    recall.append(recall_mean)
    precisionvar.append(precision_var)

# NO xticks possible!!! What could be used, rank?     
## all for xticks
        #if scoresn[sortvalues][0] < 0.01:
            #factor10 = round(-np.log10(scoresn[sortvalues][-1]))
            #factor10 = 10.**factor10
        #else:
            #factor10 = 1.
        #xlabel = [round(factor10 * scoresn[sortvalues][0],2)]
        #for i in range(1,4):
            #steps = int(float(len(sortvalues))/4.)
            #xlabel.append(round(factor10 * scoresn[sortvalues][steps*i],2))
        #xlabel.append(round(factor10 * scoresn[sortvalues][-1],2))
        #print 'Xlabel', xlabel
        #xlabels.append(xlabel)

colors = cm.rainbow(np.arange(numsets)/float(numsets))

fig = plt.figure()
ax = fig.add_subplot(111)
for n in range(numsets):
    ax.scatter(recall[n], precision[n], label = scorefiles[n], c = colors[n])
    ax.plot(recall[n], precision[n], label = scorefiles[n], c = colors[n])
#ax.set_xticks(np.arange(0, 1.25, 0.25))
#ax2 = ax.twiny()
#ax2.set_xlim(ax.get_xlim())
#ylim = ax.get_ylim()
#ax2.set_xticks(np.arange(0, 1.25, 0.25))
#ax.set_xticklabels(np.array(xlabels[0]).astype(str), rotation = 90, color =colors[0])
#ax2.set_xticklabels(np.array(xlabels[1]).astype(str), rotation = 90, color =colors[1])
ax.grid(which = 'both', linestyle = '--')
percentage = ['0', '25', '50', '75', '100']
for p, perc in enumerate(percentage):
    ax.annotate(perc, (float(p)/float(len(percentage)-1), ylim[0]+0.01))

plt.show()

'''



