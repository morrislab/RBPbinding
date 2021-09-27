# precision_recall_single_prediction.py
import numpy as np
import sys, os 
from functools import reduce
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
import seaborn as sns
from scipy.stats import iqr, pearsonr, kde
import pandas as pd

def devi(xval, devtype):
    if devtype == 'iqr':
        q75, q25 = np.percentile(xval, [75 ,25])
        q75, q25 = q75-np.median(xval), np.median(xval)-q25
    elif devtype == 'std':
        q75 = q25 = np.std(xval)
    elif devtype == 'sterr':
        q75 = q25 = np.std(xval)/float(len(xval)-1.)
    return q75, q25

def avi(xval, avtype):
    if avtype == 'mean':
        return np.mean(xval)
    elif avtype == 'median':
        return np.median(xval)



def meanfunc(x, y, steps, steptype, overlap, avtype, stdtype):
    if steptype == 'x':
        stepvals = np.copy(x)
    elif steptype == 'y':
        stepvals = np.copy(y)
    else:
        stepvals = np.argsort(np.argsort(x))
    lim = [np.amin(stepvals), np.amax(stepvals)]
    boundaries = np.linspace(lim[0], lim[1], steps + overlap)
    xmean = []
    ymean = []
    xerr = []
    yerr = []
    j = 0
    for i in range(steps):
        #print j
        stepmask = (stepvals<boundaries[i+overlap]) & (stepvals>= boundaries[j])
        if len(np.where(stepmask)[0]) > 1:
            j = j+1
            xmean.append(avi(x[stepmask], avtype))
            xerr.append(devi(x[stepmask], stdtype))
            ymean.append(avi(y[stepmask], avtype))
            yerr.append(devi(y[stepmask], stdtype))
    return xmean, ymean, np.array(xerr).T, np.array(yerr).T


def predict(xin, params):
    xtransform = []
    for r in range(len(params)):
        xtransform.append(xin**(len(params)-1-r))
    xtransform = np.array(xtransform).T
    return np.sum(xtransform*params, axis = 1)
    
def interpolate(xmean, ymean, yerr, interpoltype, overnumber, interbounds):
    xmean = np.array(xmean)
    ymean = np.array(ymean)
    yerr = np.array(yerr)
    if interbounds is None:
        xlim = [np.amin(xmean), np.amax(xmean)]
        minloc = 0
        maxloc = len(xmean)
    else:
        minloc = np.where(ymean<=interbounds[0])[0]
        maxloc = np.where(ymean>=interbounds[1])[0]
        if len(minloc) > 0:
            minloc = minloc[-1]
        else:
            minloc = 0
        if len(maxloc) > 0:
            maxloc = maxloc[0]
        else:
            maxloc = len(xmean)-1
        xlim = [xmean[minloc], xmean[maxloc]]
    if maxloc - minloc < interpoltype:
        print maxloc, minloc
        return [],[],[[],[]], [np.zeros(interpoltype+1)]
    fitparams = []
    if overnumber is None:
        params = np.polyfit(xmean[minloc:maxloc+1],ymean[minloc:maxloc+1], interpoltype)
        fitparams.append(params)
        xinter = np.linspace(xmean[minloc],xmean[maxloc], 200)
        yinter = predict(xinter, params)
        yerrparams0 = np.polyfit(xmean[minloc:maxloc+1],yerr[0][minloc:maxloc+1], interpoltype)
        yerrparams1 = np.polyfit(xmean[minloc:maxloc+1],yerr[1][minloc:maxloc+1], interpoltype)
        yintererr = np.array([predict(xinter, yerrparams0), predict(xinter, yerrparams1)]).T
        
    else:
        xinterl = []
        yinterl = []
        yintererr0 = []    
        yintererr1 = []    
        for i in range(maxloc-minloc-overnumber+1):
            xmeanloc = xmean[i:i+overnumber]
            ymeanloc = ymean[i:i+overnumber]
            yerrloc0 = yerr[0][i:i+overnumber]
            yerrloc1 = yerr[1][i:i+overnumber]
            params = np.polyfit(xmeanloc, ymeanloc, interpoltype)
            fitparams.append(params)
            xinter = np.concatenate([np.linspace(xmeanloc[t],xmeanloc[t+1], 30) for t in range(overnumber-1)])
            xinterl.append(xinter)
            yinterl.append(predict(xinter, params))
            yerrparams0 = np.polyfit(xmeanloc, ymeanloc+yerrloc0, interpoltype)
            yerrparams1 = np.polyfit(xmeanloc, ymeanloc-yerrloc1, interpoltype)
            yintererr0.append(predict(xinter, yerrparams0)-yinterl[-1])
            yintererr1.append(yinterl[-1]-predict(xinter, yerrparams1))
        xinterl = np.array(xinterl)
        yinterl = np.array(yinterl)
        yintererr0 = np.array(yintererr0)
        yintererr1 = np.array(yintererr1)
        xinter = np.unique(np.concatenate(xinterl))
        yinter = []
        yintererr = []
        for i, xi in enumerate(xinter):
            yinter.append(np.mean(yinterl[xinterl == xi]))
            yintererr.append([np.mean(yintererr0[xinterl == xi]), np.mean(yintererr1[xinterl == xi])])
    return xinter, yinter, np.array(yintererr).T, fitparams


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


def zeropoint(slope, intercept):
    zpoint = -intercept/slope
    return zpoint

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
            outnames.append(tetp)
            outtops.append(mfac*tscore[t, trainmask])
    else:
        for nr in range(numcv):
            testp, trpn, tscore = opennpz(fstart+str(nr)+fend, dtype, row)
            for t, tetp in enumerate(testp):
                trainmask = trpn != tetp
                trainnames.append(trpn[trainmask])
                outnames.append(tetp)
                outtops.append(mfac*tscore[t, trainmask])
    outnames = np.array(outnames)
    #print outtops
    return outnames, np.array(trainnames), outtops
    



if __name__ == '__main__':
    
    
    numcv = int(sys.argv[1])
    scorename = sys.argv[2]
    scotype = sys.argv[3]
    scorerow = sys.argv[4]
    
    if scorename[0] == '-':
        scorename = scorename[1:]
        mfac = -1.
    else:
        mfac = 1.
    
    
    scorenamestart, scorenameend = scorename.split(',')
    
    scorenames1, scorenames2, scores = readinnpz(scorenamestart, scorenameend, numcv, scorerow, mfac, scotype)
    
    
    if '--training2test' in sys.argv:
        scorenamestart += '-train2test'
        nsconame2 = []
        nscores = []
        for s, sname in enumerate(scorenames1):
            sns = []
            ns = []
            for t, tname in enumerate(scorenames1):
                if t != s:
                    if sname in scorenames2[t]:
                        sns.append(tname)
                        ns.append(scores[t][list(scorenames2[t]).index(sname)])
            nsconame2.append(np.array(sns))
            nscores.append(np.array(ns))
        scorenames2 = nsconame2
        scores = nscores
            
    
    #scorenames1 = scorenames1[276:]
    #scorenames2 = scorenames2[276:]
    #scores = scores[276:]
    
    identityfile= sys.argv[5]
    idnames, idmat = readinclustermat(identityfile)
    if '>' in idnames[0]:
        for i, idname in enumerate(idnames):
            idnames[i] = idname[1:]

    simfile = sys.argv[6]
    realnames, realsim = readinclustermat(simfile)
    
    identities = []
    similarities = []
    
    for s, sname in enumerate(scorenames1):
        ids = []
        sims = []
        for t, sname2 in enumerate(scorenames2[s]):
            ids.append(idmat[list(idnames).index(sname), list(idnames).index(sname2)])
            sims.append(realsim[list(realnames).index(sname), list(realnames).index(sname2)])
        identities.append(ids)
        similarities.append(sims)
        

    '''
    xdata = np.linspace(0,1.,1000)
    ydata = np.zeros(1000)
    ydata[400:] += np.arange(600)/1000.
    ydata[600:] -= 3.*np.arange(400)**2/100000.
    ydata += (np.arange(1000)/2500.)*np.random.random(1000)
    #ydata[50:] += 0.1*np.sin(np.pi*np.arange(50)/12.)


    xmean, ymean, xerr, yerr = meanfunc(xdata, ydata, 10, 'x', 1, 'mean', 'std')
    xinter, yinter, yintererr = interpolate(xmean, ymean, yerr, 4, 5, None)
    #xinter2, yinter2, yintererr2 = interpolate(xmean, ymean, yerr, 1, None, [np.amin(ymean)+0.1,np.amax(ymean)-0.05])
    '''
    
    if '--outdir' in sys.argv:
        scorenamestart = sys.argv[sys.argv.index('--outdir')+1] + os.path.split(scorenamestart)[1]
    
    
    if '--individual' not in sys.argv:
        scorefull = np.concatenate(scores)
        identityfull = np.concatenate(identities)
        similarfull = np.concatenate(similarities)
        
        # change to mean and std, before was median and iqr
        
        xmean, ymean, xerr, yerr = meanfunc(scorefull, similarfull, 50, 'x', 1, 'mean', 'std')
        xinter, yinter, yintererr, yinterpars = interpolate(xmean, ymean, yerr, 2, 5, None)
        #xinter, yinter, yintererr = interpolate(xmean, ymean, yerr, 1, None, [np.amin(ymean), np.amax(ymean)])
        
        xmeanid, ymeanid, xerrid, yerrid = meanfunc(scorefull, identityfull, 50, 'x', 1, 'mean', 'std')
        xinterid, yinterid, yintererrid, yinterparsid = interpolate(xmeanid, ymeanid, yerrid, 2, 5, None)
        
        if '--reconstruction_cut' in sys.argv:
            cut = float(sys.argv[sys.argv.index('--reconstruction_cut')+1])
            xline = [np.amin(scorefull), np.amax(scorefull)]
            simline = [cut, cut]
            simcut = xinter[np.where(np.array(yinter) >= cut)[0][0]]
            print simcut
            idcut = yinterid[np.where(np.array(xinterid) >= simcut)[0][0]]
            print idcut
            idline = [idcut, idcut]
            
        # plot seqid and similarity against cosine distance
        fig = plt.figure(figsize = (4,4), dpi = 200)
        ax = fig.add_subplot(111)
        ax.scatter(scorefull, similarfull,c='limegreen', alpha = 0.05)
        ax.errorbar(xmean, ymean, yerr=[yerr[1], yerr[0]], xerr=[xerr[1], xerr[0]], fmt='o', c = 'black')
        ax.fill_between(xinter, yinter-yintererr[1], yinter+yintererr[0], color = 'grey', alpha = 0.5, linewidth = 0.)
        ax.plot(xinter, yinter, 'k-', linewidth = 3.)
        if '--reconstruction_cut' in sys.argv:
            ax.plot(xline, simline, 'grey', ls= '--')
            ax.text(simcut, cut+0.2, str(round(-simcut,2)), ha = 'right', va = 'bottom')
            ax.plot([simcut, simcut], [cut, cut+0.2], 'red', ls= '-')
            ax.text(xline[1]+0.02, cut, str(cut), ha = 'left', va = 'center')
        ax.set_xlabel('JPLE latent distance')
        ax.set_ylabel('RNA binding similarity R')
        ax.set_xticks([-1.5, -1, -0.5, 0.])
        ax.set_xticklabels([1.5, 1.0, 0.5, 0.])
        ax.set_yticks([-0.5,0., 0.5, 1.])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        fig2 = plt.figure(figsize = (4,4), dpi = 200)
        ax2 = fig2.add_subplot(111)
        ax2.scatter(scorefull, identityfull, c='cadetblue', alpha = 0.05)
        ax2.errorbar(xmeanid, ymeanid, yerr=[yerrid[1], yerrid[0]], xerr=[xerrid[1], xerrid[0]], fmt='o', c = 'black')
        ax2.fill_between(xinterid, yinterid-yintererrid[1], yinterid+yintererrid[0], color = 'grey', alpha = 0.5, linewidth = 0.)
        ax2.plot(xinterid, yinterid, 'k-', linewidth = 3.)
        if '--reconstruction_cut' in sys.argv:
            ax2.plot(xline, idline, 'grey', ls = '--')
            ax2.text(xline[1]+0.02, idcut, str(int(idcut)), ha = 'left', va = 'center')
        ax2.set_xlabel('JPLE latent distance')
        ax2.set_ylabel('Sequence identity')
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)    
        ax2.set_xticks([-1.5, -1, -0.5, 0.])
        ax2.set_xticklabels([1.5, 1.0, 0.5, 0.])
        
        fig3 = plt.figure(figsize = (5,4), dpi = 200)
        ax3 = fig3.add_subplot(111)
        ax4 = ax3.twinx()
        ax3.errorbar(xmeanid, ymeanid, yerr=[yerrid[1], yerrid[0]], xerr=[xerrid[1], xerrid[0]], fmt='o', c = 'cadetblue')
        ax3.fill_between(xinterid, yinterid-yintererrid[1], yinterid+yintererrid[0], color = 'cadetblue', alpha = 0.5, linewidth = 0.)
        ax3.plot(xinterid, yinterid, 'k-', linewidth = 3.)
        ax3.set_xlabel('JPLE latent distance')
        ax3.set_ylabel('Sequence identity')
        ax3.spines['right'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax3.spines['left'].set_color('cadetblue')
        ax3.spines['left'].set_linewidth(2)
        ax3.set_xticks([-1.5, -1, -0.5, 0.])
        ax3.set_xticklabels([1.5, 1.0, 0.5, 0.])
        ax4.set_yticks([-0.5,0., 0.5, 1.])
        
        if '--reconstruction_cut' in sys.argv:
            ax3.plot([xline[0], simcut], idline, 'red', ls = '--')
            ax3.text(xline[0], idcut, str(int(idcut)), ha = 'center', va = 'bottom')

        ax4.errorbar(xmean, ymean, yerr=[yerr[1], yerr[0]], xerr=[xerr[1], xerr[0]], fmt='o', c = 'limegreen')
        ax4.fill_between(xinter, yinter-yintererr[1], yinter+yintererr[0], color = 'limegreen', alpha = 0.5, linewidth = 0.)
        ax4.plot(xinter, yinter, 'k-', linewidth = 3.)
        ax4.set_ylabel('RNA binding similarity R')
        ax4.spines['left'].set_visible(False)
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_color('limegreen')
        ax4.spines['right'].set_linewidth(2.)
        ax4.spines['bottom'].set_linewidth(2.)
        ax4.set_xticks([-1.5, -1, -0.5, 0.])
        ax4.set_xticklabels([1.5, 1.0, 0.5, 0.])

        if '--reconstruction_cut' in sys.argv:
            ax4.plot([simcut, xline[1]], simline, 'red', ls= '--')
            ax4.text(simcut, cut+0.25, str(round(-simcut,2)), ha = 'right', va = 'bottom')
            ax4.plot([simcut, simcut], [np.amin(ymean), cut+0.25], 'red', ls= '-')
            ax4.text(xline[1], cut, str(cut), ha = 'center', va = 'bottom')        

        fig.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        if '--savefig' in sys.argv:
            print scorenamestart+'cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]
            fig.savefig(scorenamestart+'cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-motifsim-calibration.jpg', dpi = 300, bbox_inches = 'tight')
            fig2.savefig(scorenamestart+'cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-seqid-calibration.jpg', dpi = 300, bbox_inches = 'tight')
            fig3.savefig(scorenamestart+'cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-interpolation-calibration.jpg', dpi = 300, bbox_inches = 'tight')
        else:
            plt.show()


    if '--individual' in sys.argv:
        # import motif, name and motif clusters
        if '--savefig' in sys.argv or '--showindividual' in sys.argv:
            motifs = np.genfromtxt(sys.argv[sys.argv.index('--individual')+1], dtype = str)
            protnames = np.genfromtxt(sys.argv[sys.argv.index('--individual')+2], dtype = str)
            domaintypes = np.genfromtxt(sys.argv[sys.argv.index('--individual')+3], dtype = str)
            
            #print motifs[0], protnames[0], domaintypes[0]
            
            titles = []
            for s, scorename in enumerate(scorenames1):
                titles.append(protnames[list(protnames[:,1]).index(scorename), 2]+'('+ protnames[list(protnames[:,1]).index(scorename), 3][0]+protnames[list(protnames[:,1]).index(scorename), 3].split('_')[1][0]+', '+ motifs[list(motifs[:,0]).index(scorename), 1] +', '+domaintypes[list(domaintypes[:,1]).index(scorename), 2]+')')
        
        # one figure with all in there
        
        # plot 3d distribution of slopes for id and sim
            # add single distribution of slopes for identity and similarity
            # find clusters in this 2d distribution and determine motif clusters
        
        # one figure with clusters sorted by slope
            # for identities and similarity developments
        
        xmeans = []
        ymeans = []
        xerrs = []
        yerrs = []
        
        xmeansid = []
        ymeansid = []
        xerrsid = []
        yerrsid = []
        
        xinters = []
        yinters = []
        yintererrs = []
        slopes = []
        fits = []
        
        xintersid = []
        yintersid = []
        yintererrsid = []
        slopesid = []
        fitsid = []
        
        for s, sname in enumerate(scorenames1):
            print s, sname
            scores[s] = np.array(scores[s])
            similarities[s] = np.array(similarities[s])
            identities[s] = np.array(identities[s])

            xmean, ymean, xerr, yerr = meanfunc(scores[s], similarities[s], 20, 'x', 1, 'mean', 'iqr')
            xmeans.append(xmean)
            ymeans.append(ymean)
            xerrs.append(xerr)
            yerrs.append(yerr)
            
            xinter, yinter, yintererr, yinterpars = interpolate(xmean, ymean, yerr, 1, None, [0, 1.])
            xinters.append([-2.,0.])
            yinters.append(predict(np.array([-2,0.]), yinterpars[0]))
            #yintererrs.append(yintererr)
            slopes.append([yinterpars[0][0], yinterpars[0][1], zeropoint(yinterpars[0][0],yinterpars[0][1])])
            fits.append(pearsonr(scores[s], similarities[s])[0])
            
            
            xmeanid, ymeanid, xerrid, yerrid = meanfunc(scores[s], identities[s], 20, 'x', 1, 'mean', 'iqr')
            xmeansid.append(xmeanid)
            ymeansid.append(ymeanid)
            xerrsid.append(xerrid)
            yerrsid.append(yerrid)
            
            xinterid, yinterid, yintererrid, yinterparsid = interpolate(xmeanid, ymeanid, yerrid, 1, None, [20,100.])
            xintersid.append([-2,0])
            yintersid.append(predict(np.array([-2,0.]), yinterparsid[0]))
            #yintererrsid.append(yintererrid)
            slopesid.append([yinterparsid[0][0], yinterparsid[0][1], zeropoint(yinterparsid[0][0],yinterparsid[0][1])])            
            fitsid.append(pearsonr(scores[s], identities[s])[0])
            
            if '--savefig' in sys.argv or '--showindividual' in sys.argv:
                
                fig3 = plt.figure(figsize = (5,4))
                ax3 = fig3.add_subplot(111)
                ax3.set_title(titles[s])
                ax4 = ax3.twinx()
                
                ax3.scatter(scores[s], identities[s], c = 'cadetblue', alpha = 0.4, label = 'R = '+np.around(fitsid[s],2).astype(str))
                ax3.errorbar(xmeanid, ymeanid, yerr=[yerrid[1], yerrid[0]], xerr=[xerrid[1], xerrid[0]], fmt='o', c = 'cadetblue', mec = 'k', ecolor = 'k')
                if len(xinterid) > 0:
                    ax3.fill_between(xinterid, yinterid-yintererrid[1], yinterid+yintererrid[0], color = 'cadetblue', alpha = 0.5, linewidth = 0. ,label='ID='+str(np.around(slopesid[s][0],1))+'X+'+str(np.around(slopesid[s][1],1)))
                    ax3.plot(xinterid, yinterid, 'k-', linewidth = 3.) 
                ax3.set_xlabel('JPLE latent distance')
                ax3.set_ylabel('Sequence identity')
                ax3.spines['right'].set_visible(False)
                ax3.spines['top'].set_visible(False)
                ax3.spines['left'].set_color('cadetblue')
                ax3.spines['left'].set_linewidth(2)
                ax3.set_ylim([0.,100.])
                ax3.set_xlim([np.amin(scores),0.])
                ax4.set_xlim([np.amin(scores),0.])
                ax3.set_xticks([-1.5, -1, -0.5, 0.])
                ax3.set_xticklabels([1.5, 1.0, 0.5, 0.])
                ax4.set_yticks([-0.5,0., 0.5, 1.])
                ax4.scatter(scores[s], similarities[s], c = 'limegreen', alpha = 0.4,label = 'R = '+np.around(fits[s],2).astype(str))
                ax4.errorbar(xmean, ymean, yerr=[yerr[1], yerr[0]], xerr=[xerr[1], xerr[0]], fmt='o', c = 'limegreen', mec = 'k', ecolor = 'k')
                if len(xinter) > 0:
                    ax4.fill_between(xinter, yinter-yintererr[1], yinter+yintererr[0], color = 'limegreen', alpha = 0.5, linewidth = 0.,label='SIM='+str(np.around(slopes[s][0],2))+'X+'+str(np.around(slopes[s][1],2)))
                    ax4.plot(xinter, yinter, 'k-', linewidth = 3.)
                ax4.set_ylabel('Specificity similarity R')
                ax4.spines['left'].set_visible(False)
                ax4.spines['top'].set_visible(False)
                ax4.spines['right'].set_color('limegreen')
                ax4.spines['right'].set_linewidth(2.)
                ax4.spines['bottom'].set_linewidth(2.)
                ax4.set_ylim([-0.5,1.])
                ax4.legend(loc = 'upper left', bbox_to_anchor=(0.02, 0.98))
                ax3.legend(loc = 'upper left', bbox_to_anchor=(0.02, 0.78))
                #fig3.tight_layout()
            if '--savefig' in sys.argv:
                fig3.savefig(scorenamestart+'_'+sname+'cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-calibration.jpg', dpi = 200, bbox_inches = 'tight')
            elif '--showindividual' in sys.argv:
                plt.show()
            plt.close()

        
        
        np.savez_compressed(scorenamestart+'cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'_calibration-data.npz', names = scorenames1, slopes = slopes, slopesid = slopesid, xinters = xinters, yinters = yinters, yintersid = yintersid, xintersid = xintersid, fits = fits, fitsid = fitsid)
        
            
        figsim = plt.figure(figsize = (4,3))
        axsim = figsim.add_subplot(111)
        xs = []
        ys = []
        for a in range(len(scorenames1)):
            axsim.plot(xinters[a], yinters[a], color='limegreen', alpha = 0.1, linewidth = 4)
            #axsim.plot(xinters[a], yinters[a], color='k', alpha = 0.05, linewidth = 4)
        axsim.set_ylim([-1,1.])
        axsim.spines['right'].set_visible(False)
        axsim.spines['top'].set_visible(False)
        axsim.set_xlabel('JPLE latent distance')
        axsim.set_ylabel('Specificity similarity R')
        axsim.set_xticks([-1.5, -1, -0.5, 0.])
        axsim.set_xticklabels([1.5, 1.0, 0.5, 0.])
        axsim.set_yticks([-0.5,0., 0.5, 1.])
            #xs.append(np.random.uniform(-2,0,60))
            #ys.append(xs[-1]*(yinters[a][1]-yinters[a][0])/(xinters[a][1]-xinters[a][0]) + yinters[a][1])
        # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
        #xs = np.concatenate(xs)
        #ys = np.concatenate(ys)
        #k = kde.gaussian_kde(np.array([xs, ys]))
        #nbins = 40
        #xi, yi = np.mgrid[xs.min():xs.max():nbins*1j, ys.min():ys.max():nbins*1j]
        #zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        #axsim.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.Greys)
        #axsim.contour(xi, yi, zi.reshape(xi.shape))
        

        
        
        figid = plt.figure(figsize = (4,3))
        axid = figid.add_subplot(111)
        for a in range(len(scorenames1)):
            axid.plot(xintersid[a], yintersid[a], color='cadetblue', alpha = 0.1, linewidth = 4.)
        axid.set_ylim([0,100.])
        axid.spines['right'].set_visible(False)
        axid.spines['top'].set_visible(False)
        axid.set_xlabel('JPLE latent distance')
        axid.set_ylabel('Sequence identity')
        axid.set_xticks([-1.5, -1, -0.5, 0.])
        axid.set_xticklabels([1.5, 1.0, 0.5, 0.])

        
        figslopes = plt.figure(figsize = (5,5))
        axslope = figslopes.add_subplot(111)
        axslope.scatter(np.array(slopes)[:,0], np.array(slopesid)[:,0], c = 'dimgrey')
        axslope.set_xlabel('Slopes motif similarity')
        axslope.set_ylabel('Slopes Sequence identity')
        k = kde.gaussian_kde(np.array([np.array(slopes)[:,0], np.array(slopesid)[:,0]]))
        nbins = 20
        xi, yi = np.mgrid[np.amin(np.array(slopes)[:,0]):np.amax(np.array(slopes)[:,0]):nbins*1j, np.amin(np.array(slopesid)[:,0]):np.amax(np.array(slopesid)[:,0]):nbins*1j]
        zi = np.sqrt(k(np.vstack([xi.flatten(), yi.flatten()])))
        axslope.contour(xi, yi, zi.reshape(xi.shape), colors = 'k')
 
        
        figr = plt.figure(figsize = (6,6))
        axr = figr.add_subplot(111)
        axr.scatter(fits,fitsid, c = 'dimgrey')
        axr.set_xlabel('R motif similarity')
        axr.set_ylabel('R Sequence identity')
        fits = np.nan_to_num(fits)
        fitsid = np.nan_to_num(fitsid)
        k = kde.gaussian_kde(np.array([fits, fitsid]))
        nbins = 20
        xi, yi = np.mgrid[np.amin(fits):np.amax(fits):nbins*1j, np.amin(fitsid):np.amax(fitsid):nbins*1j]
        zi = np.sqrt(k(np.vstack([xi.flatten(), yi.flatten()])))
        axr.contour(xi, yi, zi.reshape(xi.shape), colors = 'k')
        
    
        print scorenamestart+'_cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-calibration-fitindividual_seqid.jpg'
        figid.savefig(scorenamestart+'_cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-calibration-fitindividual_seqid.jpg', dpi = 300, bbox_inches = 'tight')
        figsim.savefig(scorenamestart+'_cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-calibration-fitindividual_motifsim.jpg', dpi = 300, bbox_inches = 'tight')
        figr.savefig(scorenamestart+'_cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-calibration-scatterRseqvsmotif.jpg', dpi = 300, bbox_inches = 'tight')
        figslopes.savefig(scorenamestart+'_cv'+str(numcv)+'-'+scotype+'-'+scorerow+'_'+os.path.splitext(os.path.split(simfile)[1])[0]+'-calibration-scatterslopesseqvsmotif.jpg', dpi = 300, bbox_inches = 'tight')
        plt.show()
 










