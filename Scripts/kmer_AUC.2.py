#plot2d.py
import numpy as np
import sys,os
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62
from Bio import pairwise2
import operator
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
import matplotlib
if "--savefig" in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import glob
import seaborn as sns
from scipy.stats import wilcoxon
import pandas as pd
import scipy.stats
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import gzip
from matplotlib import gridspec

def plotsecondary(secanno, seqdistance, scorenames, seqscoreA, seqscoreB, sequ, distcut):
    old_sec = 'Z'
    old_pos = 0
    x = np.arange(len(secanno))
    xlim = [0, len(secanno)]

    patches = []
    colors = []
    locations = []
    types = []
    for s, seca in enumerate(np.append(secanno, ['Z'])):
        if seca != old_sec:
            if old_sec == 'E':
                #Strand
                arrow = mpatches.FancyArrow(old_pos-0.5, 0.5, s-old_pos, 0, width=0.7, length_includes_head=True, head_width=1, head_length=min(3,s-old_pos), shape='full', overhang=0, head_starts_at_zero=False)
                patches.append(arrow)
                colors.append('gold')
            if old_sec == 'B':
                # Isolated beta-bridge residue
                arrow = mpatches.FancyArrow(old_pos-0.5, 0.5, s-old_pos, 0, width=1, length_includes_head=True, head_width=1, head_length=1, shape='full', overhang=0, head_starts_at_zero=False)
                patches.append(arrow)
                colors.append('orange')
            if old_sec == 'C':
                sp = np.linspace(old_pos-.5, s-.5, (s-old_pos)*10)
                tp = np.ones(len(sp))*0.5
                locations.append([sp, tp])
                types.append('k')
                #p = ax.plot(sp, tp, color = 'k', linewidth = 3.)
            if old_sec == 'X':
                sp = np.linspace(old_pos-.5, s-.5, (s-old_pos)*10)
                tp = np.ones(len(sp))*0.5
                locations.append([sp, tp])
                types.append('grey')
                #p = ax.plot(sp, tp, color = 'grey', linewidth = 1.)
            if old_sec == 'S':
                # Bend
                #sp = np.linspace(old_pos-.5, s-.5, (s-old_pos)*10)
                #tp = np.ones(len(sp))*0.5
                sp = np.linspace(old_pos-.5, s-.5, (s-old_pos)*10)
                tp = 0.5-np.sin(np.pi*(sp-old_pos+0.5)/np.amax(sp-old_pos+0.5))/8.
                locations.append([sp, tp])
                types.append('maroon')
                #p = ax.plot(sp, tp, color = 'cyan', linewidth = 2.)
            if old_sec == 'H':
                #Alpha helix (4-12)
                sp = np.linspace(old_pos-.5, s-.5, (s-old_pos)*10)
                tp = 0.5+np.sin(np.pi*(sp-old_pos+0.5))/2.
                locations.append([sp, tp])
                types.append('steelblue')                
                #p = ax.plot(sp, tp, color = 'blue', linewidth = 2.)
            if old_sec == 'G':
                #3-10 helix
                sp = np.linspace(old_pos-.5, s-.5, (s-old_pos)*10)
                tp = 0.5+np.sin(np.pi*(sp-old_pos+0.5))/2.
                locations.append([sp, tp])
                types.append('turquoise')                
                #p = ax.plot(sp, tp, color = 'red', linewidth = 2.)            
            if old_sec == 'I':
                #Pi helix 
                sp = np.linspace(old_pos-.5, s-.5, (s-old_pos)*10)
                tp = 0.5+np.sin(np.pi*(sp-old_pos+0.5))/2.
                locations.append([sp, tp])
                types.append(['midnightblue'])
                #p = ax.plot(sp, tp, color = 'orange', linewidth = 2.)
            if old_sec == 'T':
                # Turn
                sp = np.linspace(old_pos-.5, s-.5, (s-old_pos)*10)
                tp = 0.5+np.sin(np.pi*(sp-old_pos+0.5)/np.amax(sp-old_pos+0.5))/4.
                locations.append([sp, tp])
                types.append('mediumvioletred')
                #p = ax.plot(sp, tp, color = 'green', linewidth = 2.)
            old_sec = seca
            old_pos = s
    
    fig = plt.figure(figsize=(8.7,3.5), dpi = 200)
    gs = gridspec.GridSpec(3, 1, height_ratios= [3,3,1]) #width_ratios=[3, 1]) 
    if len(patches) != 0:
        ax = fig.add_subplot(gs[-1])
        for s in range(len(locations)):
            sp = locations[s][0]
            tp = locations[s][1]
            ckol = types[s]
            p = ax.plot(sp, tp, color = ckol, linewidth = 2.)
        ax.set_xticks(x)
        print len(sequ),6./(float(len(sequ))/85.)
        ax.set_xticklabels(list(sequ), fontsize=max(2,6./(float(len(sequ))/85.)))
        ax.tick_params(axis='y', right=False, left = False, labelleft = False, labelright = False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='x', top=False, bottom = False, labeltop = False, labelbottom = True)
        collection = PatchCollection(patches, color = colors, alpha=0.7)
        ax.add_collection(collection)
        ax.set_xlim(xlim[0], xlim[1])
        ax2 = fig.add_subplot(gs[0])
        ax2b = fig.add_subplot(gs[1])
    else:
        ax2 = fig.add_subplot(gs[0])
        ax2b = fig.add_subplot(gs[1])    
    
    
    #ax2.set_ylim(ylims[0], ylims[1]*1.05)
    ax2.set_xlim(xlim[0], xlim[1])
    ax2b.set_xlim(xlim[0], xlim[1])
    
    ax2.set_ylim([0,np.amax(seqscoreA)*1.3])
    ax2b.set_ylim([0,np.amax(seqscoreB)*1.3])
    ax2.set_yticks([np.around(np.amin(seqscoreA),1),np.around(np.amax(seqscoreA),1)])
    ax2b.set_yticks([np.around(np.amin(seqscoreB),1),np.around(np.amax(seqscoreB),1)])
    
    distcolor = 'firebrick'
    scorecolor = 'darkslategrey'
    scorecolor2 = 'slateblue'
    
    ax2.bar(x, seqscoreA, width = 1, color = scorecolor, alpha = 0.6)
    ax2b.bar(x, seqscoreB, width = 1, color = scorecolor2, alpha = 0.8)
    
    if len(patches) == 0:
        ax3 = ax2b.twiny()
        ax3.set_xticks(x)
        ax3.set_xticklabels(list(sequ), fontsize=9)
        ax3.set_xlim(ax2b.get_xlim())
        ax3.tick_params(axis='x', bottom = True, top = False, labelbottom = True, labeltop = False)
    ax2.tick_params(axis='x', bottom = False, top = False, labelbottom = False, labeltop = True)
    ax2b.tick_params(axis='x', bottom = False, top = False, labelbottom = False, labeltop = False)
    
    ax4 = ax2.twinx()
    ax4.plot(x, -seqdistance, color = distcolor)
    
    ax4b = ax2b.twinx()
    ax4b.plot(x, -seqdistance, color = distcolor)
    
    ax4.scatter(x[seqdistance < distcut], -seqdistance[seqdistance < distcut], color = distcolor)
    ax4.set_yticks(-np.arange(0,20,5))
    ax4.set_yticklabels(np.arange(0,20,5))
    ax4b.scatter(x[seqdistance < distcut], -seqdistance[seqdistance < distcut], color = distcolor)
    ax4b.set_yticks(-np.arange(0,20,5))
    ax4b.set_yticklabels(np.arange(0,20,5))
    
    ax4.spines['bottom'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.spines['left'].set_visible(False)
    ax4.spines['right'].set_color(distcolor)
    ax4.spines['right'].set_linewidth(2)
    
    ax4b.spines['bottom'].set_visible(False)
    ax4b.spines['top'].set_visible(False)
    ax4b.spines['left'].set_visible(False)
    ax4b.spines['right'].set_color(distcolor)
    ax4b.spines['right'].set_linewidth(2)
    
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_color(scorecolor)
    ax2.spines['left'].set_linewidth(2)
    
    ax2b.spines['bottom'].set_visible(False)
    ax2b.spines['top'].set_visible(False)
    ax2b.spines['left'].set_color(scorecolor2)
    ax2b.spines['left'].set_linewidth(2)
    
    ax2.set_ylabel(scorenames[0])
    ax4.set_ylabel('RNA distance (A)')
    
    ax2b.set_ylabel(scorenames[1])
    #ax4b.set_ylabel('RNA distance (A)')
    
    ax4.yaxis.label.set_color(distcolor)
    ax2.yaxis.label.set_color(scorecolor)
    
    ax4.tick_params(axis='y', colors=distcolor)
    ax2.tick_params(axis='y', colors=scorecolor)
    
    ax4b.yaxis.label.set_color(distcolor)
    ax2b.yaxis.label.set_color(scorecolor2)
    
    ax4b.tick_params(axis='y', colors=distcolor)
    ax2b.tick_params(axis='y', colors=scorecolor2)
    return fig







# AUC comparison adapted from
# https://github.com/Netflix/vmaf/

def compute_midrank(x):
    """Computes midranks.
    Args:
    x - a 1D numpy array
    Returns:
    array of midranks
    """
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=np.float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5*(i + j - 1)
        i = j
    T2 = np.empty(N, dtype=np.float)
    # Note(kazeevn) +1 is due to Python using 0-based indexing
    # instead of 1-based in the AUC formula in the paper
    T2[J] = T + 1
    return T2


def compute_midrank_weight(x, sample_weight):
    """Computes midranks.
    Args:
    x - a 1D numpy array
    Returns:
    array of midranks
    """
    J = np.argsort(x)
    Z = x[J]
    cumulative_weight = np.cumsum(sample_weight[J])
    N = len(x)
    T = np.zeros(N, dtype=np.float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = cumulative_weight[i:j].mean()
        i = j
    T2 = np.empty(N, dtype=np.float)
    T2[J] = T
    return T2


def fastDeLong(predictions_sorted_transposed, label_1_count):
    """
    The fast version of DeLong's method for computing the covariance of
    unadjusted AUC.
    Args:
    predictions_sorted_transposed: a 2D numpy.array[n_classifiers, n_examples]
        sorted such as the examples with label "1" are first
    Returns:
    (AUC value, DeLong covariance)
    Reference:
    @article{sun2014fast,
    title={Fast Implementation of DeLong's Algorithm for
            Comparing the Areas Under Correlated Receiver Oerating Characteristic Curves},
    author={Xu Sun and Weichao Xu},
    journal={IEEE Signal Processing Letters},
    volume={21},
    number={11},
    pages={1389--1393},
    year={2014},
    publisher={IEEE}
    }
    """
    # Short variables are named as they are in the paper
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    positive_examples = predictions_sorted_transposed[:, :m]
    negative_examples = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    tx = np.empty([k, m], dtype=np.float)
    ty = np.empty([k, n], dtype=np.float)
    tz = np.empty([k, m + n], dtype=np.float)
    for r in range(k):
        tx[r, :] = compute_midrank(positive_examples[r, :])
        ty[r, :] = compute_midrank(negative_examples[r, :])
        tz[r, :] = compute_midrank(predictions_sorted_transposed[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - float(m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def calc_pvalue(aucs, sigma):
    """Computes log(10) of p-values.
    Args:
    aucs: 1D array of AUCs
    sigma: AUC DeLong covariances
    Returns:
    log10(pvalue)
    """
    l = np.array([[1, -1]])
    z = np.abs(np.diff(aucs)) / np.sqrt(np.dot(np.dot(l, sigma), l.T))
    return np.log10(2) + scipy.stats.norm.logsf(z, loc=0, scale=1) / np.log(10)


def compute_ground_truth_statistics(ground_truth, sample_weight=None):
    assert np.array_equal(np.unique(ground_truth), [0, 1])
    order = (-ground_truth).argsort()
    label_1_count = int(ground_truth.sum())
    if sample_weight is None:
        ordered_sample_weight = None
    else:
        ordered_sample_weight = sample_weight[order]

    return order, label_1_count, ordered_sample_weight


def delong_roc_variance(ground_truth, predictions):
    """
    Computes ROC AUC variance for a single set of predictions
    Args:
    ground_truth: np.array of 0 and 1
    predictions: np.array of floats of the probability of being class 1
    """
    sample_weight = None
    order, label_1_count, ordered_sample_weight = compute_ground_truth_statistics(
        ground_truth, sample_weight)
    predictions_sorted_transposed = predictions[np.newaxis, order]
    aucs, delongcov = fastDeLong(predictions_sorted_transposed, label_1_count, ordered_sample_weight)
    assert len(aucs) == 1, "There is a bug in the code, please forward this to the developers"
    return aucs[0], delongcov


def delong_roc_test(ground_truth, predictions_one, predictions_two):
    """
    Computes log(p-value) for hypothesis that two ROC AUCs are different
    Args:
    ground_truth: np.array of 0 and 1
    predictions_one: predictions of the first model,
        np.array of floats of the probability of being class 1
    predictions_two: predictions of the second model,
        np.array of floats of the probability of being class 1
    """
    sample_weight = None
    order, label_1_count, ordered_sample_weight = compute_ground_truth_statistics(ground_truth)
    predictions_sorted_transposed = np.vstack((predictions_one, predictions_two))[:, order]
    aucs, delongcov = fastDeLong(predictions_sorted_transposed, label_1_count)
    return calc_pvalue(aucs, delongcov)





def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False



def readinscore(ifile, sline, info = False):
    ifile = open(ifile, 'r').readlines()
    names = []
    scores = []
    sequences = []
    intinfo = []
    for i, line in enumerate(ifile):
        if line[0] == '>':
            names.append(line[1:].split(':')[0].strip())
            if info:
                intinfo.append(line[1:].strip().split(':'))
            sequence = ifile[i+1].strip()
            score = ifile[i+sline].strip()
            if ',' in score:
                score = np.array(score.split(','))
            else:
                score = np.array(list(score))
            if 'X' in score:
                sequence = ''.join(np.array(list(sequence))[score != 'X'])
                score = score[score != 'X']
            if '?' in score:
                score[score == '?'] = '0'
            if isfloat(score[0]):
                scores.append(np.array(score,dtype = float))
            else:
                scores.append(np.array(score))
            sequences.append(sequence)
    if info:
        return np.array(names), np.array(sequences), np.array(scores), np.array(intinfo)
    else:
        return np.array(names), np.array(sequences), np.array(scores)
    
def norm01(scores):
    return (scores-np.amin(scores))/(np.amax(scores)-np.amin(scores))



if __name__ == '__main__':

    maskfile = sys.argv[1]
    scorefile = sys.argv[2]
    scorefile2 = sys.argv[3]
    scorenames = sys.argv[4].split(',')
    
    maskname, maskseq, interface, interfaceinfo = readinscore(maskfile, 2, info = True)
    inames1, iseqs1, interscore1 = readinscore(scorefile, 2)
    inames2, iseqs2, interscore2 = readinscore(scorefile2, 2)
    
    fpr = []
    tpr = []
    roc_auc = []
    rocstats = []
    totbeat = 0
    
    if '--combine' in sys.argv and len(maskname) > len(np.unique(maskname)):
        print 'Combine'
    ### align sequences and take mean 
        mylen = np.vectorize(len)
        uniquename, unicn = np.unique(maskname, return_index = True)
        uniinfo = []
        uniseqs = []
        uniinterface = []
        uniscore1 = []
        uniscore2 = []
        combinemas = []
        for u, uniseq in enumerate(uniquename):
            cns = np.where(maskname == uniseq)[0]
            if len(cns) > 1:
                print uniseq
                seqlens = mylen(maskseq[cns])
                print seqlens
                cnm = cns[np.argmax(seqlens)]
                print cnm, len(interface[cnm])
                combinemas.append(cnm)
                uniseqs.append(maskseq[cnm])
                uniinfo.append(interfaceinfo[cnm])
                uniinfo[-1][2] = ','.join(interfaceinfo[cns,2])
                uniinterface.append(interface[cnm])
                iscore1 = np.zeros(len(interface[cnm]))
                iscore2 = np.zeros(len(interface[cnm]))
                normscore = np.zeros(len(interface[cnm]))
                for cn in cns:
                    # align
                    
                    alignment = pairwise2.align.globalds(maskseq[cnm],maskseq[cn], matrix, -22,-2)
                    mask0 = np.array(list(alignment[0][0])) != '-'
                    mask1 = np.array(list(alignment[0][1])) != '-'
                    cscore1 = np.zeros(len(mask0))
                    cscore1[mask1] = interscore1[cn]
                    iscore1 += cscore1[mask0]
                    cscore2 = np.zeros(len(mask0))
                    cscore2[mask1] = interscore2[cn]
                    iscore2 += cscore2[mask0]
                    nscore = np.zeros(len(mask0))
                    nscore[mask1] += 1
                    normscore += nscore[mask0]
                # generate count array that counts non-zero entries
                uniscore1.append(iscore1/normscore)
                uniscore2.append(iscore2/normscore)
                
            else:
                cns = cns[0]
                combinemas.append(cns)
                uniinfo.append(interfaceinfo[cns])
                uniseqs.append(maskseq[cns])
                uniinterface.append(interface[cns])
                uniscore1.append(interscore1[cns])
                uniscore2.append(interscore2[cns])
                
        maskname = np.array(uniquename)
        maskseq = iseqs1 = iseqs2 = np.array(uniseqs)
        interscore1 = np.array(uniscore1)
        interscore2 = np.array(uniscore2)
        interface = np.array(uniinterface)
        interfaceinfo = np.array(uniinfo)
        print len(inames1), 'combined to ', len(maskname)
    else:
        combinemas = np.arange(len(maskname), dtype = int)
    
    datalabel = []
    pdblen = []
    for cn, cnam in enumerate(maskname):
        # at this point the scores have been aligned to the mask sequences and have the same order
        if iseqs1[cn] == iseqs2[cn] and iseqs1[cn] == maskseq[cn]:
            fp1, tp1, _ = roc_curve(interface[cn], interscore1[cn])
            fp2, tp2, _ = roc_curve(interface[cn], interscore2[cn])
            fpr.append([fp1, fp2])
            tpr.append([tp1, tp2])
            roc_auc.append([auc(fp1, tp1), auc(fp2, tp2)])
            pval = delong_roc_test(interface[cn], norm01(interscore1[cn]), norm01(interscore2[cn]))[0][0]
            pval = np.around(10.**pval,3)
            if np.isnan(pval):
                pval = 1.
            beat = int(roc_auc[cn][0] > roc_auc[cn][1])
            classbeat = int(pval<=0.1)*(beat+(beat-1))
            totbeat += classbeat
            print pval, beat, roc_auc[cn][0], roc_auc[cn][1], classbeat, ' '.join(interfaceinfo[cn])
            datalabel.append(interfaceinfo[cn][0].replace('_',':')+' '+interfaceinfo[cn][4].split(',')[0]+'('+interfaceinfo[cn][4].split(',')[1][0]+'.'+interfaceinfo[cn][4].split(',')[1].split('_')[1][0]+'.) '+interfaceinfo[cn][3]+'\n'+interfaceinfo[cn][6].split('.')[0]+'% '+interfaceinfo[cn][8].replace(' ',''))
            pdblen.append(len(interface[cn]))
            rocstats.append([pval, beat, roc_auc[cn][0], roc_auc[cn][1], classbeat])
        else:
            print iseqs1[cn], iseqs2[cn], maskseq[cn]
            print 'Not equal. Align sequences to mask'
            sys.exit()
        print totbeat
    
    
    
    savefig = False
    if "--savefig" in sys.argv:
        outname = sys.argv[sys.argv.index("--savefig")+1].rsplit('.')[0]
        print "Figurename", outname
        savefig = True



    if '--interface_visualization' in sys.argv:
        intdistcut = float(sys.argv[sys.argv.index('--interface_visualization')+1])
        intname, intseq, interdistance = readinscore(maskfile, 3)
        secname, secseq, secstruc = readinscore(maskfile, 4)
        print len(maskname)
        print maskname[1]
        print len(interscore1[1])
        for cn, cnam in enumerate(maskname):
            score1 = np.array(interscore1[cn])
            score2 = np.array(interscore2[cn])
            print cn, combinemas[cn], len(secstruc[combinemas[cn]]), len(interdistance[combinemas[cn]]), len(score1), len(score2), len(maskseq[cn])
            
            fig = plotsecondary(secstruc[combinemas[cn]], interdistance[combinemas[cn]], scorenames, score1, score2, maskseq[cn], intdistcut)
            print ' '.join(interfaceinfo[cn])
            print ' '.join(np.array(rocstats[cn] ,dtype = str))
            if '--savefig' in sys.argv:
                fig.savefig(outname+'_'+cnam+'-int'+str(intdistcut)+'_scorealign.jpg', bbox_inches = 'tight', dpi = 300)
            else:
                plt.show()
            plt.close(fig)
        

    if '--rocplots' in sys.argv:
        # row and column sharing
        ntemps = len(maskname)
        ncols = min(4,int(np.sqrt(float(ntemps))))
        nrows = int(float(ntemps)/ncols) + (ntemps % ncols > 0)
        fig = plt.figure(figsize=(ncols*3, nrows *3), dpi = 30)
        
        for cn, cnam in enumerate(maskname):
            ax = fig.add_subplot(nrows,ncols,cn+1)
            ax.plot(fpr[cn][0], tpr[cn][0], color='limegreen',
            lw=2)
            ax.set_xlim([0.,1.])
            ax.set_ylim([0.,1.])
            ax.plot(fpr[cn][1], tpr[cn][1], color='blueviolet', lw=2)
            ax.annotate(str(round(roc_auc[cn][0],2)),(0.8,0.15),  color = 'limegreen' )
            ax.annotate(str(round(roc_auc[cn][1],2)),(0.6,0.15),  color = 'blueviolet' )
            ax.annotate('p:'+str(rocstats[cn][0]),(0.6,0.05),  color = np.array(['blueviolet','grey','limegreen'])[rocstats[cn][-1]+1] )
            
            dtypeinfo = datalabel[cn]
            
            ax.set_xlabel(dtypeinfo)
            ax.plot([0, 1], [0, 1], color='grey', lw=2, linestyle='--')
            if cn % ncols == 0:
                print 'y', cn
                ax.tick_params(axis='y', which='both', left=True, right=False, labelleft = True, labelright = False)
                ax.set_yticks([0,0.5,1])
            else:
                ax.tick_params(axis='y', which='both', left=True, right=False, labelleft = False, labelright = False)
            if int(cn/ncols) == 0:
                print 'x', cn
                ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom = False, labeltop = True)
                ax.set_xticks([0,0.5,1])
            else:
                ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom = False, labeltop = False)
                
        plt.subplots_adjust(hspace = 0.4, wspace = 0.3)
        if savefig:
            fig.savefig(outname+'_rocauc.jpg', bbox_inches = 'tight', dpi = 300)   # save the figure to file
            plt.close(fig)

    
    pr = []
    rc = []
    pr_auc = []
    pr_exp = []
    for cn, cnam in enumerate(maskname):
        if iseqs1[cn] == iseqs2[cn] and iseqs1[cn] == maskseq[cn]:
            
            prc, rcc, _ = precision_recall_curve(interface[cn], interscore1[cn])
            prc2, rcc2, _ = precision_recall_curve(interface[cn], interscore2[cn])
            pr.append([prc, prc2])
            rc.append([rcc, rcc2])
            pr_auc.append([average_precision_score(interface[cn], interscore1[cn]), average_precision_score(interface[cn], interscore2[cn])])
            pr_exp.append(np.sum(interface[cn])/float(len(interface[cn])))            

        
    if '--prplots' in sys.argv:
        # row and column sharing
        ntemps = len(maskname)
        ncols = min(4,int(np.sqrt(float(ntemps))))
        nrows = int(float(ntemps)/ncols) + (ntemps % ncols > 0)
        figr = plt.figure(figsize=(ncols*3, nrows *3), dpi =30)
        #f, ax = plt.subplots(nrows, 5, sharex='col', sharey='row')
        #
        for cn, cnam in enumerate(maskname):
            ax = figr.add_subplot(nrows,ncols,cn+1)
            expected = pr_exp[cn]
            ax.step(rc[cn][0], pr[cn][0], color='limegreen', alpha = 0.9, where='post', lw = 2)
            if scorefile2 is not None:
                ax.step(rc[cn][1], pr[cn][1], color='blueviolet', alpha = 0.9, where='post', lw = 2)
                ax.annotate(str(round(pr_auc[cn][1],2)),(0.6,0.8),  color = 'blueviolet' )
                ax.fill_between(rc[cn][1], pr[cn][1], step = 'post', color='blueviolet', alpha = 0.1)
            ax.annotate(str(round(expected,2)),(0.9,expected-0.2),  color = 'r' )    
            ax.annotate(str(round(pr_auc[cn][0],2)),(0.8,0.8),  color = 'limegreen' )
    
            ax.plot([0.,1], [expected, expected], color='red', alpha = 0.5, linestyle = '--')
            ax.fill_between(rc[cn][0], pr[cn][0], step = 'post', color='limegreen', alpha = 0.2)
            ax.set_xlim([0.,1.])
            ax.set_ylim([0.,1.])
            ax.set_xlabel(datalabel[cn])
            
            if cn % ncols == 0:
                ax.tick_params(axis='y', which='both', left=True, right=False, labelleft = True, labelright = False)
                ax.set_yticks([0,0.5,1])
            else:
                ax.tick_params(axis='y', which='both', left=True, right=False, labelleft = False, labelright = False)
            if int(cn/ncols) == 0:
                ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom = False, labeltop = True)
                ax.set_xticks([0,0.5,1])
            else:
                ax.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom = False, labeltop = False)
                
        plt.subplots_adjust(hspace = 0.4, wspace = 0.3)
        if savefig:
            figr.savefig(outname+'_prauc.jpg', bbox_inches = 'tight', dpi = 300)   # save the figure to file
            plt.close(figr)
    
    print np.median(np.array(roc_auc), axis = 0), np.mean(np.array(roc_auc), axis = 0)
    if '--boxroc' in sys.argv:
        fig2 = plt.figure(figsize = (3.5,3.5), dpi = 200)
        ax2 = fig2.add_subplot(111)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        #sns.violinplot(data = np.array(roc_auc), split=True, show_median = True, ax = ax2)
        sns.boxplot(data=np.array(roc_auc), palette = ['limegreen', 'blueviolet'], ax = ax2)
        sns.swarmplot(data = np.array(roc_auc), ax = ax2, color = 'k')
        ax2.set_xticklabels(scorenames, rotation = 60)
        if savefig:
            fig2.savefig(outname+'_rocauc-box.jpg', dpi = 300, bbox_inches = 'tight')   # save the figure to file
            plt.close(fig2)
        
    print np.median(np.array(pr_auc), axis = 0), np.mean(np.array(pr_auc), axis = 0)    
    if '--boxpr' in sys.argv:
        fig3 = plt.figure(figsize = (3.5,3.5), dpi = 200)
        ax3 = fig3.add_subplot(111) 
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        #sns.violinplot(data = np.append(np.array(pr_auc).T, [pr_exp], axis = 0).T, split=True, show_median = True, ax = ax3)
        sns.boxplot(data=np.append(np.array(pr_auc).T, [pr_exp], axis = 0).T,palette = ['limegreen', 'blueviolet', 'grey'], ax = ax3)
        sns.swarmplot(data = np.append(np.array(pr_auc).T, [pr_exp], axis = 0).T, color='k', ax = ax3)
        ax3.set_xticklabels(np.append(scorenames,['Random']), rotation = 60)
        if savefig:
            fig3.savefig(outname+'_prauc-box.jpg', bbox_inches = 'tight', dpi = 300)   # save the figure to file
            plt.close(fig3)


    if '--diffscatter' in sys.argv:
        fig5 = plt.figure(figsize = (6.5,3.), dpi = 200)
        ax5 = fig5.add_subplot(121)
        
        sort = np.argsort(np.absolute(np.array(rocstats)[:,4]))
        ax5.scatter(np.array(roc_auc)[:,1][sort],np.array(roc_auc)[:,0][sort], c = np.array(['blueviolet','grey','limegreen'])[np.array(rocstats, dtype = int)[sort,-1]+1], marker = 'o', s = np.array(pdblen)[sort]*1.5, edgecolor = 'k')
        
        lims = [np.amin(roc_auc)-0.03, 1.03]
        ax5.plot([lims[0],lims[1]], [lims[0],lims[1]], c = 'dimgrey')
        ax5.set_xlim(lims)
        ax5.set_ylim(lims)
        ax5.spines['top'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        
        ax6 = fig5.add_subplot(122)
        ax6.scatter(np.array(pr_auc)[:,1][sort],np.array(pr_auc)[:,0][sort], c = np.array(['blueviolet','grey','limegreen'])[np.array(rocstats, dtype = int)[sort,-1]+1], marker = 'o', s = np.array(pdblen)[sort]*1.5, edgecolor = 'k')
        lims = [np.amin(pr_auc)-0.03, 1.03]
        ax6.plot([lims[0], lims[1]], [lims[0], lims[1]], c = 'dimgrey')
        ax6.set_xlim(lims)
        ax6.set_ylim(lims)
        ax6.spines['top'].set_visible(False)
        ax6.spines['right'].set_visible(False)
        ax5.set_title('AUROC')
        ax6.set_title('PRAUC')
        
        ax5.set_ylabel(scorenames[0])
        ax5.set_xlabel(scorenames[1])
        ax6.set_ylabel(scorenames[0])
        ax6.set_xlabel(scorenames[1])        
        
        
        
        if savefig:
            fig5.savefig(outname+'_auc-scatter.jpg', dpi = 300, bbox_inches = 'tight')   # save the figure to file
            plt.close(fig5)    
            


    if "--wilcoxon" in sys.argv:
        print 'Wilcoxon correlation'
        print np.shape(roc_auc)
        roc_auc = np.array(roc_auc)
        pr_auc = np.array(pr_auc)
        wilroc =  wilcoxon(roc_auc[:,0], roc_auc[:,1])[-1]
        wilpr = wilcoxon(pr_auc[:,0], pr_auc[:,1])[-1]
        print wilroc
        print wilpr
        
    
    if '--saveresults' in sys.argv:
        resfile = open(os.path.splitext(scorefile)[0]+'_intAUC.dat', 'w')
        for cn, cnam in enumerate(maskname):
            resfile.write(' '.join(interfaceinfo[cn])+' '+ str(roc_auc[cn][0])+' '+ str(pr_auc[cn][0])+' '+ str(pr_exp[cn])+'\n')
        resfile.close()
        if scorefile2 is not None:
            resfile2 = open(os.path.splitext(scorefile2)[0]+'_intAUC.dat', 'w')
            for cn, cnam in enumerate(maskname):
                resfile2.write(' '.join(interfaceinfo[cn])+' '+ str(roc_auc[cn][1])+' '+ str(pr_auc[cn][1])+' '+ str(pr_exp[cn])+'\n')
            resfile2.close()
    
    if not savefig:
        plt.show()







