#make_newick_phylo.py
import sys, os
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from functools import reduce
import logomaker as lm
import pandas as pd
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62
from scipy.stats import pearsonr
from matplotlib.colors import LogNorm
from matplotlib import cm


def plotdistribution(x, y, xticks = None, xticklabels = None, yticks = None, yticklabels = None, xlabel = None, ylabel = None, plotdist = True, xlim = None, ylim = None, style = 'scatter', countscale = 'linear', color = None, title = None, xtresh = None, ytresh = None):
    bottom = 0.1
    top = 0.9
    left = 0.1
    right = 0.9
    if plotdist:
        mainwidth = 0.75*(right - left)
        mainheight = 0.75*(top - bottom)
        vwidth = 0.2*(right - left)
        hheight = 0.2*(top - bottom)
    else:
        mainwidth = right - left
        mainheight = top - bottom
    if plotdist:    
        fig = plt.figure(figsize = (7,7))
    else:
        fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot(221)
    ax.set_position([left,bottom,mainwidth,mainheight])

    if xlim is not None:
        xlimdist = 0.025*(xlim[1] - xlim[0])
        ixlim = [xlim[0]-xlimdist, xlim[1] +xlimdist]
    else:
        ixlim = ax.get_xlim()
    if ylim is not None:
        ylimdist = 0.025*(ylim[1] - ylim[0])
        iylim = [ylim[0]-ylimdist, ylim[1] +xlimdist]
    else:
        iylim = ax.get_ylim()
    
    if xtresh is not None:
        ax.axvspan(xtresh[0], xtresh[1], facecolor='purple', alpha=0.2)
        ax.plot([xtresh[1], xtresh[1]], iylim, ls = '--', c = 'purple')
        ax.text(xtresh[1], 1.005*iylim[1], str(int(np.sum(x<=xtresh[1]))), va = 'bottom', ha= 'right')
    if ytresh is not None:
        ax.axhspan(ytresh[0], ytresh[1], facecolor='steelblue', alpha=0.2)        
        ax.plot(ixlim, [ytresh[0], ytresh[0]], ls = '--', c = 'steelblue') 
        ax.text(1.005*ixlim[1],ytresh[0], str(int(np.sum(y>=ytresh[0]))), va = 'bottom', ha= 'left', rotation = 270)
    
    if style == 'scatter':
        if color is None:
            ax.scatter(x,y, alpha = 0.5, s = 100)
        else:
            ax.scatter(x,y, alpha = 0.9, cmap= cm.viridis_r, c = color, vmin = ylim[0], vmax = ylim[1], s = 100)
            axc = fig.add_subplot(224)
            axc.set_position([left+0.8*mainheight,bottom+0.775*mainheight,0.075*mainwidth,0.15*mainheight])
            cmat = np.linspace(ylim[0], ylim[1], 100)
            axc.imshow(cmat.reshape(-1,1), cmap = cm.viridis_r, aspect = 'auto', origin = 'lower')
            axc.spines['top'].set_visible(False)
            axc.spines['bottom'].set_visible(False)
            axc.spines['right'].set_visible(False)
            axc.spines['left'].set_visible(False)
            axc.tick_params(left = False, labelleft = False, right = False, labelright = True, bottom = False, labelbottom = False, labeltop = False, top = False)
            axc.set_title('Highest SID')
            #if xtresh is not None:
                #axc.plot([-0.5,0.5], [np.where(cmat > np.amin(color[x<=xtresh[1]]))[0][0], np.where(cmat > np.amin(color[x<=xtresh[1]]))[0][0]], c = 'r')
                #axc.set_yticks([np.where(cmat > np.amin(color[x<=xtresh[1]]))[0][0],99])
                #axc.set_yticklabels([np.around(np.amin(color[x<=xtresh[1]]),0), ylim[1]])
            #else:
            axc.plot([-0.5,0.5], [np.where(cmat > np.amin(color))[0][0], np.where(cmat > np.amin(color))[0][0]], c = 'r')
            axc.set_yticks([np.where(cmat > np.amin(color))[0][0],99])
            axc.set_yticklabels([np.around(np.amin(color),0), ylim[1]])
            
    elif style == 'hexbin':
        cscale = 'log'
        ax.hexbin(x,y, gridsize = (100,100), bins = cscale, )
    elif style == 'hist2d':
        cscale =LogNorm()
        ax.hist2d(x,y, bins = [100,100], norm = cscale)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels)
    if xlim is not None:
        ax.set_xlim(ixlim)
    if ylim is not None:
        ax.set_ylim(iylim)
    
    if plotdist:
        axh = fig.add_subplot(222)
        axh.set_position([left,bottom+0.05+mainheight,mainwidth,hheight])
        axh.hist(x, bins = 100, color = cm.viridis(0))
        if countscale == 'log':
            axh.set_yscale('log')
        axh.spines['top'].set_visible(False)
        axh.spines['bottom'].set_visible(False)
        axh.spines['right'].set_visible(False)
        axh.tick_params(left = False, labelleft = True, right = False, labelright = False, bottom = False, labelbottom = False)
        
        axv = fig.add_subplot(223)
        axv.set_position([left+0.05+mainwidth,bottom,vwidth,mainheight])
        axv.hist(y, bins = 100, orientation = 'horizontal', color = cm.viridis(0))
        if countscale == 'log':
            axv.set_xscale('log')        
        axv.spines['top'].set_visible(False)
        axv.spines['left'].set_visible(False)
        axv.spines['right'].set_visible(False)
        axv.tick_params(left = False, labelleft = False, right = False, labelright = False, bottom = False, labelbottom = True)
        
        if xtresh is not None:
            axh.axvspan(xtresh[0], xtresh[1], facecolor='purple', alpha=0.2)
        if ytresh is not None:
            axv.axhspan(ytresh[0], ytresh[1], facecolor='steelblue', alpha=0.2) 
            
        print ixlim, iylim    
        axh.set_xlim([ixlim[0], ixlim[1]])        
        axv.set_ylim([iylim[0], iylim[1]]) 
        
            
    if title is not None:
        fig.text(right,top + 0.04,title, ha = 'right', va = 'bottom', fontsize = 13)
    return fig







jpledistance = np.genfromtxt(sys.argv[1], dtype = str)
sequence_maxseqid = np.genfromtxt(sys.argv[2], dtype = str)

jplesort = np.argsort(jpledistance[:,0])
sidsort = np.argsort(sequence_maxseqid[:,0])
jpledistance = jpledistance[jplesort]
sequence_maxseqid = sequence_maxseqid[sidsort]
if not np.array_equal(jpledistance[:,0], sequence_maxseqid[:,0]):
    jpledistance = jpledistance[np.isin(jpledistance[:,0], sequence_maxseqid[:,0])]
    sequence_maxseqid = sequence_maxseqid[np.isin(sequence_maxseqid[:,0],jpledistance[:,0])]


xthreshold = float(sys.argv[3])
ythreshold = float(sys.argv[4])




fig2 = plotdistribution(jpledistance[:,-1].astype(float), sequence_maxseqid[:,-1].astype(float), xticks = None, xticklabels = None, yticks = None, yticklabels = None, xlabel = 'JPLE latent distance', ylabel = 'Closest Seq ID', plotdist = True, style = 'hist2d', countscale = 'linear', title = 'Reconstruction 692 eukaryotes',xlim = [0.,1.0], ylim = [0,100], xtresh = [-.1, xthreshold], ytresh = [ythreshold, 110] )


if '--outdir' in sys.argv:
    outdir = sys.argv[sys.argv.index('--outdir')+1]+os.path.splitext(os.path.split(sys.argv[1])[1])[0]
else:
    outdir = os.path.splitext(sys.argv[1])[0]

if '--savefig' in sys.argv:
    fig2.savefig(outdir+'cosine-closeid.jpg', dpi = 300, bbox_inches = 'tight', transparent=True)
else:
    plt.show()





    

