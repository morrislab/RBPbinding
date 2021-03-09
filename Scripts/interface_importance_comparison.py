import numpy as np
import sys, os
import matplotlib
if "--savefig" in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
#plt.rcdefaults()
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
matrix = matlist.blosum62
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import pandas as pd
import logomaker
               
      


def dompos(totseq, singdseqs):
    alignment = pairwise2.align.globalds(totseq, singdseqs, matlist.blosum62, -11, -.5)
    alignment = alignment[0]
    print alignment[0]
    print alignment[1]
    posmask = (np.array(list(alignment[1])) != '-')
    seqmask = (np.array(list(alignment[0])) != '-')
    return seqmask, posmask

def identity(seqa, seqb):
    alignment = pairwise2.align.globalds(totseq, singdseqs, matlist.blosum62, -11, -1)
    alignment = alignment[0]
    ide = 0.
    for a in range(len(alignment[0])):
        ide += float(alignment[0][a] == alignment[1][a])
    return ide/float(len(alignment[0]))

def domainlocation(fseq, dseq, start):
    for dl in range(start, len(fseq)):
        if dseq == fseq[dl :dl+ len(dseq)]:
            return dl, dl+len(dseq)

    return start, start+len(dseq)


def checkname(name, namelist):
    for na in namelist:
        if na in name:
            return na


## add confidence score to not make jpg for all

# make choice to visualize on whole sequence or per domains
# plot on top of each other... 



if __name__ == '__main__':
    # domain secondary structure file 
    intersequences = []
    inames = []
    interscores = []
    interface_file = open(sys.argv[1], 'r').readlines()
    for s, sline in enumerate(interface_file):
        if sline[0] == '>':
            inames.append(sline[1:].strip())
            intersequences.append(interface_file[s+1].strip())
            interscores.append(interface_file[s+2].strip().split(','))
  
    inames = np.array(inames)
    interscores = np.array(interscores)
    intersequences = np.array(intersequences)
    
    confident = sys.argv[2].split(',')
    imask = np.isin(inames, confident)
    inames = inames[imask]
    interscores = interscores[imask]
    intersequences = intersequences[imask]

    
    if '--savefig' in sys.argv:
        plotname = sys.argv[sys.argv.index('--savefig')+1]
        plotfmt = plotname.rsplit('.',1)[-1]
        plotname = plotname.rsplit('.',1)[0]
    
    seqscore = np.zeros((len(inames), len(intersequences[0])))
    seqalignment = np.chararray((len(inames),len(intersequences[0]) ))
    seqalignment[:] = '-'
    for u, uname in enumerate(inames):
        #ds = np.where(inames == uname)
        print uname
        dscores = np.array(interscores[u], dtype = float)
        dsseqs = intersequences[u]
        print dsseqs
            
        if '--yzscale' in sys.argv:
            dscores = dscores - np.mean(dscores)
            dscores = dscores/np.std(dscores)
            ylims = [0, np.amax(dscores)]
        elif '--yzscaletotal' in sys.argv:
            allscores = np.concatenate(interscores).astye(float)
            dscores = dscores - np.mean(allscores)
            dscores = dscores/np.std(allscores)
            ylims = [0, np.amax(dscores)]
        else:
            dscores = dscores - np.amin(dscores)
            ylims = [np.amin(dscores), np.amax(dscores)]        
        
        if u == 0:
            seqscore[u] = dscores
            seqalignment[u] = np.array(list(dsseqs))
        else:
            seqmask, positions = dompos(intersequences[0].replace('*', 'X'), intersequences[u].replace('*', 'X'))
            alignscore = np.zeros(len(seqmask))
            alignscore[positions] = dscores
            seqscore[u] = alignscore[seqmask]
            alignme = np.chararray(len(seqmask))
            alignme[:] = '-'
            alignme[positions] = np.array(list(intersequences[u]))
            seqalignment[u] = alignme[seqmask]
    
    aminoacid = list(np.sort(['A', 'R', 'N', 'D', 'P', 'V', 'I', 'C', 'Y', 'H', 'T', 'G', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S'])) 
    protframes = []
    for u, uname in enumerate(inames):
        dscores = seqscore[u]
        dsseqs = seqalignment[u]
        pwmdataf = pd.DataFrame()
        for a, amino in enumerate(aminoacid):
            occscore = np.zeros(len(dscores))
            occscore[dsseqs == amino] = dscores[dsseqs == amino]
            pwmdataf[amino] = occscore
        protframes.append(pwmdataf) 
    fig = plt.figure(figsize = (18, len(inames)), dpi = 150)
    for u, uname in enumerate(inames):
        ax = fig.add_subplot(len(inames),1,u+1)
        logomaker.Logo(protframes[u], ax = ax, color_scheme='chemistry')
        ax.set_ylim(ylims)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelleft = False, left = False)
        if u != len(inames)-1:
            #ax.spines['bottom'].set_visible(False)
            ax.tick_params(labelbottom =False, bottom = False,labelleft = False, left = False)
    dpi = 300
    if '--savefig' in sys.argv:
        fig.savefig(plotname+'.'+plotfmt, dpi = dpi, bbox_inches = 'tight')
        print plotname+'.'+plotfmt
        plt.close()
    else:
        plt.show()

    
