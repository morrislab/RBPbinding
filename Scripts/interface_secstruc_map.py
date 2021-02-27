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


                
 


def plotsecondary(secannot, seqscoret, sequt, ylims, xlims, title):
    
    fig = plt.figure(figsize=(17.4,len(secannot)*2.), dpi = 100)
    
    for ai in range(len(secannot)):
        
        secanno = secannot[ai]
        seqscore = seqscoret[ai]
        sequ = sequt[ai]
        xticklabels = xlims[ai]
        
        old_sec = 'Z'
        old_pos = 0
        x = np.arange(len(secanno))
        if '--cutedges' in sys.argv:
            lim = np.where((seqscore != 0) | (np.array(secanno) != 'X'))[0]
            #print secanno, seqscore
            #print lim[1], lim[-2]
            xlim = [lim[1]-1, lim[-2]+1]
        else:
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
        
        
        if len(patches) != 0:
            print 2*len(secannot),1,ai*2+2
            ax = fig.add_subplot(2*len(secannot),1,ai*2+2)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            for s in range(len(locations)):
                sp = locations[s][0]
                tp = locations[s][1]
                ckol = types[s]
                p = ax.plot(sp, tp, color = ckol, linewidth = 4.)
            ax.set_xticks(x)
            ax.set_xticklabels(list(sequ), fontsize=11)
            ax.tick_params(right=False, left = False, labelleft = False, labelright = False, bottom = False)
            collection = PatchCollection(patches, color = colors, alpha=0.7)
            ax.add_collection(collection)
            ax.set_xlim(xlim[0], xlim[1])
            ax2 = fig.add_subplot(2*len(secannot), 1, ai*2+1)
        
        else:
            ax2 = fig.add_subplot(len(secannot), 1, ai+1)
        ax2.set_ylim(ylims[0], ylims[1]*1.05)
        ax2.set_xlim(xlim[0], xlim[1])
        
        
        ax2.bar(x, seqscore, width = 1, color = 'darkslategrey', alpha = 0.8)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        if len(patches) == 0:
            ax3 = ax2.twiny()
            ax3.set_xticks(x)
            ax3.tick_params(top = False, labeltop = False, labelbottom = True)
            ax3.set_xticklabels(list(sequ), fontsize=11, direction = 'in')
            ax3.set_xlim(ax2.get_xlim())
            ax3.tick_params(axis='x', bottom = True, top = False, labelbottom = True, labeltop = False)
        ax2.tick_params(axis='x', bottom = True, top = False, labelbottom = True, labeltop = False)
        
        xlabels = np.arange(0,xticklabels[1], 20)
        xlabels = xlabels[xlabels >= xticklabels[0]]
        xticks = xlabels - xticklabels[0]
        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xlabels, fontsize = 11)
        
        ax2.set_yticks(ylims)
        ax2.set_yticklabels(np.around(ylims,0), fontsize = 11)
        
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
    return fig
        

      


def dompos(totseq, singdseqs):
    alignment = pairwise2.align.globalds(totseq, singdseqs, matlist.blosum62, -11, -.5)
    alignment = alignment[0]
    print alignment[0]
    print alignment[1]
    posmask = np.array(list(alignment[1])) != '-'
    return posmask

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
    ssfile = open(sys.argv[1], 'r').readlines()
    secsequences = []
    secstr = []
    secnames = []
    secpronames = []
    for s, sline in enumerate(ssfile):
        if sline[0] == '>':
            secpronames.append(sline[1:].strip().split('_')[0])
            secnames.append(sline[1:].strip())
            secsequences.append(ssfile[s+1].strip())
            secstr.append(ssfile[s+2].strip())
    secsequences = np.array(secsequences)
    secstr = np.array(secstr)
    secnames = np.array(secnames)
    secpronames = np.array(secpronames)
    
    intersequences = []
    inames = []
    interscores = []
    interface_file = open(sys.argv[2], 'r').readlines()
    for s, sline in enumerate(interface_file):
        if sline[0] == '>':
            inames.append(checkname(sline, secpronames))
            intersequences.append(interface_file[s+1].strip())
            interscores.append(interface_file[s+2].strip().split(','))
  
    inames = np.array(inames)
    interscores = np.array(interscores)
    intersequences = np.array(intersequences)
    
    # provide a confidence value to filter some sequences with bad confidence
    if '--confidence' in sys.argv:
        confident = np.genfromtxt(sys.argv[sys.argv.index('--confidence')+1] , dtype = str)
        confcut = float(sys.argv[sys.argv.index('--confidence')+2])
        confident = confident[confident[:, -1].astype(float) <= confcut]
        for c, conf in enumerate(confident[:,0]):
            confident[c,0] = conf.split('_')[0]
        smask = np.isin(secpronames, confident[:,0])
        secnames = secnames[smask]
        secpronames = secpronames[smask]
        secstr= secstr[smask]
        secsequences = secsequences[smask]
        imask = np.isin(inames, confident[:,0])
        inames = inames[imask]
        interscores = interscores[imask]
        intersequences = intersequences[imask]
    elif '--proteinset' in sys.argv:
        confident = sys.argv[sys.argv.index('--proteinset')+1]
        if os.path.isfile(confident):
            confident = np.genfromtxt(confident, dtype = str)
        else:
            confident = confident.split(',')
        smask = np.isin(secpronames, confident)
        secnames = secnames[smask]
        secpronames = secpronames[smask]
        secstr= secstr[smask]
        secsequences = secsequences[smask]
        imask = np.isin(inames, confident)
        inames = inames[imask]
        interscores = interscores[imask]
        intersequences = intersequences[imask]
        

    # map the map onto the full sequences, and show the domain in context of the sequence
    fullseq = False
    if '--fullsequences' in sys.argv:
        fullsequences = []
        fullnames = []
        wholeseq_file = open(sys.argv[sys.argv.index('--fullsequences')+1], 'r').readlines()
        for s, sline in enumerate(wholeseq_file):
            if sline[0] == '>':
                fullnames.append(checkname(sline[1:],inames))
                fullsequences.append(wholeseq_file[s+1].strip())
        fullsequences = np.array(fullsequences)
        fullnames = np.array(fullnames)
        fullseq = True
    
    # show only the domain sequences but provide location of the domain in the sequence
    fullseqloc = False
    if '--fullseqloc' in sys.argv:
        fullsequences = []
        fullnames = []
        wholeseq_file = open(sys.argv[sys.argv.index('--fullseqloc')+1], 'r').readlines()
        for s, sline in enumerate(wholeseq_file):
            if sline[0] == '>':
                fullnames.append(checkname(sline[1:],inames))
                fullsequences.append(wholeseq_file[s+1].strip())
        fullsequences = np.array(fullsequences)
        fullnames = np.array(fullnames)
        fullseqloc = True
    
    #sort them for the sequences with interface information
    #unames = np.unique(inames)
    
    if '--savefig' in sys.argv:
        plotname = sys.argv[sys.argv.index('--savefig')+1]
        plotfmt = plotname.rsplit('.',1)[-1]
        plotname = plotname.rsplit('.',1)[0]
    
    for u, uname in enumerate(inames):
        #ds = np.where(inames == uname)
        print uname
        
        dscores = np.array(interscores[u], dtype = float)
        dsseqs = intersequences[u]
        print dsseqs
        if '--correctedgescores' in sys.argv:
            if '*' in dsseqs:
                loc = np.concatenate([[-1], np.where(dsseqs == '*')[0] ,[len(dsseqs)]])
            else:
                loc = [-1, len(dsseqs)]
            for l in range(len(loc)-1):
                lowd = np.amin(dscores[loc[l]+6:loc[l+1]-6])
                lowscore = dscores[loc[l]+1:loc[l]+6]
                lowscore[lowscore < lowd] = lowd
                dscores[loc[l]+1:loc[l]+6] = lowscore

                highscore = dscores[loc[l+1]-6:loc[l+1]]
                highscore[highscore < lowd] = lowd
                dscores[loc[l+1]-6:loc[l+1]] = highscore
            
        if '--ylimmean' in sys.argv:
            #dscores = dscores - np.mean(dscores)
            ylims = [np.mean(dscores), np.amax(dscores)]
        elif '--yzscale' in sys.argv:
            dscores = dscores - np.mean(dscores)
            dscores = dscores/np.std(dscores)
            ylims = [0, np.amax(dscores)]
        else:
            dscores = dscores - np.amin(dscores)
            ylims = [np.amin(dscores), np.amax(dscores)]        
        
        
        secid = np.where(secpronames == uname)[0]
        
        if fullseq:
            dsmask = np.array(list(dsseqs)) != '*'
            dsseqs = ''.join(np.array(list(dsseqs))[dsmask])
            dscores = dscores[dsmask]
            
            fseq = fullsequences[list(fullnames).index(uname)]
            dompositions = dompos(fseq, dsseqs)
            
            fseqsec = np.array(['X' for x in range(len(fseq))])
            fseqscore = np.zeros(len(fseq))
            
            fseqscore[dompositions] = dscores
            
            
            #print secid
            if len(secid) > 0:
                for d in secid:
                    dompositions = dompos(fseq, secsequences[d])
                    fseqsec[dompositions] = np.array(list(secstr[d]))
       
            fseq = [fseq]
            fseqscore = [fseqscore]
            fseqsec = [fseqsec]
            xlims = [[0, len(fseq[0])]]
            
        else:
            if '*' in dsseqs:
                dsmask = np.concatenate([[-1],np.where(np.array(list(dsseqs)) == '*')[0],[len(dsseqs)]])
                fseq = []
                fseqscore = []
                dsseqs = dsseqs.split('*')
                for ds in range(len(dsmask)-1):
                    fseq.append(dsseqs[ds])
                    fseqscore.append(dscores[dsmask[ds]+1:dsmask[ds+1]])
                
            else:
                fseq = [dsseqs]
                fseqscore = [dscores]
            fseqsec = []
            xlims=[]
            xstart = 0
            for fs in fseq:
                found = False
                if fullseqloc:
                    checkseq = fullsequences[list(fullnames).index(uname)]
                    dompositions = domainlocation(checkseq, fs, xstart)
                    xlims.append(dompositions)
                    xstart = dompositions[1]
                else:
                    xlims.append([xstart, xstart + len(fs)])
                    xstart += len(fs) + 1
                for si in secid:
                    idid = fs == secsequences[si]
                    if idid:
                        fseqsec.append(np.array(list(secstr[si])))
                        found = True
                        break
                if found == False:
                    print fs, secsequences[secid], 'has no structure'
                    sys.exit()
        
        #print fseqsec, fseqscore, fseq, ylims, xlims 
        fig = plotsecondary(fseqsec, fseqscore, fseq, ylims, xlims, uname)
        #plt.tight_layout(pad=0.2)
        dpi = 100 # 50
        if '--savefig' in sys.argv:
            fig.savefig(plotname+'-'+inames[u]+'.'+plotfmt, format=plotfmt, dpi = dpi, bbox_inches = 'tight')
            print plotname+'-'+inames[u]+'.'+plotfmt
            plt.close()
        else:
            plt.show()

    
