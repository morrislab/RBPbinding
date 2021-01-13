import numpy as np
import sys, os 
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62



class protein:
    def __init__(self):
        self.proteins = []
        self.sequences = []
        self.domainlocation = []
        self.domainseq = []
        self.domaintype = []
        
    def readdata(self, mfile):
        sobj = open(mfile, 'r').readlines()
        headers = sobj[0].split('\t')
        data = []
        for l,line in enumerate(sobj):
            if l > 0:
                data.append(line.split('\t'))
                if ',' in data[-1][1]:
                    rncmpt= data[-1][1].split(',')
                else:
                    rncmpt = [data[-1][1]]
                for rnc in rncmpt:
                    self.proteins.append(data[-1][0]+'||'+rnc)
                    self.sequences.append(data[-1][headers.index('Construct_AA_seq')])
                    self.domaintype.append(data[-1][headers.index('Domains')].split(','))
                    self.domainlocation.append(np.array(data[-1][headers.index('Domain_Boundaries')].split(','),dtype = int))
                    self.domainseq.append(data[-1][headers.index('RBD_or_RBR_AA_Seq')].split(','))
  

def ident(align1, align2):
    redid = False
    for i in range(len(align1)):
        if align2[i] != '-':
            start = i
            break
    for i in range(len(align1)):
        if align2[-i-1] != '-':
            end = len(align1)-i
            break
    id = 0.
    for a in range(start, end):
        if align1[a] == align2[a]:
            id += 1
    return float(id)/float(end-start), start, end


def findpos(dseq, pseq):
    ldom = len(dseq)
    location = []
    for f in range(len(pseq)-ldom+1):
        if pseq[f:f+ldom] == dseq:
            location.append([f, f+ldom])
    if len(location) > 0:
        return location
    else:
        print 'Not found'
        # Search for high id alternative with biopython
        alignments = pairwise2.align.globalds(pseq, dseq, matrix, -11, -1)
        print 'Best match:'
        print alignments[0][0]
        print alignments[0][1]
        idalign, sta, en = ident(alignments[0][0], alignments[0][1])
        if idalign > 0.75:
            print idalign
            return [[sta, en]]
        else:
            print 'Could not find', dseq, 'in', pseq
            sys.exit()





def extenddomain(protseq, domainseqs, seqextend):
    dompos = []
    for d, dom in enumerate(domainseqs):
        dpos = findpos(dom, protseq)
        dompos.append(dpos)

    # sometimes the same domain sequennce can occur twice
    ifset = []
    ifnum = []
    for p, pos in enumerate(dompos):
        if len([pos]) > 1:
            if '-'.join(np.array(pos)[:,0].astype(str)) in ifset:
                ifn = np.copy(ifnum[ifset.index('-'.join(np.array(dompos[pos])[:,0].astype(str)))])
                ifnum[ifset.index('-'.join(np.array(dompos[pos])[:,0].astype(str)))] += 1
            else:
                ifset.append('-'.join(np.array(pos)[:,0].astype(str)))
                ifnum.append(1)
                ifn = 0
            print pos
            dompos[p] = pos[ifn]
            print dompos[p]
        else:
            dompos[p] = pos[0]

    
    lpseq = len(protseq)
    ndompos = []
    if seqextend > 0:
        for p in range(len(dompos)):
            if p == 0:
                ns = max(0, dompos[p][0]-seqextend)
            else:
                ns = max(dompos[p][0] - seqextend, ndompos[p-1][1])
            if p == len(dompos) -1:
                ne = min(lpseq, dompos[p][1]+seqextend)
            else:
                ne = min(dompos[p][1]+seqextend, max(dompos[p][1],dompos[p+1][0]))
            ndompos.append([ns, ne])
    else:
        ndompos = dompos


    undomseq = []
    for dp in ndompos:
        if dp[0] < dp[1]:
            undomseq.append(protseq[dp[0]: dp[1]])
        else:
            print dp, 'something wronng'
            sys.exit()

    return undomseq, ndompos


def fuse_domains(edom, eloc, domtype):
    fusedomains =[]
    fusedomnames = []
    fudomain = edom[0]
    fudomtype = domtype[0]
    for e in range(len(edom)):
        if e < len(edom)-1:
            if eloc[e][1] == eloc[e+1][0]:
                fudomain += edom[e+1]
                fudomtype += '-'+domtype[e+1]
            else:
                fusedomains.append(fudomain)
                fusedomnames.append(fudomtype)
                fudomain = edom[e+1]
                fudomtype = domtype[e+1]
    fusedomains.append(fudomain)
    fusedomnames.append(fudomtype)
    return fusedomains, fusedomnames

if __name__ == '__main__':

    masterfile = sys.argv[1]
    exten = int(sys.argv[2])
    outname = sys.argv[3]
    rbps = protein()
    rbps.readdata(masterfile)
    edoms = []
    faobj = open(outname+'_domain.fasta', 'w')
    gobj = open(outname+'_domain_fused.fasta', 'w')
    cobj = open(outname+'_combined.fasta', 'w')
    fobj = open(outname+'_pseq.fasta', 'w')
    singfa = False
    for r in range(len(rbps.proteins)):
        edom, eloc = extenddomain(rbps.sequences[r], rbps.domainseq[r], exten)
        edomfuse, fusenames = fuse_domains(edom, eloc, rbps.domaintype[r])
        celoc = np.array(eloc)
        celoc[:,0] = celoc[:,0] +1
        if not np.array_equal(rbps.domainlocation[r], celoc.flatten()):
            print rbps.proteins[r],'before', rbps.domainlocation[r], 'extension', celoc.flatten()
        for e in range(len(edom)):
            faobj.write('>'+rbps.proteins[r]+'__'+rbps.domaintype[r][e]+'_'+str(int(e))+'\n'+edom[e]+'\n')
        for e in range(len(edomfuse)):
            gobj.write('>'+rbps.proteins[r]+'__'+fusenames[e]+'_'+str(int(e))+'\n'+edomfuse[e]+'\n')
        
        cobj.write('>'+rbps.proteins[r]+'\n'+''.join(np.array(edom))+'\n')
        fobj.write('>'+rbps.proteins[r]+'\n'+rbps.sequences[r]+'\n')

