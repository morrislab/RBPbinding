import sys,os
import numpy as np
import scipy.stats as st

alphabet = {list('ACGT')[k]:k for k in range(4)}
def readinpwm(pwmfile):
    obj = open(pwmfile, 'r').readlines()
    pwms = []
    pwm = []
    pwmname = []
    for l, line in enumerate(obj):
        line = line.strip()
        if len(line) != 0:
            line = line.split()
            if line[0] == 'Motif':
                if len(pwm) > 0:
                    pwms.append(np.array(pwm, dtype = float))
                pwmname.append(line[1])
                pwm = []
            elif line[0].isdigit():
                pwm.append(line[1:])
    pwms.append(np.array(pwm, dtype = float))
    return pwmname, pwms

def ic(pw):
    icout = pw*np.log2(pw/0.25)
    icout[icout < -2] = -2
    return icout



if __name__ == "__main__":
    # Read in fasta file
    # Example:
    #>Chr1:5630-5680 Parent=AT1G01010.1 block1len50
    #GAGGTCAAATCGGATTCTTGCTCAAAATTTGTATTTCTTAGAATGTGTGT
    #>Chr1:5680-5730 Parent=AT1G01010.1 block2len50
    #TTTTTTTTGTTTTTTTTTCTTTGCTCTGTTTTCTCGCTCCGGAAAAGTTT

    fastafile = np.load(sys.argv[1])
    ids = fastafile['ids']
    sequences = fastafile['sequences']    
    

    for i, idn in enumerate(ids):
        ids[i] = idn.split('Parent=transcript:')[1].split('.')[0]    

    pwmname, pwms = readinpwm(sys.argv[2])
    if '--summax' in sys.argv:
        maxs = float(sys.argv[sys.argv.index('--summax')+1]) 

    if '--maxscore' in sys.argv:
        add = '-maxscore'
    elif '--summax' in sys.argv:
        add = '-summaxgt'+str(maxs)
    else:
        add = ''


    if '--geneset' in sys.argv:
        geneset = np.genfromtxt(sys.argv[sys.argv.index('--geneset')+1], dtype = str)
        print len(geneset), len(ids)
        gmask = np.isin(ids, geneset)
        ids = ids[gmask]
        sequences = sequences[gmask]
        print len(sequences), len(ids)
    
 
    if '--random_shuffle' in sys.argv:
        np.random.seed(int(sys.argv[sys.argv.index('--random_shuffle')+1]))
        add+='_random_shuffle'+sys.argv[sys.argv.index('--random_shuffle')+1]
        for s, seq in enumerate(sequences):
            sequences[s] = seq[np.random.permutation(len(seq))]
    

    # print 'Topfasta', np.prod(np.amax(PSAM, axis = 1))
    # For each sequence in the fasta get the PSAM for each k-mer of len PSAM
    # in the sequence
    scoremat = np.zeros((len(pwmname),len(ids)), dtype = np.float32)
    for p, pname in enumerate(pwmname):
        PSAM = (pwms[p].T/np.amax(pwms[p], axis = 1)).T
    	for i in range(len(sequences)):
            sequence = sequences[i]
            # if sequence is shorter than PSAM just output an empty list
            if len(sequence) >= len(PSAM):
                seqscores = [0]
                # Iterate through k-mers and score
                for j in range(1+len(sequence)-len(PSAM)):
                    seq = sequence[j:j+len(PSAM)]
                    # Score is the product of the base scores at each position in the PSAM
                    seqscores.append(np.prod(PSAM[seq==1]))
                if '--maxscore' in sys.argv:
                    scoremat[p,i] = np.amax(seqscores)
	        elif '--summax' in sys.argv:
		    seqscores = np.array(seqscores)
		    scoremat[p,i] = np.sum(seqscores[seqscores >= maxs])
                else:
                    scoremat[p,i] = np.sum(seqscores)

    outfile = os.path.splitext(sys.argv[3])[0]+add+'.npz'
    np.savez_compressed(outfile, bindingmat = scoremat, rbpnames = pwmname, genenames = ids) 
    
    for p, pprot in enumerate(pwmname):
        print pprot, np.sum(scoremat[p]>0), len(ids)
    print np.percentile(np.sum(scoremat>0, axis =0), [1,10,20,30,50,70,80,90,99])


    

