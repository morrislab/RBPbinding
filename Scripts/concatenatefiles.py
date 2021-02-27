# concatenatefiles.py
import numpy as np
import glob 
import sys


# reads in all files with defined start and end and splits at 'separ'

startfile = sys.argv[1]
separ = sys.argv[2]
endfile = sys.argv[3]
outname = sys.argv[4]

listfiles = np.sort(glob.glob(startfile+"*"+endfile))
if len(listfiles) == 0:
    print 'Not found', startfile+"*"+endfile

def checknumloc(string):
    for c, s in enumerate(string):
        if not s.isdigit():
            return string[c:], string[:c]
            
            
    

numbers = []
mainfiles = []
for f, filel in enumerate(listfiles):
    if separ in filel:
        filel = filel.split(separ, 1)
        fname = filel[0]+separ+":"
        fnameleft, fnum = checknumloc(filel[1])
        fname += fnameleft

        # if files are different from each other besides from the set number they get separated 
        if fname in mainfiles and fnum != '': 
            ind = mainfiles.index(fname)
            numbers[ind].append(fnum)
        elif outname not in fnum and fnum != '':
            mainfiles.append(fname)
            numbers.append([fnum])

# define the minimum number of sets needed
if '--minnum' in sys.argv:
    mn = int(sys.argv[sys.argv.index('--minnum')+1])
    keep = []
    for n, numb in enumerate(numbers):
        if len(numb) >= mn:
            keep.append(n)
        else:
            print mainfiles[n], 'does not contain enough files', len(numb)
            print 'Missing:'
            for i in range(mn):
                if str(i) not in numb:
                    print i,
    numbers = np.array(numbers)[keep]
    mainfiles = np.array(mainfiles)[keep]


if '--definesets' in sys.argv:
    dsets = sys.argv[sys.argv.index('--definesets')+1].split(',')
    for n, numb in enumerate(numbers):
        keep = []
        for nn, fnum in enumerate(numb):
            if fnum in dsets:
                keep.append(nn)
        numbers[n] = np.array(numb)[keep]

print 'Created'
for i, mfile in enumerate(mainfiles):
    mfile = mfile.split(":")
    wobj = open(mfile[0]+'_'+outname+"_"+mfile[1], 'w')
    print mfile[0]+'_'+outname+"_"+mfile[1]
    
    for numm in numbers[i]:
        obj = open(mfile[0]+numm+mfile[1], 'r')
        for line in obj:
            wobj.write(line)
        obj.close()
    wobj.close()
    

        
