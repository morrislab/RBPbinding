import numpy as np
import sys, os
import logomaker as lm
import pandas as pd
import matplotlib.pyplot as plt
#plt.ion()

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
    for p, pwm in enumerate(pwms):
        pwms[p] = pd.DataFrame({'A':pwm[:,0],'C':pwm[:,1], 'G':pwm[:,2], 'U':pwm[:,3]})
    return pwmname, pwms

def ic(pw):
    icout = pw*np.log2(pw/0.25)
    icout[icout < -2] = -2
    return icout

pnames, pwwm = readinpwm(sys.argv[1])
#print pnames, pwwm

if '--outdir' in sys.argv:
    outname = sys.argv[sys.argv.index('--outdir')+1]+os.path.splitext(os.path.split(sys.argv[1])[1])[0]
else:
    outname = os.path.splitext(sys.argv[1])[0]
if '--infocont' in sys.argv:
    outname += '-infocont'

for p, pname in enumerate(pnames):
    pwm = pwwm[p]
    print pwm
    print len(pwm)
    fig = plt.figure('Logo'+pname, (len(pwm)*0.4, 1.))
    ax = fig.add_subplot(111)
    if '--infocont' in sys.argv:
        pwm = ic(pwm)
    lm.Logo(pwm, ax = ax)
    if '--infocont' in sys.argv:
        ax.set_yticks([0,1,2])
        ax.set_yticklabels([0,1,2])
        ax.set_ylim([0,2])
    else:
        ax.set_yticks([0,1])
        ax.set_yticklabels([0,1])
        ax.set_ylim([0,1])
    if '--removeaxis' in sys.argv:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.tick_params('both', left=False, bottom = False, labelleft = False, labelbottom = False)
    if '--removeframe' in sys.argv:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params('both', left=True, bottom = False, labelleft = True, labelbottom = False)
    
    #fig.tight_layout()
    
    dpi = 300
    
    if '--format' in sys.argv:
        formt = '.'+sys.argv[sys.argv.index('--format')+1]
    else:
        formt = '.jpg'
    
    print 'Saved as', outname+'_'+pname+formt
    if '--transparent' in sys.argv:
        fig.savefig(outname+'_'+pname+formt, dpi = dpi, bbox_inches = 'tight', transparent = True)
    else:
        fig.savefig(outname+'_'+pname+formt, dpi = dpi, bbox_inches = 'tight')
    plt.close()
    
    
    
