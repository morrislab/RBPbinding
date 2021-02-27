#matrix_factorization.py 
#Identity_regression.py
import numpy as np
import scipy
import sys, os
import matplotlib
if "--savefig" in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import svd
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
import time

def readin(proteinfeatures, profiles, proteinset = 'none'):
    if os.path.isfile(profiles):
        Y = np.genfromtxt(profiles)[:,1:]
        Zfeatures = np.genfromtxt(profiles, dtype = str)[:,0]
        ynames = np.array(open(profiles, 'r').readline().strip().split()[1:])
        Y = np.nan_to_num(Y).T
        pfile = True
    else:
        pfile = False
        Y = None
        ynames = None
        Zfeatures = None

    if os.path.splitext(proteinfeatures)[1] == '.npz':
        pobj = np.load(proteinfeatures)
        P = pobj['features']
        sequencefeatures = pobj['kmers']
        protnames = pobj['protnames']
        pnames = pobj['expnames']
        protnames = np.append([pnames], [protnames], axis = 0)
    else:
        pobj = open(proteinfeatures, 'r').readlines()
        protnames = []
        sequencefeatures = []
        ph = 0
        for line in pobj:
            if line[0] == '#':
                ph +=1
                protnames.append(line.strip().split()[1:])
            else:
                sequencefeatures.append(line.split()[0])
        P = np.genfromtxt(proteinfeatures, skip_header = ph)[:,1:]

    seqfeatsort = np.argsort(sequencefeatures)
    sequencefeatures = sequencefeatures[seqfeatsort]
    P = P[seqfeatsort]
    P = P.T

    if pfile:
        #sort after pnames
        if len(np.intersect1d(ynames, protnames[0])) != 0:
            protnames = np.array(protnames[0])
        else:
            protnames = np.array(protnames[1])
            
        indexes = []
        for i, ina in enumerate(protnames):
            if ina in ynames:
                indexes.append([i, list(ynames).index(ina)])
        indexes = np.array(indexes).T
        Y = Y[indexes[1]]
        P = P[indexes[0]]
        protnames = protnames[indexes[0]]
        ynames = ynames[indexes[1]]
        mprot = np.shape(P)[0]
        nprobes = np.shape(Y)[1]
    
    if os.path.isfile(proteinset):
        pset = np.genfromtxt(proteinset, dtype = str)
        pmask = np.isin(protnames, pset)
        Y = Y[pmask]
        P = P[pmask]
        protnames = protnames[pmask]
        ynames = ynames[pmask]
    
    print np.shape(Y), np.shape(P)
    return Y, ynames, Zfeatures, P, protnames, sequencefeatures



profiles = sys.argv[1]
proteinfeatures = sys.argv[2]
proteinset =  sys.argv[3]
pearsonlow = float(sys.argv[4])
if len(sys.argv) > 5:
    outname = sys.argv[5]
else:
    outname = ''

Y, datnames, Zfeatures, P, pnames, sequencefeatures = readin(proteinfeatures, profiles, proteinset = proteinset)
XY = np.append(P,Y, axis = 1)
np.shape(XY)

if '--normP2' in sys.argv:
    P = P/np.sqrt(np.sum(P*P, axis = 1))[:, None]
if '--normY2' in sys.argv:
    Y = Y/np.sqrt(np.sum(Y*Y, axis = 1))[:, None]

# SVD for concatenated specificity and sequence k-mers
u,s,v = np.linalg.svd(XY, full_matrices=False)
# Variance explained for joint features
s2 = s**2/np.sum(s**2)
# eigenvector parts of joint eigenvector for z-scores
vr = v[:, -len(Y[0]):]
# eigenvector parts of joint eigenvector for protein features
vp = v[:,: -len(Y[0])]

# Rescaled explained variance captured for joint matrices
sr2 = s**2 * np.sum(v[:, -len(Y[0]):]**2, axis = 1)/ np.sum(s**2 * np.sum(v[:, -len(Y[0]):]**2, axis = 1))
sp2 = s**2 * np.sum(v[:,: -len(Y[0])]**2, axis = 1)/ np.sum(s**2 * np.sum(v[:, :-len(Y[0])]**2, axis = 1))

# SVD only specificities
uy, sy, vy = np.linalg.svd(Y, full_matrices=False)
sy2 = sy**2/np.sum(sy**2)

# SVD only protein features
ux, sx, vx = np.linalg.svd(P, full_matrices=False)
sx2 = sx**2/np.sum(sx**2)



# measure top 100 overlap
def topx(s1, s2, top):
    tout = float(len(np.intersect1d(np.argsort(-s1)[:top],np.argsort(-s2)[:top])))/float(top)
    return tout


## resort joint data according to variance explained for specificity
rearrange = np.argsort(-sr2)


## calculate pearson correlation of reconstructed specificities from subsets of its eigenvectors
hightop = []
lowtop =[]
mediantop = []

high = []
low =[]
median = []
for i in range(1, len(sx)):
    yrec = np.dot(uy[:, :i]*sy[:i],vy[:i])
    corr = np.nan_to_num([pearsonr(yrec[y], Y[y])[0] for y in range(len(Y))])
    #top100 = np.nan_to_num([topx(yrec[y], Y[y], 100) for y in range(len(Y))])
    high.append(np.amax(corr))
    low.append(np.amin(corr))
    median.append(np.median(corr))
    #hightop.append(np.amax(top100))
    #lowtop.append(np.amin(top100))
    #mediantop.append(np.median(top100))


high = np.array(high)
low = np.array(low)
median = np.array(median)

hightop = np.array(hightop)
lowtop = np.array(lowtop)
mediantop = np.array(mediantop)


# Figure that captures reconstruction quality of specificities from subset of eigenvectore
figrec = plt.figure(figsize = (4,4), dpi = 150)
axr = figrec.add_subplot(111)
axr.spines['right'].set_visible(False)
axr.spines['top'].set_visible(False)
axr.plot(np.arange(len(sy2)), np.cumsum(sy2), label = 'Sum(Variance(R))', c = 'steelblue')
axr.scatter(np.arange(len(sx2)-1), median, label = 'Pearson median', c = 'k', s = 30)
axr.fill_between(np.arange(len(sx2)-1), low, high, color = 'grey', alpha = 0.3, linewidth = 0)
axr.scatter(np.arange(len(sx2)-1), low, label = 'Pearson min', marker = '_', c = 'k')
axr.scatter(np.arange(len(sx2)-1), high, label = 'Pearson max', marker = '_', c = 'grey')
axr.legend(prop={'size':8})
axr.set_xticks([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])
axr.set_xticklabels([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350], rotation = 60)
axr.set_xlabel('Number maintained eigenvectors')
axr.set_ylabel('Variance explained/\nDistribution over Pearson reconstructed R')
low09 = np.where(low>=pearsonlow)[0][0]
axr.plot([low09, low09],[0, np.cumsum(sy2)[low09]], c = 'r', ls = '--')
axr.plot([0, low09], [np.cumsum(sy2)[low09], np.cumsum(sy2)[low09]], c = 'r', ls = '--')
axr.text(low09+5, 0.85, 'Eigenvalue '+str(low09)+'\nat min Pearson='+str(pearsonlow)+'\nSum(Var(R))='+str(np.around(np.cumsum(sy2)[low09], 2)), va = 'top', ha = 'left' )


# Variance explained from eigenvectors from 3 types of SVD at defined variance/reconstruction cut-off
fig = plt.figure(figsize = (4,4), dpi = 150)
ax = fig.add_subplot(111)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.plot(np.arange(len(s2)), np.cumsum(s2), label = 'Sum(Variance(R,P))', c = 'grey')
ax.plot(np.arange(len(sx2)), np.cumsum(sx2), label = 'Sum(Variance(P))', ls = 'dotted', c = 'green')
ax.plot(np.arange(len(sy2)), np.cumsum(sy2), label = 'Sum(Variance(R))', ls = 'dotted', c = 'steelblue')

ax.plot(np.arange(len(sr2)), np.cumsum(sr2[rearrange]), label = 'Sum(Variance(R,P)_R)', c = 'red')
ax.plot(np.arange(len(sp2)), np.cumsum(sp2[rearrange]), label = 'Sum(Variance(R,P)_P)', c = 'goldenrod')

ax.set_xticks([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])
ax.set_xticklabels([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350], rotation = 60)

ax.legend(prop={'size':8})
ax.set_xlabel('Number maintained eigenvectors')
ax.set_ylabel('Variance explained\n(R of reconstructed specificities)')

low09 = np.where(low>=pearsonlow)[0][0]
ax.plot([low09, low09],[0, np.cumsum(sy2)[low09]], c = 'r', ls = '--')
ax.plot([0, low09], [np.cumsum(sy2)[low09], np.cumsum(sy2)[low09]], c = 'r', ls = '--')
ax.plot([0, low09], [np.cumsum(sp2)[low09], np.cumsum(sp2)[low09]], c = 'r', ls = '--')
ax.set_title('Eigenvalue '+str(low09)+' at min Pearson='+str(pearsonlow)+'\nCumsum(E(R))='+str(np.around(np.cumsum(sy2)[low09], 2))+'\nCumsum(E(P))='+str(np.around(np.cumsum(sp2)[low09], 2)))





# dot product between the individual and joint eigenvectors for P and R
fig2 = plt.figure(figsize = (8,4), dpi = 100)
ax2 = fig2.add_subplot(121)
ax2.imshow(np.absolute(np.dot(vr[rearrange],vy.T)), vmin = 0, vmax = 1., cmap = 'gist_heat_r', origin = 'lower')
ax2.set_xlabel('Eigenvectors R only')
ax2.set_ylabel('Eigenvector_R from P,R')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.set_xticks([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])
ax2.set_xticklabels([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350], rotation = 60)
ax2.set_yticks([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])
ax2.set_yticklabels([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])

axb = fig2.add_subplot(122)
imb = axb.imshow(np.absolute(np.dot(vp[rearrange], vx.T)), vmin = 0, vmax = 1., cmap = 'gist_heat_r', origin = 'lower')
axb.set_xlabel('Eigenvector P only')
axb.set_ylabel('Eigenvector_P from P,R ')
axb.spines['right'].set_visible(False)
axb.spines['top'].set_visible(False)
axb.set_xticks([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])
axb.set_xticklabels([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350], rotation = 60)
axb.set_yticks([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])
axb.set_yticklabels([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])

colorbar = fig2.add_subplot(311)
colorbar.set_position([0.935, 0.75, 0.03, 0.1])
colorbar.imshow(np.linspace(0,1,100).reshape(-1,1), cmap = 'gist_heat_r', vmin = 0, vmax = 1., aspect = 'auto', origin = 'lower')
colorbar.tick_params(left = False, labelleft = False, bottom = False, labelbottom = False, right = True, labelright = True)
colorbar.spines['right'].set_visible(False)
colorbar.spines['top'].set_visible(False)
colorbar.spines['left'].set_visible(False)
colorbar.spines['bottom'].set_visible(False)
colorbar.set_yticks([0, 50,100])
colorbar.set_yticklabels([0, 0.5, 1])


# sum of squared dot-product between the individual and joint eigenvectors for P and R in first X eigenvecotrs and last eigenvectors
fig3 = plt.figure(figsize = (9,3), dpi = 100)
ax3 = fig3.add_subplot(121)
print np.shape(vr[rearrange]), np.shape(vy)
ax3.bar(np.arange(len(rearrange)), np.sum(np.dot(vr[rearrange], vy.T)**2, axis = 0), width = 1., color = 'grey', alpha = 0.5)
ax3.bar(np.arange(len(rearrange)), np.sum((np.dot(vr[rearrange], vy.T)**2)[:low09], axis = 0), width = 1., color = 'brown')
ax3.set_xlabel('Eigenvectors R only')
ax3.set_ylabel('VarianceR in \nEigenvector_R(P,R) <'+str(low09))
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.set_xticks([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])
ax3.set_xticklabels([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350], rotation = 60)

axb3 = fig3.add_subplot(122)
imb3 = axb3.bar(np.arange(len(rearrange)), np.sum(np.dot(vp[rearrange], vx.T)**2, axis = 0), width = 1., color = 'grey', alpha = 0.5)
imb3 = axb3.bar(np.arange(len(rearrange)), np.sum((np.dot(vp[rearrange], vx.T)**2)[:low09], axis = 0), width = 1., color = 'brown')
axb3.set_xlabel('Eigenvector P only')
axb3.set_ylabel('VarianceP in \nEigenvector_P(P,R) <'+str(low09))
axb3.spines['right'].set_visible(False)
axb3.spines['top'].set_visible(False)
axb3.set_xticks([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350])
axb3.set_xticklabels([0,25,50, 75, 100,125, 150, 175, 200, 225, 250, 275, 300, 325, 350], rotation = 60)


figrec.savefig(outname+'EigenvectorR_variance_reconstructions.jpg', dpi = 300, bbox_inches = 'tight')
fig.savefig(outname+'EigenvectorPR_variance_explanation.jpg', dpi = 300, bbox_inches = 'tight')
fig2.savefig(outname+'Eigenvector_variance_composition.jpg', dpi = 300, bbox_inches = 'tight')
fig3.savefig(outname+'Eigenvector_variance_composition_sum.jpg', dpi = 300, bbox_inches = 'tight')

plt.show()








