import numpy as np
import scipy
import sys, os
from scipy.optimize import nnls
from numpy.linalg import lstsq
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso, Lars
from sklearn.linear_model import HuberRegressor
from sklearn.linear_model import LinearRegression
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.linalg import orth
from numpy.linalg import inv
from scipy.linalg import svd
from sklearn.decomposition import NMF, SparsePCA
import sklearn.linear_model as lm 
from sklearn import preprocessing
from scipy.spatial.distance import euclidean, cosine
from sklearn import decomposition
from sklearn.externals import joblib
#from scipy.stats import logistic
from scipy.spatial.distance import cdist
import time
import multiprocessing
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.metrics.pairwise import polynomial_kernel, rbf_kernel, laplacian_kernel, chi2_kernel
from scipy.optimize import NonlinearConstraint
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt

# Explicit solution of regression. 
# Seems to be much faster than scikitlearn regression for high dimensional output space
# No latent representation 
class regression():
    def __init__(self, alpha = 0., fit_intercept = False, lasso = False):
        self.alpha = alpha
        self.coef_ = None
        self.fit_intercept = fit_intercept
        self.X = None
        self.Y = None
        
        
    def predict(self, X):
        if self.fit_intercept:
            X = np.append(X, np.ones((len(X),1)), axis = 1)
        return np.dot(X,self.coef_)
    
    def predict_backward(self, Y):
        out = np.dot(Y, self.coef_reverse)
        if self.fit_intercept:
            out = out[:,:-1]
        return out
    
    def predict_similarity_totrain(self, X):
        yp = self.predict(X)
        return np.dot(yp,self.Y.T)
    
    def fit(self, X, Y):
        shapeX = np.shape(X)
        shapeY = np.shape(Y)
        X = np.array(X, dtype = float)
        Y = np.array(Y, dtype = float)
        self.X = X
        self.Y = Y
        
        if self.fit_intercept:
            X = np.append(X, np.ones((len(X),1)), axis = 1)
        
        if lasso:
            las = Lasso(fit_intercept = False, alpha = self.alpha)
            las.fit(X, Y)
            self.coef_ = las.coef_.T
            
        
        elif shapeX[0] > shapeX[1]:
            invcovar = np.linalg.pinv(np.dot(X.T, X) + np.eye(len(X.T))*self.alpha)
            self.coef_ = np.dot(np.dot(invcovar,X.T), self.Y)
        else:
            invcovar = np.linalg.pinv(np.dot(X, X.T) + np.eye(len(X))* self.alpha)
            self.coef_ = np.dot(np.dot(X.T, invcovar), self.Y)
        
        if shapeY[0] > shapeY[1]:
            invcovar = np.linalg.pinv(np.dot(Y.T, Y) + np.eye(len(Y.T))*self.alpha)
            self.coef_reverse = np.dot(np.dot(invcovar,Y.T), self.X)
        else:
            invcovar = np.linalg.pinv(np.dot(Y, Y.T) + np.eye(len(Y))* self.alpha)
            self.coef_reverse = np.dot(np.dot(Y.T, invcovar), self.X)
        
        

## K-nearest neighbor predictor (latent representation is input protein representation)
class knearest():
    def __init__(self, k = 1, cutoff = False, mixturefunc = 'equal', similarityfunc = 'cosine', **kwargs):
        self.k = k # k nearest neighbors to combine
        self.mixfunc = mixturefunc # equal, weighted, squared
        self.similarityfunc = similarityfunc # euclidean, cosine, correlation, precomputed
        self.X = None
        self.Y = None
        self.coef_ = None
        self.kwargs = kwargs 
        self.cutoff = cutoff
    
    def choosetop(self, dist):
        # can also use k as maximum distance threshold
        if self.cutoff:
            fullfill = np.where(dist >= float(self.k))
            fullmin = np.argmax(dist, axis = 1)
            combine = [fullfill[1][fullfill[0] == i] for i in range(len(dist))]
            for c, com in enumerate(combine):
                if len(com) == 0:
                    combine[c] = np.array([fullmin[c]])
            return combine
        # otherwise take k-nearest neighbors
        else:
            return np.argsort(-dist, axis = 1)[:, :int(self.k)]
        
    def mixture(self, weights, profiles):
        if self.mixfunc == 'equal':
            return np.mean(profiles, axis = 0)
        if self.mixfunc == 'weighted':
            return np.sum(weights[:, None] * profiles, axis = 0)/np.sum(weights)
        if self.mixfunc == 'squared':
            return np.sum(weights[:, None]**2 * profiles, axis = 0)/np.sum(weights**2)
    
    def predict(self, X):
        dista = self.simfunc(X)
        chosensets = self.choosetop(dista)
        targetmat = self.Y
        
        outprofile = []
        outscore = []
        for t, tes in enumerate(X):
            if len(chosensets[t]) > 1:
                outscore.append(np.sum(self.coef_[chosensets[t]][:,chosensets[t]][np.triu_indices(len(chosensets[t]),1)]*np.dot(dista[t][chosensets[t]].reshape(-1,1), dista[t][chosensets[t]].reshape(1,-1))[np.triu_indices(len(chosensets[t]),1)])/np.sum(np.dot(dista[t][chosensets[t]].reshape(-1,1), dista[t][chosensets[t]].reshape(1,-1))[np.triu_indices(len(chosensets[t]),1)])) # weighted similarity of training set examples
                outprofile.append(self.mixture(dista[t][chosensets[t]], targetmat[chosensets[t]]))
            elif len(chosensets[t]) == 1:
                outscore.append(dista[t][chosensets[t][0]])
                outprofile.append(targetmat[chosensets[t][0]])
            else:
                outprofile.append(np.zeros(len(targetmat[0])))
                outscore.append(0.)
        self.knnscores = np.array(outscore)
        return np.array(outprofile)

    def latent(self, X, reverse = False):
        return X
    
    def simfunc(self,x):
        if self.similarityfunc == 'rbf_cosine':
            return np.exp(-cdist(x, self.X, 'cosine')**2/self.kwargs['rbf_std']**2)
        if self.similarityfunc == 'rbf':
            return np.exp(-cdist(x, self.X, 'euclidean')**2/self.kwargs['rbf_std']**2)
        else:
            return 1.-cdist(x, self.X, self.similarityfunc)
    
    def fit(self, X, Y):
        self.X = X
        self.Y = Y
        if 'rbf_std' in self.kwargs.keys():
            self.kwargs['rbf_std'] = float(self.kwargs['rbf_std'])
        else:
            self.kwargs['rbf_std'] = 0.25
        self.coef_ = self.simfunc(self.X)
        self.coefmax_ = np.amax(self.coef_)
        self.coefmin_ = np.amin(self.coef_)
        

## Non negative Regression for multiple outputs
class NonNegativeRegression():
    def __init__(self, fit_intercept = True):
        self.fit_intercept = fit_intercept
        self.coef_ = None
    def predict(self, X):
        if self.fit_intercept == True:
            X = np.append(X, np.ones((len(X),1)), axis = 1)
        return np.dot(X, self.coef_.T)
    def fit(self, X, Y):
        if self.fit_intercept:
            Xt = np.append(X, np.ones(len(X)).reshape(-1, 1), axis = 1)
        else:
            Xt = np.copy(X)
            
        if len(np.shape(Y)) == 1:
            self.coef_ = nnls(Xt,Y)[0]
            
        else:
            coef_ = []
            for i in range(np.shape(Y)[1]):
                coef_.append(nnls(Xt, Y[:, i])[0])
            coef_ = np.array(coef_)
        self.coef_ = coef_
        self.score = self.predict(X)


        

# embedding or decomposition function
class decomposition():
    def __init__(self, embedmethod = 'svd', perP = 0.9, varrem = 'significant', **kwargs):
        self.embedmethod = embedmethod ### nmf, svd, kernelpca
        # svd: singular value decomposition
        # nmf: non negative matrix factorization
        # kernelpca: svd on sum of covariance kernels of X and Y, *kwargs define Xkerneltype, Ykerneltype, Xkernelgamma, Ykernelgamma
        
        self.perP = perP # number of maintained eigenvectors, or variance explained to keep
        self.varrem = varrem # method to calculate the eigenvalues/variance for each eigenvector: significant, significant_response, reconstruction, threshold, threshold_per_response
        # significant: overall variance explained in both
        # significant_response: variance explained in response 
        # reconstruction: each specificity is reconstructed and only the largest eigenvectors that are required to get a reconstruction correlation until perP are kept
        # threshold: minimum variance that an individual eigenvector can explain of the response

        self.kwargs = kwargs
    
    def reducemat(self, X, Y):
        XY = np.append(X,Y, axis = 1)
        dimsplit = len(Y[1])
        
        if self.embedmethod == 'svd' or self.embedmethod == 'SVD':
            print 'SVD'
            u,s,v = np.linalg.svd(XY, full_matrices=False)
        
        if self.embedmethod == 'sparsepca' or self.embedmethod == 'sparsePCA':
            # put sparse also in kernel method
            if 'l1' in self.kwargs:
                self.kwargs['l1'] = float(self.kwargs['l1'])
            else:
                self.kwargs['l1'] = 1./float(len(XY))
            
            print 'Sparse PCA', self.kwargs['l1']
            spca = SparsePCA(n_components=4,alpha=self.kwargs['l1'], ridge_alpha=0.0001, max_iter=1000, tol=1e-08, method='lars', n_jobs=-1)
            spca.fit(np.dot(XY,XY.T))
            u = spca.components_
            u = u.T
            su = np.sqrt(np.sum(u*u, axis = 0))
            u = u/su
            vt = LinearRegression(fit_intercept = False).fit(u, XY).coef_ 
            s = np.sqrt(np.sum(vt*vt,axis = 0))
            v = (vt/s).T
            
        elif self.embedmethod == 'kernelpca' or self.embedmethod == 'kpca':
            print 'Kernelpca'
            if 'Xkernelgamma' not in self.kwargs.keys():
                self.kwargs['Xkernelgamma'] = 1./np.mean(np.sqrt(np.sum(X*X,axis = 0)))
            else:
                self.kwargs['Xkernelgamma'] = float(self.kwargs['Xkernelgamma'])
            if 'Ykernelgamma' not in self.kwargs.keys():
                self.kwargs['Ykernelgamma'] = 1./np.mean(np.sqrt(np.sum(Y*Y,axis = 0)))
            else:
                self.kwargs['Ykernelgamma'] = float(self.kwargs['Ykernelgamma'])
            if 'Xkerneldegree' not in self.kwargs.keys():
                self.kwargs['Xkerneldegree'] = 2
            else:
                self.kwargs['Xkerneldegree'] = int(self.kwargs['Xkerneldegree'])
            if 'Ykerneldegree' not in self.kwargs.keys():
                self.kwargs['Ykerneldegree'] = 2
            else:
                self.kwargs['Ykerneldegree'] = int(self.kwargs['Ykerneldegree'])
            
                
            
            
            if self.kwargs['Xkerneltype'] == 'linear':
                Xk = np.dot(X, X.T)*self.kwargs['Xkernelgamma']
            elif self.kwargs['Xkerneltype'] == 'polynomial':
                Xk = polynomial_kernel(X, degree = self.kwargs['Xkerneldegree'], gamma = self.kwargs['Xkernelgamma'], coef0 = 0)
                # gamma represents normalization of dot-product between X
                # default = 1/nfeatures
                # maybe use 1/mean(sum(X))
            elif self.kwargs['Xkerneltype'] == 'rbf':
                Xk = rbf_kernel(X, gamma = self.kwargs['Xkernelgamma'])
            elif self.kwargs['Xkerneltype'] == 'laplacian':
                Xk = laplacian_kernel(X, gamma = self.kwargs['Xkernelgamma'])
            elif self.kwargs['Xkerneltype'] == 'chi2':
                Xk = chi2_kernel(X, gamma = self.kwargs['Xkernelgamma'])
            elif self.kwargs['Xkerneltype'] == 'rbf_cosine':
                Xk = np.exp(-cdist(X,X,'cosine')**2/self.kwargs['Xkernelgamma'])
            
            if self.kwargs['Ykerneltype'] == 'linear':
                Yk = np.dot(Y, Y.T)*self.kwargs['Ykernelgamma']
            elif self.kwargs['Ykerneltype'] == 'polynomial':
                Yk = polynomial_kernel(Y, degree = self.kwargs['Ykerneldegree'], gamma = self.kwargs['Ykernelgamma'], coef0 = 0)
            elif self.kwargs['Ykerneltype'] == 'rbf':
                Yk = rbf_kernel(Y, gamma = self.kwargs['Ykernelgamma'])
            elif self.kwargs['Ykerneltype'] == 'laplacian':
                Yk = laplacian_kernel(Y, gamma = self.kwargs['Ykernelgamma'])
            elif self.kwargs['Ykerneltype'] == 'chi2':
                Yk = chi2_kernel(Y, gamma = self.kwargs['Ykernelgamma'])
            elif self.kwargs['Ykerneltype'] == 'rbf_cosine':
                Yk = np.exp(-cdist(Y,Y,'cosine')**2/self.kwargs['Ykernelgamma'])

            if 'Multiplykernel' in self.kwargs.keys():
                XYk = Xk*Yk
            else:
                XYk = Xk + Yk
        
            u,s2,ut = np.linalg.svd(XYk, full_matrices=False)
            s = np.sqrt(s2)
            # Don't have to do this, could do the predictions only using u and s
            vt = LinearRegression(fit_intercept = False).fit(u*s, XY).coef_ 
            v = vt.T
            
        if self.varrem == 'significant_response':
            # determines the variance of the eigenvector in Y direction
            sig = np.sum(v[:, -dimsplit:]*v[:, -dimsplit:], axis = 1)*s**2
            
        elif self.varrem == 'significant':
            # represents the overall variance represented by the eigenvector
            sig = s**2
        
        if 'significant' in self.varrem:
            sortsig = np.argsort(-sig)
            if self.perP > 1:
                perP = int(self.perP)
            else:
                perP = max(1, len(np.where(np.cumsum(sig[sortsig])/np.sum(sig)<self.perP)[0]))
            print 'maintained eigenvectors', perP
            sortsig = sortsig[:perP]
            print 'deleted', np.delete(np.arange(perP), sortsig), 'maint', sortsig[sortsig > perP]
        
        elif self.varrem == 'reconstruction':
            #uses the most important eigenvalues for each training example
            sortsig = []
            reconst = []
            #pst = 0
            for jj, uv in enumerate(u):
                sortu = np.argsort(-np.absolute(uv*s))
                reconstz = np.zeros(len(Y[0]))
                for su in sortu:
                    sortsig.append(su)
                    reconstz += uv[su]*s[su]*v[su, -dimsplit:]
                    checkp = pearsonr(reconstz, Y[jj])[0]
                    #reconst.append(checkp)
                    if checkp >= self.perP:
                        break
                #print np.sort(sortsig[pst:])
                #print np.sort(reconst[pst:])
                #pst = len(sortsig)
            sortsig = np.unique(sortsig)
            print 'maintained eigenvectors', sortsig
            print 'removed eigenvectors', np.delete(np.arange(len(u)), sortsig)
        
        elif self.varrem == 'threshold':
            sig = np.sum(v[:, -dimsplit:]*v[:, -dimsplit:], axis = 1)*s**2
            sortsig = np.where(sig/np.sum(sig) >= float(self.perP))[0]
            print 'maintained eigenvectors',sortsig, len(sortsig), len(s)
 
        
        
        removed = np.delete(np.arange(len(s), dtype = int), sortsig)
        removed = removed[s[removed]>1e-10]
        
        print len(s)
        ur = u[:,removed[:-1]]
        vr = v[removed[:-1]]
        sr = s[removed[:-1]]
        
        u = u[:,sortsig]
        v = v[sortsig]
        s = s[sortsig]
        
        print len(s), len(sr)
        return u,s,v, ur, sr, vr


class jple():
    def __init__(self, latentmap = 'SVD', perP=0.9, varrem='significant_response', mapalg='lsq', decoder = 'global', coef_ = None, xloc = [0,0], yloc = [0,0], **kwargs):
        self.perP = perP # number of eigenvalues, or fraction of variance explained by eigenvectors
        self.varrem = varrem 
        # significant, significant_response, reconstruction, method to determine the final set of eigenvectors
        self.coef_ = coef_ # variables with embedding and eigenvectors
        if decoder in ['local', 'global']:
            self.decoder = decoder
        else:
            self.decoder = 'global'
            print 'Decoder not specifified, set to global'
         
        if self.coef_ is not None:
            self.Wtrain = self.coef_[0]*self.coef_[1]
        
        self.latentalg = decomposition(embedmethod = latentmap, perP = perP, varrem = varrem, **kwargs) # parameters for embedding method: svd, kernelpca, jointkernelpca, jkpca
        self.latentmap = latentmap
        if mapalg == 'lsq' or mapalg == "LinearRegression" or mapalg == 'lr':
            self.mapalg = LinearRegression(fit_intercept = False)
        elif mapalg == 'nnls':
            self.mapalg = NonNegativeRegression(fit_intercept = False)
        elif mapalg == 'lasso':
            self.mapalg = Lasso(fit_intercept = False, alpha = float(kwargs['regalpha']))
        elif mapalg == 'ridge':
            self.mapalg = Ridge(fit_intercept = False, alpha = float(kwargs['regalpha']))
        elif mapalg == 'lars':
            self.mapalg = Lars(fit_intercept = False, n_nonzero_coefs = int(kwargs['regalpha']))
        else:
            print mapalg, 'not understood'
        self.xloc = xloc # will be defined at fit: start and end of X variable
        self.yloc = yloc # will be defined at fit: start and end of Y variable
        
        if 'xmeans' in kwargs.keys():
            self.xmean = kwargs['xmeans']
        if 'ymeans' in kwargs.keys():
            self.ymean = kwargs['ymeans']
        
        ## parameters for local reconstruction if wanted
        # default values
        self.k = 0.01 # only considers close proteins when similarty over 0.01
        self.kcut = True # works with cutoff rather than a fixed number of surrounding proteins
        self.mfunc = 'weighted' # mixture of profiles is weighted by their similarity
        self.sfunc = 'rbf_cosine' # similarity function is an gaussian on the cosine similarity 
        if 'k' in kwargs.keys():
            self.k = kwargs['k']
            del kwargs['k']
        if 'cutoff' in kwargs.keys():
            self.kcut = check(kwargs['cutoff'])
            del kwargs['cutoff']
        if 'mixturefunc' in kwargs.keys():
            self.mfunc = kwargs['mixturefunc']
            del kwargs['mixturefunc']
        if 'similarityfunc' in kwargs.keys():
            self.sfunc = kwargs['similarityfunc']
            del kwargs['similarityfunc']
        
        self.kwargs = kwargs
        if not 'rbf_std' in kwargs.keys():
            self.kwargs['rbf_std'] = 0.2   # standard deviation for similarity function, cosine similar ~0.2
        self.latentdistances = None
        self.bestlatentdistances = None
        self.confidence = None
        
    def fit(self, X, Y):
        self.X = np.copy(X)
        self.Y = np.copy(Y)
        if self.coef_ is None:
            self.xloc = [0, len(X[0])]
            self.yloc = [len(X[0]), len(X[0]) + len(Y[0])]
            
            # better approximation for mean can be provided, for example from all available sequences
            if 'xmeans' not in self.kwargs.keys():
                self.xmean = np.mean(X, axis = 0)
            
            if 'ymeans' not in self.kwargs.keys():
                self.ymean = np.mean(Y, axis = 0)

            X = X - self.xmean
            Y = Y - self.ymean
            
            ue, se, ve, ur, sr, vr = self.latentalg.reducemat(X, Y)
            self.coef_ =[ue, se, ve]
            self.rcoef_ =[ur, sr, vr]
        self.Wtrain = self.coef_[0]*self.coef_[1]
        
    def predict(self, X, reverse= False):
        X = np.array(X)
        wembed = self.latent(X, reverse = reverse)
        
        if reverse:
            ypred = self.predict_from_latent(wembed, direction = 'backward', add_mean = True)
        else:
            ypred = self.predict_from_latent(wembed, direction = 'forward', add_mean = True)
        return ypred
    
    def local_reconstruction(self,latrep):
        fr = knearest(k = self.k, cutoff = self.kcut, mixturefunc = self.mfunc, similarityfunc = self.sfunc, **self.kwargs)
        fr.fit(self.Wtrain, self.Y)
        outp = fr.predict(latrep)
        self.confidence = fr.knnscores
        return outp
        
    def global_reconstruction(self, latrep, decoder, add_mean, tmean):
        outp = np.dot(latrep, decoder)
        if add_mean:
            outp += tmean
        self.confidence = self.bestlatentdistances
        return outp
        

    def predict_from_latent(self, latrep, direction = 'forward', add_mean = False):
        if direction == 'forward':
            if self.decoder == 'global':
                outp = self.global_reconstruction(latrep, self.coef_[2][:, self.yloc[0]:self.yloc[1]], add_mean, self.ymean)
            elif self.decoder == 'local':
                outp = self.local_reconstruction(latrep)
                
        elif direction == 'backward':
            if self.decoder == 'global':
                outp = self.global_reconstruction(latrep, self.coef_[2][:, self.xloc[0]:self.xloc[1]], add_mean, self.xmean)
        return outp
            
     
    def latent(self, X, reverse = False):
        
        if reverse:
            X = X - self.ymean
            vx = self.coef_[2][:, self.yloc[0]:self.yloc[1]]
        else:
            X = X - self.xmean
            vx = self.coef_[2][:, self.xloc[0]:self.xloc[1]]        
        self.mapalg.fit(vx.T, X.T)
        wembed = self.mapalg.coef_
        if self.latentdistances is None:
            self.latentdistances = cdist(wembed, self.Wtrain, 'cosine')
            self.bestlatentdistances = np.amin(self.latentdistances , axis = 1)
        #print wembed
        return wembed

def check(val):
    return val == 'True' or val == 'true'
        

# load the protein features 

proteinfeatures = sys.argv[1]
pobj = np.load(proteinfeatures)
P = pobj['features']
sequencefeatcompare = pobj['kmers']
protnames = pobj['protnames']
pnames = pobj['expnames']
outname = os.path.splitext(proteinfeatures)[0]

# load the coefficients of the trained model to embed protein features
recoefile = sys.argv[2]
Loadcoef = np.load(recoefile, allow_pickle = True) # coefficients of trained model can be loaded
coef_= Loadcoef['coef_']
xmean= Loadcoef['xmean']
ymean= Loadcoef['ymean']
sequencefeatures = Loadcoef['sequencefeatures']
outname += os.path.splitext(os.path.split(recoefile)[1])[0]


# align protein features to features in trainingset
if not np.array_equal(sequencefeatcompare, sequencefeatures):
    print 'Protein features differ', len(sequencefeatcompare), len(sequencefeatures)
    protnumbers = np.sum(P, axis = 1)
    if np.array_equal(sequencefeatcompare[np.isin(sequencefeatcompare, sequencefeatures)], sequencefeatures[np.isin(sequencefeatures, sequencefeatcompare)]):
        print np.shape(P)
        Pnew = np.zeros((len(sequencefeatures), len(P[0])))
        Pnew[np.isin(sequencefeatures, sequencefeatcompare)] = P[np.isin(sequencefeatcompare, sequencefeatures)]
        P = Pnew
    else:
        print 'Protein features could not be aligned'
        sys.exit()

P = P.T
print np.shape(P)

# initiate jple model with loaded parameters
params = {'xmeans':xmean, 'ymeans':ymean}
pr = jple(coef_ = coef_, xloc = [0,len(sequencefeatures)], **params)

Wtrain = pr.Wtrain

permuts = 100

distlat = []
lengthlat = []
varlat = []

if '--fraction' in sys.argv:
    outname += '-fraction'
    fcuts = np.array(sys.argv[sys.argv.index('--fraction')+1].split(','), dtype = float)
    outname += '-'.join(fcuts.astype(str))

# iterate over all proteins in protein set
for p, pna in enumerate(pnames):
    porig = P[p]
    lockmer = np.where(porig>0)[0]
    numkmers = len(lockmer)
    print p, '/', len(pnames), numkmers
    # choose random percentage between 30 and 70% of k-mers that should be deleted
    randomperc = (np.random.uniform(0.3,0.7, permuts)*float(numkmers)).astype(int)
    Ptest = [porig/np.sqrt(np.sum(porig**2))]
    # generate a hundred permutations of that protein sequence by setting a random permutation to 0
    for i in range(permuts):
        ptest = np.copy(porig)
        ptest[np.random.permutation(lockmer)[:randomperc[i]]] = 0
        ptest = ptest/np.sqrt(np.sum(ptest**2))
        Ptest.append(ptest)
    Ptest = np.array(Ptest)
    # embed into latent space
    Wtest = pr.latent(Ptest)
    
    wlength = np.sqrt(np.sum(Wtest[0]**2))
    lengthlat.append(wlength)
    # Measure distance of all permuted latent vectors to original one / assumption is that they should still be the close to each other. 
    bestrain = np.amin(cdist(Wtest[[0]],Wtrain, 'cosine'))
    distlat.append(bestrain)
    
    cosines = cdist(Wtest[[0]], Wtest[1:], 'cosine')
    if '--fraction' in sys.argv:
        var = []
        for fcut in fcuts:
            var.append(1.-np.mean(cosines<=fcut))
        varlat.append(var)
    else:
        varlat.append(np.mean(cosines))
    
    print wlength, bestrain, varlat[-1]

varlat = np.array(varlat)
lengthlat = np.array(lengthlat)
distlat = np.array(distlat)

'''
if '--fraction' in sys.argv:
    np.savetxt(outname + '_variance.txt', np.concatenate([[pnames[:p]], [distlat], [lengthlat], varlat.T], axis = 0).T, fmt = '%s')
    varlat = varlat.T
else:
    np.savetxt(outname + '_variance.txt', np.array([pnames[:p], distlat, lengthlat, varlat]).T, fmt = '%s')
'''

varlat = varlat.T
fig = plt.figure(figsize = (8,3), dpi = 200)
ax = fig.add_subplot(121)

if '--fraction' in sys.argv:
    for v, var in enumerate(varlat):
        if '--meanx' in sys.argv:
            x = []
            y = []
            yerr = []
            bins = 8
            ol = 2
            sorting = np.argsort(distlat)
            wz = len(distlat)/(ol*bins) + int(len(distlat)%(ol*bins)>0)
            bins = np.arange(bins*ol, dtype = int)
            bins = bins[bins*wz<len(distlat)]
            print bins, len(distlat)
            for i in bins:
                print i*wz,(i+ol)*wz
                x.append(np.mean(distlat[sorting][i*wz:(i+ol)*wz]))
                y.append(np.mean(var[sorting][i*wz:(i+ol)*wz]))
                yerr.append(np.percentile(var[sorting][i*wz:(i+ol)*wz], [25,75]))
            yerr = np.array(yerr)
            ax.plot(x,y, label=str(fcuts[v]), marker = 'o')
            ax.fill_between(x, y, yerr[:,0], color = 'silver', alpha = 0.3, lw = 0)
            ax.fill_between(x, yerr[:,1],y, color = 'silver', alpha = 0.3, lw = 0)
        else:
            ax.scatter(distlat, var, label=str(fcuts[v]))
else:
    ax.scatter(distlat, varlat)

ax.set_ylim([0,np.amax(varlat)])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('min e-dist')
if '--fraction' in sys.argv:
    ax.set_ylabel('Fraction outside cutoff')
else:
    ax.set_ylabel('Mean cosine distance (100 shuffle)')
ax.legend()


ax2 = fig.add_subplot(122)
if '--fraction' in sys.argv:
    for v, var in enumerate(varlat):
        if '--meanx' in sys.argv:
            x = []
            y = []
            yerr = []
            bins = 8
            ol = 2
            sorting = np.argsort(lengthlat)
            wz = len(lengthlat)/(ol*bins) + int(len(lengthlat)%(ol*bins)>0)
            bins = np.arange(bins*ol, dtype = int)
            bins = bins[bins*wz<len(lengthlat)]
            print bins, len(lengthlat)
            for i in bins:
                print i*wz,(i+ol)*wz
                x.append(np.mean(lengthlat[sorting][i*wz:(i+ol)*wz]))
                y.append(np.mean(var[sorting][i*wz:(i+ol)*wz]))
                yerr.append(np.percentile(var[sorting][i*wz:(i+ol)*wz], [25,75]))
            yerr = np.array(yerr)
            ax2.plot(x,y, label=str(fcuts[v]), marker = 'o')
            ax2.fill_between(x, y, yerr[:,0], color = 'silver', alpha = 0.3, lw = 0)
            ax2.fill_between(x, yerr[:,1],y, color = 'silver', alpha = 0.3, lw = 0)
        else:
            ax2.scatter(lengthlat, var, label=str(fcuts[v]))
else:
    ax2.scatter(lengthlat, varlat)

ax2.set_xscale('log')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlabel('Length latent')
ax2.set_ylim([0,np.amax(varlat)])
#ax2.set_ylabel('Mean cosine distance (100 shuffle)')
ax2.legend()

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4) 
fig.tight_layout()

if '--savefig' in sys.argv:
    print 'Saved as', outname+'.svg'
    fig.savefig(outname+'.svg', dpi = 250, bbox_inches = 'tight')
plt.show()


# Use clusters from real latent space, with additional proteins in those clusters and exclude them from training. Determine then false positive and negative rate after traningn without cluster






        
