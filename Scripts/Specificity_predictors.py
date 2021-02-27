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
from joblib import Parallel, delayed
import multiprocessing
from scipy.optimize import minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.metrics.pairwise import polynomial_kernel, rbf_kernel, laplacian_kernel, chi2_kernel
from scipy.optimize import NonlinearConstraint
from sklearn.ensemble import RandomForestRegressor

## Affintiy regression determines W from YDWP^T = YY^T or DWP^T = Y^T 
class affinity_regression():
    def __init__(self, alpha = 1., center = True, redX = 0.9, redY = 0.9, ridge = True, D = None, kronecker = False, xmean = 0., ymean = 0., coef_ = None):
        self.alpha = alpha # regularization 
        self.coef_ = coef_ # predictor
        # If mapping of 7-mers to lower dimensional space not give just keep in high space
        if D is None:
            self.D = np.eye(len(Y[0])) # This is 1-6mer specificity mapping to 7-mer mapping
        else:
            self.D = D
        self.center = center # center X and Y 
        self.redX = redX # maintained variance for X
        self.redY = redY # maintained variance for Y
        self.kronecker = kronecker # fitting similarities YY^T from kronecker product between pca X and pca Y.
        self.ridge = ridge
        self.xmean = xmean
        self.ymean = ymean
    
    def fit(self, X, Y):
        
        shapeX = np.shape(X)
        shapeY = np.shape(Y)
        X = np.array(X, dtype = float)
        Y = np.array(Y, dtype = float)

        # center X and Y before processing
        if self.center:
            if self.xmean == 0.:
                self.xmean = np.mean(X, axis = 0)
            if self.ymean == 0.:
                self.ymean = np.mean(Y, axis = 0)
            X = X - self.xmean
            Y = Y - self.ymean
        else:
            self.xmean = np.zeros(len(X[0]))
            self.ymean = np.zeros(len(Y[0]))
        
        # store training specificities for later
        self.Y = Y
        # deconvolve all binding specificities into eigenvector space
        uy, sy, ey = np.linalg.svd(self.Y, full_matrices = False)
        # coefficient of training examples and eigenvectors for reconstruction
        self.c = uy*sy 
        self.ey = ey
        
        # left multiply both sides with Y
        YD = np.dot(Y,self.D)
        Y = np.dot(Y, Y.T)
        
        if self.redY is not None:
            # use pca to make problem solvable on standard computer
            uy,sy,vy = np.linalg.svd(YD, full_matrices=False)
            scy = np.cumsum(sy)/np.sum(sy)
            Sky = scy < self.redY
            uy =uy[:, Sky]
            sy =sy[Sky]
            vy = vy[Sky]
            YD = uy
            
        if self.redX is not None:
            # use pca to make problem solvable on standard computer
            ux,sx,vx = np.linalg.svd(X, full_matrices=False)
            scx = np.cumsum(sx)/np.sum(sx)
            Skx = scx < self.redX
            ux =ux[:, Skx]
            sx =sx[Skx]
            vx = vx[Skx]
            X = ux
    
        # vectorize YDWP by kron(YD,P): result is features f1xf2 for all protein pairs n1xn2
        # This approach is performed in the original code but memory usage grow exponentially with larger k-mers and larger training set.
        if self.kronecker:
            print np.shape(X), np.shape(YD)
            Xin = np.kron(X, YD)
            print np.shape(X)
            # vectorize target: target is specificity similarity represented by dot product
            Y = Y.reshape(len(Xin))
            # adjust alpha to length of data
            self.alpha = self.alpha/float(len(X))
        else: 
            # If not run on a computer with more than 128GB of ram (limited memory) memory can be saved by directly fitting the matrix W by transforming Y with the pseudoinverse of YD: WP = D^-1Y^-1YY^T = D^-1Y^T
            self.alpha = self.alpha/float(len(X)**2)
            Y = np.dot(np.linalg.pinv(YD),Y).T
            Xin = np.copy(X)
        
    
        if self.ridge:
            # Lasso will take forever
            print 'Ridge',self.alpha 
            las = Ridge(fit_intercept = False, alpha = self.alpha)
            las.fit(Xin, Y)
            if self.kronecker == False:
                self.coef_ = las.coef_.T
        # Explicitly computing the solution of regression seems to be way faster!
        elif shapeX[0] > shapeX[1]:
            # Regular solution of linear regression 
            invcovar = np.dot(np.linalg.pinv(np.dot(Xin.T, X) + np.eye(len(Xin.T))*self.alpha), Xin.T)
            self.coef_ = np.dot(invcovar, Y)
        else:
            # Dual form solution of regression if more features than data point available. 
            invcovar = np.dot(Xin.T, np.linalg.pinv(np.dot(Xin, Xin.T) + np.eye(len(Xin))* self.alpha))
            self.coef_ = np.dot(invcovar, Y)
        

        # if knonecker was used, reshape to original form
        if self.kronecker:
            self.coef_ = self.coef_.reshape(len(X[0]), len(YD[0]))
        
        # reconstruction of full W as described in supplementary of AR
        if self.redY is not None:
            self.coef_ = np.dot(self.coef_/sy,vy)
            print np.shape(vy)
        if self.redY is not None:  
            self.coef_ = np.dot(vx.T*sx**-1, self.coef_)
            print np.shape(vx)
        
    def predict_similarity_totrain(self, X):
        if self.center:
            X = X - self.xmean
        spec = np.dot(np.dot(X, self.coef_), self.D.T)
        print np.shape(spec)
        return np.dot(spec, self.Y.T)
    
    def predict(self, X):
        A = self.predict_similarity_totrain(X)
        # find coefficients for test examples that produce A if multiplied by coefficients from test examples
        ct = np.dot(A, np.linalg.pinv(self.c.T))
        # reconstruct specificity with these test coefficients
        return np.dot(ct, self.ey) + self.ymean
    
    def predict_backward(self, Y):
        if self.center:
            Y = Y - self.ymean
        return np.dot(np.dot(Y, self.D), self.coef_.T)



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
    else:
        protnames = protnames[0]
    
    if os.path.isfile(proteinset) or isinstance(proteinset, (np.ndarray, list)):
        if os.path.isfile(proteinset):
            pset = np.genfromtxt(proteinset, dtype = str, delimiter = '\n')
        else:
            pset = proteinset
        if len(np.shape(pset)) == 0:
            pset = [pset]
        pmask = []
        for p, ps in enumerate(pset):
            if ps in protnames:
                pmask.append(list(protnames).index(ps))
            else:
                print 'Attention', ps, 'not in protnames'
        if Y is not None:
            Y = Y[pmask]
            ynames = ynames[pmask]
        P = P[pmask]
        protnames = protnames[pmask]

    return Y, ynames, Zfeatures, P, protnames, sequencefeatures



coefgiven = False
compute = False
# provide outname for run
if "--outname" in sys.argv:
    Wname = sys.argv[sys.argv.index("--outname")+1]
    print "Matrix outname is", Wname
    

# model runs fast enough to rerun on training data before reconstructions
elif "--recompute" in sys.argv:
    # to impute z-score profiles for unmeasured protein just provide a protein profile and make sure that recprofiles and proteinset are None
    # To perform test on additional measurements/proteins, for example from different assay add recprofiles
    # If protein features are stored in big file and only a few need to be accessed, provide proteinset
    
    recoefile = sys.argv[sys.argv.index("--recompute")+1]
    if os.path.isfile(recoefile):
        Loadcoef = np.load(recoefile, allow_pickle = True) # coefficients of trained model can be loaded
        coef_= Loadcoef['coef_']
        xmean= Loadcoef['xmean']
        ymean= Loadcoef['ymean']
        trainprots = Loadcoef['trainprots']
        testprots = Loadcoef['testprots']
        coefgiven = True
        Zfeatures, sequencefeatures = Loadcoef['Zfeatures'], Loadcoef['sequencefeatures']
    recproteinfeatures = sys.argv[sys.argv.index("--recompute")+2] ## Protein sequence features
    recprofiles = sys.argv[sys.argv.index("--recompute")+3] # Test specificities
    proteinset =  sys.argv[sys.argv.index("--recompute")+4] # Set of proteins to be tested or reconstructed, use None if not reconstructed
    if proteinset == 'testprots':
        proteinset = testprots
    Ytest, tpnames, Zfeatcompare, Ptest, testprots, sequencefeatcompare = readin(recproteinfeatures, recprofiles, proteinset)
    readtrain = False
    if Ytest is not None or recprofiles == 'local':
        readtrain = True
    if Ytest is None:
        Ytest = np.ones((2,len(Zfeatures))) 
    Wname = sys.argv[sys.argv.index("--recompute")+5] # Use '' if don't want to define
    if os.path.isfile(recoefile):
        Wname += os.path.splitext(os.path.split(recoefile)[1])[0]
    compute = True
else: 
    print "Please define whether to recompute results or save the generated matrix" 
    sys.exit()



if not coefgiven:
    profiles = sys.argv[1] 
    proteinfeatures = sys.argv[2]
    Y, datnames, Zfeatures, P, pnames, sequencefeatures = readin(proteinfeatures, profiles)
elif readtrain: # read in training profiles if we want to test something with Ytest
    profiles = sys.argv[1]
    proteinfeatures = sys.argv[2]
    Y, datnames, Zfeatures2, P, pnames, sequencefeatures2 = readin(proteinfeatures, profiles, trainprots)
    if not np.array_equal(sequencefeatures, sequencefeatures2):
        P = P[:,np.isin(sequencefeatures2,sequencefeatures)]
        if not np.array_equal(sequencefeatures, sequencefeatures2[np.isin(sequencefeatures2, sequencefeatures)]):
            print 'Training protein feautures not equal to model features'
            sys.exit()
else:
    # dummy data that is not used but too lazy to put if everywhere
    Y = np.ones((2,len(Zfeatures)))
    P = np.ones((2, len(sequencefeatures)))
    datnames = np.array(['f1', 'f2'])
    
if compute:
    #compare sequence and Zscore features between training and recompute set
    if not np.array_equal(sequencefeatcompare, sequencefeatures):
        print 'Protein features differ', len(sequencefeatcompare), len(sequencefeatures)
        protnumbers = np.sum(Ptest, axis = 1)
        if np.array_equal(sequencefeatcompare[np.isin(sequencefeatcompare, sequencefeatures)], sequencefeatures[np.isin(sequencefeatures, sequencefeatcompare)]):
            Pnew = np.zeros((len(sequencefeatures), len(Ptest)))
            Pnew[np.isin(sequencefeatures, sequencefeatcompare)] = Ptest.T[np.isin(sequencefeatcompare, sequencefeatures)]
            Ptest = Pnew.T
        else:
            print 'Protein features could not be aligned'
            sys.exit()
        
        print 'Num protfeatures diff', protnumbers - np.sum(Ptest, axis = 1)
    if Zfeatcompare is not None:    
        if not np.array_equal(Zfeatcompare, Zfeatures):
            print 'Zscore features differ', len(Zfeatcompare), len(Zfeatures)
            Ytest = Ytest[:, np.isin(Zfeatcompare, Zfeatures)]
            


        


# split data into training and test set
if compute:
    Ptrain = P[:]
    Ytrain = Y[:]
    trainprots = datnames

elif "--trainingset" in sys.argv:
    numset = sys.argv[sys.argv.index("--trainingset")+1]
    trainfile = sys.argv[sys.argv.index("--trainingset")+2]
    print "trainset:", trainfile, numset
    tobj = open(trainfile, 'r')
    tlines = tobj.readlines()
    for i, tline in enumerate(tlines):
        if tline[:6]=="###Set" and tline.strip().split()[-1] == numset:
            trainprots = np.array([tlines[i+2].strip().split()])
            if tlines[i+3][:6] == '##Test':
                trainprots = np.append(trainprots, [tlines[i+2].strip().split()], axis = 0).T
                testprots = np.array([tlines[i+4].strip().split(), tlines[i+4].strip().split()]).T
            else:
                trainprots = np.append(trainprots, [tlines[i+3].strip().split()], axis = 0).T
                testprots = np.array([tlines[i+5].strip().split(), tlines[i+6].strip().split()]).T
            break
    ### sort into training and test set
    if len(np.intersect1d(trainprots[:,0], datnames)) !=0:
        trainprots = trainprots[:,0]
        testprots = testprots[:,0]
    else: 
        trainprots = trainprots[:,1]
        testprots = testprots[:,1]
    
    if compute:
        trainprots = trainprotsin[:]

    indexestrain = []
    indexestest = []
    for i, ina in enumerate(datnames):
        if ina in trainprots:
            indexestrain.append(i)
        if ina in testprots:
            indexestest.append(i)
    indexestrain = np.array(indexestrain)
    indexestest = np.array(indexestest)    
    
    Ytrain = Y[indexestrain]
    Ptrain = P[indexestrain]
    Ytest = Y[indexestest]
    Ptest = P[indexestest]
    trainprots = datnames[indexestrain]
    testprots = datnames[indexestest]

elif '--proteinset' in sys.argv:
    psetred = np.genfromtxt(sys.argv[sys.argv.index('--proteinset')+1], dtype = str)
    mask = np.isin(datnames, psetred)
    Ptrain = Ptest = P[mask]
    Ytrain = Ytest = Y[mask]
    trainprots = testprots = datnames[mask]

elif "--crossvalidation" in sys.argv:
    ncross = int(mprot/float(sys.argv[sys.argv.index("--crossvalidation")+1]))
    print sys.argv[sys.argv.index("--crossvalidation")+1]+"-fold crossvalidation" 
    Ptest = P[-ncross:]
    Ptrain = P[:mprot-ncross]
    Ytest = Y[-ncross:]
    Ytrain = Y[:mprot-ncross]
    if sim == False:
        trainprots = datnames[-ncross:]
        testprots = datnames[:mprot-ncross]
else:
    Ptrain = Ptest = P[:]
    Ytrain = Ytest = Y[:]
    trainprots = testprots = datnames



    

# remove Pfeatures that don't occur in training set
ptmask = np.sum(Ptrain, axis = 0) > 0
Ptrain = Ptrain[:,ptmask].astype(np.float)
Ptest = Ptest[:,ptmask].astype(np.float)
sequencefeatures = np.array(sequencefeatures)[ptmask]



# check if any protein has zero features
pmask = np.sum(Ptrain*Ptrain, axis = 1) != 0
if np.sum(pmask) != len(Ptrain):
    print 'Missing', trainprots[~pmask]
    Ptrain = Ptrain[pmask]
    trainprots = trainprots[pmask]
    Ytrain = Ytrain[pmask]
pmask = np.sum(Ptest*Ptest, axis = 1) != 0
if np.sum(pmask) != len(Ptest):
    print 'Missing', testprots[~pmask]
    Ptest = Ptest[pmask]
    testprots = testprots[pmask]
    Ytest = Ytest[pmask]




### Options to adjust protein features
# if you have a list of k-mers that you want to keep and remove the rest
if '--filterpfeatures' in sys.argv:
    featurefile = np.genfromtxt(sys.argv[sys.argv.index('--filterpfeatures')+1], dtype = str)
    featuremask = np.isin(sequencefeatures, featurefile)
    sequencefeatures = sequencefeatures[featuremask]
    print 'Pfilter', len(P[0])
    Ptrain = Ptrain[:,featuremask]
    Ptest = Ptest[:,featuremask]
    print 'to', len(P[0])
        


if '--atleastkfeaturessum' in sys.argv:
    kmini = float(sys.argv[sys.argv.index('--atleastkfeaturessum')+1])
    kminmask = np.sum(Ptrain, axis = 0) >= kmini
    sequencefeatures = np.array(sequencefeatures)[kminmask]
    Ptrain = Ptrain[:,kminmask]
    Ptest = Ptest[:,kminmask]
    print 'All features removed that appeared less than', kmini, 'times!', len(sequencefeatures), 'remain'

if "--cutP" in sys.argv:
    print "P contains binary values"
    Ptrain[Ptrain>0] = 1
    Ptest[Ptest>0] = 1

if '--TfilterP' in sys.argv:
    zmin = float(sys.argv[sys.argv.index('--TfilterP')+1])
    print 'Use only significant features'
    zstats = np.zeros(len(sequencefeatures))
    for s, sfeat in enumerate(sequencefeatures):
        if s % 1000 == 1:
            print s, np.percentile(zstats[:s], [1,10,20,50,80,90,99])
        Yp = Ytrain[Ptrain[:,s]>0]
        Yn = Ytrain[Ptrain[:,s]<=0]
        npo = float(len(Yp))
        nn = float(len(Yn))
        meanp = np.mean(Yp, axis = 0)
        meann = np.mean(Yn, axis = 0)
        stdp = np.var(Yp, axis = 0)
        stdn = np.var(Yn, axis = 0)
        tstat = np.amax((meanp-meann)/(np.sqrt((npo+nn)/(npo*nn))*np.sqrt(((npo-1.)*stdp+(nn-1.)*stdn)/(npo+nn-2.))))
        zstats[s] = np.amax(tstat)
    keep = []
    minp = np.percentile(zstats, int(zmin))
    for p, ptrain in enumerate(Ptrain):
        perszmin = np.sort(zstats[ptrain > 0])[-int(np.sum(ptrain>0)*(100.-zmin)/200.)-1]
        keep.append(np.where((ptrain>0)*(zstats>=min(perszmin, minp)))[0])
    keep = np.unique(np.concatenate(keep))
    
    sequencefeatures = sequencefeatures[keep]
    Ptrain = Ptrain[:, keep]
    Ptest = Ptest[:, keep]

    
if '--normP2' in sys.argv:
    print "normP2"
    Ptrain = (Ptrain.T/np.sqrt(np.sum(Ptrain*Ptrain, axis = 1))).T
    Ptest = (Ptest.T/np.sqrt(np.sum(Ptest*Ptest, axis = 1))).T
    
if '--normP1' in sys.argv:
    print "normP1"
    Ptrain = (Ptrain.T/np.sum(np.absolute(Ptrain), axis = 1)).T
    Ptest = (Ptest.T/np.sum(np.absolute(Ptest), axis = 1)).T
    print Ptrain[np.nonzero(Ptrain)]







# manipulate Y, so that order at protein level is not changed/specificity not changed
if '--normY2' in sys.argv:
    print "normY2"
    Ytrain = (Ytrain.T/np.sqrt(np.sum(Ytrain*Ytrain, axis = 1))).T
    Ytest = (Ytest.T/np.sqrt(np.sum(Ytest*Ytest, axis = 1))).T
    
if '--quartileY' in sys.argv:
    Y = np.append(Ytrain, Ytest, axis = 0)
    Ysort = np.argsort(Y, axis = 1)
    Ymean = Y
    for i, y in enumerate(Ymean):
        Ymean[i] = y[Ysort[i]]
    Ymean = np.mean(Ymean, axis = 0)
    Ytrain = Ymean[np.argsort(np.argsort(Ytrain, axis = 1), axis = 1)]
    Ytest = Ymean[np.argsort(np.argsort(Ytest))]
    








    
if "--JPLE" in sys.argv:
    latent = sys.argv[sys.argv.index("--JPLE")+1] # nmf, svd, kernelpca
    perP = float(sys.argv[sys.argv.index("--JPLE")+2]) # number of maintained features or variance of eigenvetors
    reduct = sys.argv[sys.argv.index("--JPLE")+3] # significant, significant_response, reconstruction
    mapping = sys.argv[sys.argv.index("--JPLE")+4] # lsq, lasso, ridge, nnls
    decoding = sys.argv[sys.argv.index("--JPLE")+5] # global or local, if local also give 'cutoff:True', 'k:1', 'mixturefunc:weighted', and 'similarityfunc:rbf_cosine', additionally one can specify the standard deviation of the similarityfunc, rbf_std:0.02 and whether we want to look at a markov field in the latent space, network:2 or network:11 
    # (for details see knearestneighbor class
    ### additional parameter sets for kernelpca: Xkerneltype, Xkernelgamma, Ykerneltype, Ykernelgamma
            
    parmapping = []
    params = {}
    if len(sys.argv) > sys.argv.index("--JPLE")+6:
        if ':' in sys.argv[sys.argv.index("--JPLE")+6]:
            parmapping = sys.argv[sys.argv.index("--JPLE")+6]
            if ',' in parmapping:
                parmapping = parmapping.split(',') 
            else:
                parmapping = [parmapping]
            params = {}
            for pp in parmapping:
                pp = pp.split(':')
                params[pp[0]] = pp[1]
        else:
            params = {}
            parmapping = ['']
        if 'xmeans' in params.keys():
            if check(params['xmeans']):
                params['xmeans'] = np.mean(np.append(Ptrain, Ptest, axis = 0), axis = 0)
            else:
                del params['xmeans']
     
    if '--recompute' not in sys.argv:
        Wname += '_'+latent+reduct+str(perP)+'_map'+mapping+'_dec'+decoding
        if len(parmapping) > 0:
            Wname += '_'+'_'.join(np.array(parmapping)).replace(':','-')
    else:
        Wname+='_recdec'+decoding
    print Wname
    model = latent+reduct+str(perP)+'_map'+mapping+'_'+'_'.join(np.array(parmapping)).replace(':','-')

    if coefgiven:
        params = {'xmeans':xmean, 'ymeans':ymean}
        pr = jple(decoder = decoding, coef_ = coef_, xloc = [0,len(sequencefeatures)], yloc = [len(sequencefeatures),len(sequencefeatures)+len(Zfeatures)], **params)
    else:
        pr = jple(latentmap = latent, perP=perP, varrem=reduct, mapalg=mapping, decoder = decoding, **params)
        


elif '--Affinity_regression' in sys.argv:
    reducez = sys.argv[sys.argv.index('--Affinity_regression')+1] # largest sub k-mer in D
    
    print 'making subkmers', reducez
    if reducez.isdigit():
        reducez = int(reducez)
        subZfeaturelist = []
        for zfeatc in Zfeatures:
            subz = []
            for zk in range(1, reducez+1):
                for l in range(len(zfeatc)-zk+1):
                    subz.append(zfeatc[l:l+zk])
            subZfeaturelist.append(np.unique(subz))
        subZfeatures = np.unique(np.concatenate(subZfeaturelist))
        Zfeaturesubmat = np.zeros((len(Zfeatures), len(subZfeatures)))
        for z, zfeatc in enumerate(Zfeatures):
            Zfeaturesubmat[z, np.isin(subZfeatures, subZfeaturelist[z])] = 1.
        # Normalization of D, same as in AR matlab code
        Zfeaturesubmat = Zfeaturesubmat/np.sqrt(np.sum(Zfeaturesubmat**2, axis = 0))
    else:
        Zfeaturesubmat = None
    print 'done'
    
    redX = float(sys.argv[sys.argv.index('--Affinity_regression')+2]) # PCA contained variance of P
    redY = float(sys.argv[sys.argv.index('--Affinity_regression')+3]) # PCA contained variance of YD
    center = check(sys.argv[sys.argv.index('--Affinity_regression')+4]) # center X and Y before fitting
    alpha = float(sys.argv[sys.argv.index('--Affinity_regression')+5]) # regularization
    ridge = check(sys.argv[sys.argv.index('--Affinity_regression')+6]) # use scikitlearn ridge
    kronecker = check(sys.argv[sys.argv.index('--Affinity_regression')+7]) # use kronecker product
    
    Wname += 'AR-P'+str(redX)+'_YD'+str(redY)+'_cent'+str(center)+'_alph'+str(alpha)+'_ridge'+str(ridge)+'_D'+str(reducez)
    pr = affinity_regression(redX= redX, redY = redY, center = center, alpha = alpha, ridge = ridge, kronecker = kronecker, D = Zfeaturesubmat)
    
elif '--LinearRegression' in sys.argv:
    ralpha = float(sys.argv[sys.argv.index('--LinearRegression')+1])
    fit_intercept = check(sys.argv[sys.argv.index('--LinearRegression')+2])
    lasso = check(sys.argv[sys.argv.index('--LinearRegression')+3])
    if lasso:
        lasadd = 'lasso'
    else:
        lasadd = ''
    Wname += '_LR'+lasadd+str(ralpha)+'fi'+str(fit_intercept)
    model = 'LR'+lasadd+str(ralpha)+'fi'+str(fit_intercept)
    pr = regression(fit_intercept = fit_intercept, alpha = ralpha, lasso = lasso)
    

elif '--knearestneighbor' in sys.argv:
    kcut = check(sys.argv[sys.argv.index('--knearestneighbor')+1]) # If True: use value as minimum similarity
    k = float(sys.argv[sys.argv.index('--knearestneighbor')+2]) # either number of cosindered training profiles or lower cutoff for similarity
    distfunc = sys.argv[sys.argv.index('--knearestneighbor')+3] # euclidean, correlation, cosine
    mixfunc = sys.argv[sys.argv.index('--knearestneighbor')+4] # equal, weighted, squared
    
    parmapping = []
    params = {}
    if len(sys.argv) > sys.argv.index("--knearestneighbor")+5:
        if ':' in sys.argv[sys.argv.index("--knearestneighbor")+5]:
            parmapping = sys.argv[sys.argv.index("--knearestneighbor")+5]
            if ',' in parmapping:
                parmapping = parmapping.split(',') 
            else:
                parmapping = [parmapping]
            params = {}
            for pp in parmapping:
                pp = pp.split(':')
                params[pp[0]] = pp[1]
        else:
            params = {}
            parmapping = ['']
        
    Wname += '_KNN'
    if kcut:
        Wname += '_cutoff'
    Wname+=str(k)+'-'+distfunc+'-'+mixfunc
    if len(parmapping) > 0:
        Wname += '_'+'_'.join(np.array(parmapping)).replace(':','-')
    print Wname
    pr = knearest(k = k, cutoff = kcut, mixturefunc = mixfunc, similarityfunc = distfunc, **params)
    
    

#train model
print 'fitmodel'

pr.fit(Ptrain, Ytrain)

if '--knearestneighbor' in sys.argv or '--JPLE' in sys.argv:
    print 'transform'
    Wtest = pr.latent(Ptest)
    nonzerow = np.sum(Wtest!= 0,axis = 1)
    if (nonzerow != len(Wtest[0])).any():
        for t, nt in enumerate(np.sum(Wtest!= 0,axis = 1)):
            print testprots[t], nt
    Wtrain = pr.Wtrain
    print 'predict'
    Ypred = pr.predict_from_latent(Wtest, add_mean = True)
else:
    print 'predict'
    Ypred = pr.predict(Ptest)

if '--savemodel' in sys.argv and '--JPLE' in sys.argv:
    # saving the model parameters only possible for JPLE, interaction matrices for Affinity regression and others are too big
    np.savez_compressed(Wname+'_coef.npz',coef_=pr.coef_, Zfeatures = Zfeatures, sequencefeatures = sequencefeatures, xmean = pr.xmean, ymean = pr.ymean, trainprots = trainprots, testprots = testprots)
    print 'Model saved under', Wname+'_coef.npz'





### Compute and save accuracy of test reconstruction performances

if '--savetestcorrelation' in sys.argv:
    pearsont = np.zeros((len(testprots),2))
    for i in range(len(testprots)):
        pearsont[i,0] = pearsonr(Ypred[i], Ytest[i])[0]
        pearsont[i,1] = spearmanr(Ypred[i], Ytest[i])[0]
        print testprots[i], pearsont[i,0], pearsont[i,1]
    print "Mean", np.mean(pearsont, axis = 0)
    print 'Median', np.median(pearsont, axis = 0)
    np.savetxt(Wname+"-testset_profile_pcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')
    
if '--savetopintersection' in sys.argv:
    pearsont = np.zeros((len(testprots),3))
    perctop = int(len(Ypred[0])*0.01)
    for i in range(len(testprots)):
        pearsont[i,0] = len(np.intersect1d(np.argsort(Ypred[i])[-perctop:], np.argsort(Ytest[i])[-perctop:]))/float(perctop)
        pearsont[i,1] = len(np.intersect1d(np.argsort(Ypred[i])[-100:], np.argsort(Ytest[i])[-100:]))/(100.)
        pearsont[i,2] = len(np.intersect1d(np.argsort(Ypred[i])[-10:], np.argsort(Ytest[i])[-10:]))/(10.)       
        print testprots[i], pearsont[i,0], pearsont[i,1], pearsont[i,2]
    print "Mean", np.mean(pearsont, axis = 0)
    print 'Median', np.median(pearsont, axis = 0)
    np.savetxt(Wname+"-testset_profile_topintersection.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')

if '--savetopcorrelation' in sys.argv:
    pearsont = np.zeros((len(testprots),2))
    perctop = int(len(Ypred[0])*0.01)
    for i in range(len(testprots)):
        weight100 = np.zeros(len(Ypred[i]))
        weight1p = np.zeros(len(Ypred[i]))
        w1 = np.argsort(Ypred[i])[-perctop:]
        w2 = np.argsort(Ytest[i])[-perctop:]
        weight1p[w1] += 0.5
        weight1p[w2] += 0.5
        weight100[w1[-100:]] +=0.5
        weight100[w2[-100:]] +=0.5
        mask1p = weight1p > 0
        mask100 = weight100 > 0
        weight1p = weight1p[mask1p]
        weight100 = weight100[mask100]
        pearsont[i,0] = np.sum((Ypred[i][mask1p]-np.mean(Ypred[i][mask1p]))*weight1p*(Ytest[i][mask1p]-np.mean(Ytest[i][mask1p])))/(np.sum(weight1p)*np.std(Ypred[i][mask1p])*np.std(Ytest[i][mask1p]))
        pearsont[i,1] = np.sum((Ypred[i][mask100]-np.mean(Ypred[i][mask100]))*weight100*(Ytest[i][mask100]-np.mean(Ytest[i][mask100])))/(np.sum(weight100)*np.std(Ypred[i][mask100])*np.std(Ytest[i][mask100]))
    np.savetxt(Wname+"-testset_profile_topcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), pearsont.astype(str), axis=1), fmt = '%s')    





# compute clostest measured profile to testprofile    
if '--saveclosestcorrelation' in sys.argv:
    pearsontesttrain = 1. - cdist(Ypred, Ytrain, 'correlation')
    pearsonclose = np.argmax(pearsontesttrain, axis = 1)
    pearsont = np.array([trainprots[pearsonclose], pearsontesttrain[np.arange(len(pearsonclose)), pearsonclose]]).T
    np.savetxt(Wname+"-directfit_highcorrelation.dat", np.append(testprots.reshape((len(testprots),1)), np.array(pearsont).astype(str), axis=1), fmt = '%s')

# save reconstructed specificities
if "--savetestprofiles" in sys.argv:
        print "Test Z Profiles saved"
        np.savez_compressed(Wname+"-testset_profiles.npz", names=testprots, profiles=Ypred.T, kmers = Zfeatures)

# save a light version of top 100 of reconstructed specificities that use less space
# can be given to 'extract_motifs_zscores.py' to constuct PWM
if "--savetestprofileslight" in sys.argv:
       tk = []
       testprotnames = []
       for t, tprot in enumerate(testprots):
           ksort = np.argsort(-Ypred[t])[:100]
           tk.append(Zfeatures[ksort])
           testprotnames.append(tprot)
           tk.append(Ypred[t][ksort])
           testprotnames.append(tprot+'.Z')
       np.savetxt(Wname+"-testset_profiles_light.txt", np.array(tk).T, fmt= '%s', header=' '.join(np.array(testprotnames)))






# latent representation of test protein
if '--savelatenttest' in sys.argv:
    np.savez_compressed(Wname+"-testset_Platent.npz", names=testprots, profiles=Wtest)
    
# latent representation of training protein
if '--savelatenttrain' in sys.argv:
    np.savez_compressed(Wname+"-traingset_Platent.npz", names=trainprots, profiles=Wtrain)


    


# distances in latent space
if '--savelatentstats' in sys.argv:
    testcosine = cdist(Wtest, Wtrain, metric='cosine') # test to train distances
    testtestcosine = cdist(Wtest, Wtest, metric='cosine') # test to test distances
    tstats = np.amin(testcosine, axis = 1) # minimum distance to trainset
    
    np.savez_compressed(Wname+"-testset_latentdists.npz", testcosine = testcosine, testtestcosine = testtestcosine, testprots = testprots, trainprots = trainprots)

    np.savetxt(Wname+"-testset_latentstats.dat", np.append(testprots.reshape((len(testprots),1)), tstats.reshape((len(tstats),1)).astype(str), axis=1), fmt = '%s')

if '--save_confidence' in sys.argv: # for knn reconstruction
    np.savetxt(Wname+"-confidencescore.dat", np.append(testprots.reshape((len(testprots),1)), np.array(pr.confidence).reshape(-1,1), axis=1), fmt = '%s')
    




if '--save_predicted_similarity' in sys.argv:# like latent stats for affinity regression
    # this is for affinity regression only 
    testcosine = 1.-pr.predict_similarity_totrain(Ptest)
    #traincosine = 1.-pr.predict_similarity_totrain(Ptrain)
    testtestcosine = 1. - np.diag(cdist(Ypred, Ypred, 'correlation'))
    tstats = np.amin(testcosine, axis = 1)
    
    np.savez_compressed(Wname+"-testset_latentdists.npz", testcosine = testcosine, testtestcosine = testtestcosine, testprots = testprots, trainprots = trainprots)

    np.savetxt(Wname+"-testset_latentstats.dat", np.append(testprots.reshape((len(testprots),1)), tstats.reshape((len(tstats),1)).astype(str), axis=1), fmt = '%s')

# measure distance of protein vector to latent space, also like reconstruction error
if '--distance2latentspace' in sys.argv:
    YtestPtoP = pr.predict_from_latent(Wtest, direction = 'backward', add_mean = True)
    cosineback = []
    euclideanback = []
    for t, tprot in enumerate(testprots):
        cosineback.append(cosine(Ptest[t], YtestPtoP[t]))
        euclideanback.append(euclidean(Ptest[t], YtestPtoP[t]))
    np.savetxt(Wname+"-testset_distance2latentspace.txt", np.concatenate([testprots.reshape((len(testprots),1)), np.array(cosineback).reshape(-1,1).astype(str), np.array(euclideanback).reshape(-1,1).astype(str)], axis=1), fmt = '%s')





# reconstructed protein profiles from Ytest
if '--save_reconstructionP' in sys.argv:
    WtestP = pr.latent(Ytest, reverse = True)
    YtestZtoP = pr.predict_from_latent(WtestP, direction = 'backward')
    np.savez_compressed(Wname+"-testset_PreconstructionZ.npz", names=testprots, profiles=YtestZtoP, kmers = sequencefeatures)

# sparse reconstruction for Ytest with only the most significant eigenvectors
if '--savesparse_reconstructionP' in sys.argv:
    WtestP = pr.latent(Ytest, reverse = True)
    wimportance = WtestP**2*np.sum(pr.coef_[2][:,-len(Ytest[0]):]**2, axis = 1)
    YtestP =[]
    # this might be right thing to do for reconstruction option as well
    for p, wtest in enumerate(WtestP):
        wimp = wimportance[p]
        yrec = np.zeros(len(Ytest[p]))
        impset = []
        for wi in np.argsort(-wimp):
            impset.append(wi)
            yrec += WtestP[p][wi] * pr.coef_[2][wi,-len(Ytest[0]):]
            if pearsonr(yrec, Ytest[p])[0] >= 0.85:
                break
        YtestP.append(np.dot(WtestP[p,impset], pr.coef_[2][impset,:-len(Ytest[0])]))
    np.savez_compressed(Wname+"-testset_PreconstructionSparse.npz", names=testprots, profiles=YtestP, kmers = sequencefeatures)
   
# reconstructed protein profiles from Ptest
if '--save_reconstructionPtoP' in sys.argv:
    YtestPtoP = pr.predict_from_latent(Wtest, direction = 'backward')
    np.savez_compressed(Wname+"-testset_PreconstructionP.npz", names=testprots, profiles= YtestPtoP, kmers = sequencefeatures)
# reconstructed protein profiles for training proteins
if '--save_reconstructionPZtoP' in sys.argv:
    YtrainPZtoP = pr.predict_from_latent(Wtrain, direction = 'backward')
    np.savez_compressed(Wname+"-testset_PreconstructionPZ.npz", names=trainprots, profiles = YtrainPZtoP, kmers = sequencefeatures)
    

if '--determine_sequence' in sys.argv:
    dummyproteinfeatures = sys.argv[sys.argv.index('--determine_sequence')+1]
    if os.path.isfile(dummyproteinfeatures):
        Yd, Ydnames, Zd, Pdummy, pdumnames, sequencefeaturesdummy= readin(dummyproteinfeatures, 'none', 'none')
    else:
        print dummyproteinfeatures, 'not a file'
        sys.exit()
    if not np.array_equal(sequencefeatures, sequencefeaturesdummy):
        Pdum = np.zeros((len(Pdummy), len(sequencefeatures)))
        if np.array_equal(sequencefeatures[np.isin(sequencefeatures, sequencefeaturesdummy)], sequencefeaturesdummy[np.isin(sequencefeaturesdummy, sequencefeatures)]):
            Pdum[:, np.isin(sequencefeatures, sequencefeaturesdummy)] = Pdummy[:,np.isin(sequencefeaturesdummy,sequencefeatures)]
            Pdummy = Pdum
        else:
            print 'Dummy protein feautures not equal to model features'
            sys.exit()
    
    if '--JPLE' in sys.argv:
        WtestY = pr.latent(Ytest, reverse = True)
        Wdummy = pr.latent(Pdummy)
        Wtrain = pr.latent(Ptrain)
    elif '--LinearRegression' in sys.argv:
        Wtest = Ptest
        WtestY = pr.predict_backward(Ytest)
        Wdummy = Pdummy
        Wtrain = Ptrain
    rdist = cdist(Wtest, WtestY, 'cosine')
    ddist = cdist(np.append(Wtrain,Wdummy, axis = 0), WtestY, 'cosine')
    rank = np.argsort(np.argsort(np.append(rdist, ddist, axis = 0),axis = 0), axis = 0)
    rank = rank[[i for i in range(len(WtestY))],[i for i in range(len(WtestY))]]
    np.savetxt(Wname+os.path.splitext(os.path.split(dummyproteinfeatures)[1])[0]+'_rank.txt', np.concatenate([testprots.reshape(-1,1), rank.reshape(-1,1), np.diag(rdist).reshape(-1,1)], axis = 1), header = str(len(rdist)+len(ddist)), fmt = '%s')                                                                                                                




if '--predict_backward' in sys.argv: ## used by affinity regression and regression to predict important peptides from output to input
    YtestZtoP = pr.predict_backward(Ytest)
    np.savez_compressed(Wname+"-testset_PreconstructionZ.npz", names=testprots, profiles=YtestZtoP, kmers = sequencefeatures)
    YtrainPZtoP = pr.predict_backward(Ytrain)
    np.savez_compressed(Wname+"-testset_PreconstructionPZ.npz", names=trainprots, profiles = YtrainPZtoP, kmers = sequencefeatures)
    
    


    










    
    
    
    

