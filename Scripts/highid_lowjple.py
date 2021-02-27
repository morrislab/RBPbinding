import numpy as np
import sys, os
import glob
jplenpzs = sys.argv[1]
idnpzfile = sys.argv[2]
species = np.sort(glob.glob('*/'))

obj = open('CisBP_highid_jplemin_species.stats', 'w')
for s, speci in enumerate(species):
    if os.path.isfile(speci+jplenpzs) and os.path.isfile(speci+speci.strip('/')+idnpzfile):
        jpdist = np.load(speci+jplenpzs)
        trainprots = jpdist['trainprots'].astype(str)
        testprots = jpdist['testprots'].astype(str)
        testcosine = jpdist['testcosine']
    
        idnpzs = np.load(speci+speci.strip('/')+idnpzfile)
        identmat = idnpzs['identmat']
        itestprots = idnpzs['names2'].astype(str)
        itrainprots = idnpzs['names'].astype(str)
        for it, itr in enumerate(itestprots):
            itestprots[it] = itr.split('__')[0]
        for it, itr in enumerate(itrainprots):
            itrainprots[it] = itr.split('||')[1]
    
        if not np.array_equal(itestprots,testprots):
            itestsort = np.argsort(itestprots)
            testsort = np.argsort(testprots)
            testcosine = testcosine[testsort]
            identmat = identmat[itestsort]
        if not np.array_equal(itrainprots, trainprots):
            itrainsort = np.argsort(itrainprots)
            trainsort = np.argsort(trainprots)
            testcosine = testcosine[:,trainsort]
            identmat = identmat[:, itrainsort]
                
        print np.shape(identmat), np.shape(testcosine)
        jplemin = np.argmin(testcosine, axis = 0)
        idmax = np.argmax(identmat, axis = 0)
        
        jplemin = zip(jplemin, np.arange(len(jplemin), dtype = int))
        idmax = zip(idmax,np.arange(len(idmax), dtype = int))
 
        jplecut = zip(np.where(testcosine < 0.3)[0], np.where(testcosine < 0.15)[1])
        idcut = zip(np.where(identmat > 30)[0],np.where(identmat > 50)[1])
        if len(jplecut) == 0:
             jplecut = jplemin
        if len(idcut) == 0:
            idcut = idmax
        tuples = np.unique(np.concatenate([jplemin, idmax, jplecut, idcut], axis = 0), axis = 0)
        
        for t in tuples:
            obj.write(testprots[t[0]]+' '+trainprots[t[1]]+' '+str(identmat[t[0], t[1]])+' '+str(testcosine[t[0],t[1]])+'\n')


