# -*- coding: utf-8 -*-
"""
@author: 
"""
from myPDB import *
import cPickle
import numpy as np
from more_itertools import sort_together

def scale_between0and1(X):
    X = ( X - X.min())/(X.max() - X.min() )
    return X

def fill_list(my_list, length, fill=None):
    if len(my_list) >= length:
        return my_list
    else:
        return my_list + (length - len(my_list)) * [fill]

def mergeFeaturesfromReceptorandLigand(L_name, L_file, R_file):
    featDict = {}
    featDict["complex_code"] = L_name 
    
    L=myPDB.loader(L_file)
    R=myPDB.loader(R_file)
    N=len(L.R)
    temp_S = L.S[0]
    for i in range(N):
        temp_S[i] = sort_together( [L.S[1][i], L.S[0][i]] ,reverse=True)[1]
        temp_S[i] = list(temp_S[i])
        temp_S[i].pop(0)
        temp_S[i] = fill_list(temp_S[i], 30, N-1)
    featDict["l_hood_indices"] = np.zeros((N,30,1),dtype=np.int64)
    for i in range(N):
        for j in range(30):
            featDict["l_hood_indices"][i][j] = temp_S[i][j]

    featDict["l_edge"] = np.zeros((N,25,2),dtype=np.float64)
    featDict["l_vertex"] = np.zeros((N,64),dtype=np.float64)
    featDict["l_vertex"][:,0:21] =  L.oneHotFeats   # represent one hot encoding of the residue. 
    featDict["l_vertex"][:,21:41] =  L.pssm.T  # represent PSSM entry
    temp_ASA = scale_between0and1(L.ASA)
    featDict["l_vertex"][:,41] = temp_ASA.T # represent solvent accessibility
    featDict["l_vertex"][:,42] =  L.rASA.T # represent solvent accessibility
    
    temp_HSAAC = (L.UC + L.DC).T   #  Half Sphere Amino Acid Composition
    for i in range(N):
        hsaas = temp_HSAAC[i] 
        hsaas[hsaas==0] = 999
        hsaas[hsaas==999] = 1.0e-05 * hsaas.min()
        temp_HSAAC[i] = hsaas
        
    featDict["l_vertex"][:,43:64] =  temp_HSAAC # represent neighborhood composition
    featDict["l_vertex"] = np.delete(featDict["l_vertex"],63,axis=1)

    featDict["label"] = np.zeros((N,2),dtype=np.int64)
    featDict["label"][:,0] = range(N)
    featDict["label"][:,1] = 1

    N_r=len(R.R)
    featDict["label_r"] = np.zeros((N_r,2),dtype=np.int64)
    featDict["label_r"][:,0] = range(N_r)
    featDict["label_r"][:,1] = 1
    featDict["r_edge"] = np.zeros((N_r,25,2),dtype=np.float64)
    
    featDict["r_hood_indices"] = np.zeros((N_r,30,1),dtype=np.int64)
    temp_S_r = R.S[0]
    for i in range(N_r):
        temp_S_r[i] = sort_together( [R.S[1][i], R.S[0][i]] ,reverse=True)[1]
        temp_S_r[i] = list(temp_S_r[i])
        temp_S_r[i].pop(0)
        temp_S_r[i] = fill_list(temp_S_r[i], 30, N_r-1)
    featDict["r_hood_indices"] = np.zeros((N_r,30,1),dtype=np.int64)
    for i in range(N_r):
        for j in range(30):
            featDict["r_hood_indices"][i][j] = temp_S_r[i][j]
    
    featDict["r_vertex"] = np.zeros((N_r,64),dtype=np.float64)
    featDict["r_vertex"][:,0:21] =  R.oneHotFeats   # represent one hot encoding of the residue. 
    featDict["r_vertex"][:,21:41] =  R.pssm.T  # represent PSSM entry
    temp_ASA_r = scale_between0and1(R.ASA)
    featDict["r_vertex"][:,41] =  temp_ASA_r.T # represent solvent accessibility
    featDict["r_vertex"][:,42] =  R.rASA.T # represent solvent accessibility
    
    temp_HSAAC_r = (R.UC + R.DC).T   #  Half Sphere Amino Acid Composition
    for i in range(N_r):
        hsaas_r = temp_HSAAC_r[i] 
        hsaas_r[hsaas_r==0] = 999
        hsaas_r[hsaas_r==999] = 1.0e-05 * hsaas_r.min()
        temp_HSAAC_r[i] = hsaas_r
    
    featDict["r_vertex"][:,43:64] =  temp_HSAAC_r # represent neighborhood composition
    featDict["r_vertex"] = np.delete(featDict["r_vertex"],63,axis=1)

    return featDict


if __name__=="__main__":

    featDictList = []
    pklpath_ag = '..\\..\\example\\pkl\\ag\\'
    fname_ag = pklpath_ag + 'ag.pkl'
    print fname_ag
    pklpath_ab = '..\\..\\example\\pkl\\ab\\'
    fnames_ab=glob.glob(pklpath_ab + '*.pkl')
    for (idx,f) in enumerate(fnames_ab):
        print f
        dirStr, ext = os.path.splitext(f)
        sampleName = dirStr.split("\\")[-1]
        print sampleName
        newDict = mergeFeaturesfromReceptorandLigand(sampleName,f,fname_ag)
        featDictList.append(newDict)
        featDictList.append(newDict)
    
    ofname= "..\\..\\example\\pkl\\input.cpkl"
    output = open(ofname, 'wb')
    cPickle.dump(featDictList, output,-1)
    output.close() 
    
    print "Done."

