# -*- coding: utf-8 -*-
"""
# Original code by Afsar and Fayyaz
# Modified by Zhengshan Chen
# Changes made: Added new function __getOneHotFeats and global variable PSIBLASTBASE
#               modified function __getPSSM
#               deleted function __getPSAIA, __calcRD, showSurface and related dependent code
# Original code can be found at: https://combi.cs.colostate.edu/supplements/pairpred/PAIRPred.zip
# This code is licensed under GNU Lesser General Public License.
# Last modified on 2024-09-01
"""
from BISEPutils import *
from scipy.spatial.distance import *
from scipy import linalg as LA
from psiblast import *
import stride
import numpy
import os
import glob

PSIBLASTBASE="/amax/meta/db/NCBI_nr/nr" 
class myPDB:
    def __init__(self,fname,auxdir=None):
        self.HWS_seqFV=10.0
        self.HWS_PSSM=5.0
        
        (dpath,name,ftype)=getFileParts(fname)
        if auxdir is None:
            auxdir=dpath
        self.ftype=ftype.lower()
        self.name=name
        self.ifname=fname
        self.resi=dict()
        if self.ftype == '.pdb':
            (self.stx,self.R,self.pp,self.seq,self.S2Ri)=readPDB(fname)        
            self.FV=getSeqFV(self.stx,self.R,self.HWS_seqFV)
            self.Coords=getCoords(self.R)
            self.__getStride()
            self.getSimMtx()
            self.__getBvalues()
            self.__getOneHotFeats()            
            for (idx,r) in enumerate(self.R):               
                self.resi[getResiId(r.get_full_id())]=idx
            self.__assignAux(auxdir)
        elif self.ftype == '.fasta':
            self.seq=readFASTA(fname)
            self.R=self.seq
            self.S2Ri=range(len(self.R))
            self.__getPSSM(auxdir)
            for (idx,r) in enumerate(self.R):               
                self.resi[('A',str(idx))]=idx
        else:
            estring = 'Error: Unknown File Type. Input file must be PDB or FASTA.'
            print estring
            raise Exception(estring)
             
        
    def __assignAux(self,auxpath):
        try:
            self.__getPSSM(auxpath)
        except Exception as e:
            print "Error getting profile features for",self.name,":",e," Is the .mat file present?"
    def save2FASTA(self,fpath,saveHeader=False):        
        f = open(fpath, 'w+')
        if saveHeader:
            f.write('>'+self.name+'\n')
        f.write(self.seq)
        f.close()
    
    def __getOneHotFeats(self):  
        """
        Constructs the 21 one-hot representation
        """
        AA_dic = "ACDEFGHIKLMNPQRSTVWY*"
        letter_to_id = dict((let, id) for (id, let)
                                      in enumerate(AA_dic))
        N = len(self.seq)
        self.oneHotFeats = numpy.zeros((N,21), dtype=float)
        for (id, letter) in enumerate(self.seq):
            self.oneHotFeats[id,letter_to_id.get(letter)] = 1
            
    def __getPSSM(self,auxpath):
        #       Assigning PSSM Profile and using spinex
        
        ppath=os.path.join(auxpath,self.name+'.mat')        
        if not (os.path.exists(ppath) and os.path.isfile(ppath)): 
            print 'Extracting PSI BLAST profile'
            if self.ftype == '.pdb':
                fpath=os.path.join(auxpath,self.name+'.fasta')
                self.save2FASTA(fpath)
            elif self.ftype == '.fasta':
                fpath=self.ifname
            runPSIBLAST(fpath,db=PSIBLASTBASE,ofile=auxpath+self.name,niter=3)
        else:
            print 'Using existing PSI-BLAST Profile',ppath
        pfile=parsePSSMfile(ppath)#parsed file
        N=len(self.R)
        if pfile is not None:
            (pssm,psfm,info)=pfile
            self.pssm=np.zeros((20,N))
            self.psfm=np.zeros((20,N))
            self.info=np.zeros(N);self.info.fill(np.nan)
            for k in range(len(self.seq)):
                idx=self.S2Ri[k]
                self.pssm[:,idx]=pssm[:,k]
                self.psfm[:,idx]=psfm[:,k]
                self.info[idx]=info[k]
        else:
            print 'Encountered some error processing spinex or psi-blast',self.name            
            raise Exception(('Encountered some error processing spinex or psi-blast',self.name))
            
    def getSpectralVector(self):
        C=[]
        for r in self.R:
            C.append(r['CA'].get_coord())
        C=np.array(C)
        D=squareform(pdist(C))
        S=np.exp(-0.05*D)
        e_val,e_vec = LA.eig(S)
        v=(e_vec[:,0])
        dv=np.abs(v-v[:,np.newaxis])
        return (S,dv)
    def save(self,ofname=None,bdir='./'):
        """
        Save the object
        """
        if ofname is None:
            ofname=bdir+self.name+'.pdb.pkl'
        output = open(ofname, 'wb')
        cPickle.dump(self, output,-1)
        output.close()
    
    @classmethod   
    def loader(self,pklfile):
        """
        Load the class object from a pickel file
        """
        return cPickle.load(open(pklfile, "rb" ) )
                      
    def getSimMtx(self,sg=2.0,thr=1e-3):
        self.sg=sg
        self.thr=thr
        sg=2*(sg**2)
        I=[[] for k in range(len(self.Coords))]
        V=[[] for k in range(len(self.Coords))]
        for i in range(len(self.Coords)):
            for j in range(i,len(self.Coords)):
                d=spatial.distance.cdist(self.Coords[i], self.Coords[j]).min()
                s=np.exp(-(d**2)/sg)
                if(s>thr):# and not np.isnan(self.Phi[i]) and not np.isnan(self.Phi[j])
                    I[i].append(j)
                    V[i].append(s)
                    if i!=j:
                        I[j].append(i)
                        V[j].append(s)        
        self.S=(I,V)
        self.CN=np.array([len(a) for a in self.S[0]])
        self.__getHSE()        
        
    def __getHSE(self):
        """
        Compute the Half sphere exposure statistics
        The up direction is defined as the direction of the side chain and is 
        calculated by taking average of the unit vectors to different side chain
        atoms from the C-alpha atom
        Anything within the up half sphere is counted as up and the rest as 
        down    
        """        
        N=len(self.R)
        Na=len(AA)
        self.UN=np.zeros(N)
        self.DN=np.zeros(N)
        self.UC=np.zeros((Na,N))
        self.DC=np.zeros((Na,N))
        for (i,r) in enumerate(self.R):
            u=getSideChainV(r)
            if u is None:
                self.UN[i]=np.nan
                self.DN[i]=np.nan
                self.UC[:,i]=np.nan
                self.DC[:,i]=np.nan
            else:                
                idx=aaidx[getResLetter(r)]
                self.UC[idx,i]=self.UC[idx,i]+1
                self.DC[idx,i]=self.DC[idx,i]+1
                n=self.S[0][i]
                for j in n:
                    r2=self.R[j]
                    if is_aa(r2) and r2.has_id('CA'):
                        v2=r2['CA'].get_vector()
                        scode=getResLetter(r2)
                        idx=aaidx[scode]
                        if u[1].angle((v2-u[0])) < np.pi/2.0:
                            self.UN[i]=self.UN[i]+1
                            self.UC[idx,i]=self.UC[idx,i]+1
                        else:
                            self.DN[i]=self.DN[i]+1
                            self.DC[idx,i]=self.DC[idx,i]+1
        self.UC=self.UC/(1.0+self.UN)
        self.DC=self.DC/(1.0+self.DN)
        
    def __getDSSP(self):
        """
        Use DSSP to compute SS,ASA,rASA,Phi,Psi
        """
        ssk="HBEGITS-" #list of possible secondary structures
        sskey=dict(zip(ssk,range(len(ssk))))
        dssp=DSSP(self.stx[0],self.ifname)
        N=len(self.R)
        self.SS=np.zeros((len(ssk),N))
        self.ASA=np.zeros(N)
        self.rASA=np.zeros(N)
        self.Phi=np.zeros(N)
        self.Psi=np.zeros(N)
        for (idx,r) in enumerate(self.R):
            (_,_,cid,rid)=r.get_full_id()
            key=(cid,rid)
            if dssp.has_key(key):
                (_,s,self.ASA[idx],self.rASA[idx],self.Phi[idx],self.Psi[idx])=dssp[key]
                self.SS[sskey[s],idx]=1
            else:
                self.SS[:,idx]=np.nan
                self.ASA[idx]=np.nan
                self.rASA[idx]=np.nan
                self.Phi[idx]=np.nan
                self.Psi[idx]=np.nan
            
    def __getStride(self):        
        #ssk="HBEGITS-" #list of possible secondary structures
        ssk="HGIEBTC"
        sskey=dict(zip(ssk,range(len(ssk))))
        stridex=stride.stride_dict_from_pdb_file(self.ifname)
        N=len(self.R)
        self.SS=np.zeros((len(ssk),N))
        self.ASA=np.zeros(N)
        self.rASA=np.zeros(N)
        self.Phi=np.zeros(N)
        self.Psi=np.zeros(N)
        for (idx,r) in enumerate(self.R):
            key=getResiId(r.get_full_id())
            #pdb.set_trace()
            if stridex.has_key(key):
                # (aa,ss,phi,psi,asa,rasa)
                (_,s,self.Phi[idx],self.Psi[idx],self.ASA[idx],self.rASA[idx])=stridex[key]
                if not sskey.has_key(s):
                    print "unknown scondary structure! Add to dictionary!"
                    pdb.set_trace()
                self.SS[sskey[s],idx]=1
            else:
                print('key not found in stride processing!')
                #pdb.set_trace()
                self.SS[:,idx]=np.nan
                self.ASA[idx]=np.nan
                self.rASA[idx]=np.nan
                self.Phi[idx]=np.nan
                self.Psi[idx]=np.nan  
                
    def __getBvalues(self):
        """
        Save the b value for each residue        
        """
        self.B=np.zeros(len(self.R))
        for (idx,r) in enumerate(self.R):
            #pdb.set_trace()
            self.B[idx]=max([ak.get_bfactor() for ak in r])
    
    def plotRamachandran(self):        
        """
        Plot the ramachandran plot (Psi vs. Phi)
        """
        plt.plot(self.Phi,self.Psi,'.')
        plt.xlabel("Phi")
        plt.ylabel("Psi")
        plt.show()

        
    def toFile(self,ofname=None):
        """
        saves the structure stx to the pdb file path (ofname)
        """
        if ofname is None:
            ofname=self.ifname
        io=PDBIO()
        io.set_structure(self.stx)
        io.save(ofname)

    def __str__( self ):
        s="myPDB object with Original File Name: "+self.name
        return s

if __name__=="__main__":
    usage="python myPDB.py ifile profile_path [ofile.pdb.pkl]"
    ofile=None
    alist=sys.argv
    if len(alist)==4:
        ofile=alist[3]
    if len(alist)<3:
        print 'Error: Improper arguments. Proper Usage is: ', usage
        exit(1)
    ifile=alist[1]        

    ppath=alist[2]    
    L=myPDB(ifile,ppath)
    L.save(ofile)
