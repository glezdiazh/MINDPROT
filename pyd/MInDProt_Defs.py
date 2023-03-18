## Functions for MIProt software = TIs for PDBs
## import modulos python

from numpy import *
from Bio.PDB import *
from Bio import SeqIO
import os # from os.path import join, getsize
import sys, urllib, time, datetime
from os.path import join, getsize
from os.path import dirname
from string import *
import tempfile

############################
#   function definitions   #
############################

## ---------------------------------------------------------------------
def GetSeqsFromFASTA(sInput):
    # take a FAST file as input and return a dictionary of proteins
    dSeqs={}
    for seq_record in SeqIO.parse(sInput, "fasta"):
        sID=str(seq_record.id)
        EntireID=sID.split("|")
        dSeqs[EntireID[1]]=str(seq_record.seq)
    return dSeqs # dictionary of GI numbers and Sequences
## ---------------------------------------------------------------------
def GetFastaCM(Seq):
    # take a sequence as input and return the connectivity maxtrix
    l=len(Seq)
    CM=zeros((l,l)) # create empty CM Seq*Seq dimension
    minL=0    # minimal position in the seq
    maxL=l-1  # last position in the seq
    for i in range(l):
        # seq links    
        # at least the second AA in the Seq
        if i>minL:CM[i][i-1]=1.0
        # at least one AA before the end of the seq
        if i<maxL:CM[i][i+1]=1.0

        # recurrence
        CM[i][i]=1.0
        for j in range(l):
            if Seq[i]==Seq[j]: CM[i][j]=1.0
    return CM
## ------------------------------------------------------------------------
def CheckPDBDown(pathPDB,PDB):
    e=0
    # test if the file exists in the path
    PDBfile= os.path.join(pathPDB,PDB+".pdb")
    if os.path.isfile(PDBfile) == False : # if the file is not inside the path
        # Note the %s where the PDB code should go:
        src = "http://www.pdb.org/pdb/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s"
 
        print "-> Downloading " + PDB+".pdb ...",
        try:
            urllib.urlretrieve(src % PDB, PDBfile)
            if int(getsize(PDBfile))<int(1000): # if is not a good PDB
                e=1
                print "\nERROR: "+ PDB + " is not in the Web database! Please check the name of the PDB file in the list."
                os.remove(PDBfile)
            else: #if the download was fine print the file lenght
                print "%(lung)s bytes" %\
                      {"lung": getsize(PDBfile)}
        except: # if thee server connection has problems
            print "ERROR: The Web server is down! Please try later!"
            e=1
    return e
## ------------------------------------------------------------------------
def MPower(M,np):
    # matrix M to the power np
    nM=M.shape[0]
    if np==0:
        Mfin=identity(nM)
    if np==1:
        Mfin=M.copy()
    if np>1:
        Mfin=M.copy()
        for npi in range(1,np):
            Mfin=dot(Mfin,M)
    return Mfin
# ..............................................................
def printM(M):
    # friendly column print of a matrix M
    strM=''
    l=M.shape[0]
    c=M.shape[1]
    for ll in range(l):
        for cc in range(c):
            if cc==c-1:
                strM+=" "+str(M[ll][cc])+"\n"
            else:
                strM+=" "+str(M[ll][cc])
    return strM
# ----------------------------------------------------------------------------------------------------------------------
## functions for AA info

def AAinfo():
    # take the entire AA info (TAB delimited) from a text file and return a list
    # AAFullName[0],A3[1],A1[2],vdWradius[3],NetCh[4],AmberCh[5],Polar_KJ[6]
    # AtContrib2P[7],AtRefr[8],vdWArea[9],hardness_I-A[10],Electrophilicity[11],ElectroMulliken[12]
    AAp=[]
    AAFile = open("AAconstants.txt","r")
    for line in AAFile.readlines():
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        AAline=line.split("\t")
        AAp.append(AAline)
    AAFile.close()
    return AAp
# ..............................................................
def AAinfoCol(AAname, AAp,col):
    # return the values for one AAA (AAname 3 letter) and one column property (col)
    # AAFullName[0],A3[1],A1[2],vdWradius[3],NetCh[4],AmberCh[5],Polar_KJ[6]
    # AtContrib2P[7],AtRefr[8],vdWArea[9],hardness_I-A[10],Electrophilicity[11],ElectroMulliken[12]
    r=-1 # bad value
    for AA in AAp:
        AA3=AA[1]
        if AAname.find(AA3)!=-1:
            r=float(AA[col])
            return r
    # no found, use the Unknown value
    Unknown=AAp[len(AAp)-1]
    r=float(Unknown[col])
    return r

# ..............................................................
def AAinfoColA1(AAname, AAp,col):
    # return the values for one AAA (AAname 3 letter) and one column property (col)
    r=0
    for AA in AAp:
        AA1=AA[2]
        if AAname==AA1:
            r=float(AA[col-1]) # -1 because col is the real column no.
            return r
    # no found, use the Unknown value
    Unknown=AAp[len(AAp)-1]
    r=float(Unknown[col-1])
    return r
# ..............................................................
def PDBchainAtoms(pathPDB,PDBchain,atom):
    # get the list of a specific type of atom with coordinates and labels for one PDB protein chain
    # {model, chain, AAA, residue position, X, Y, Z}
    
    PDB=PDBchain[:4]   # take the PDB name
    sChain=PDBchain[4] # take the chain

    PDBfile= os.path.join(pathPDB,PDB+".pdb")
    
    # extract the alpha-Carbon list with the coords
    p=PDBParser(PERMISSIVE=1) # permissive read for errors
    s=p.get_structure("ProtChain",PDBfile) # read the PDB structure
    model=s[0] # choose the first model if there are many

    ChainList=[] # chain list to analyse

    # gets the chain list of the PDB
    PDBchainList=[] #chain list from PDB
    # checking if the chain exists in the PDB
    for chain in model.get_iterator():
        ch=str(chain.get_id())
        PDBchainList.append(ch)

    # if no chain is founded use the entire protein
    ccc=0 # flag to check if a chain exists
    if sChain!="*": # if you have a chain in the PDB list 4+1
        # checking if the chain exists in the PDB
        for iCh in PDBchainList:
            if iCh==sChain: # if the chain exists inside the PDB
                ccc=1 # at least one chain exists
        if ccc==1:
            ChainList.append(sChain) # use the chain from the list
        else:
            ChainList.append(PDBchainList[0]) # if the PDB chain is wrong we are taking the first chain
        
    else:          # if there is no chain, use all the chains
        print "--> CHAIN missing! Using the entire proteins!"
        # get the chain list from PDF
        for chain in model.get_iterator():
            ch=str(chain.get_id())              
            ChainList.append(ch) # add to process all the chains

    # print "ChainList=",ChainList,"-"
    # process the protein
    AtomList=[] # alpha C list
    for chain in model.get_iterator(): # for each chain in the PDB
        # verify a list of chains (1 if the chain is in the PDBlist; all chains if the chain is not specified)
        for iChain in ChainList: # compare with each chain to be processed
            # if there is an error in the chain lists, avoid it
            try:
                # if there is the chain that should be processes
                if chain==model[iChain]: # only one chain
                    for residue in chain.get_iterator():
                        if is_aa(residue): # if the residue is an AA
                            if residue.has_id(atom): # if it's an alpha C
                                ca=residue[atom]
                                coords=ca.get_coord()           
                                resseq=residue.get_id()[1]
                                resname=residue.get_resname()
                                model_id=model.get_id()
                                chain_id=chain.get_id()
                                AtomList.append((model_id, chain_id, resname, resseq, float(coords[0]),float(coords[1]),float(coords[2])))
            except:
                pass
            
    return AtomList # atom list with coords and AA info for one chain5
# ..............................................................
def AAcentroid(AtomList):
    # calculate the atom coordinate centroid from the AtomList of PDBchainAtoms(pathPDB,PDBchain,atom)
    # input: atoms list containing coordinates
    # {model, chain, AAA, residue position, X, Y, Z}
    # output: the coordinates of the centroid {Xc,Yc,Zc}

    Xc=float(0) # coordinates of the AA centroid
    Yc=float(0)
    Zc=float(0)
    nCa=len(AtomList) # number of AA in the list
    for i in range(nCa):
        (model, chain, resname, resseq, X, Y, Z)=AtomList[i]
        Xc=Xc+float(X) # sum of all the X of Ca
        Yc=Yc+float(Y) # sum of all the Y of Ca
        Zc=Zc+float(Z) # sum of all the Z of Ca
    # CENTROID coords
    Xc=float(Xc)/float(nCa) # x of Ca centroid
    Yc=float(Yc)/float(nCa) # y of Ca centroid
    Zc=float(Zc)/float(nCa) # z of Ca centroid
    return (Xc,Yc,Zc) # return the coords of the centroid as (Xc,Yc,Zc)
# ..............................................................
def CentroidDistVect(AtomList,(Xc,Yc,Zc)):
    # calculate the vector with the distances to a centroid (Xc,Yc,Zc) of a list of atoms (AtomList)
    # input: list of atoms and the centroid coords
    nCa=len(AtomList)
    D0j=zeros((nCa))  # centriod distance vector
    for j in range(nCa):
        (model1, chain1, resname1, resseq1, x1, y1, z1)=AtomList[j]
        D0j[j]=sqrt(float(Xc-x1)*float(Xc-x1)+float(Yc-y1)*float(Yc-y1)+float(Zc-z1)*float(Zc-z1))
    return D0j # return the vector with the distances to a centroid (Xc,Yc,Zc) of a list of atoms (AtomList)
# ..............................................................
def DistMat(AtomList):
    # calculate the (alpha-Carbon or other) distance matrix for a PDB(chain) = Dij
    # input: atom/AA list (AtomList)
    nCa=len(AtomList)
    Dij=zeros((nCa,nCa)) # Ca-Ca distance matrix initialization
    AAp=AAinfo() # list with 22 items/AA info
    for i in range(nCa):
        (model1, chain1, resname1, resseq1, x1, y1, z1)=AtomList[i]
        r1=float(AAinfoCol(resname1,AAp,3)) # vdW radius of the 1st AA
        for j in range(i,nCa):
            (model2, chain2, resname2, resseq2, x2, y2, z2)=AtomList[j]
            dCa1Ca2=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1) # distance between the Ca1 and Ca2
            if dCa1Ca2>0.0:  # different Ca atoms
                dCa1Ca2=sqrt(dCa1Ca2)
            else:            # the same Ca => i=j
                dCa1Ca2=r1   # vdW radius of the 1st AA
            Dij[i][j]=dCa1Ca2 # elements of the distance matrix
            Dij[j][i]=dCa1Ca2
    return Dij # return the distance matrix of an atom list
# ..............................................................
def ConnectMat(AtomList,Dij,iCut):
    # calculate the connectivity matrix for a list of atoms
    # input: atom list (AtomList) and atom distance matrix (Dij)
    # dependence: PDBchainAtoms(pathPDB,PDBchain,atom), DistMat(AtomList), AAinfo(), AAradius()
    nCa=len(AtomList)
    CM=zeros((nCa,nCa))   # Ca-Ca connectivity matrix initialization ("Float" to avoid divide problems)
    AAp=AAinfo() # list with 22 items/AA info
    for i in range(nCa):
        (model1, chain1, resname1, resseq1, x1, y1, z1)=AtomList[i]
        r1=float(AAinfoCol(resname1,AAp,3)) # vdW radius of the 1st AA
        for j in range(i,nCa):
            (model2, chain2, resname2, resseq2, x2, y2, z2)=AtomList[j]
            r2=float(AAinfoCol(resname2,AAp,3)) # vdW radius of the 2nd AA
            Fact=float(0.0)
            # the type of the cutoff and the value come from outside in iCut as (type[Int],val/Roff[Float],Ron[float])
            Roff=iCut[1]
            Ron =iCut[2]
            r   =Dij[i][j]
            if i==j: r=r2 # if i=j, the distance is the vdW radius
            if iCut[0]==1:   # Direct cutoff iCut[1]=7 A default
                if r <= iCut[1]:
                    Fact=1.0
                else:
                    Fact=0.0
            elif iCut[0]==2: # Abrupt truncation iCut[1]=0.5 default
                if r <= iCut[1]*(r1+r2):
                    Fact=1.0
                else:
                    Fact=0.0
            elif iCut[0]==3: # Shifting function
                if r <= Roff:
                    Fact=(1-(float(r)/float(Roff))**2)**2
                else:
                    Fact=0.0
            elif iCut[0]==4: # Force Shifting
                if r <= Roff:
                    Fact=(1-float(r)/float(Roff))**2
                else:
                    Fact=0.0
            elif iCut[0]==5: # Switching function
                if r <= Ron:
                    Fact=1.0
                else:
                    if (r > Ron) and (r < Roff):
                        Fact=(((Roff*Roff-r*r)**2)*(Roff*Roff+2*r+r-3*Ron*Ron))/((Roff*Roff-Ron-Ron)**3)
                    else:
                        Fact=0.0
            CM[i][j]=Fact
            CM[j][i]=Fact  # the matrix is symetric
    return CM # CM= connectivity matrix
# ..............................................................
def WMatFasta(Seq,CM,colW):
    # calculate the weighted matrix for the AA of a Fasta seq
    # input: Seq,CM, colW = column with the AA property
    nAt=len(Seq)
    W=zeros((nAt,nAt))  # electro matrix initialization
    AAp=AAinfo() # list with AA info
    for i in range(nAt): # all i and j values because is wj!
        for j in range(nAt):W[i][j]=CM[i][j]*float(AAinfoColA1(Seq[j],AAp,colW))
    return W # return the weighted matrix
# ..............................................................
def WMat(AtomList,CM,Dij,iW):
    # calculate the weighted matrix for the AA PDB chain
    # input: list of atoms/AAs (AtomList),CM,Dij,iW = type of weight, col= column with the AA property
    nCa=len(AtomList)
    W=zeros((nCa,nCa))  # electro matrix initialization
    AAp=AAinfo() # list with 22 items/AA info
    for i in range(nCa):
        for j in range(nCa): # all i and j values because is wj!
            (model2, chain2, resname2, resseq2, x2, y2, z2)=AtomList[j]
            W[i][j]=CM[i][j]*float(AAinfoCol(resname2,AAp,iW))
    return W # return the weighted matrix
# ..............................................................
def ProbMatrix(W,CM):
    # calculate the interaction probability matrix for W (weigthed matrix) using CM (contact matrix)
    # input: W (or other matrix such as vdW), CM
    nCa=CM.shape[0]
    PiW=zeros((nCa,nCa)) # moment matrix
    for i in range(nCa):
        sumWij=float(0.0)
        for j in range(nCa):sumWij+=W[i][j]*CM[i][j]
        for j in range(nCa):
            if (sumWij!=0) and (CM[i][j]!=0):PiW[i][j]=float(W[i][j])/float(sumWij)
    return PiW
# ..............................................................
def PropMeansFasta(Seq,AAp,colW,PI2n,OrbitLim):
    # OrbitLim=[float(0.0),float(120.0),float(240),float(360),float(len(Seq))] # list with the limits of the orbitals
    # calculate mean properties
    # colW = column for prop of AA!!
    LenSeq=len(Seq)
    (PI0,PI1,PI2,PI3,PI4,PI5)=PI2n   # take the mom matrices to power
    P0j=ones((LenSeq))/float(LenSeq) # P0j = absolut initial probabilities aprox. 1/n
    vPot=zeros((LenSeq))             # potential vector 
    for j in range(LenSeq):vPot[j]=float(AAinfoColA1(Seq[j],AAp,colW))
    # first multiplication from the end of the formula
    Temp0=dot(PI0,vPot)
    Temp1=dot(PI1,vPot)
    Temp2=dot(PI2,vPot)
    Temp3=dot(PI3,vPot)
    Temp4=dot(PI4,vPot)
    Temp5=dot(PI5,vPot)
    
    AEN0O0=float(0)
    AEN1O0=float(0)
    AEN2O0=float(0)
    AEN3O0=float(0)
    AEN4O0=float(0)
    AEN5O0=float(0)
    AEN0O1=float(0)
    AEN1O1=float(0)
    AEN2O1=float(0)
    AEN3O1=float(0)
    AEN4O1=float(0)
    AEN5O1=float(0)
    AEN0O2=float(0)
    AEN1O2=float(0)
    AEN2O2=float(0)
    AEN3O2=float(0)
    AEN4O2=float(0)
    AEN5O2=float(0)
    AEN0O3=float(0)
    AEN1O3=float(0)
    AEN2O3=float(0)
    AEN3O3=float(0)
    AEN4O3=float(0)
    AEN5O3=float(0)
    AEN0O4=float(0)
    AEN1O4=float(0)
    AEN2O4=float(0)
    AEN3O4=float(0)
    AEN4O4=float(0)
    AEN5O4=float(0)
    
    # FASTA pseudo ORBITALS
    # ex. for window 120:
    # O0: 0 <= c <= 120 | O1: 121 < i <= 240 | O2: 241 < m <= 360 | O3: 361 < s <= end | O4: all = t
    for j in range(LenSeq):
        # calculation of each power terms
        t0=P0j[j]*Temp0[j]
        t1=P0j[j]*Temp1[j]
        t2=P0j[j]*Temp2[j]
        t3=P0j[j]*Temp3[j]
        t4=P0j[j]*Temp4[j]
        t5=P0j[j]*Temp5[j]
        # orbit total (O4)
        AEN0O4+=t0
        AEN1O4+=t1
        AEN2O4+=t2
        AEN3O4+=t3
        AEN4O4+=t4
        AEN5O4+=t5
        if j>=OrbitLim[0] and j<=OrbitLim[1]: # only orbit c (0) for all 0-5 powers
            AEN0O0+=t0
            AEN1O0+=t1
            AEN2O0+=t2
            AEN3O0+=t3
            AEN4O0+=t4
            AEN5O0+=t5
        else:
            if j>OrbitLim[1] and j<=OrbitLim[2]: # only orbit i (1) for all 0-5 powers
                AEN0O1+=t0
                AEN1O1+=t1
                AEN2O1+=t2
                AEN3O1+=t3
                AEN4O1+=t4
                AEN5O1+=t5
            else:
                if j>OrbitLim[2] and j<=OrbitLim[3]: # only orbit m (2) for all 0-5 powers
                    AEN0O2+=t0
                    AEN1O2+=t1
                    AEN2O2+=t2
                    AEN3O2+=t3
                    AEN4O2+=t4
                    AEN5O2+=t5
                else:
                    if j>OrbitLim[3] and j<=OrbitLim[4]: # only orbit m (3) for all 0-5 powers
                        AEN0O3+=t0
                        AEN1O3+=t1
                        AEN2O3+=t2
                        AEN3O3+=t3
                        AEN4O3+=t4
                        AEN5O3+=t5     
##    return [AEN0O0,AEN1O0,AEN2O0,AEN3O0,AEN4O0,AEN5O0,AEN0O1,AEN1O1,AEN2O1,AEN3O1,AEN4O1,AEN5O1,AEN0O2,AEN1O2,AEN2O2,AEN3O2,AEN4O2,AEN5O2,AEN0O3,AEN1O3,AEN2O3,AEN3O3,AEN4O3,AEN5O3,AEN0O4,AEN1O4,AEN2O4,AEN3O4,AEN4O4,AEN5O4]
##    # return the average electrostatic numer
    return (float(AEN0O0+AEN1O0+AEN2O0+AEN3O0+AEN4O0+AEN5O0)/float(6),float(AEN0O1+AEN1O1+AEN2O1+AEN3O1+AEN4O1+AEN5O1)/float(6),float(AEN0O2+AEN1O2+AEN2O2+AEN3O2+AEN4O2+AEN5O2)/float(6),float(AEN0O3+AEN1O3+AEN2O3+AEN3O3+AEN4O3+AEN5O3)/float(6),float(AEN0O4+AEN1O4+AEN2O4+AEN3O4+AEN4O4+AEN5O4)/float(6))   # return the average by orbital for all ks
# ..............................................................
def PropMeans(AtomList,AAp,iW,PI2n,D0j,OrbitLim):
    # calculate the average electrostatic numbers = mean properties
    # input: list with the atom/AA list, mom matrices at power PI2n
    # iW= column for property
    nCa=len(AtomList)
    (PI0,PI1,PI2,PI3,PI4,PI5)=PI2n # take the mom matrices to power
    P0j=ones((nCa))/float(nCa) # P0j = absolut initial probabilities aprox. 1/n
    vPot=zeros((nCa))          # potential vector 
    for j in range(nCa):
        (model2, chain2, resname2, resseq2, x2, y2, z2)=AtomList[j]
        q2=AAinfoCol(resname2,AAp,iW) # any property of the second AA
        #print "AA property:  ",resname2,shift+iW, q2
        if D0j[j]!=0:vPot[j]=float(q2)
    # first multiplication from the end of the formula
    Temp0=dot(PI0,vPot)
    Temp1=dot(PI1,vPot)
    Temp2=dot(PI2,vPot)
    Temp3=dot(PI3,vPot)
    Temp4=dot(PI4,vPot)
    Temp5=dot(PI5,vPot)
    
    AEN0O0=float(0)
    AEN1O0=float(0)
    AEN2O0=float(0)
    AEN3O0=float(0)
    AEN4O0=float(0)
    AEN5O0=float(0)
    AEN0O1=float(0)
    AEN1O1=float(0)
    AEN2O1=float(0)
    AEN3O1=float(0)
    AEN4O1=float(0)
    AEN5O1=float(0)
    AEN0O2=float(0)
    AEN1O2=float(0)
    AEN2O2=float(0)
    AEN3O2=float(0)
    AEN4O2=float(0)
    AEN5O2=float(0)
    AEN0O3=float(0)
    AEN1O3=float(0)
    AEN2O3=float(0)
    AEN3O3=float(0)
    AEN4O3=float(0)
    AEN5O3=float(0)
    AEN0O4=float(0)
    AEN1O4=float(0)
    AEN2O4=float(0)
    AEN3O4=float(0)
    AEN4O4=float(0)
    AEN5O4=float(0)
    
    maxDist=max(D0j) # maxim distancia Ca to centroid = origin
    #orbits
    # O0: 0 <= c <= 25 | O1: 25 < i <= 50 | O2: 50 < m <= 75 | O3: 75 < s <= 100 | O4: all = t
    for j in range(nCa):
        # calculation of each power terms
        t0=P0j[j]*Temp0[j]
        t1=P0j[j]*Temp1[j]
        t2=P0j[j]*Temp2[j]
        t3=P0j[j]*Temp3[j]
        t4=P0j[j]*Temp4[j]
        t5=P0j[j]*Temp5[j]
        # orbit total (O4)
        AEN0O4+=t0
        AEN1O4+=t1
        AEN2O4+=t2
        AEN3O4+=t3
        AEN4O4+=t4
        AEN5O4+=t5
        d_rel=float(D0j[j])*float(100)/float(maxDist)
        if d_rel>=OrbitLim[0] and d_rel<=OrbitLim[1]: # only orbit c (0) for all 0-5 powers
            AEN0O0+=t0
            AEN1O0+=t1
            AEN2O0+=t2
            AEN3O0+=t3
            AEN4O0+=t4
            AEN5O0+=t5
        else:
            if d_rel>OrbitLim[1] and d_rel<=OrbitLim[2]: # only orbit i (1) for all 0-5 powers
                AEN0O1+=t0
                AEN1O1+=t1
                AEN2O1+=t2
                AEN3O1+=t3
                AEN4O1+=t4
                AEN5O1+=t5
            else:
                if d_rel>OrbitLim[2] and d_rel<=OrbitLim[3]: # only orbit m (2) for all 0-5 powers
                    AEN0O2+=t0
                    AEN1O2+=t1
                    AEN2O2+=t2
                    AEN3O2+=t3
                    AEN4O2+=t4
                    AEN5O2+=t5
                else:
                    if d_rel>OrbitLim[3] and d_rel<=OrbitLim[4]: # only orbit m (3) for all 0-5 powers
                        AEN0O3+=t0
                        AEN1O3+=t1
                        AEN2O3+=t2
                        AEN3O3+=t3
                        AEN4O3+=t4
                        AEN5O3+=t5
                        
    ## return [AEN0O0,AEN1O0,AEN2O0,AEN3O0,AEN4O0,AEN5O0,AEN0O1,AEN1O1,AEN2O1,AEN3O1,AEN4O1,AEN5O1,AEN0O2,AEN1O2,AEN2O2,AEN3O2,AEN4O2,AEN5O2,AEN0O3,AEN1O3,AEN2O3,AEN3O3,AEN4O3,AEN5O3,AEN0O4,AEN1O4,AEN2O4,AEN3O4,AEN4O4,AEN5O4]
    # return the average electrostatic numer
    return (float(AEN0O0+AEN1O0+AEN2O0+AEN3O0+AEN4O0+AEN5O0)/float(6),float(AEN0O1+AEN1O1+AEN2O1+AEN3O1+AEN4O1+AEN5O1)/float(6),float(AEN0O2+AEN1O2+AEN2O2+AEN3O2+AEN4O2+AEN5O2)/float(6),float(AEN0O3+AEN1O3+AEN2O3+AEN3O3+AEN4O3+AEN5O3)/float(6),float(AEN0O4+AEN1O4+AEN2O4+AEN3O4+AEN4O4+AEN5O4)/float(6))   # return the average by orbital for all ks
    

### testing
##
##AAp=AAinfo()
##for i in range(3,13):
##    print AAinfoCol("GLX",AAp,i)
