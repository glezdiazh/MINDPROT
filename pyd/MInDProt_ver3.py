# calculate the mean properties for proteins and drugs using PDB 3D info and the SMILES
# creates pairs of prot-prot, drug-drug and prot-drug

from Bio.PDB import *
import os, sys, wx, random, time, datetime
from wx.lib.buttons import GenBitmapTextButton
from numpy import *
from math import *
import operator
import datetime

import MInDProt_Defs         # my function file
import TIs2Pairs             # transform the results PDBchain-TIs in PDBChain Pairs-Mixed TIs
import ParseSelecPDBHeader   # parse selective fields of the PDB header
import PDBHeaderAvg          # add averages by classes

import MInDProt_HTMLDoc      # html help supportorig_dir

orig_dir=os.getcwd()+"\\" # get the current path to be used for all functions

def TIs4Fasta(sPDBlistFile,sRezFile,OrbitLim):
    # for FASTA input there is no average calculation or input classes!!!!!
    pathRes="" # if empty -> the current folder

    # AAFullName[0],A3[1],A1[2],vdWradius[3],NetCh[4]
    # iWn: ["ElectroMulliken", "Polar_KJ", "vdWArea", "AtomContrib2P"]
    # iWn=1 # the iWn th column of mean props

    # open Protein Simple result file
    fRez=open(sRezFile,"w")

    sOut="" # print the header
    sOut+=".GenBankID"

    # generate the labels of TIs by orbital 
    IndTypes=["ElectroMulliken", "Polar_KJ", "vdWArea", "AtomContrib2P"]
    for Indice in IndTypes:
        for iOrb in range(5):            # for each orbital
            sOut+="\t"+Indice+"Orb"+str(iOrb) # write TAB delimited [Label]Orb[iOrb]

    # print sOut in Protein simple result
    fRez.write(sOut+"\n")

    # take the list of PDB chains
    sInput=sPDBlistFile
    dSeqs=MInDProt_Defs.GetSeqsFromFASTA(sInput) # dictionary with the Fasta PDBs
    i=0
    iShift=6
    # process each line of PDBchains
    for PDB, Seq in dSeqs.iteritems():
##        try:
        i=i+1
        print "======================================================================="
        print "Processing: "+PDB+" ("+str(i)+" of "+str(len(dSeqs))+")"
        sOut=PDB
        # calculate TIs for each class
        for x in range(len(IndTypes)):
            colW=iShift+x+1
            sOut+="\t"+TIcalcFromCM(PDB,Seq,MInDProt_Defs.GetFastaCM(Seq),colW,OrbitLim)
        fRez.write(sOut+"\n")
##        except:
##            print "--> !!! ERROR for PDB: "+PDB
    fRez.close()
    print "Done!"
    return
## ----------------------------------------------------------------------------------------------
def TIcalcFromCM(PDB,Seq,CM,colW,OrbitLim):
    # CM = contact matrix, Dij = CA distance matrix
    # D0j = CA distance matrix to the centroid, W = weighted matrix

    # calculate connectivity matrix D0j and W of AAs in the sequence
    #         iW=the type of the weigths (int) linked with the column in the AAconsts.txt (1-8)
    # output: TIs
    LenSeq=len(Seq)
    sOut="" # output string
    AAp=MInDProt_Defs.AAinfo() # list with 22 items/AA info

    W=zeros((LenSeq,LenSeq))               # weighted matrix initialization
    W=MInDProt_Defs.WMatFasta(Seq,CM,colW) # sequence, CM, column weights inside AAconstant.txt

    PiW=zeros((LenSeq,LenSeq))             # moment matrix
    PiW=MInDProt_Defs.ProbMatrix(W,CM)  # calculate the interaction probability matrix = moment matrix for W
    
    # power of the moment spectral matrices PI[power]=PI0,PI1,PI2,PI3,PI4,PI5
    PI2n=(MInDProt_Defs.MPower(PiW,0),MInDProt_Defs.MPower(PiW,1),MInDProt_Defs.MPower(PiW,2),MInDProt_Defs.MPower(PiW,3),MInDProt_Defs.MPower(PiW,4),MInDProt_Defs.MPower(PiW,5))
    
    # type of calculations
    AENlist=MInDProt_Defs.PropMeansFasta(Seq,AAp,colW,PI2n,OrbitLim) # any prop means using the column in the Q column in the const file
    for esm in AENlist:
        sOut+="\t"+str(esm)

    return sOut[1:] # return the string with all the information TAB delimited (Prot, Chain, TIs) (exclude the first TAB)
## ----------------------------------------------------------------------------------------------
def TIcalc(pathPDB,PDBchain,iCut,OrbitLim,iW):
    # iCut=1 # the type of the cut chosen from the interface
    # get the descriptors for the calculations
    # CM = contact matrix, Dij = CA distance matrix
    # D0j = CA distance matrix to the centroid, W = weighted matrix

    # calculate connectivity matrix CM, D0j and W of the alpha carbons in one PDB protein chain
    # input : path of the PDB (+ final saparator) and the PDBchain label (PDB=4c + chain=1c)
    #         iCut = (type, val/Roff,Ron) the type of cutoff and the values
    #         OrbitLim=[float(0.0),float(25.0),float(50.0),float(75.0),float(100.0)] the limists of the orbitals
    #         iW=the type of the weigths (int) linked with the column in the AAconsts.txt (1-8)
    # output: connectivity map CM (dimension n*n, n = AA number in the chain)
    #         and the AA label list {model, chain, AAA, residue position}
    
    PDB=PDBchain[:4]     # take the PDB name
    sChain=PDBchain[4:5] # take the chain

    sOut="" # output string
    e=0
    e=MInDProt_Defs.CheckPDBDown(pathPDB,PDB)
    if e==1: # if the is not PDB inside the path or PDB web
        sOut=PDB+" missing!"
        return sOut # return the error message as the outpu string for this PDBchain
    
    AAp=MInDProt_Defs.AAinfo() # list with 22 items/AA info

    AtomList=MInDProt_Defs.PDBchainAtoms(pathPDB,PDBchain,"CA") # get alpha-Carbon (CA) info list from a PDBchain {model, chain, AAA, residue position, X, Y, Z}
    nCa=len(AtomList)

    #if no CA founded
    if nCa==0:
        return "" # nothing

    (Xc,Yc,Zc)=MInDProt_Defs.AAcentroid(AtomList)           # calculate the atom centroid coords
    D0j=zeros((nCa))                              # centriod distance vector D0j for all atoms from AtomList
    D0j=MInDProt_Defs.CentroidDistVect(AtomList,(Xc,Yc,Zc))

    Dij=zeros((nCa,nCa))             # Ca-Ca distance matrix initialization
    Dij=MInDProt_Defs.DistMat(AtomList)        # calculate the distance matrix 

    CM=zeros((nCa,nCa))              # Ca-Ca connectivity matrix initialization
    CM=MInDProt_Defs.ConnectMat(AtomList,Dij,iCut)  # calculate connectivity matrix
    #print "CM\n",CM

    W=zeros((nCa,nCa))               # weighted matrix initialization
    W=MInDProt_Defs.WMat(AtomList,CM,Dij,iW)
    #print "W\n",W
    
    PiW=zeros((nCa,nCa))             # moment matrix
    PiW=MInDProt_Defs.ProbMatrix(W,CM)  # calculate the interaction probability matrix = moment matrix for W
    
    # power of the moment spectral matrices PI[power]=PI0,PI1,PI2,PI3,PI4,PI5
    PI2n=(MInDProt_Defs.MPower(PiW,0),MInDProt_Defs.MPower(PiW,1),MInDProt_Defs.MPower(PiW,2),MInDProt_Defs.MPower(PiW,3),MInDProt_Defs.MPower(PiW,4),MInDProt_Defs.MPower(PiW,5))

    # output the results
    # print the values
    
    # sOut+=PDB+sChain
        
    ##################################
    # only Means

    AENlist=MInDProt_Defs.PropMeans(AtomList,AAp,iW,PI2n,D0j,OrbitLim) # any prop means using the column in the const file
    for esm in AENlist:
        sOut+="\t"+str(esm)
  
    return sOut[1:] # return the string with all the information TAB delimited (Prot, Chain, TIs) (exclude the first TAB)

#################################################
# main function
#################################################

def TIs4PDBs(sPDBlistFile,sRezFile,spathPDB,iCut,OrbitLim,iByChain,iHeader,iInClass): # main function to calculate the TIs for PDBs
    # calculates all the TIs depending on the flag iByChain
    # if iByChain=1: all the chains in the PDBs, even if you are giving a chain in the input
    # if iByChain=0: only the chains given in the input and if no chain founded, the full protein will be as one chain

    ## PDBlist="PDBlist.txt"
    ## pathPDB="PDB" # default sub-folder
    pathRes="" # if empty -> the current folder

    # type of the cut [0], direct cutoff / Roff [1], Ron / 0.0 [2](interval in the future/list)
    # iCut[0] can be
    # 1-Direct cutoff      iCut[1]=7 A default
    # 2-Abrupt truncation  iCut[1]=0.5 default
    # 3-Shifting function  iCut[1]=Roff
    # 4-Force Shifting     iCut[1]=Roff
    # 5-Switching function iCut[1]=Roff, iCut[2]=Ron
    # Ron must be > 4.8 (Gly Radius) and Roff > Ron

    # iCut=(int(1),float(7),float(6))
    # OrbitLim=[float(0.0),float(25.0),float(50.0),float(75.0),float(100.0)] # list with the limits of the orbitals

    # AAFullName[0],A3[1],A1[2],vdWradius[3],NetCh[4]
    # iWn: 1 - AmberCh[5],2 - Polar_KJ[6], 3- AtContrib2P[7],4-AtRefr[8],5-vdWArea[9],6-hardness_I-A[10],7-Electrophilicity[11],8-ElectroMulliken[12]

    # iWn=1 # the iWn th column of mean props (1-8), 1=electro

    ########################################################################################################################
    # Header for results:
    # TIs labels : TIsLabels=[("Mean","Mean properties")]
    ########################################################################################################################

    # open Protein Simple result file
    fRez=open(sRezFile,"w")

    sOut="" # print the header
    sOut+=".PDBChain"

    # if input classes are used
    if iInClass==True:
        sOut+="\tClass"

    # header flag to take add the info in the header output before the TIs labels
    h=iHeader

    # if use asked for header info
    if h==1:
        # split the path of the simple output for getting the output folder for all the other output files!!
        sOutDir,sFile=os.path.split(sRezFile)
        # open the full header info file (full inflo + TIs)
        fHeader=open(sOutDir+"\\PROT_FullHeaderRes.txt","w")
        # write all the fields from the PDB header to be checked
        sOut2=sOut # +"\thead\texpression_system\texpression_system_taxid\tname\tchain\torganism_scientific\tmolecule\texpression_system_vector_type\tec\torganism_common\texpression_system_plasmid\tengineered\texpression_system_strain\tcell_line\tcellular_location\tgene\torganism_taxid"


    # generate the labels of TIs by orbital and power
    IndTypes=["ElectroMulliken", "Polar_KJ", "vdWArea", "AtomContrib2P"]
    for Indice in IndTypes:
        for iOrb in range(5):            # for each orbital
            sOut+="\t"+Indice+"_Orb"+str(iOrb) # write TAB delimited [Label]Orb[iOrb]
##            if h==1:
##                sOut2+="\t"+Indice+"Orb"+str(iOrb)
  
    # print sOut in Protein simple result
    fRez.write(sOut+"\n")
##    # if header enabled, print out2 in full header info file
##    if h==1:
##        fHeader.write(sOut2+"\n")    

    # take the list of PDB chains
    input_file = open(sPDBlistFile,"r")
    lines = input_file.readlines()
    input_file.close()
    
    i=0
    hhFlag=0 # flag to write one time only the header of all possible PDB headers
    # process each line of PDBchains
    for line in lines :
##        try:
        i=i+1
        
        # removing strange characters
        if line[-1]=="\n" : line = line[:-1]
        if line[-1]=="\r" : line = line[:-1]

        # if line info is shorter than 4 = no PDB name, skip it and go to the next line
        if line<4:
            print "Input error for "+line+". You need at least 4 character PDB name."
            continue

        # verify the chains if exist; if not use "*" that means all the chains for the next software
        PDB=line[:4]     # take the PDB name
        CurrR=line.split("\t")
        
        # if input classes are used take the class info
        if iInClass==True:
            sClass=CurrR[1]

        # if each PDBchain line has more than 4 letters and last character is not TAB
        
        # if chain exist, take it; if not, chain is "*"
        if len(line)>4 and line[4]!="\t":
            # if By Chain is enabled take only the PDB and chain=*
            if iByChain==1:
                chain="*"
            else:
                chain=line[4] # take the input chain
        else:
            chain="*"

        print "======================================================================="
        print "Processing: "+PDB+chain+" ("+str(i)+" of "+str(len(lines))+")"
        
        PDBfile= os.path.join(spathPDB,PDB+".pdb")

        # check if a PDB exists into the local folder
        e=0
        e=MInDProt_Defs.CheckPDBDown(spathPDB,PDB)
        # if the is not PDB inside the path or PDB web, print an error and go to next line
        if e==1: 
            print PDB+" missing!"
            # fRez.write(PDB+"\t"+chain+"ERROR!\n")
            continue

        # read all the possible chains in the PDB
        p=PDBParser(PERMISSIVE=1) # permissive read for errors
        s=p.get_structure("ProtChain",PDBfile) # read the PDB structure
        model=s[0] # choose the first model if there are many

        # chain list in one protein to be processed
        # (if we need by chain will take all the PDB chains, if not only the declared chain)
        ChainList=[]

        # if "By Chain" enabled, take all the chain of the protein
        # no matther if there is a chain in the input!
        if iByChain==1:
            # get all the chains in one protein
            for chain in model.get_iterator():
                ch=str(chain.get_id())
                ChainList.append(ch)
        else:
            # get only the chain from the input
            ChainList.append(chain)

        # print "ChainList=",ChainList,"-"

        # process each chain
        # if chain is "*" inside TIcalc all the chains will be generated
        
        for sChain in ChainList:
            #  calculations for each chair and each weigth
            # in the case of iW=1 electrostatic potentials all the 3 types of indices are calculated
            # else only the mean properties are calculated
            #for PDBchain in PDBchainList: # for each PDBchain
            
            PDBchain=PDB+sChain
            if sChain=="*":
                PDBchain2=PDB # only for output reason if entire protein when chain is "*"
            else:
                PDBchain2=PDB+sChain

            sOut=""
            # if use wants header info, the items are added before the TIs
            if h==1:
                # d=ParseSelecPDBHeader.Get1DictPDBHeader(PDBfile) # all possible header
                # only the contol header fields
                d=ParseSelecPDBHeader.GetSelectivePDBHeader(AllPDBHeader=ParseSelecPDBHeader.Get1DictPDBHeader(PDBfile))
                # print d
                if hhFlag==0: # if no header wrote for the PDB headers
                    hhFlag=1
                    # add the PDB header fields
                    for k, y in d.iteritems():
                        sOut2+="\t"+str(k)
                    # generate the labels of TIs by orbital
                    for Indice in IndTypes:
                        for iOrb in range(5):# for each orbital
                            sOut2+="\t"+Indice+"_Orb"+str(iOrb)
                    fHeader.write(sOut2+"\n")
                    
                sOut2=PDBchain2
                if iInClass==True:
                    sOut2+="\t"+sClass
                
                # d=ParseSelecPDBHeader.GetSelectivePDBHeader(AllPDBHeader=ParseSelecPDBHeader.Get1DictPDBHeader(PDBfile))
                for k, y in d.iteritems():
                    sOut2+="\t"+str(y)

            sOut1=PDBchain2
            if iInClass==True:
                sOut1+="\t"+sClass   

            #======================================
            # calculate TIs
            
            sOut=""
            for iW in [7,8,9,10]: # the columns for the properties in AAconstants.txt
                sOut+="\t"+TIcalc(spathPDB,PDBchain,iCut,OrbitLim,int(iW)-1)

            # if no CAs founded in the PDB 
            if sOut=="":
                continue # go to the next PDBchain
            
            # print "Processing PDB: "+PDBfile+" ("+str(i)+" of "+str(len(lines))
            fRez.write(sOut1+sOut+"\n")

            if h==1:
                fHeader.write(sOut2+sOut+"\n") # removing the PDBchain from the sOut
##        except:
##            # fRez.write(CurrR[0]+"ERROR!\n") #the labels of PDBchains in the output file is will not be recognised by future processes
##            print "--> !!! ERROR for PDB: "+line
            
    fRez.close()
    
    if h==1:
        fHeader.close()

    print "Done!"
    return

#################################################
# GUI
#################################################

class MyMenu(wx.Frame):
    def __init__(self, parent, id, title):
        colourGUI='#FFCC66'
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(1200, 620),style=wx.CAPTION | wx.SYSTEM_MENU | wx.MINIMIZE_BOX | wx.CLOSE_BOX)
        self.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        self.SetBackgroundColour(colourGUI)
        self.statusbar = self.CreateStatusBar()        
        self.statusbar.SetStatusText("Markov Indices for Drugs and Proteins")

        #########################################################################################################
        # PROTEIN box
        #########################################################################################################
        wx.StaticBox(self, -1, 'PROTEINS', (5, 5), size=(600, 455))

        self.nProteins=wx.CheckBox(self, -1 ,'Protein calculation', (20, 30))
        self.nProteins.SetValue(True)
        self.PDB = wx.RadioButton(self, -1, 'PDB', (160,30), style=wx.RB_GROUP)
        self.FASTA = wx.RadioButton(self, -1, 'FASTA',(220, 30))

        ##################################################
        # File parameters
        ##################################################
        dy=30
        x=10
        y=50 # 20
        
        wx.StaticBox(self, -1, 'Files', (x+5, y+5), size=(340, 120))
        
        wx.StaticText(self, -1, 'Input', (x+20, y+dy*1))
        self.PDBListFile = wx.TextCtrl(self, -1, orig_dir+'PDBlist.txt',(x+90,y+dy*1) , (120,-1) ,  style=wx.TE_LEFT)
        
        br1=wx.Button(self, 20, '...', (x+215,y+dy*1-5), (30,-1)) # browse 1
        edit1=wx.Button(self, 30, 'Edit', (x+250,y+dy*1-5)) # edit 1
        
        wx.StaticText(self, -1, 'Results', (x+20, y+dy*2))
        self.results = wx.TextCtrl(self, -1, orig_dir+'PROT_SimpleRes.txt',(x+90,y+dy*2-5) , (120,-1) ,  style=wx.TE_LEFT)
        br2=wx.Button(self, 21, '...', (x+215,y+dy*2-5), (30,-1)) # browse 2
        edit2=wx.Button(self, 31, 'Edit', (x+250,y+dy*2-5)) # edit 2

        # PDB folder
        br1=wx.Button(self, 32, 'PDB Folder', (x+10,y+dy*3))
        self.PDBFolder = wx.TextCtrl(self, -1, "",(x+20+70,y+dy*3) , (250,-1) ,  style=wx.TE_LEFT)
        self.PDBFolder.SetValue(os.getcwd()+"\\PDB\\") # get the current path to be used for all functions

        # for my calculations
        #self.PDBFolder.SetValue("F:\\ProteinDataBank\\PDB2\\") # get the current path to be used for all functions

        
        ##################################################
        # Protein alpha Carbon Network parameters
        ##################################################
        x=10
        y=180
        dy=30
        
        wx.StaticBox(self, -1, 'Networks', (x+5, y+5), size=(340, 150))
        wx.StaticBox(self, -1, 'PDB CA net', (x+10, y+20), size=(320, 70))
        wx.StaticText(self, -1, 'Cutoff type', (x+20, y+dy*1+10))
        self.CutoffType=wx.SpinCtrl(self, -1, '1', (x+20, y+dy*2), (50, -1), min=1, max=5)
        wx.StaticText(self, -1, 'Roff/Cutoff', (x+90, y+dy*1+10))
        self.Roff=wx.SpinCtrl(self, -1, '7', (x+90, y+dy*2), (50, -1), min=4, max=15)
        wx.StaticText(self, -1, 'Ron', (x+160, y+dy*1+10))
        self.Ron=wx.SpinCtrl(self, -1, '6', (x+160, y+dy*2), (50, -1), min=1, max=100)

        self.fByChain=wx.CheckBox(self, -1 ,'By Chain', (x+250, y+dy*2))
        self.fByChain.SetValue(False)

        wx.StaticText(self, -1, 'PDB/FASTA Orbital Region Limits (%/AAs)', (x+20,  y+dy*3+5))
        wx.StaticText(self, -1, 'c-i-m-o', (x+20,  y+dy*4))
        self.zOrbital=wx.SpinCtrl(self, -1, '0', (x+60, y+dy*4-5), (50, -1),min=0, max=999)
        self.cOrbital=wx.SpinCtrl(self, -1, '25', (x+115, y+dy*4-5), (50, -1),min=0, max=999)
        self.iOrbital=wx.SpinCtrl(self, -1, '50', (x+170, y+dy*4-5), (50, -1),min=0, max=999)
        self.mOrbital=wx.SpinCtrl(self, -1, '75', (x+225, y+dy*4-5), (50, -1),min=0, max=999)
        self.oOrbital=wx.SpinCtrl(self, -1, '100', (x+280, y+dy*4-5), (50, -1),min=0, max=999)

      
        ##################################################
        # full header
        # pair result
        ##################################################
        # y=y+dy
        
        self.nFullHeader=wx.CheckBox(self, -1 ,'Full header results => PROT_FullHeaderRes.txt', (x+10, y+dy*6-10))
        self.nFullHeader.SetValue(False)

        y=y+dy*7
        self.nPairOut=wx.CheckBox(self, -1 ,'PROT pairs => PROT_Pairs.txt: ', (x+10, y-10))
        self.nPairOut.SetValue(False)
        wx.StaticText(self, -1, 'Until ', (x+190,  y-10))
        self.NegPairs=wx.SpinCtrl(self, -1, '1', (x+220, y-10-5), (50, -1), min=1, max=10)
        wx.StaticText(self, -1, ' time(s) of random disconnected protein pairs', (x+270, y-10))

        # activity/property of protein-chain pairs = other input file: [prot1] [prot2] [activity]
        self.nProtPairActiv=wx.CheckBox(self, -1 ,'Use Activity File', (x+30, y+dy*1-10))
        self.nProtPairActiv.SetValue(False)
        self.ProtPairActiv = wx.TextCtrl(self, -1, orig_dir+'ProtPairActivity.txt',(x+130,y+dy*1-10-5) , (330,-1) ,  style=wx.TE_LEFT)
        
        br7=wx.Button(self, 24, '...', (x+465,y+dy*1-10-5), (30,-1)) # browse 7
        edit7=wx.Button(self, 35, 'Edit', (x+500,y+dy*1-10-5)) # edit 7
        
        ##################################################
        # Class Average parameters
        ##################################################
        # Averaged results by header class or by classes inside the input file (PDBchain\tCLASS_NAME)
        x=360
        y=20
        dy=30
        
        wx.StaticBox(self, -1, 'Averaged Indices', (x, y+5), size=(230, 310))

        self.nHeader=wx.CheckBox(self, -1 ,'Class Averages => PROT_ClassAvgs.txt', (x+10, y+dy*1))
        self.nHeader.SetValue(False)
        wx.StaticText(self, -1, 'By:', (x+10, y+dy*2-10))
        self.AvgTypeInput = wx.RadioButton(self, -1, 'Input classes', (x+30,y+dy*2), style=wx.RB_GROUP)
        self.AvgTypeHeader = wx.RadioButton(self, -1, 'PDB Header',(x+120, y+dy*2))

        #################################################
        # PDB header class averaging of the TIs
        #
        # head expression_system expression_system_taxid
        # name chain organism_scientific molecule expression_system_vector_type
        # ec organism_common expression_system_plasmid engineered expression_system_strain
        # cell_line cellular_location gene organism_taxid

        dy=20

        wx.StaticText(self, -1, 'Header fields:', (x+15, y+dy*4))
        self.nCellular_location=wx.CheckBox(self, -1 ,'cellular_location', (x+15, y+dy*5))
        self.nTissue=wx.CheckBox(self, -1 ,'tissue', (x+160, y+dy*5))
        self.nOrgan=wx.CheckBox(self, -1 ,'organ', (x+15, y+dy*6))
        self.nOrganism_scientific=wx.CheckBox(self, -1 ,'organism_scientific', (x+100, y+dy*6))
        self.nOrganism_common=wx.CheckBox(self, -1 ,'organism_common', (x+15, y+dy*7))
        self.nEC=wx.CheckBox(self, -1 ,'ec', (x+160, y+dy*7))
        self.nExpression_system_vector_type=wx.CheckBox(self, -1 ,'expression_system_vector_type', (x+15, y+dy*8))
        self.nExpression_system_taxid=wx.CheckBox(self, -1 ,'expression_system_taxid', (x+15, y+dy*9))
        self.nExpression_system=wx.CheckBox(self, -1 ,'expression_system', (x+15, y+dy*10))
        self.nEngineered=wx.CheckBox(self, -1 ,'engineered', (x+15, y+dy*11))
        self.nHead=wx.CheckBox(self, -1 ,'head', (x+160, y+dy*11))
        wx.StaticText(self, -1, 'Negative cases:', (x+15, y+dy*13))
        wx.StaticText(self, -1, 'Until ', (x+15, y+dy*14))
        self.ProtAvgNegPairs=wx.SpinCtrl(self, -1, '1', (x+45, y+dy*14-5), (50, -1), min=1, max=10)
        wx.StaticText(self, -1, ' time(s) positives', (x+100, y+dy*14))
        
        ########################################################################################
        # DRUGS box
        ########################################################################################
        wx.StaticBox(self, -1, 'DRUGS', (620, 5), size=(565, 280))

        x=625
        y=20
        dy=30

        self.nDrugs=wx.CheckBox(self, -1 ,'Drug calculation', (x+5, 30))
        self.nDrugs.SetValue(False)
        
        ########################
        # File box
        ########################

        dy=30
        x=625
        y=50 # 20
        
        wx.StaticBox(self, -1, 'Files', (x+5, y+5), size=(340, 90))
        wx.StaticText(self, -1, 'SMILE list', (x+20, y+dy*1))
        self.SMILEfile = wx.TextCtrl(self, -1, orig_dir+'SMILEs.txt',(x+90,y+dy*1-5) , (120,-1) ,  style=wx.TE_LEFT)
        
        br3=wx.Button(self, 22, '...', (x+215,y+dy*1-5), (30,-1)) # browse 3
        edit3=wx.Button(self, 33, 'Edit', (x+250,y+dy*1-5)) # edit 3
        
        wx.StaticText(self, -1, 'Results', (x+20, y+dy*2))
        self.DrugResults = wx.TextCtrl(self, -1, orig_dir+'DRUG_SimpleRes.txt',(x+90,y+dy*2-5) , (120,-1) ,  style=wx.TE_LEFT)
        br4=wx.Button(self, 23, '...', (x+215,y+dy*2-5), (30,-1)) # browse 4
        edit4=wx.Button(self, 34, 'Edit', (x+250,y+dy*2-5)) # edit 4

        self.nDrugClass=wx.CheckBox(self, -1 ,'Averages by Input Classes => DRUG_ClassAvgs.txt', (x+10, y+dy*4-10))
        self.nDrugClass.SetValue(False)

        y=y+dy*5
        self.nDrugPairs=wx.CheckBox(self, -1 ,'Node pairs => DRUG_Pairs.txt', (x+10, y-10))
        self.nDrugPairs.SetValue(False)

        # activity/property of drug-drug pairs = other input file: [drug1] [drug2] [activity]
        wx.StaticText(self, -1, 'Using activity pair file: ', (x+20,  y+dy-10))
        self.DrugPairActiv = wx.TextCtrl(self, -1, orig_dir+'DrugPairActivity.txt',(x+20,y+dy*2-10-5) , (205,-1) ,  style=wx.TE_LEFT)
        br8=wx.Button(self, 36, '...', (x+235,y+dy*2-10-5), (30,-1)) # browse 8
        edit8=wx.Button(self, 37, 'Edit', (x+270,y+dy*2-10-5)) # edit 8

        # CHEMBL options
        wx.StaticBox(self, -1, 'CHEMBL Averaged Indices', (x+350, 20+5), size=(200, 250))
        # FIELD to use in averaging: STANDARD_TYPE, ASSAY_CHEMBLID, ASSAY_TYPE, TARGET_CHEMBLID
        # ORGANISM, TARGET_TYPE, TARGET_MAPPING
        self.nCHEMBL=wx.CheckBox(self, -1 ,'CHEMBL Fields', (x+350+15, 25+20))
        
        self.nSTANDARD_TYPE=wx.CheckBox(self, -1 ,'STANDARD_TYPE', (x+350+35, 25+40))
        self.nASSAY_CHEMBLID=wx.CheckBox(self, -1 ,'ASSAY_CHEMBLID', (x+350+35, 25+60))
        self.nASSAY_TYPE=wx.CheckBox(self, -1 ,'ASSAY_TYPE', (x+350+35, 25+80))
        self.nTARGET_CHEMBLID=wx.CheckBox(self, -1 ,'TARGET_CHEMBLID', (x+350+35, 25+100))
        self.nORGANISM=wx.CheckBox(self, -1 ,'ORGANISM', (x+350+35, 25+120))
        self.nTARGET_TYPE=wx.CheckBox(self, -1 ,'TARGET_TYPE', (x+350+35, 25+140))
        self.nTARGET_MAPPING=wx.CheckBox(self, -1 ,'TARGET_MAPPING', (x+350+35, 25+160))

        #self.nSTANDARD_TYPE.Enable(False)

        wx.StaticText(self, -1, 'NEG vs. POS cases:', (x+350+15,  25+200))
        # cutoff for Zij -> C=1/0
        self.CutOffZij=wx.SpinCtrl(self, -1, '0', (x+350+15, 25+220), (50, -1), min=-99999, max=99999)
        wx.StaticText(self, -1, ' Z-score cut-off', (x+350+70, 25+225))
        
        
        #self.WeightList = ['AmberCh', 'Polar_KJ', 'AtContrib2P', 'AtRefr', 'vdWArea','hardness_I-A','Electrophilicity','ElectroMulliken']

        #################################
        # Prot-Drug pairs
        #################################

        wx.StaticBox(self, -1, 'PROTEIN-DRUG', (620, 300), size=(565, 130))

##        x=10
        y=330
        dy=20
        
        self.nProtDrugPairs=wx.CheckBox(self, -1 ,'PROTEIN-DRUG PAIRS (->PROT_DRUG_Pairs.txt) using activity file:', (x+10, y))
        self.nProtDrugPairs.SetValue(False)
        
        # activity/property of drug-drug pairs = other input file: [drug1] [drug2] [activity]
        self.ProtDrugPairActiv = wx.TextCtrl(self, -1, orig_dir+'ProtDrugPairActivity.txt',(x+20,y+dy*1) , (420,-1) ,  style=wx.TE_LEFT)
        br9=wx.Button(self, 38, '...', (x+445,y+dy*1), (30,-1)) # browse 9
        edit9=wx.Button(self, 39, 'Edit', (x+480,y+dy*1)) # edit 9

        wx.StaticText(self, -1, 'If only positives cases: ', (x+20, y+dy*2+10))
        self.PDPairs=wx.SpinCtrl(self, -1, '1', (x+20, y+dy*3-5+10), (50, -1), min=1, max=10)
        wx.StaticText(self, -1, ' time(s) of negative random generated PDB-drug pairs', (x+70, y+dy*3+10))

        self.ProtDrugPairActiv.Enable(False)
        self.PDPairs.Enable(False)

        ####################################################
        # New Drugs against the pre-calculated averages
        ####################################################

        x=615
        y=330+135
        dy=30
        
        wx.StaticBox(self, -1, 'NEW DRUG deviations from pre-calculated averages', (x+5, y+5-30), size=(565, 120))
        
        wx.StaticText(self, -1, 'New SMILES', (x+20, y+5))
        self.NewSmileList = wx.TextCtrl(self, -1, orig_dir+'SMILEs.txt',(x+90,y) , (250,-1) ,  style=wx.TE_LEFT)
        br45=wx.Button(self, 45, '...', (x+350,y), (30,-1)) # browse 45
        edit46=wx.Button(self, 46, 'Edit', (x+385,y))       # edit 46

        wx.StaticText(self, -1, 'Pre-calc AVGs', (x+20, y+dy*1+5))
        self.PreAverages = wx.TextCtrl(self, -1, orig_dir+'DRUG_SimpleRes.txt',(x+90,y+dy*1) , (250,-1) ,  style=wx.TE_LEFT)
        br40=wx.Button(self, 40, '...', (x+350,y+dy*1), (30,-1)) # browse 40
        edit41=wx.Button(self, 41, 'Edit', (x+385,y+dy*1))       # edit 41
        
        wx.StaticText(self, -1, 'Results', (x+20, y+dy*2))
        self.PreResults = wx.TextCtrl(self, -1, orig_dir+'DRUG_New_DevsByCHEMBL_Avgs.txt',(x+90,y+dy*2-5) , (250,-1) ,  style=wx.TE_LEFT)
        br42=wx.Button(self, 42, '...', (x+350,y+dy*2-5), (30,-1)) # browse 42
        edit43=wx.Button(self, 43, 'Edit', (x+385,y+dy*2-5))       # edit 43

        self.bNewDrugs=GenBitmapTextButton(self, 44, wx.Bitmap(orig_dir+'images/gtk-execute.png'), 'Exec', (x+480, 500))
            
        ##################################################################################
        ## BUTTONS
        ################################################################################## 630
        
        self.bSubmit=GenBitmapTextButton(self, 25, wx.Bitmap(orig_dir+'images/gtk-execute.png'), 'Run', (20, 500))
        self.bHelp=GenBitmapTextButton(self, 26, wx.Bitmap(orig_dir+'images/gtk-help.png'), 'Help', (100, 500))
        self.bAbout=GenBitmapTextButton(self, 27, wx.Bitmap(orig_dir+'images/gtk-about.png'), 'About', (180, 500))
        self.bQuit=GenBitmapTextButton(self, 28, wx.Bitmap(orig_dir+'images/gtk-quit.png'), 'Quit', (260, 500))

        wx.EVT_BUTTON(self, 20, self.OnBrowsePDBList)
        wx.EVT_BUTTON(self, 21, self.OnBrowseRes)
        wx.EVT_BUTTON(self, 22, self.OnBrowseSMILEList)
        wx.EVT_BUTTON(self, 23, self.OnBrowseDrugRes)
        wx.EVT_BUTTON(self, 24, self.OnBrowseProtPairActiv)
        wx.EVT_BUTTON(self, 36, self.OnBrowseDrugPairActiv)
        wx.EVT_BUTTON(self, 38, self.OnBrowseProtDrugPairActiv)
        
        wx.EVT_BUTTON(self, 25, self.OnSubmit)
        wx.EVT_BUTTON(self, 26, self.OnHelp)
        wx.EVT_BUTTON(self, 27, self.OnAbout)
        wx.EVT_BUTTON(self, 28, self.OnQuit)
##        wx.EVT_BUTTON(self, 29, self.OnPlotGraph)
        
        wx.EVT_BUTTON(self, 30, self.OnEditPDBList)
        wx.EVT_BUTTON(self, 31, self.OnEditRes)
        wx.EVT_BUTTON(self, 32, self.OnBrowsePDBFolder)
        wx.EVT_BUTTON(self, 33, self.OnEditSMILEList)
        wx.EVT_BUTTON(self, 34, self.OnEditDrugRes)
        wx.EVT_BUTTON(self, 35, self.OnEditProtPairActiv)
        wx.EVT_BUTTON(self, 37, self.OnEditDrugPairActiv)
        wx.EVT_BUTTON(self, 39, self.OnEditProtDrugPairActiv)

        wx.EVT_BUTTON(self, 45, self.OnBrowseNewSMILE)
        wx.EVT_BUTTON(self, 46, self.OnEditNewSMILE)
        
        wx.EVT_BUTTON(self, 40, self.OnBrowsePreAvgIn)
        wx.EVT_BUTTON(self, 41, self.OnEditPreAvgIn)

        wx.EVT_BUTTON(self, 42, self.OnBrowsePreAvgOut)
        wx.EVT_BUTTON(self, 43, self.OnEditPreAvgOut)

        wx.EVT_BUTTON(self, 44, self.OnExec)

        # self.Bind(wx.EVT_COMBOBOX, self.OnSelect2S, )
        wx.EVT_CHECKBOX(self,self.nDrugs.GetId(),self.EventDrugs)
        wx.EVT_CHECKBOX(self,self.nProteins.GetId(),self.EventProteins)
        wx.EVT_CHECKBOX(self,self.nHeader.GetId(),self.EventProtAveraged)
        wx.EVT_CHECKBOX(self,self.nProtDrugPairs.GetId(),self.EventProtDrugPairs)
        
        self.Bind(wx.EVT_RADIOBUTTON, self.EventProtInput, id=self.AvgTypeHeader.GetId())

        # header
        print "\n***************************************************************************"
        print "MInD-Prot\nMarkov Indices for Drugs and Proteins\n(ver. 3.0, 2012)\n\nby\nCristian Robert Munteanu (muntisa@gmail.com)\nHumberto Gonzalez-Diaz (gonzalezdiazh@yahoo.es)"
        print "***************************************************************************\n"
        
        self.orig_dir = globals()["orig_dir"] #takes the original folder for the un modified files in the browse controls
        return

    def EventProtDrugPairs(self, event):
        # disable all controls
        if self.nProtDrugPairs.GetValue()==False:
            self.ProtDrugPairActiv.Enable(False)
            self.PDPairs.Enable(False)

        # enable all controls
        if self.nProtDrugPairs.GetValue()==True:
            self.ProtDrugPairActiv.Enable(True)
            self.PDPairs.Enable(True)
            # enable proteins
            if self.nDrugs.GetValue()==False:
                self.nDrugs.Enable(True)
            # enable drugs
            if self.nProteins.GetValue()==False:
                self.nProteins.Enable(True)            
        return
    
    def EventDrugs(self, event):
        # disable all controls
        if self.nDrugs.GetValue()==False:
            self.SMILEfile.Enable(False)
            self.DrugResults.Enable(False)
            self.nDrugClass.Enable(False)
            self.nDrugPairs.Enable(False)
            self.DrugPairActiv.Enable(False)

            self.nProtDrugPairs.Enable(False)
            self.nProtDrugPairs.SetValue(False)
            self.ProtDrugPairActiv.Enable(False)
            self.PDPairs.Enable(False)
            
        # enable all controls
        if self.nDrugs.GetValue()==True:
            self.SMILEfile.Enable(True)
            self.DrugResults.Enable(True)
            self.nDrugClass.Enable(True)
            self.nDrugPairs.Enable(True)
            self.DrugPairActiv.Enable(True) 

            # enable prot-drug pairs only if protein is enabled
            if self.nProteins.GetValue()==True:
                self.nProtDrugPairs.Enable(True)
                self.ProtDrugPairActiv.Enable(True)
                self.PDPairs.Enable(True)
            
        return

    def EventProteins(self, event):
        # disable all controls
        if self.nProteins.GetValue()==False:
            self.PDBListFile.Enable(False)
            self.results.Enable(False)
            self.PDBFolder.Enable(False)
            self.CutoffType.Enable(False)
            self.Roff.Enable(False)
            self.Ron.Enable(False)
            self.zOrbital.Enable(False)
            self.cOrbital.Enable(False)
            self.iOrbital.Enable(False)
            self.mOrbital.Enable(False)
            self.oOrbital.Enable(False)
            self.fByChain.Enable(False)
            self.nFullHeader.Enable(False)
            self.nPairOut.Enable(False)
            self.NegPairs.Enable(False)
            self.nProtPairActiv.Enable(False)
            self.ProtPairActiv.Enable(False)
            self.nHeader.Enable(False)
            self.AvgTypeHeader.Enable(False)
            self.AvgTypeInput.Enable(False)
            
            self.nCellular_location.Enable(False)
            self.nTissue.Enable(False)
            self.nOrgan.Enable(False)
            self.nOrganism_scientific.Enable(False)
            self.nOrganism_common.Enable(False)
            self.nEC.Enable(False)
            self.nExpression_system_vector_type.Enable(False)
            self.nExpression_system_taxid.Enable(False)
            self.nExpression_system.Enable(False)
            self.nEngineered.Enable(False)
            self.nHead.Enable(False)
    
            self.nProtDrugPairs.Enable(False)
            self.nProtDrugPairs.SetValue(False)
            self.ProtDrugPairActiv.Enable(False)
            self.PDPairs.Enable(False)
                      
        # enable all controls
        if self.nProteins.GetValue()==True:
            self.PDBListFile.Enable(True)
            self.results.Enable(True)
            self.PDBFolder.Enable(True)
            self.CutoffType.Enable(True)
            self.Roff.Enable(True)
            self.Ron.Enable(True)
            self.zOrbital.Enable(True)
            self.cOrbital.Enable(True)
            self.iOrbital.Enable(True)
            self.mOrbital.Enable(True)
            self.oOrbital.Enable(True)
            self.fByChain.Enable(True)
            self.nFullHeader.Enable(True)
            self.nPairOut.Enable(True)
            self.NegPairs.Enable(True)
            self.nProtPairActiv.Enable(True)
            self.ProtPairActiv.Enable(True)
            self.nHeader.Enable(True)
            self.AvgTypeHeader.Enable(True)
            self.AvgTypeInput.Enable(True)

            self.nCellular_location.Enable(True)
            self.nTissue.Enable(True)
            self.nOrgan.Enable(True)
            self.nOrganism_scientific.Enable(True)
            self.nOrganism_common.Enable(True)
            self.nEC.Enable(True)
            self.nExpression_system_vector_type.Enable(True)
            self.nExpression_system_taxid.Enable(True)
            self.nExpression_system.Enable(True)
            self.nEngineered.Enable(True)
            self.nHead.Enable(True)

            if self.nDrugs.GetValue()==True:
                self.nProtDrugPairs.Enable(True)
                self.ProtDrugPairActiv.Enable(True)
                self.PDPairs.Enable(True)
            
        return

    def EventProtAveraged(self, event):
        # disable all controls
        if self.nHeader.GetValue()==False:
            self.AvgTypeHeader.Enable(False)
            self.AvgTypeInput.Enable(False)

            self.nCellular_location.Enable(False)
            self.nTissue.Enable(False)
            self.nOrgan.Enable(False)
            self.nOrganism_scientific.Enable(False)
            self.nOrganism_common.Enable(False)
            self.nEC.Enable(False)
            self.nExpression_system_vector_type.Enable(False)
            self.nExpression_system_taxid.Enable(False)
            self.nExpression_system.Enable(False)
            self.nEngineered.Enable(False)
            self.nHead.Enable(False)
                      
        # enable all controls
        if self.nHeader.GetValue()==True:
            self.nFullHeader.SetValue(True) # if you need averages by PDB you need full header file info
            self.AvgTypeHeader.Enable(True)
            self.AvgTypeInput.Enable(True)
            self.nCellular_location.Enable(True)
            self.nTissue.Enable(True)
            self.nOrgan.Enable(True)
            self.nOrganism_scientific.Enable(True)
            self.nOrganism_common.Enable(True)
            self.nEC.Enable(True)
            self.nExpression_system_vector_type.Enable(True)
            self.nExpression_system_taxid.Enable(True)
            self.nExpression_system.Enable(True)
            self.nEngineered.Enable(True)
            self.nHead.Enable(True)
        return
    
    def EventProtInput(self, event):
        # disable all controls
        if self.AvgTypeInput.GetValue()==True:
            self.nCellular_location.Enable(False)
            self.nTissue.Enable(False)
            self.nOrgan.Enable(False)
            self.nOrganism_scientific.Enable(False)
            self.nOrganism_common.Enable(False)
            self.nEC.Enable(False)
            self.nExpression_system_vector_type.Enable(False)
            self.nExpression_system_taxid.Enable(False)
            self.nExpression_system.Enable(False)
            self.nEngineered.Enable(False)
            self.nHead.Enable(False)
                      
        # enable all controls
        if self.AvgTypeInput.GetValue()==False:
            self.nCellular_location.Enable(True)
            self.nTissue.Enable(True)
            self.nOrgan.Enable(True)
            self.nOrganism_scientific.Enable(True)
            self.nOrganism_common.Enable(True)
            self.nEC.Enable(True)
            self.nExpression_system_vector_type.Enable(True)
            self.nExpression_system_taxid.Enable(True)
            self.nExpression_system.Enable(True)
            self.nEngineered.Enable(True)
            self.nHead.Enable(True)
        return

    def OnEditProtDrugPairActiv(self, event):
        Net=self.ProtDrugPairActiv.GetValue()
        if Net.find("\\")==-1:
            Net=self.orig_dir+Net 
        iFile=Net        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnBrowseProtDrugPairActiv(self, event):
        dir = os.getcwd()
        wildcard = "All files (*.*)|*.*"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.ProtDrugPairActiv.SetValue(path)
            
    def OnEditProtPairActiv(self, event):
        Net=self.ProtPairActiv.GetValue()
        if Net.find("\\")==-1:
            Net=self.orig_dir+Net 
        iFile=Net        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnBrowseProtPairActiv(self, event):
        dir = os.getcwd()
        wildcard = "All files (*.*)|*.*"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.ProtPairActiv.SetValue(path)

    def OnEditDrugPairActiv(self, event):
        Net=self.DrugPairActiv.GetValue()
        if Net.find("\\")==-1:
            Net=self.orig_dir+Net 
        iFile=Net        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnBrowseDrugPairActiv(self, event):
        dir = os.getcwd()
        wildcard = "All files (*.*)|*.*"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.DrugPairActiv.SetValue(path)

    def OnBrowsePDBList(self, event):
        dir = os.getcwd()
        wildcard = "All files (*.*)|*.*"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.PDBListFile.SetValue(path)
    def OnBrowseRes(self, event):
        dir = os.getcwd()
        wildcard = "Text files (*.txt)|*.txt"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.results.SetValue(path)

    def OnEditPDBList(self, event):
        Net=self.PDBListFile.GetValue()
        if Net.find("\\")==-1:
            Net=self.orig_dir+Net 
        iFile=Net        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            
    def OnEditRes(self, event):
        results=self.results.GetValue()
        if results.find("\\")==-1:
            results=self.orig_dir+results
        iFile=results
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            
    def OnBrowsePDBFolder(self,event):
        dialog = wx.DirDialog(None, "Choose a directory:",os.getcwd(), style=wx.DD_DEFAULT_STYLE|wx.DD_DIR_MUST_EXIST|wx.DD_CHANGE_DIR)
        if dialog.ShowModal() == wx.ID_OK:
            self.SetStatusText('You selected: %s\n' %dialog.GetPath())
            self.dirname=dialog.GetPath() # variable with the folder
            self.PDBFolder.SetValue(self.dirname)
            os.chdir(orig_dir)
        dialog.Destroy
        return

    def OnBrowseSMILEList(self, event):
        dir = os.getcwd()
        wildcard = "All files (*.*)|*.*"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.SMILEfile.SetValue(path)
            
    def OnBrowseDrugRes(self, event):
        dir = os.getcwd()
        wildcard = "Text files (*.txt)|*.txt"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.results.SetValue(path)

    def OnEditSMILEList(self, event):
        Net=self.SMILEfile.GetValue()
        if Net.find("\\")==-1:
            Net=self.orig_dir+Net 
        iFile=Net        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            
    def OnEditDrugRes(self, event):
        results=self.DrugResults.GetValue()
        if results.find("\\")==-1:
            results=self.orig_dir+results
        iFile=results
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnHelp(self, event):
        frame3 = wx.Frame(None, -1, "MInD-Prot Help", size=(720, 560),style=wx.CAPTION | wx.SYSTEM_MENU | wx.MINIMIZE_BOX | wx.CLOSE_BOX)
        frame3.SetIcon(wx.Icon(orig_dir+'faviconMC.ico', wx.BITMAP_TYPE_ICO))
        frame3.SetPosition((5,5))
        MInDProt_HTMLDoc.MyHtmlPanel(frame3,-1)
        frame3.Show(True)

    def OnBrowsePreAvgIn(self, event):
        dir = os.getcwd()
        wildcard = "All files (*.*)|*.*"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.OnBrowsePreAvgIn.SetValue(path)
            
    def OnBrowsePreAvgOut(self, event):
        dir = os.getcwd()
        wildcard = "All files (*.*)|*.*"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.OnBrowsePreAvgOut.SetValue(path)
# ----------------------------------------------------------
    def OnBrowseNewSMILE(self, event):
        dir = os.getcwd()
        wildcard = "All files (*.*)|*.*"
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', wildcard=wildcard, style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.NewSmileList.SetValue(path)

    def OnEditNewSMILE(self, event):
        Net=self.NewSmileList.GetValue()
        if Net.find("\\")==-1:
            Net=self.orig_dir+Net 
        iFile=Net        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            
    def OnEditPreAvgIn(self, event):
        Net=self.PreAverages.GetValue()
        if Net.find("\\")==-1:
            Net=self.orig_dir+Net 
        iFile=Net        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnEditPreAvgOut(self, event):
        Net=self.PreResults.GetValue()
        if Net.find("\\")==-1:
            Net=self.orig_dir+Net 
        iFile=Net        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','MInD-Prot: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'MInD-Prot: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    # CHEMBL calculation
    def OnExec(self, event):
        # inputs for the CHEMBL calc
        sFileInSmiles=self.NewSmileList.GetValue()  # SMILEs.txt
        sFileInAvg=self.PreAverages.GetValue()      # DRUG_SimpleRes.txt
        sFileOut=self.PreResults.GetValue()         # DRUG_New_DevsByCHEMBL_Avgs.txt
        
        sFieldNames=["STANDARD_TYPE","ASSAY_CHEMBLID","ASSAY_TYPE","TARGET_CHEMBLID","ORGANISM","TARGET_TYPE","TARGET_MAPPING"]

        # header indices
        IndLabels=[]
        AtomTypes=['MP_tot','MP_Csat','MP_Cinst','MP_Halog','MP_Hetero','MP_Hx']
        IndTypes=["ElectroMulliken", "Polar_KJ", "vdWArea", "AtomContrib2P"]
        for Indice in IndTypes:
            for atom in AtomTypes:                     # for each atom type
                IndLabels.append(Indice+"_"+str(atom)) # [Propriety][AtomType]

        TIsHeader=IndLabels # list with TI names

        # Calculate deviations of New Drugs against pre-calculated averages
        NewDrugDevs_PreCalcAvgs(sFileInSmiles,sFileInAvg,sFileOut,sFieldNames,TIsHeader)

        return
        
    def OnAbout(self, event):
        description = """
MInD-Prot is a free Windows application that calculates the Markov indices
for drugs (given as SMILE codes) and proteins (given as PDB file) based on the Mean Properties indices.

You can calculate separately the indices for PROTEINs and DRUGs but you can get mixed indices too.
The tool can generate pair outputs for PDBchain-PDBchain pairs, drug-drug pairs and PDBchain-drug pairs.
In addition, MInD-Prot can calculate averaged indices by classes and get the PDB header full info for the proteins.

All these indices are using the node probabilities resulted as a Markov normalization of
the classical connectivity matrix.
It is a Python application, with wxPython GUI.
"""
        licence = """
MInD-Prot ver. 3.0
Copyright MInD-Prot 2012
"""
        info = wx.AboutDialogInfo()
        info.SetIcon(wx.Icon(self.orig_dir+'images/logo.png', wx.BITMAP_TYPE_PNG))
        info.SetName('MInD-Prot: Markov Indices for Drugs and Proteins')
        info.SetVersion('3.0')
        info.SetDescription(description)
        info.SetCopyright('2012 MInD-Prot')
        info.SetLicence(licence)
        info.AddDeveloper('Cristian Robert Munteanu, Spain (muntisa@gmail.com)')
        info.AddDeveloper('Humberto Gonzalez-Diaz, Spain (gonzalezdiazh@yahoo.es)')
        wx.AboutBox(info)
        return
        
    def OnQuit(self, event):
        self.Close()
        return

##    def OnPlotGraph(self, event):
##        pass
##        spathPDB=self.PDBFolder.GetValue()
##        iCut=(int(self.CutoffType.GetValue()),float(self.Roff.GetValue()),float(self.Ron.GetValue()))
##        # iCut=(int(1),float(7),float(6))
##        OrbitLim=[float(self.zOrbital.GetValue()),float(self.cOrbital.GetValue()),float(self.iOrbital.GetValue()),float(self.mOrbital.GetValue()),float(self.oOrbital.GetValue())]
##        try:
##            Plot1PDB_CAnet(spathPDB,self.PDB4Graph.GetValue(),iCut,OrbitLim)
##            print "\nThe correspondent graph plot was created in the same folder with the outputs files as PNG."
##        except:
##            print "This PDBchain don't exist. Please try a real one."
##            
##        return
    
    def OnSubmit(self, event):
        ###############################################################################
        # MAIN calculation

        tt0=time.clock()
        self.statusbar.SetStatusText("Running the calculations ...")

        sPDBlistFile=self.PDBListFile.GetValue()
        sRezFile=self.results.GetValue()
        spathPDB=self.PDBFolder.GetValue()

        iInClass=self.AvgTypeInput.GetValue()
        if iInClass==True:
            iShift=2
        else:
            iShift=1

        # split the path of the simple output for getting the output folder for all the other output files!!
        sOutDir,sFileX=os.path.split(sRezFile)

        ################################################################################################
        # PROTEINS

        # if Proteins enabled
        if self.nProteins.GetValue()==True:
            print "-> PROTEIN calculations enabled ..."

            if self.PDB.GetValue()==True: # if PDB inputs
                print "--> PDB input calculations enabled ..."
                ##########################
                # Protein simple results: TIs
                
                # type of the cut [0], direct cutoff / Roff [1], Ron / 0.0 [2](interval in the future/list)
                # iCut[0] can be
                # 1-Direct cutoff      iCut[1]=7 A default
                # 2-Abrupt truncation  iCut[1]=0.5 default
                # 3-Shifting function  iCut[1]=Roff
                # 4-Force Shifting     iCut[1]=Roff
                # 5-Switching function iCut[1]=Roff, iCut[2]=Ron
                # Ron must be > 4.8 (Gly Radius) and Roff > Ron

                iCut=(int(self.CutoffType.GetValue()),float(self.Roff.GetValue()),float(self.Ron.GetValue()))
                # iCut=(int(1),float(7),float(6))

                OrbitLim=[float(self.zOrbital.GetValue()),float(self.cOrbital.GetValue()),float(self.iOrbital.GetValue()),float(self.mOrbital.GetValue()),float(self.oOrbital.GetValue())]
                # OrbitLim=[float(0.0),float(25.0),float(50.0),float(75.0),float(100.0)] # list with the limits of the orbitals

                # AAFullName[0],A3[1],A1[2],vdWradius[3],NetCh[4]
                # iWn: 1 - AmberCh[5],2 - Polar_KJ[6], 3- AtContrib2P[7],4-AtRefr[8],5-vdWArea[9],6-hardness_I-A[10],7-Electrophilicity[11],8-ElectroMulliken[12]

             
                # iWn=7-10 # the iWn th column of mean props

                if self.fByChain.GetValue()==True:
                    iByChain=1
                else:
                    iByChain=0

                # if user need full header result file
                iHeader=0 # no header result file
                if self.nFullHeader.GetValue()==True:
                    iHeader=1
        
                TIs4PDBs(sPDBlistFile,sRezFile,spathPDB,iCut,OrbitLim,iByChain,iHeader,iInClass)

            if self.FASTA.GetValue()==True: # if FASTA inputs
                print "--> FASTA input calculations enabled ..."

                # to modified the limits!!!
                iWin=round(float(10)/float(4)) # example of Windows
                OrbitLim=[float(0.0),float(iWin),float(iWin*2),float(iWin*3),float(10-1)]
                TIs4Fasta(sPDBlistFile,sRezFile,OrbitLim)

            ###################################################################
            # Protein Pairs

            # if Protein pairs enabled
            if self.nPairOut.GetValue()==True:
                print "--> PROTEIN pairs enabled ..."
                outFile=sOutDir+"\\PROT_Pairs.txt"

                # if protein pair activity enabled
                if self.nProtPairActiv.GetValue()==True:
                    print "---> PROTEIN pair activity enabled ..."
                    # use the activity file [ProtPairActiv]
                    TIs2Pairs.TIs2PairsActivFile(self.ProtPairActiv.GetValue(),sRezFile,outFile,iShift)
                
                else:
                    # by considering the same PDB and generating negative pairs until X time
                    print "---> PROTEIN pairing using the same PDB and random negative pairs ..."
                    
                    m=int(self.NegPairs.GetValue()) # times the linked
                    outFile=sOutDir+"\\PROT_Pairs.txt"
                    # external function to transform results in pairs
                    TIs2Pairs.TIs2PairFile(sRezFile,outFile,m,iShift)

            ################################
            # Protein Averaged Indices

            # if CLASS AVERAGES enabled
            if self.nHeader.GetValue()==True and self.FASTA.GetValue()==False: # not avg by input for FASTA
                print "--> PROTEIN class averages enabled ..."
                
                # if Input classes enabled
                if self.AvgTypeInput.GetValue()==True: 
                    print "---> PROTEIN averages by input classes enabled ..."
                    sInn=self.results.GetValue() #input the simple results
                    sOutt=sOutDir+"\\PROT_ClassAvgs.txt"

                    sActiv=self.PDBListFile.GetValue() # classes from input file
                    Col=2 # column in the input file of the class (PDBchain\tClass)
                    # averages by input file
                    # check the input classes
                    InputProtClasses=PDBHeaderAvg.GetInputClass(sActiv,Col)
                    # if there is at least one class
                    if len(InputProtClasses)!=0:
                        PDBHeaderAvg.AverageTIs4Col(sInn,sOutt,InputProtClasses,Col)
                    else:
                        print "!!! Error for PROTEIN Averages !!! The class column in the input is missing!"
                
                # if PDB header enabled
                if self.AvgTypeHeader.GetValue()==True:
                    print "---> PROTEIN PDB header fields enabled ..."
                    # getting the choise of header fields for average
                    # HeaderList=['cellular_location','tissue','organ','organism_scientific','organism_common','ec','expression_system_vector_type','expression_system_taxid','expression_system','engineered','head']
                    # flag_Header= zeros((len(HeaderList)))  # matrix with flags for each header type, default=0
##                    if self.nCellular_location.GetValue()==True:flag_Header[0]=1
##                    if self.nTissue.GetValue()==True:flag_Header[1]=1
##                    if self.nOrgan.GetValue()==True:flag_Header[2]=1
##                    if self.nOrganism_scientific.GetValue()==True:flag_Header[3]=1
##                    if self.nOrganism_common.GetValue()==True:flag_Header[4]=1
##                    if self.nEC.GetValue()==True:flag_Header[5]=1
##                    if self.nExpression_system_vector_type.GetValue()==True:flag_Header[6]=1
##                    if self.nExpression_system_taxid.GetValue()==True:flag_Header[7]=1
##                    if self.nExpression_system.GetValue()==True:flag_Header[8]=1
##                    if self.nEngineered.GetValue()==True:flag_Header[9]=1
##                    if self.nHead.GetValue()==True:flag_Header[10]=1
                    
                    HeaderList=[] # user list with fields
                    if self.nCellular_location.GetValue()==True:HeaderList.append("cellular_location")
                    if self.nTissue.GetValue()==True:HeaderList.append("tissue")
                    if self.nOrgan.GetValue()==True:HeaderList.append("organ")
                    if self.nOrganism_scientific.GetValue()==True:HeaderList.append("organism_scientific")
                    if self.nOrganism_common.GetValue()==True:HeaderList.append("organism_common")
                    if self.nEC.GetValue()==True:HeaderList.append("ec")
                    if self.nExpression_system_vector_type.GetValue()==True:HeaderList.append("expression_system_vector_type")
                    if self.nExpression_system_taxid.GetValue()==True:HeaderList.append("expression_system_taxid")
                    if self.nExpression_system.GetValue()==True:HeaderList.append("expression_system")
                    if self.nEngineered.GetValue()==True:HeaderList.append("engineered")
                    if self.nHead.GetValue()==True:HeaderList.append("head")

                    # if it is at least one choise of header item
                    # if sum(flag_Header)!=0:
                    if len(HeaderList)!=0:
                        sIn=sOutDir+"\\PROT_FullHeaderRes.txt"       # full header file
                        sOut=sOutDir+"\\PROT_ClassAvgs.txt"          # output averaged indices by class of specific header
                        # run CLASS AVERAGES by Header idems
                        # PDBHeaderAvg.AverageTIs_PDBHeader(sIn,sOut,HeaderList,flag_Header,iShift) # old
                        PDBHeaderAvg.AverageTIs_PDBHeader(sIn,sOut,HeaderList,iShift)
                        
        #######################################################################
	# DRUGS

        # Drugs enabled
        if self.nDrugs.GetValue()==True:
            print "\n-> DRUG calculation enabled ..."
            sSMILEfile=self.SMILEfile.GetValue()
            sRezFile=self.DrugResults.GetValue()
            sOutDir,sFileX=os.path.split(sRezFile)

            # GET pairs LABEL-SMILE If Simple Calculation or Input class Averages
            if self.nCHEMBL.GetValue()==False:
                # get the SMILES from the input file
                input_file = open(sSMILEfile,"r")
                lines = input_file.readlines()
                input_file.close()

                # from the SMILE file: list of pairs : Label - SMILE for each drug
                ListSmiles=[]
                for line in lines :
                    # removing strange characters
                    if line[-1]=="\n" : line = line[:-1]
                    if line[-1]=="\r" : line = line[:-1]
                    DataLine=line.split("\t")
                    if len(DataLine)>1:
                        ListSmiles.append(DataLine)

            # GET pairs LABEL-SMILE if CHEMBL calculation
            # Orig SMILE file from GUI -> Intermediate TI simple calc result
            # -> use both to generate -> Intermediate EXTRA input for CHEMBL calc -> output file from GUI
            if self.nCHEMBL.GetValue()==True:
                #------------------------------------------------------------
                # Read the INPUT
                ListSmiles=GetDrugID_SMILEPairs(sSMILEfile)

                #------------------------------------------------------------
                # CHANGE simple OUTPUT with temp file and take the val for the real output
                # the last output file = from the interface!
                sRezFile3=sRezFile
                # intermediate output = simple TI calculation output to be used for CHEMBL calc
                sRezFile=sRezFile[:-4]+"_temp.txt"

            # header indices
            IndLabels=[]
            AtomTypes=['MP_tot','MP_Csat','MP_Cinst','MP_Halog','MP_Hetero','MP_Hx']
            IndTypes=["ElectroMulliken", "Polar_KJ", "vdWArea", "AtomContrib2P"]
            for Indice in IndTypes:
                for atom in AtomTypes:                     # for each atom type
                    IndLabels.append(Indice+"_"+str(atom)) # [Propriety][AtomType]
            
            fDrugTI=open(sRezFile,"w") # open to write the output file

            # write the header
            sOut="DrugID\tSMILE_formula"

            # if input averages enabled add header in the simple results
            if self.nDrugClass.GetValue()==True:
                sOut+="\t"+"Activity"

            # add TI labels
            for Indice in IndLabels:
                sOut+="\t"+Indice
            fDrugTI.write(sOut+"\n")
            
            i=0
            for Drug in ListSmiles:
                try:
                    Smile=Drug[1] # get the SMILE for one drug
                    Label=Drug[0] # get the label for one drug
                    i+=1
                    print "\n\nDrug No. "+str(i)+" from "+str(len(ListSmiles))+":", Label
                    fDrugTI.write(Label+"\t"+Smile)

                    # if input averages enabled
                    if self.nDrugClass.GetValue()==True:
                        Activity=Drug[2] # get the activity drug
                        fDrugTI.write("\t"+Activity)
                    
                    for col in [1,2,3,4]: # 4 types of weigths
                        drugTIs=CalcDrugMeanProps(Label,Smile,col) # calculate drug mean properties for one type of weigths
                        for TI in drugTIs:
                            fDrugTI.write("\t"+str(TI))
                    fDrugTI.write("\n")
                except:
                    # fDrugTI.write(str(Drug[0])+"ERROR!\n")
                    print "--> !!! ERROR for Drug: "+str(Drug[0])+" ("+str(i)+" of "+str(len(ListSmiles))+")"
     
            fDrugTI.close()


            # CHEMBL calculation
            # Orig SMILE file from GUI -> Intermediate TI simple calc result
            # -> use both to generate -> Intermediate EXTRA input for CHEMBL calc -> output file from GUI
            if self.nCHEMBL.GetValue()==True:
                #-------------------------------------------------------------
                # take the TIs and Add to the input file

                # get original input file
                fInOrig = open(sSMILEfile,"r")
                LinesInOrig = fInOrig.readlines()
                fInOrig.close()

                for i in range(len(LinesInOrig)):
                    CurrLine=LinesInOrig[i]
                    # removing strange characters
                    if CurrLine[-1]=="\n" : CurrLine = CurrLine[:-1]
                    if CurrLine[-1]=="\r" : CurrLine = CurrLine[:-1]
                    LinesInOrig[i]=CurrLine


                # get intermediate output file with the TIs
                fIntermTIs = open(sRezFile,"r")
                LinesIntermTIs = fIntermTIs.readlines()
                fIntermTIs.close()
                
                for i in range(len(LinesIntermTIs)):
                    CurrLine=LinesIntermTIs[i]
                    # removing strange characters
                    if CurrLine[-1]=="\n" : CurrLine = CurrLine[:-1]
                    if CurrLine[-1]=="\r" : CurrLine = CurrLine[:-1]
                    LinesIntermTIs[i]=CurrLine

                #------------------------------------------------------------
                # Write the input file for the CHEMBL calculation
                sRezFile2=sRezFile3[:-4]+"_temp_extra.txt"
                fIntermTI_extra = open(sRezFile2,"w")

                # write original input + intermediate TI simple calc line by line
                for i in range(len(LinesInOrig)):
                    fIntermTI_extra.write(LinesInOrig[i]+"\t"+LinesIntermTIs[i]+"\n")
                fIntermTI_extra.close()              
                
                # -----------------------------------------------------------
                # use the new function for CHEMBL

                ## PARAMETERS
                ##
                ## sFileIn="InputTEST.txt"
                ## sFileOut="OutputTEST.txt"
                ## cutoff=0.0
                ##
                ### lists for FIELD to be use for TI averaging
                ## sFieldNames=["STANDARD_TYPE","ASSAY_CHEMBLID","ASSAY_TYPE","TARGET_CHEMBLID","ORGANISM","TARGET_TYPE","TARGET_MAPPING"]
                ## TIsHeader=["TI1","TI2"] # list with TI names
                ### -----------------------------------------------------------------

                sFileIn=sRezFile2
                sFileOut=sRezFile3
                cutoff=float(self.CutOffZij.GetValue())

                sFieldNames=[] # user list with fields
                if self.nSTANDARD_TYPE.GetValue()==True  :sFieldNames.append("STANDARD_TYPE")
                if self.nASSAY_CHEMBLID.GetValue()==True :sFieldNames.append("ASSAY_CHEMBLID")
                if self.nASSAY_TYPE.GetValue()==True     :sFieldNames.append("ASSAY_TYPE")
                if self.nTARGET_CHEMBLID.GetValue()==True:sFieldNames.append("TARGET_CHEMBLID")
                if self.nORGANISM.GetValue()==True       :sFieldNames.append("ORGANISM")
                if self.nTARGET_TYPE.GetValue()==True    :sFieldNames.append("TARGET_TYPE")
                if self.nTARGET_MAPPING.GetValue()==True :sFieldNames.append("TARGET_MAPPING")
                
                # lists for FIELD to be use for TI averaging                
                TIsHeader=IndLabels # list with TI names

                # CHEMBL calculation
                CHEMBL_AvgIndice(sFileIn,sFileOut,cutoff,sFieldNames,TIsHeader)

                print "\n\n*****************************************************"
                print "\nResult files for CHEMBL calculation:"
                print "\n\n*****************************************************"
                print "\n-> Final output:\n"+sFileOut
                print "\n-> Simple drug output:\n"+sRezFile3[:-4]+"_temp.txt"
                print "\n-> Extra drug output:\n"+sRezFile3[:-4]+"_temp_extra.txt"
                
            #########################################
            # Drug Avgs by input classes

            # if Avg by input classes enabled
            if self.nDrugClass.GetValue()==True:
                print "--> DRUG averages by input classes enabled ..."
                
                sInn=sRezFile #input the simple results
                sOutt=sOutDir+"\\DRUG_ClassAvgs.txt"
                sActiv=self.SMILEfile.GetValue() # classes from input file
                # averages by input file
                # check if there is any class
                InputDrugClasses=PDBHeaderAvg.GetInputClass(sActiv,3)
                if len(InputDrugClasses)!=0:
                    PDBHeaderAvg.DrugAverageTIs4Col(sInn,sOutt,InputDrugClasses,3)
                else:
                    print "!!! Error for DRUG Averages !!! The class column in the input is missing!"

            #########################################
            # Drug Pairs by activity input only

            # if Drug pairs enabled
            if self.nDrugPairs.GetValue()==True:
                print "--> DRUG pairs by activity input enabled ..."
                sOutFile=sOutDir+"\\DRUG_Pairs.txt"
                if self.nDrugClass.GetValue()==True:
                    Offset=3
                else:
                    Offset=2
                TIs2Pairs.DrugTIs2PairsActivFile(self.DrugPairActiv.GetValue(),sRezFile,sOutFile,Offset)
            

        #######################################################################
	# PROTEIN - DRUG PAIRS
	#######################################################################

        # Pairs Prot-Drug enabled
        if self.nProtDrugPairs.GetValue()==True:
            print "\n-> PROTEIN - DRUG PAIRS enabled ..."
            sOutDir,sFileX=os.path.split(self.results.GetValue())

            actFile=self.ProtDrugPairActiv.GetValue() # activities from input file
            inFile1=self.results.GetValue()           # input PROT simple results
            inFile2=self.DrugResults.GetValue()       # input DRUG simple results
            outFile=sOutDir+"\\PROT_DRUG_Pairs.txt"

            # if we found only positive pairs (one class of pairs)
            # create extra negative pairs in other activity file and use it as input
            TypesOfClasses=PDBHeaderAvg.GetInputClass(actFile,3)
            if len(TypesOfClasses)==1:
                print "--> PROTEIN - DRUG pairing using random negative pairs ..."
                m=int(self.PDPairs.GetValue()) # times the linked

                sDir,sfX=os.path.split(self.ProtDrugPairActiv.GetValue())
                
                actFile2=sOutDir+"\\PROT_DRUG_ExtraRndPairs.txt"
                
                # generate random pairs not included in the original input as secondary activity file!
                # if you remove it you can use the existent ExtraPair file!!!
                
                TIs2Pairs.ExtraNegRndPairs(actFile,inFile1,inFile2,actFile2,self.PDPairs.GetValue())
                actFile=actFile2 # change the pair file!
        
            # columns before TIs
            iShift1=1 # for PROTs without input classes
            iShift2=2 # for DRUGS without input classes
            
            # for PROTs with input classes
            if self.AvgTypeInput.GetValue()==True:
                iShift1=2

            # for DRUGS with input classes
            if self.nDrugClass.GetValue()==True:
                iShift2=3
                
            TIs2Pairs.PD_TIs2PairsActivFile(actFile,inFile1,inFile2,outFile,iShift1,iShift2)

        # final message
        if self.nProteins.GetValue()==True or self.nDrugs.GetValue()==True:
            print "\nCalculations finished! Please find the protein output files in the folder: "+sOutDir+"\\"
        else:
            print "\nPlease enable at least one of the types of calculations (PROTEINS, DRUGS) and retry to run!"

        # list the total execution time
        ttn=time.clock()
        print "\nTotal Execution Time: %(ddiffsec).2f min\n\n\n" %\
              {"ddiffsec": (ttn-tt0)/60}

        self.statusbar.SetStatusText("Done!")
        return

# DRUG FUNCTIONS

# ..............................................................
def ProbMatrix(W,CM):
    # calculate the interaction probability matrix for W (weigthed matrix) using CM (contact matrix)
    # input: W (or other matrix such as vdW), CM
    nCa=CM.shape[0]
    PiW=zeros((nCa,nCa)) # moment matrix
    for i in range(nCa):
        sumWij=float(0.0)
        for j in range(nCa):
            sumWij+=W[i][j]*CM[i][j]
        for j in range(nCa):
            if (sumWij!=0) and (CM[i][j]!=0):
                PiW[i][j]=float(W[i][j])/float(sumWij)
    return PiW
# .........................................................
def FileData2Lists(sFile):
    # take the entire AA info (TAB delimited) from a text file and return a list
    # AAFullName[0],A3[1],A1[2],vdWradius[3],NetCh[4],AmberCh[5],Polar_KJ[6]
    # AtContrib2P[7],AtRefr[8],vdWArea[9],hardness_I-A[10],Electrophilicity[11],ElectroMulliken[12]
    AAp=[]
    AAFile = open(sFile,"r")
    for line in AAFile.readlines():
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        AAline=line.split("\t")
        AAp.append(AAline)
    return AAp
# ..............................................................
def AtomInfoCol(atom, AAp,col):
    # return the value that correspond with a column for a specific atom
    r=-1 # bad value
    for A in AAp:
        AA=A[0]
        if atom==AA:
            r=float(A[col])
            return r
    # no found, use the Unknown value
    Unknown=AAp[len(AAp)-1]
    r=float(Unknown[col])
    return r
# ........................................................
def GetWeights(atList,colW):
    DataList=FileData2Lists("ATOMconstants.txt")
    # print DataList
    #print "atList=", atList
    n=len(atList)
    vW=zeros((n))

    for iAt in range(len(atList)):vW[iAt]=AtomInfoCol(atList[iAt],DataList,colW)
    return vW
# ---------------------------------------------------------------------------------------------------------
def WMat(CM,vW):
    # add vector weights to a contact map
    # inputs: connectivity matrix and weight vector
    n=CM.shape[0]
    W=zeros((n,n))  # electro matrix initialization
    for i in range(n):
        for j in range(n):
            if CM[i][j]!=0: # only if nodes are connected
                W[i][j]=vW[j]
    return W # return the weighted matrix

# ..............................................................
def ProbMatrix(W):
    # calculate the interaction probability matrix for W (weigthed matrix)
    # input: W (or other matrix such as vdW)
    n=W.shape[0]
    P=zeros((n,n)) # moment matrix = probability matrix
    for i in range(n):
        sumWij=float(0.0)
        for j in range(n):
            sumWij+=W[i][j]
        for j in range(n):
            if (sumWij!=0) and (W[i][j]!=0):
                P[i][j]=float(W[i][j])/float(sumWij)
    return P

# --------------------------------------------------------------------
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

# ........................................................
def GetAtomTypeVectors(AtomList,CM):
    n=len(AtomList)
    vCsat=zeros((n))    # Saturated carbon vector Csat # 0 or 1
    vCinst=zeros((n))   # Insaturated carbon vector Cinst # 0 or 1
    vHalog=zeros((n))   # halogens
    vHetero=zeros((n))  # heteroatoms
    vHx=zeros((n))      # H linked to the Heteroatoms

    testCsat=sum(CM,1) # sum by line
    # vCsat = if sum by line =4 and it is a C atom, it's Csat=1, else Cinst=1
    for j in range(n):
        atom=AtomList[j]
        if atom=='C': # C atoms
            if testCsat[j]==4: # C saturated
                vCsat[j]=1
            else:              # C insaturated
                vCinst[j]=1
        else:
            if atom=='F' or atom.find("Cl")!=-1 or atom.find("Br")!=-1 or atom=='I': # Halogens
                vHalog[j]=1
            else:
                if atom!='H': # Hetero atoms
                    vHetero[j]=1
    # Hx atoms = H linked to Hetero atoms
    for j in range(n):
        atom=AtomList[j]
        if atom=='H': # if H atom
            for k in range(n):
                if vHetero[k]==1: # if found a Hetero
                    for l in range(n):
                        if CM[k][l]==1 : # if is linked with a Hetero atom
                            vHx[j]=1
   
    return (vCsat,vCinst,vHalog,vHetero,vHx) # return the vectors for each type of atom

# ........................................................
def DrugMeanProps(AtomList,CM,vW,P2n):
    # calculate the mean properties
    # input: list with the atoms, contact matrix, weight vector, powered probability matrices after weight
    n=len(AtomList)
    (PI0,PI1,PI2,PI3,PI4,PI5)=P2n # take the mom matrices to power
    P0j=ones((n))/float(n) # P0j = absolut initial probabilities aprox. 1/n

    # first multiplication from the end of the formula
    Temp0=dot(PI0,vW)
    Temp1=dot(PI1,vW)
    Temp2=dot(PI2,vW)
    Temp3=dot(PI3,vW)
    Temp4=dot(PI4,vW)
    Temp5=dot(PI5,vW)

    MP_tot00=MP_tot01=MP_tot02=MP_tot03=MP_tot04=MP_tot05=float(0)
    MP_Csat00=MP_Csat01=MP_Csat02=MP_Csat03=MP_Csat04=MP_Csat05=float(0)
    MP_Cinst00=MP_Cinst01=MP_Cinst02=MP_Cinst03=MP_Cinst04=MP_Cinst05=float(0)
    MP_Halog00=MP_Halog01=MP_Halog02=MP_Halog03=MP_Halog04=MP_Halog05=float(0)
    MP_Hetero00=MP_Hetero01=MP_Hetero02=MP_Hetero03=MP_Hetero04=MP_Hetero05=float(0)
    MP_Hx00=MP_Hx01=MP_Hx02=MP_Hx03=MP_Hx04=MP_Hx05=float(0)

    vCsat=zeros((n))    # Saturated carbon vector Csat # 0 or 1
    vCinst=zeros((n))   # Insaturated carbon vector Cinst # 0 or 1
    vHalog=zeros((n))   # halogens
    vHetero=zeros((n))  # heteroatoms
    vHx=zeros((n))      # H linked to the Heteroatoms
    # calculate these values
    (vCsat,vCinst,vHalog,vHetero,vHx)=GetAtomTypeVectors(AtomList,CM)
    
    for j in range(n):
        # calculation of each power terms
        t0=P0j[j]*Temp0[j]
        t1=P0j[j]*Temp1[j]
        t2=P0j[j]*Temp2[j]
        t3=P0j[j]*Temp3[j]
        t4=P0j[j]*Temp4[j]
        t5=P0j[j]*Temp5[j]
        # total MP (MP_tot)
        MP_tot00+=t0
        MP_tot01+=t1
        MP_tot02+=t2
        MP_tot03+=t3
        MP_tot04+=t4
        MP_tot05+=t5
        if vCsat[j]==1: # only if Csat
            MP_Csat00+=t0
            MP_Csat01+=t1
            MP_Csat02+=t2
            MP_Csat03+=t3
            MP_Csat04+=t4
            MP_Csat05+=t5
        else:
            if vCinst[j]==1: # only if Cinst
                MP_Cinst00+=t0
                MP_Cinst01+=t1
                MP_Cinst02+=t2
                MP_Cinst03+=t3
                MP_Cinst04+=t4
                MP_Cinst05+=t5
            else:
                if vHalog[j]==1: # only if Halog
                    MP_Halog00+=t0
                    MP_Halog01+=t1
                    MP_Halog02+=t2
                    MP_Halog03+=t3
                    MP_Halog04+=t4
                    MP_Halog05+=t5
                else:
                    if vHetero[j]==1: # only if Halog
                        MP_Hetero00+=t0
                        MP_Hetero01+=t1
                        MP_Hetero02+=t2
                        MP_Hetero03+=t3
                        MP_Hetero04+=t4
                        MP_Hetero05+=t5
                    else:
                        if vHx[j]==1: # only if Halog
                            MP_Hx00+=t0
                            MP_Hx01+=t1
                            MP_Hx02+=t2
                            MP_Hx03+=t3
                            MP_Hx04+=t4
                            MP_Hx05+=t5
                            
    return [float(MP_tot00+MP_tot01+MP_tot02+MP_tot03+MP_tot04+MP_tot05)/float(6),float(MP_Csat00+MP_Csat01+MP_Csat02+MP_Csat03+MP_Csat04+MP_Csat05)/float(6),float(MP_Cinst00+MP_Cinst01+MP_Cinst02+MP_Cinst03+MP_Cinst04+MP_Cinst05)/float(6),float(MP_Halog00+MP_Halog01+MP_Halog02+MP_Halog03+MP_Halog04+MP_Halog05)/float(6),float(MP_Hetero00+MP_Hetero01+MP_Hetero02+MP_Hetero03+MP_Hetero04+MP_Hetero05)/float(6),float(MP_Hx00+MP_Hx01+MP_Hx02+MP_Hx03+MP_Hx04+MP_Hx05)/float(6)]
    # return the mean properties

def SmileToConnectMatrix(sSmile): # reads a smaile and gives the connectiviy matrix
    # write SMILE in a file
    sFileSmile="smile.txt"
    sFileMol="mol.txt"
    print sSmile

    # write Smile file
    fSmile= open(sFileSmile,"w")
    fSmile.write(sSmile)
    fSmile.close()
    
    cmd = 'babel.exe -ismi %s -omol %s' % (sFileSmile, sFileMol)
    os.system(cmd)
    
    fmol = open("mol.txt","r")
    lines = fmol.readlines()
    fmol.close()  
    
    line=lines[3] # take the 4th line with the number of atoms
    sHeader=line.split() # take that header
    atoms=int(line[:3]) # number of atoms
    CM=zeros((atoms,atoms))
    # gets the atom list
    atList=[]
    for a in range(atoms):
        line=lines[4+a]
        detLine=line.split()
        atom_name=detLine[3]
        atList.append(atom_name)

    for a in range(atoms):
        CM[a-1][a-1]=1 # diagonal 1 for a-a pairs (self atom links)
    for r in range(4+atoms,len(lines)-1):
        line=lines[r]
        if line[0] != 'M':
            detLine=line.split()
            a1=int(line[:3])
            a2=int(line[3:6])
            CM[a1-1][a2-1]=1 # -1 because the index in matrix begin with 0, not 1
            CM[a2-1][a1-1]=1 # the invers case
    return (CM,atList) # return the connectivity matrix and the atom list

# ----------------------------------------------------------------------------------------------------------------
def CalcDrugMeanProps(sDrug,sSMILE,col):
    drugTIs=[]
    # CM = contact map
    (CM,atList)=SmileToConnectMatrix(sSMILE)
    n=len(atList) # all the following matrices will have dimension n by n

    # W = weight matrix
    W=zeros((n,n))        # weighted matrix initialization
    vW=zeros((n))         # vector with the weights
    
    vW=GetWeights(atList,col)
    #print "vW=",vW
    
    # Weights matrix
    W=WMat(CM,vW)
    #print "W\n",W
    
    # P = probability/moment matrix
    P=zeros((n,n))      # moment matrix = probability matrix
    P=ProbMatrix(W)            # calculate the interaction probability matrix = moment matrix for W
    #print "P\n",P
    
    # P2n = power of the moment spectral matrices PI[power]=PI0,PI1,PI2,PI3,PI4,PI5
    P2n=(MPower(P,0),MPower(P,1),MPower(P,2),MPower(P,3),MPower(P,4),MPower(P,5))

    # calculate mean properties for drug by atom type
    drugTIs=DrugMeanProps(atList,CM,vW,P2n)
    
    return drugTIs # return the list with the drug TIs

# CHEMBL PART
# ---------------------------------------------------------------------------------------------------------------
def unique(s):
    n = len(s)
    if n == 0:
        return []
    u = {}
    try:
        for x in s:
            u[x] = 1
    except TypeError:
        del u  # move on to the next method
    else:
        return u.keys()
    try:
        t = list(s)
        t.sort()
    except TypeError:
        del t  # move on to the next method
    else:
        assert n > 0
        last = t[0]
        lasti = i = 1
        while i < n:
            if t[i] != last:
                t[lasti] = last = t[i]
                lasti += 1
            i += 1
        return t[:lasti]
    # Brute force is all that's left.
    u = []
    for x in s:
        if x not in u:
            u.append(x)
    return u

# ----------------------------------------------------------
# Calculate the Zij: Z-score for STANDARD TYPE (EC50 ..)
# (parameters: VALUESi,nj,AVGj for a STANDARD TYPE)
# ----------------------------------------------------------

def Calc_Zij(iVal,nj,AVGj):
    # AVGj= sum(VALUESi)/nj
    SDj=math.sqrt(float(((iVal-AVGj)**2))/float(nj+1))
    if SDj>0:
        Zij=float(iVal-float(AVGj))/float(SDj)
    else:
        Zij=float(0)
    return Zij

# -----------------------------------------------------------
# Get unique values for FIELD NAME in a TXT file
# -----------------------------------------------------------
def Read_CHEMBL_UniqueFieldValues(sFileIn,sFieldName):
    FieldValues=[]
    fIN = open(sFileIn,"r")
    lines=fIN.readlines()
    FirstLine=lines[0]
    if FirstLine[-1]=="\n" : FirstLine = FirstLine[:-1] # removing strange characters
    if FirstLine[-1]=="\r" : FirstLine = FirstLine[:-1]
    pFirstLine=FirstLine.split("\t")

    # localize the field in the header
    for i in range(len(pFirstLine)):
        if pFirstLine[i]==sFieldName:
            Col=i # get the column of the header field
    for i in range(1,len(lines)):
        CurrL=lines[i] # current line
        if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
        if CurrL[-1]=="\r" : CurrL = CurrL[:-1]
        CurrR=(CurrL).split("\t")  # split each line by columns tab separated
        item=CurrR[Col]
        FieldValues.append(item) # add to the non-unique list a specific value for the column Col
    fIN.close()
    
    return unique(FieldValues)

### -----------------------------------------------------------
### Get unique values for FIELD NAME in a TXT file
### -----------------------------------------------------------
##def Read_CHEMBL_UniqueFieldValues(sFileIn,sFieldName):
##    FieldValues=[]
##    fIN = open(sFileIn,"r")
##    lines=fIN.readlines()
##    FirstLine=lines[0]
##    if FirstLine[-1]=="\n" : FirstLine = FirstLine[:-1] # removing strange characters
##    if FirstLine[-1]=="\r" : FirstLine = FirstLine[:-1]
##    pFirstLine=FirstLine.split("\t")
##
##    # localize the field in the header
##    for i in range(len(pFirstLine)):
##        if pFirstLine[i]==sFieldName:
##            Col=i # get the column of the header field
##    for i in range(1,len(lines)):
##        CurrL=lines[i] # current line
##        if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
##        if CurrL[-1]=="\r" : CurrL = CurrL[:-1]
##        CurrR=(CurrL).split("\t")  # split each line by columns tab separated
##        item=CurrR[Col]
##        FieldValues.append(item) # add to the non-unique list a specific value for the column Col
##    fIN.close()
##    
##    return unique(FieldValues)

# ----------------------------------------------------------------------------
# Get all values for a column in a TXT file with FIELD NAME header row
#  (sORf = reading strings or floats)
# ----------------------------------------------------------------------------
def Read_CHEMBL_FieldValues(sFileIn,sFieldName,sORf):
    FieldValues=[]
    fIN = open(sFileIn,"r")
    lines=fIN.readlines()
    FirstLine=lines[0]
    if FirstLine[-1]=="\n" : FirstLine = FirstLine[:-1] # removing strange characters
    if FirstLine[-1]=="\r" : FirstLine = FirstLine[:-1]
    pFirstLine=FirstLine.split("\t")

    Col=-1 #default wrong column

    # localize the field in the header
    for i in range(len(pFirstLine)):
        if pFirstLine[i]==sFieldName:
            Col=i # get the column of the header field
    # if no field found, error print and return empty list
    if Col==-1:
        print "ERROR: the field was not found in the header!"
        return []
    
    for i in range(1,len(lines)):
        CurrL=lines[i] # current line
        if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
        if CurrL[-1]=="\r" : CurrL = CurrL[:-1]
        CurrR=(CurrL).split("\t")  # split each line by columns tab separated
        item=CurrR[Col]
        if sORf=="s":
            FieldValues.append(item) # add to the non-unique list a specific value for the column Col
        if sORf=="f":
            if len(item)==0:
                item="0"
            FieldValues.append(float(item)) # add to the non-unique list a specific value for the column Col 
    fIN.close()
    
    return FieldValues

# ----------------------------------------------------------------------------
# Get all values for a column 1 and 2 in a TXT file with FIELD NAME header row
# as Field1(Field2)
# ----------------------------------------------------------------------------
def Read_CHEMBL_2FieldValues(sFileIn,sFieldName1,sFieldName2):
    FieldValues=[]
    fIN = open(sFileIn,"r")
    lines=fIN.readlines()
    FirstLine=lines[0]
    if FirstLine[-1]=="\n" : FirstLine = FirstLine[:-1] # removing strange characters
    if FirstLine[-1]=="\r" : FirstLine = FirstLine[:-1]
    pFirstLine=FirstLine.split("\t")

    Col1=-1 #default wrong column
    Col2=-1 #default wrong column

    # localize the field in the header
    for i in range(len(pFirstLine)):
        if pFirstLine[i]==sFieldName1:
            Col1=i # get the column of the header field 1
        if pFirstLine[i]==sFieldName2:
            Col2=i # get the column of the header field 2
    # if no field found, error print and return empty list
    if Col1==-1 or Col2==-1:
        print "ERROR: the field was not found in the header!"
        return []
    
    for i in range(1,len(lines)):
        CurrL=lines[i] # current line
        if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
        if CurrL[-1]=="\r" : CurrL = CurrL[:-1]
        CurrR=(CurrL).split("\t")  # split each line by columns tab separated
        item1=CurrR[Col1]
        FieldValues.append(CurrR[Col1]+"("+CurrR[Col2]+")") # add to the non-unique list a specific value for the column Col

    fIN.close()
    
    return FieldValues

# -----------------------------------------------------------
# Calculate the averages values for each type of class
#
# inputs: full list with classes, list with values to average, list with unique types of classes
# output: full list with averages corresponding to each class in the FullClassesList
# -----------------------------------------------------------
def AVGsbyClasses(FullClasses,Values,ClassTypes):
    Full_Avgs=[]                    # list with averages depending of the class type
    Avgs=zeros((len(ClassTypes)))   # averages for each class
    Counts=zeros((len(ClassTypes))) # counts of values for each class

    # process each value
    for i in range(len(Values)):
        # seach for each type of class
        for j in range(len(ClassTypes)):
            # if the type of class is identified, add values
            if FullClasses[i]==ClassTypes[j]:
                Avgs[j]+=Values[i]
                Counts[j]+=1

    for j in range(len(Avgs)):
        if Counts[j]!=0:
            Avgs[j]=float(Avgs[j])/float(Counts[j])
        else:
            Avgs[j]=0.0

    # process each value to obtain Full list with averages
    for i in range(len(Values)):
        # seach for each type of class
        for j in range(len(ClassTypes)):
            # if the type of class is identified, add values
            if FullClasses[i]==ClassTypes[j]:
                Full_Avgs.append(Avgs[j])

    return Full_Avgs

# -----------------------------------------------------------
# Calculate the averages values for each type of class using C as filter (only active/positive cases)
#
# inputs: full list with classes, list with values to average, list with unique types of classes, C values (1/0)
# output: full list with averages corresponding to each class in the FullClassesList
# -----------------------------------------------------------
def AVGsbyClassesAndC(FullClasses,Values,ClassTypes,Cs):
    Full_Avgs=[]                    # list with averages depending of the class type
    Avgs=zeros((len(ClassTypes)))   # averages for each class
    Counts=zeros((len(ClassTypes))) # counts of values for each class

    # process each value
    for i in range(len(Values)):
        # seach for each type of class
        for j in range(len(ClassTypes)):
            # if the type of class is identified and C=1 (active), add values
            if FullClasses[i]==ClassTypes[j] and Cs[i]==1:
                Avgs[j]+=Values[i]
                Counts[j]+=1
    for j in range(len(Avgs)):
        if Counts[j]!=0:
            Avgs[j]=float(Avgs[j])/float(Counts[j])
        else:
            Avgs[j]=0.0
            
    # process each value to obtain Full list with averages
    for i in range(len(Values)):
        # seach for each type of class
        for j in range(len(ClassTypes)):
            # if the type of class is identified, add values
            if FullClasses[i]==ClassTypes[j]:
                Full_Avgs.append(Avgs[j])

    return Full_Avgs


# -----------------------------------------------------------
# Counts for each type of class
#
# inputs: full list with classes, list with values to average, list with unique types of classes
# output: counts for each unique class in the FullClassesList
# -----------------------------------------------------------
def COUNTbyClasses(FullClasses,Values,ClassTypes):
    Counts=zeros((len(ClassTypes))) # counts of values for each class

    # process each value
    for i in range(len(Values)):
        # seach for each type of class
        for j in range(len(ClassTypes)):
            # if the type of class is identified, add values
            if FullClasses[i]==ClassTypes[j]:
                Counts[j]+=1

    return list(Counts)

def GetDrugID_SMILEPairs(sFileIn):
    DrugIDs_SMILEs=[]
    # get the list with drug IDs
    sORf="s"
    DrugIDs=Read_CHEMBL_FieldValues(sFileIn,"CMPD_CHEMBLID",sORf)   # get all values from CMPD_CHEMBLID
    SMILEs=Read_CHEMBL_FieldValues(sFileIn,"CANONICAL_SMILES",sORf) # get all values from CANONICAL_SMILES

    # if the input info is not complete, return empty list
    if len(DrugIDs) != len(SMILEs):
        print " >>> Input file error! There are drugs without IDs or SMILES. Please check the input and retry."
        return []

    # create the pairs
    for i in range(len(DrugIDs)):
        DrugIDs_SMILEs.append((DrugIDs[i],SMILEs[i]))

    return DrugIDs_SMILEs

###############################################################
# CHEMBL Avg Indices

#-------------------------------------------------------------------
# PARAMETERS
#
# sFileIn="InputTEST.txt"
# sFileOut="OutputTEST.txt"
# cutoff=0.0
#
# lists for FIELD to be use for TI averaging
# sFieldNames=["STANDARD_TYPE","ASSAY_CHEMBLID","ASSAY_TYPE","TARGET_CHEMBLID","ORGANISM","TARGET_TYPE","TARGET_MAPPING"]
# TIsHeader=["TI1","TI2"] # list with TI names
# -----------------------------------------------------------------
def CHEMBL_AvgIndice(sFileIn,sFileOut,cutoff,sFieldNames,TIsHeader):
    # read TI DRUG SIMPLE calculation file
    fIN = open(sFileIn,"r")
    lines=fIN.readlines()
    fIN.close()

    # read HEADER
    Header=lines[0]   # header text
    if Header[-1]=="\n" : Header = Header[:-1] # removing strange characters
    if Header[-1]=="\r" : Header = Header[:-1]


    # -------------------------------------------------------------------------------
    # Verify STANDARD TYPE Z-score based on STANDARD VALUES and Activity Class (C)
    # -------------------------------------------------------------------------------

    # Merge STANDARD_TYPE + STANDARD_UNITS into one class

    # get full list with Classes
    FullClasses=Read_CHEMBL_2FieldValues(sFileIn,"STANDARD_TYPE","STANDARD_UNITS")
    
    # get the unique values of FIELD
    ClassTypes=unique(FullClasses)

    # get the list with values to be averaged
    sORf="f" # get float values
    Values= Read_CHEMBL_FieldValues(sFileIn,"STANDARD_VALUE",sORf)

    # AVGs and Counts for each class type
    AVGs=AVGsbyClasses(FullClasses,Values,ClassTypes)
    Counts=COUNTbyClasses(FullClasses,Values,ClassTypes)

    Zijs=[] # list with all the Zij values
    Cs=[]   # list with all the C values
    # calculate Zij and C
    for i in range(len(AVGs)):
        # verify the Class to calculate Zij
        for j in range(len(ClassTypes)):
            if FullClasses[i]==ClassTypes[j]:
                # calculate Zij
                Zij=Calc_Zij(Values[i],Counts[j],AVGs[i])
                # adding Zij
                Zijs.append(Zij)
                # calculate ACTIVITY CLASS using Zij and cutoff
                if Zij>=cutoff:
                    C=1
                else:
                    C=0
                # adding C
                Cs.append(C)

    nC1=Cs.count(1) # count 1s
    nC0=Cs.count(0) # count 0s

    pC1=100.0*float(nC1)/float(nC1+nC0)
    pC0=100.0*float(nC0)/float(nC1+nC0)

    print "\n\nNOTE: There are %d (%.2f %%) positive activity values (C=1) and %d(%.2f %%) negative activity values (C=0)." \
          % (nC1,pC1,nC0,pC0)

    sDecision = raw_input('\n--> Do you want to continue the calculation? (y/n)')

    if sDecision == "n" or sDecision == "N": exit(0) # EXIT the function (USE RETURN!)

    print "--> Creating the outputs, please wait ..."

    # Output HEADER

    # Zij, C, P_curate + original header + TIs + Avgs + Deviations
    Header="Z-score"+"\t"+"Activity_Class"+"\t"+"P_curate"+"\t"+Header

    # for each FIELD to be use for TI averaging write AVG
    for sFieldName in sFieldNames:
        for TI in TIsHeader:
            Header=Header+"\t"+"Avg_"+TI+"["+sFieldName+"]"

    # for each FIELD to be use for TI averaging write Deviation header
    for sFieldName in sFieldNames:
        for TI in TIsHeader:
            Header=Header+"\t"+"Dev_"+TI+"["+sFieldName+"]"
            
    Header=Header+"\n"
    lines[0]=Header

    # -------------------------------------------------------------------------------
    # Write STANDARD TYPE Z-score and Activity Class (C)
    # -------------------------------------------------------------------------------

    # get full list with CURATED_BY field and calculate P_curation values (0.5, 0.75,1.0)
    sORf="s"
    Curations=Read_CHEMBL_FieldValues(sFileIn,"CURATED_BY",sORf)# list with all Curation values
    P_curations=[] # list with all P_curation values

    # calculate P_curation: Autocuration=0.5, Intermediate=0.75, Expert =1.0
    for iCur in Curations:
        if iCur=="Autocuration":P_curations.append(0.5)
        if iCur=="Intermediate":P_curations.append(0.75)
        if iCur=="Expert":P_curations.append(1.0)

    # Zij & C
    for i in range(len(AVGs)):
        CurrL=lines[i+1] # current line
        if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
        if CurrL[-1]=="\r" : CurrL = CurrL[:-1]

        # adding Zij, C, P_curate at the begining of output lines
        lines[i+1]=str(Zijs[i])+"\t"+str(Cs[i])+"\t"+str(P_curations[i])+"\t"+CurrL # CurrL includes TIS!!!! from the input file!

    # for each FIELD to be use for TI averaging

    # AVGs
    for sFieldName in sFieldNames:
        # get the unique values of FIELD
        ClassTypes=Read_CHEMBL_UniqueFieldValues(sFileIn,sFieldName)

        # get full list with Classes
        sORf="s"
        FullClasses=Read_CHEMBL_FieldValues(sFileIn,sFieldName,sORf)

        # process all the TIs
        for TI in TIsHeader:
            # get the list with values to be averaged
            sORf="f"
            Values= Read_CHEMBL_FieldValues(sFileIn,TI,sORf)

            # average list corrsponding with the class type
            AVGs=AVGsbyClassesAndC(FullClasses,Values,ClassTypes,Cs)

            # adding new AVGs to LINES data
            for i in range(len(AVGs)):
                CurrL=lines[i+1] # current line
                if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
                if CurrL[-1]=="\r" : CurrL = CurrL[:-1]
                # adding information
                lines[i+1]=CurrL+"\t"+str(AVGs[i])+"\n"

    # DEVIATIONS
    for sFieldName in sFieldNames:
        # get the unique values of FIELD
        ClassTypes=Read_CHEMBL_UniqueFieldValues(sFileIn,sFieldName)

        # get full list with Classes
        sORf="s"
        FullClasses=Read_CHEMBL_FieldValues(sFileIn,sFieldName,sORf)

        # process all the TIs
        for TI in TIsHeader:
            # get the list with values to be averaged
            sORf="f"
            Values= Read_CHEMBL_FieldValues(sFileIn,TI,sORf)

            # average list corrsponding with the class type
            AVGs=AVGsbyClassesAndC(FullClasses,Values,ClassTypes,Cs)

            # adding new Deviations to LINES data
            for i in range(len(AVGs)):
                CurrL=lines[i+1] # current line
                if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
                if CurrL[-1]=="\r" : CurrL = CurrL[:-1]
                # adding information
                lines[i+1]=CurrL+"\t"+str((Values[i]-AVGs[i])*P_curations[i])+"\n"

    # open to write the output file
    fOUT=open(sFileOut,"w")
    # write all LINES into the output file
    for line in lines:
        fOUT.write(line)
    fOUT.close()

    print "\nDone!"

    return

###############################################################
# Calculate the deviations for new drugs against
# precalculated unique combinations of field averages
# -------------------------------------------------------------------
# PARAMETERS
#
# sFileInSmiles="SMILEs.txt"
# sFileInAvg="DRUG_NewDrugsPreCalcAverages.txt"
# sFileOut="DRUG_New_Deviations.txt"

# sFieldNames=["STANDARD_TYPE","ASSAY_CHEMBLID","ASSAY_TYPE","TARGET_CHEMBLID","ORGANISM","TARGET_TYPE","TARGET_MAPPING"]
# TIsHeader=["TI1","TI2"] # list with TI names
# -----------------------------------------------------------------
def NewDrugDevs_PreCalcAvgs(sFileInSmiles,sFileInAvg,sFileOut,sFieldNames,TIsHeader):
    # --------------------------------------------------------------
    # Read INPUT file with NEW DRUGS (DrugName[TAB]SMILES)
    # --> input file from GUI => ListSmiles

    print "*******************************************"
    print " > New Drug Deviations based on pre-calculated averages (CHEMBL database) ..."
    
    # GET pairs DrugName-SMILE
    input_file = open(sFileInSmiles,"r")
    lines = input_file.readlines()
    input_file.close()

    # from the SMILE file: list of pairs : Label - SMILE for each drug
    ListSmiles=[]
    for line in lines :
        # removing strange characters
        if line[-1]=="\n" : line = line[:-1]
        if line[-1]=="\r" : line = line[:-1]
        DataLine=line.split("\t")
        if len(DataLine)>1:
            ListSmiles.append((DataLine[0],DataLine[1]))
            
    # --------------------------------------------------------------
    # Read PRE-calculated result file
    # --> Get the unique header field combinations
    # (search for them in the pre-calculated file)

    # Read Pre-Calc AVG file
    input_file = open(sFileInAvg,"r")
    lines = input_file.readlines()
    input_file.close()

    HeaderLine=lines[0]
    # removing strange characters
    if HeaderLine[-1]=="\n" : HeaderLine = HeaderLine[:-1]
    if HeaderLine[-1]=="\r" : HeaderLine = HeaderLine[:-1]
    Header=HeaderLine.split("\t") # get the header list

    # create Avg headers
    IndLabels=[]
    for Field in sFieldNames:
        for TI in TIsHeader:
            IndLabels.append("Avg_"+TI+"["+Field+"]") # Avg_ProprietyAtomType[Field]

    # check each header
    HeaderFlags=zeros((len(IndLabels)))
    for j in range(len(Header)):
        # check each field if exist
        for i in range(len(IndLabels)):
            # if found the field in the header => it was calculated = there are all the TIs for it
            if Header[j]==IndLabels[i]:
                HeaderFlags[i]=1.
                
    # create the founded AVG headers
    NewTIsHeader=[]
    for i in range(len(IndLabels)):
        if HeaderFlags[i]==1.:
            NewTIsHeader.append(IndLabels[i])

    # create new Field list
    NewFieldNameFlags=zeros((len(sFieldNames)))
    for i in range(len(NewTIsHeader)):
        for j in range(len(sFieldNames)):
            # if found the field in the AvgTI name
            if sFieldNames[j] in NewTIsHeader[i]:
                NewFieldNameFlags[j]=1.
                
    NewFieldNames=[]
    for i in range(len(sFieldNames)):
        if NewFieldNameFlags[i]==1.:
            NewFieldNames.append(sFieldNames[i])

    # if no HEADER found, exit
    if len(NewTIsHeader)==0:
        print "\n\n --> ERROR: The Pre-calculated Average File has no Average Header. Please check this input file and retry."
        return

    # ---------------------------------------------------------------------
    # Get the unique combination of Fields with the correspondent averages

    sMixData=[] # list with all the brute info to make unique

    # get column with first FIELD
    # start a header of MixData
    sMixData.append(NewFieldNames[0])
    Field_Col=Read_CHEMBL_FieldValues(sFileInAvg,NewFieldNames[0],"s")
    # add the field to sMixData using TAB
    for i in range(len(Field_Col)): sMixData.append(Field_Col[i])
    
    # get info from the other FIELD from Pre-calc AVGs
    for i in range(1,len(NewFieldNames)):
        # header of MixData
        sMixData[0]=sMixData[0]+"\t"+NewFieldNames[i]

        # get column with header FIELD
        Field_Col=Read_CHEMBL_FieldValues(sFileInAvg,NewFieldNames[i],"s")

        # add the field to sMixData using TAB
        for j in range(len(Field_Col)): sMixData[j+1]=sMixData[j+1]+"\t"+Field_Col[j]

    # add averages TIs
    for i in range(len(NewTIsHeader)):
        # header of MixData
        sMixData[0]=sMixData[0]+"\t"+NewTIsHeader[i]

        
        # get column with header FIELD
        Field_Col=Read_CHEMBL_FieldValues(sFileInAvg,NewTIsHeader[i],"s")
        # add the field to sMixData using TAB
        for j in range(len(Field_Col)): sMixData[j+1]=sMixData[j+1]+"\t"+Field_Col[j]

    # remove header of MixData
    NoHeaderData=sMixData[1:]

    # header of Mix Data
    HeaderMixData=sMixData[0]
    # get the unique List of combination for Fields and averages to be used for the new drug
    UMixData=unique(NoHeaderData)

    # create a list from Header and Unique values for MixData
    AvgCombination=[]
    AvgCombination.append(HeaderMixData.split("\t"))
    for UData in UMixData:AvgCombination.append(UData.split("\t"))

    AvgCombinationH=AvgCombination[0] # header list
    AvgCombinationD=AvgCombination[1:] # data list
    
    # ---------------------------------------------------------------------
    # Calculate TIs for the new NEW DRUGS (DrugName[TAB]SMILES)
    # --> simple output result (DrugName, SMILES, TIs)

    sOut_NewTIs=sFileInSmiles[:-4]+"_NewTIs.txt"
    fDrugTI=open(sOut_NewTIs,"w") # open to write the output file

    # write the header

    # take the data in lists
    NewTIList=[]

    HLine=[] # each line 
    sOut="DrugID\tSMILES"
    HLine.append(sOut)
    
    # add TI labels
    for TI in TIsHeader:
        sOut+="\t"+TI
        HLine.append(TI)
    fDrugTI.write(sOut+"\n")

    # add header in the list
    NewTIList.append(HLine)

    # Calculate TIs for the New Drugs
    for Drug in ListSmiles:
        HLine=[]
        try:
            Smile=Drug[1] # get the SMILE for one drug
            Label=Drug[0] # get the label for one drug
            fDrugTI.write(Label+"\t"+Smile)
            HLine.append(Label+"\t"+Smile)
            
            for col in [1,2,3,4]: # 4 types of weigths
                drugTIs=CalcDrugMeanProps(Label,Smile,col) # calculate drug mean properties for one type of weigths
                for TI in drugTIs:
                    fDrugTI.write("\t"+str(TI))
                    HLine.append(TI)
            fDrugTI.write("\n")
            # add list of TI to the big list
            NewTIList.append(HLine)
        except:
            print "\n-->> ERROR: Drug "+Label+" is excluded from the results! Please check the SMILE code.\n"
            fDrugTI.write("\n")
    fDrugTI.close()

    print "\n ---> Created "+sFileInSmiles[:-4]+"_NewTIs.txt file = only the new calculated TIs."

    # ----------------------------
    # DEVIATIONS

    # AvgCombinationH = header list for AVGs = Field1,...,Fieldn,AvgHeader1,...,AvgHeadern
    # AvgCombinationD = data list if lists   =[[FieldVal1,...,FieldValn,AvgVal1,...,AvgValn]]

    NewTIListH=NewTIList[0]  # header list for TIs
    NewTIListD=NewTIList[1:] # data list for TIs

    Results=NewTIList # list of lists with results, initialy with the TI results only (to add the Dev)

    # create output file to write
    fFinalOutput=open(sFileOut,"w") # open to write the output file

    # write Header
    sOut=""
    for TIh in NewTIListH:
        sOut=sOut+"\t"+str(TIh)
    for Hd in AvgCombinationH:
        sOut=sOut+"\t"+str(Hd)
    for i in range(1,len(NewTIListH)): # -1 = minus the first column with Drug\tSMILE
        for j in range(len(AvgCombinationH)): # no of header in Avg - no of CEHMBL fields
            if NewTIListH[i] in AvgCombinationH[j]: # if TI exists inside the Avg label
                sOut=sOut+"\t"+"Dev_"+AvgCombinationH[j]
    fFinalOutput.write(sOut[1:]+"\n") 

    # add the deviation if TI found in AvgTICHMEBLfield

    # for each new drug TI calculation
    for w in range(len(NewTIListD)):
        # for each unique combination of Avg
        for k in range(len(AvgCombinationD)):
            sOut=""
            # add info from TI res
            CurrDataTI=NewTIListD[w]
            for cTI in CurrDataTI:
                sOut=sOut+"\t"+str(cTI)
                
            # add info from AVGs
            CurrDataAVG=AvgCombinationD[k]
            for cAvg in CurrDataAVG:
                sOut=sOut+"\t"+str(cAvg)

            # DEVIATIONS
            for i in range(1,len(NewTIListH)): # for each TI header search for AvgTI
                for j in range(len(AvgCombinationH)):
                    if NewTIListH[i] in AvgCombinationH[j]: # if TI exists inside the Avg label
                        sOut=sOut+"\t"+str((float(CurrDataTI[i])-float(CurrDataAVG[j])))
            sOut=sOut+"\n"
            fFinalOutput.write(sOut[1:])
    fFinalOutput.close()

    print "\n ---> Created the final output file: "+sFileOut+"."
    print "\nDone!"
    
    return
 
###############################################################################################
# main GUI application
###############################################################################################

class MyApp(wx.App):
    def OnInit(self):
        frame = MyMenu(None, -1, 'MInD-Prot: Markov Indices for Drugs and Proteins (ver. 3)')
        frame.SetIcon(wx.Icon(orig_dir+'faviconMC.ico', wx.BITMAP_TYPE_ICO))
        frame.Centre()
        frame.Show(True)
        return True

app = MyApp(0)
app.MainLoop()
