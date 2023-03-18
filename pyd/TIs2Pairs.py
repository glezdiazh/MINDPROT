# transform the output file PDBchains - TIs in Pair PDBchain interactions - Mixed TIs

from numpy import *
from random import *

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

# ---------------------------------------------------------------------------------------------------------------
def GetUniqueListFromColumn(sInFile,Col):
    # list with the unique classes excepted the first line (no header line processed)
    ClassList=[]
    
    fIn=open(sInFile,"r") #open file to read
    lines=fIn.readlines()
    
    l=len(lines)                 # no. of lines in the input file
    FirstLine=lines[0]           # get the header
    c=len(FirstLine.split("\t")) # numbers of header items = columns
    # only if there are more columns
    if c>1:
        for i in range(1,l): # process each line in the file
            CurrL=lines[i] # current line
            if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
            if CurrL[-1]=="\r" : CurrL = CurrL[:-1]
            CurrR=(CurrL).split("\t")  # split each line by columns tab separated
            item=CurrR[Col-1]
            ClassList.append(item) # add to the non-unique list a specific value for the column Col

    return unique(ClassList) # return the unique classes from the input

# ---------------------------------------------------------------------------------------------------------------
def PDBRndPairsList(PairUnique,n):
    # use a unique list of PDB to create random pairs (removed invers pairs)
    # not the same PDB
    
    RndListPairs={} #new dictionary that include the original pairs (positive) and the new pairs (negative)
    ProtItemList=[] # list with unique PDBchain items from the PairUnique list

    # add to dictionary the positive pairs
    for Pairx in PairUnique:
        # add all the positive chains in order to check the new ones
        # at the and all positive will be removed
        RndListPairs[Pairx[0]+"&"+Pairx[1]]= "1"
        # add all the PDBchain founded in the positive list
        ProtItemList.append(Pairx[0])
        ProtItemList.append(Pairx[1])

    # make the unique item list of PDBchains
    ProtItemList=unique(ProtItemList)
    ProtItemList.sort()
    
    # search for the negative pairs
    iP=len(ProtItemList)
    
    # if pairs*pairs < n => take all the pairs
    if n>iP*iP:
        print " ... Generate all possible protein pairs ..."
        for P1 in ProtItemList:
            for P2 in ProtItemList:
                # if there is not in the pairs and the invers
                if (P1+"&"+P2 not in RndListPairs) and (P2+"&"+P1 not in RndListPairs):
                    # add negative pairs to the dictionaty
                    RndListPairs[P1+"&"+P2]= "0"
    else:
        # if not, take a number n of random unique pairs
        print " ... Generate random protein pairs ..."
        iPairs=0
        fLimit=0 # limit number when it is not possible to generate more pairs
        # if pair no. is less then we need, find a pair
        while iPairs<n:
            # if the limit was reach stop the searching (if no pair found n times = brutal condition)
            if fLimit==n:
                break
            # if the random data cant be generated more
            P1=choice(ProtItemList)
            P2=choice(ProtItemList)
            # if there is not in the pairs
            if (P1+"&"+P2 not in RndListPairs) and (P2+"&"+P1 not in RndListPairs):
                iPairs+=1
                # add negative pairs to the dictionaty
                RndListPairs[P1+"&"+P2]= "0"
                fLimit=0
            else:
                # increase the limit number
                fLimit+=1
                
    # generate the final list of negative pairs of PDBchains
    RndListPairsFin=[]
    for k, v in RndListPairs.iteritems():
        # get only the negative pairs
        if v=="0":
            RndListPairsFin.append((k[:5],k[-5:]))

    return RndListPairsFin # return only n pairs

# ---------------------------------------------------------------------------------------------------------------
def TIs2PairFile(inFile,outFile,m,iShift): # m=multiple of linked
    # iShift = no of columns in input files to skip before TIs
    # list to take the file containt
    PDBChainTIs=[]
    f = open(inFile,"r")
    lines=f.readlines()
    lines.sort()
    
    for line in lines:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        PDBChainTIs.append(line.split("\t"))

    l=len(PDBChainTIs)  # number of lines
    r = PDBChainTIs[0]  # record = one line of data
    c = len(r)          # number of columns

    o = open(outFile,"w") # open the pair file (columns separated by space)
    s = "\t"              # column separator
    
    # write the header
    sOut=""  # output variable to write in Pair file
    # col 1 = L = link
    sOut+="Linked"
    # col 2-3 = PDBChain1 PDBChain2
    sOut+=s+r[0]+"(1)"+s+r[0]+"(2)"
    # cols TI_1
    for i in range(iShift,c):
        sOut+=s+r[i]+"(1)"
    # cols TI_2
    for i in range(iShift,c):
        sOut+=s+r[i]+"(2)"
    # cols Diffs: abs(TI_1-TI_2)
    for i in range(iShift,c):
        sOut+=s+"Diff_"+r[i]   
    # cols Prod: TI_1*TI_2
    for i in range(iShift,c):
        sOut+=s+"Prod_"+r[i] 
    # cols Avg: (TI_1+TI_2)/2
    for i in range(iShift,c):
        sOut+=s+"Avg_"+r[i]
    sOut+="\n"
    # end of heading
    
    o.write(sOut) # write the header
    sOut=""       # initialize the output var

    # LINKED chains (the same protein) in PDBChainTIs
    # (from record 1 up, first line 0 is the header)

    PairUnique=[] # unique pairs to check the nonLinked pairs
    
    Linked=0 # number of linked pairs
    for i in range(1,l): # P1 list of PDB chain
        r1=PDBChainTIs[i]
        P1_PDBch=r1[0]      # PDBchain
        P1_PDB=P1_PDBch[:4] # PDB
        for j in range(1,l): # P2 list of the PDB chain
            r2=PDBChainTIs[j]
            P2_PDBch=r2[0]      # PDBchain
            P2_PDB=P2_PDBch[:4] # PDB
            # only positive LINKED = the same PDB
            if (P1_PDB==P2_PDB) and (P1_PDBch!=P2_PDBch):
                u=1
                # check the list of unique pairs to avoid x-y and y-x
                if len(PairUnique)!=0:
                    for pair in PairUnique:
                        P1=pair[0]
                        P2=pair[1]
                        if (P1_PDBch==P1 and P2_PDBch==P2) or (P1_PDBch==P2 and P2_PDBch==P1):
                            u=0 # to not take this pair again
                    
                if u==1: # if the pair is not repeated
                    PairUnique.append((P1_PDBch,P2_PDBch))
                    Linked+=1
                    sOut+="1"+s+P1_PDBch+s+P2_PDBch # write PDBchains of the pair
                    # cols TI_1
                    for i1 in range(iShift,c):
                        sOut+=s+r1[i1]
                    # cols TI_2
                    for i2 in range(iShift,c):
                        sOut+=s+r2[i2]
                    # cols Diffs: abs(TI_1-TI_2)
                    for i12 in range(iShift,c):
                        sOut+=s+str(abs(float(r1[i12])-float(r2[i12])))  
                    # cols Prod: TI_1*TI_2
                    for i12 in range(iShift,c):
                        sOut+=s+str(float(r1[i12])*float(r2[i12])) 
                    # cols Avg: (TI_1+TI_2)/2
                    for i12 in range(iShift,c):
                        sOut+=s+str((float(r1[i12])+float(r2[i12]))/float(2))
                    sOut+="\n"
                    o.write(sOut) # write the pair
                    sOut=""       # initialize the output var

    # non-LINKED chains (the same protein) in PDBChainTIs
    # (from record 1 up, first line 0 is the header)

    # get the unique list of PDBchains
    # ItemList=GetUniqueListFromColumn(inFile,1)
    # get the random negative pairs
    RndListPairs=PDBRndPairsList(PairUnique,Linked*m)
    
    # introduce each negative pair
    for P1,P2 in RndListPairs:
        # list for TIs for each pair
        TIs1=[]
        TIs2=[]
        for i1 in range(1,l): # each line of the input
            r1=PDBChainTIs[i1]
            PDBch1=r1[0]      
            for i2 in range(1,l): # each line of the input
                r2=PDBChainTIs[i2]
                PDBch2=r2[0]
                if P1==PDBch1 and P2==PDBch2:
                    for j in range(iShift,c):
                        TIs1.append(str(r1[j]))
                    for jj in range(iShift,c):
                        TIs2.append(str(r2[jj]))
                    
                    # if we found both pairs
                    if len(TIs1)!=0 and len(TIs2)!=0:
                        sOut+="0"+s+P1+s+P2 # write PDBchains of the pair
                        # cols TI_1
                        for i1 in range(len(TIs1)):
                            sOut+=s+TIs1[i1]
                        # cols TI_2
                        for i2 in range(len(TIs2)):
                            sOut+=s+TIs2[i2]
                        # cols Diffs: abs(TI_1-TI_2)
                        for i12 in range(len(TIs1)):
                            sOut+=s+str(abs(float(TIs1[i12])-float(TIs2[i12])))
                        # cols Prod: TI_1*TI_2
                        for i12 in range(len(TIs1)):
                            sOut+=s+str(float(TIs1[i12])*float(TIs2[i12]))
                        # cols Avg: (TI_1+TI_2)/2
                        for i12 in range(len(TIs1)):
                            sOut+=s+str((float(TIs1[i12])+float(TIs2[i12]))/float(2))
                        sOut+="\n"
                        o.write(sOut) # write the pair
                        sOut=""       # initialize the output var
                        continue #break??
    o.close()
  
    return

# ---------------------------------------------------------------------------------------------------------------
def TIs2PairsActivFile(actFile,inFile,outFile,iShift):
    # create a file with pairs of PDBchains using an activity pair file
    
    # iShift = no of columns in input files to skip before TIs
    # actFile = activity file Pair1 Pair2 Activ [tab]
    # inFile = simple result file with TIs

    PairsActiv=[]  # list with pairs and activity from activity file

    # get the pairs from activity file
    fAct = open(actFile,"r")
    liness=fAct.readlines()
    for line in liness:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        PairsActiv.append(line.split("\t"))    
    
    # process the simple result file with TIs
    PDBChainTIs=[] # PDBchains and the TIs for all the simple result file 
    f = open(inFile,"r")
    lines=f.readlines()
    
    for line in lines:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        PDBChainTIs.append(line.split("\t"))

    l=len(PDBChainTIs)  # number of lines
    r = PDBChainTIs[0]  # record = one line of data
    c = len(r)          # number of columns

    o = open(outFile,"w") # open the pair file (columns separated by space)
    s = "\t"              # column separator
    
    #----------------------
    # HEADER
    sOut=""  # output variable to write in Pair file
    # col 1 = L = link
    sOut+="Activity"
    # col 2-3 = PDBChain1 PDBChain2
    sOut+=s+r[0]+"(1)"+s+r[0]+"(2)"
    # cols TI_1
    for i in range(iShift,c):
        sOut+=s+r[i]+"(1)"
    # cols TI_2
    for i in range(iShift,c):
        sOut+=s+r[i]+"(2)"
    # cols Diffs: abs(TI_1-TI_2)
    for i in range(iShift,c):
        sOut+=s+"Diff_"+r[i]   
    # cols Prod: TI_1*TI_2
    for i in range(iShift,c):
        sOut+=s+"Prod_"+r[i] 
    # cols Avg: (TI_1+TI_2)/2
    for i in range(iShift,c):
        sOut+=s+"Avg_"+r[i]
    sOut+="\n"
    # end of heading
    
    o.write(sOut) # write the header
    sOut=""       # initialize the output var

    # OUTPUT pairs
    # for each line in the activity file
    for pairs in PairsActiv:
        P1=pairs[0]
        P2=pairs[1]
        Act=pairs[2]
        # check the TIs
        for i in range(l): # TIs from simple result file     
            r1=PDBChainTIs[i]
            P1_PDBch=r1[0]      # PDBchain
            P1_PDB=P1_PDBch[:4] # PDB
            for j in range(l):
                r2=PDBChainTIs[j]
                P2_PDBch=r2[0]      # PDBchain
                P2_PDB=P2_PDBch[:4] # PDB
                # only positive LINKED = the same PDB
                if (P1==P1_PDBch) and (P2==P2_PDBch):
                    sOut+=str(Act)+s+P1+s+P2 # write PDBchains of the pair
                    # cols TI_1
                    for i1 in range(iShift,c):
                        sOut+=s+r1[i1]
                    # cols TI_2
                    for i2 in range(iShift,c):
                        sOut+=s+r2[i2]
                    # cols Diffs: abs(TI_1-TI_2)
                    for i12 in range(iShift,c):
                        sOut+=s+str(abs(float(r1[i12])-float(r2[i12])))  
                    # cols Prod: TI_1*TI_2
                    for i12 in range(iShift,c):
                        sOut+=s+str(float(r1[i12])*float(r2[i12])) 
                    # cols Avg: (TI_1+TI_2)/2
                    for i12 in range(iShift,c):
                        sOut+=s+str((float(r1[i12])+float(r2[i12]))/float(2))
                    sOut+="\n"
                    o.write(sOut) # write the pair
                    sOut=""       # initialize the output var
                    continue #break??
    o.close()
    return

# ---------------------------------------------------------------------------------------------------------------
def DrugTIs2PairsActivFile(actFile,inFile,outFile,iShift):
    # create a file with pairs of drugs using an activity pair file
    
    # iShift = no of columns in input files to skip before TIs
    # actFile = activity file Drug_name1 Drug_name2 Activ [tab]
    # inFile = simple result file with TIs

    PairsActiv=[]  # list with pairs and activity from activity file

    # get the pairs from activity file
    fAct = open(actFile,"r")
    liness=fAct.readlines()
    for line in liness:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        PairsActiv.append(line.split("\t"))    
    
    # process the simple result file with TIs
    DrugTIs=[] # PDBchains and the TIs for all the simple result file 
    f = open(inFile,"r")
    lines=f.readlines()
    
    for line in lines:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        DrugTIs.append(line.split("\t"))

    l=len(DrugTIs)  # number of lines
    r = DrugTIs[0]  # record = one line of data
    c = len(r)          # number of columns

    o = open(outFile,"w") # open the pair file (columns separated by space)
    s = "\t"              # column separator
    
    #----------------------
    # HEADER
    sOut=""  # output variable to write in Pair file
    # col 1 = L = link
    sOut+="Activity"
    # col 2-3 = PDBChain1 PDBChain2
    sOut+=s+r[0]+"(1)"+s+r[0]+"(2)"
    # cols TI_1
    for i in range(iShift,c):
        sOut+=s+r[i]+"(1)"
    # cols TI_2
    for i in range(iShift,c):
        sOut+=s+r[i]+"(2)"
    # cols Diffs: abs(TI_1-TI_2)
    for i in range(iShift,c):
        sOut+=s+"Diff_"+r[i]   
    # cols Prod: TI_1*TI_2
    for i in range(iShift,c):
        sOut+=s+"Prod_"+r[i] 
    # cols Avg: (TI_1+TI_2)/2
    for i in range(iShift,c):
        sOut+=s+"Avg_"+r[i]
    sOut+="\n"
    # end of heading
    
    o.write(sOut) # write the header
    sOut=""       # initialize the output var

    # OUTPUT pairs
    # for each line in the activity file
    for pairs in PairsActiv:
        P1=pairs[0]
        P2=pairs[1]
        Act=pairs[2]
        # check the TIs
        for i in range(l): # TIs from simple result file     
            r1=DrugTIs[i]
            Drug1=r1[0]      # Drug pair 1
            for j in range(l):
                r2=DrugTIs[j]
                Drug2=r2[0]  # Drug pair 2
                # only pairs
                if (P1==Drug1) and (P2==Drug2):
                    sOut+=str(Act)+s+P1+s+P2 # write PDBchains of the pair
                    # cols TI_1
                    for i1 in range(iShift,c):
                        sOut+=s+r1[i1]
                    # cols TI_2
                    for i2 in range(iShift,c):
                        sOut+=s+r2[i2]
                    # cols Diffs: abs(TI_1-TI_2)
                    for i12 in range(iShift,c):
                        sOut+=s+str(abs(float(r1[i12])-float(r2[i12])))  
                    # cols Prod: TI_1*TI_2
                    for i12 in range(iShift,c):
                        sOut+=s+str(float(r1[i12])*float(r2[i12])) 
                    # cols Avg: (TI_1+TI_2)/2
                    for i12 in range(iShift,c):
                        sOut+=s+str((float(r1[i12])+float(r2[i12]))/float(2))
                    sOut+="\n"
                    o.write(sOut) # write the pair
                    sOut=""       # initialize the output var
                    continue # break ?
    o.close()
    return


# ---------------------------------------------------------------------------------------------------------------
def PD_TIs2PairsActivFile(actFile,inFile1,inFile2,outFile,iShift1,iShift2):
    # create a file with pairs of protein - drug using an activity pair file
    
    # iShift1,2 = no of columns in input files to skip before TIs
    # actFile = activity file PDBChain Drug Activ [tab]
    # inFile1 = protein simple result file with TIs
    # inFile1 = drug simple result file with TIs

    PairsActiv=[]  # list with pairs and activity from activity file

    # get the pairs from activity file
    fAct = open(actFile,"r")
    liness=fAct.readlines()
    for line in liness:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        PairsActiv.append(line.split("\t"))    
    
    # process PROTEIN simple result file with TIs
    Prot_TIs=[] # PDBchains and the TIs for all the simple result file 
    f = open(inFile1,"r")
    lines=f.readlines()
    f.close()
    
    for line in lines:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        Prot_TIs.append(line.split("\t"))

    l1=len(Prot_TIs)  # number of lines
    r1 = Prot_TIs[0]  # record = one line of data
    c1 = len(r1)      # number of columns

    # process DRUG simple result file with TIs
    Drug_TIs=[] # PDBchains and the TIs for all the simple result file 
    f = open(inFile2,"r")
    lines=f.readlines()
    f.close()
    
    for line in lines:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
        Drug_TIs.append(line.split("\t"))

    l2=len(Drug_TIs)  # number of lines
    r2 = Drug_TIs[0]  # record = one line of data
    c2 = len(r2)      # number of columns

    o = open(outFile,"w") # open the pair file (columns separated by space)
    s = "\t"              # column separator
    
    #----------------------
    # HEADER

    # header indices
    ProtTItypes=[]
    DrugTItypes=[]
    IndTypes=["EM", "Pol", "vdWA", "AtCon2P"]
    for Indice in IndTypes:
        for iOrb in range(5):            # for each orbital
            ProtTItypes.append(Indice+"_Orb"+str(iOrb)) # write TAB delimited [Label]Orb[iOrb]  
    AtomTypes=['TotA','Cs','Ci','Hal','Het','Hx']
    for Indice in IndTypes:
        for atom in AtomTypes:            # for each orbital
            DrugTItypes.append(Indice+"_"+str(atom)) # write TAB delimited [Label]Orb[iOrb]
            
    sOut=""  # output variable to write in Pair file
    # col 1 = L = link
    sOut+="Activity"
    # col 2-3 = PDBChain Drug
    sOut+=s+"Protein"+s+"Drug"
    # cols TI_1 headers
    for i in range(iShift1,c1):sOut+=s+r1[i]+"(P)"
    # cols TI_2 headers
    for j in range(iShift2,c2):sOut+=s+r2[j]+"(D)"

    # mixed indices: Diff, Prod, Avg
    MixedTypes=["Diff","Prod","Avg"]
    for mixed in MixedTypes:
        for Indice in IndTypes:
            for iOrb in range(5):            # for each orbital
                for atom in AtomTypes:            # for each orbital
                    sOut+=s+mixed+"_"+Indice+"_O"+str(iOrb)+"-"+str(atom) # write TAB delimited mixed indices Diff/Prod/Avg, Weigth, Orbital, AtomType
    sOut+="\n"
    # end of heading
    
    o.write(sOut) # write the header
    sOut=""       # initialize the output var

    #--------------------------------------------------------------------------
    # OUTPUT pairs

    # averages by TI type over all k=0-5 (6)
    lP=len(ProtTItypes) # 5
    lD=len(DrugTItypes) # 6
    
    # for each line in the activity file
    for pairs in PairsActiv:
        P1=pairs[0]  # PDBchain
        P2=pairs[1]  # drug name
        Act=pairs[2] # activity

        # flags for founding the pairs
        P1founded=0
        P2founded=0
        
        # find PAIR 1 = Prot
        for i in range(l1): # TIs from simple result file
            r1=Prot_TIs[i]
            Prot=r1[0]      # Prot pair 1
            # if protein found
            if P1==Prot:
                P1founded=1
                break
        # find PAIR 2 = Drug
        for j in range(l2): # TIs from simple result file     
            r2=Drug_TIs[j]
            Drug=r2[0]      # Drug pair 2
            # if drug found
            if P2==Drug:
                P2founded=1
                break
        # if both pair founded in the activity input
        if (P1founded==1) and (P2founded==1):
            sOut+=str(Act)+s+P1+s+P2 # write PDBchains of the pair
            
            # get the indices for proteins and drugs to be used for mixed indices
            NewProtTIs=zeros((lP))
            NewDrugTIs=zeros((lD))
            # write cols TI_1 for Protein
            for i1 in range(iShift1,c1):
                sOut+=s+r1[i1]
               
            # write cols TI_2 for Drug
            for i2 in range(iShift2,c2):
                sOut+=s+r2[i2]
                
            # mixed pairs: averaged values for type of TI over all k
            # cols Diffs: abs(TI_1-TI_2)
            
            for iTI in range(lP):NewProtTIs[iTI]=float(r1[iShift1+iTI])
            for iTI in range(lD):NewDrugTIs[iTI]=float(r2[iShift2+iTI])

            for iWin in range(len(IndTypes)): # for each type of weigth
                PInd=[]
                for iTI1 in range(5): # 5 orbitals
                    PInd.append(float(NewProtTIs[iTI1*iWin]))
                # get the indices for drug for one weight
                DInd=[]
                for iTI2 in range(len(AtomTypes)): # atom types
                    DInd.append(float(NewDrugTIs[iTI2*iWin]))
                # calculate the mixed indices
                
                # cols Diff: TI_1-TI_2
                for Pr in PInd:
                    for Dr in DInd:
                        sOut+=s+str(abs(float(Pr)-float(Dr)))
                # cols Prod: TI_1*TI_2
                for Pr in PInd:
                    for Dr in DInd:
                        sOut+=s+str(float(Pr)*float(Dr))
                # cols Prod: (TI_1+TI_2)/2
                for Pr in PInd:
                    for Dr in DInd:
                        sOut+=s+str((float(Pr)+float(Dr))/float(2))                    
            sOut+="\n"
            o.write(sOut) # write the pair
            sOut=""       # initialize the output var
        else:
            print "NOT FOUNDED: "+ P1+" - "+P2
            
    o.close()
    return

# ---------------------------------------------------------------------------------------------------------------
def ExtraNegRndPairs(actFile,inProtRes,inDrugRes,actFile2,m): # m=multiple of linked
    # actFile = original activity file with positive pairs
    # inProtRes = PROTEIN simple result
    # inDrugRes = DRUG simple result
    # actFile2 = new acitivity file including negative pairs (randoms or all combinations)
    
    RndListPairs={} # new dictionary that include the original pairs (positive) and the new pairs (negative)
    s = "\t"        # column separator in the files
    
    # take the original pairs (positive)
    f = open(actFile,"r")
    lines=f.readlines()
    for line in lines:
        # removing strange characters
        if line[-1]=="\n" : line = line[:-1]
        if line[-1]=="\r" : line = line[:-1]
        PDC=line.split(s) # prot drug class
        # add to the new dictionary the original positive pairs if the line contains at least 3 columns (pair1, pair2, activ)
        if len(PDC)>=3:
            RndListPairs[PDC[0]+"&"+PDC[1]]= PDC[2]

    # search for random negative pairs
    # get the sorted unique list of PDBchains/drugs from simple result files
    ProtItemList=GetUniqueListFromColumn(inProtRes,1)
    DrugItemList=GetUniqueListFromColumn(inDrugRes,1)

    # sort the lists
    ProtItemList.sort()
    DrugItemList.sort()

    # dimension of thje lists
    iP=len(ProtItemList)
    iD=len(DrugItemList)

    o = open(actFile2,"w") # open the new pair file (columns separated by space)
    # add the original pairs
    for line in lines:
        # removing strange characters
        if line[-1]=="\n" : line = line[:-1]
        if line[-1]=="\r" : line = line[:-1]
        o.write(line+"\n")
        
    # if pairs*pairs < m => take all the pairs
    m=m*len(RndListPairs) # m times the number of positive pairs
    if m>iP*iD:
        print " ... Generate all possible pairs ..."
        for P1 in ProtItemList:
            for P2 in DrugItemList:
                # if there is not in the pairs
                if P1+"&"+P2 not in RndListPairs:
                    # add negative pairs to the dictionaty
                    RndListPairs[P1+"&"+P2]= "0"
                    # write negative pair
                    o.write(P1+s+P2+s+"0\n")

    else:
        # if not, take a number m of random unique pairs
        print " ... Generate random pairs ..."
        
        iPairs=0 # no of new pairs
        fLimit=0 # limit number when it is not possible to generate more pairs
        # if pair no. is less then we need, find a pair
        while iPairs<m:
            # if the limit was reach stop the searching (if no pair found n times = brutal condition)
            if fLimit==m:
                break
            P1=choice(ProtItemList)
            P2=choice(DrugItemList)
            # if there is not in the pairs
            if P1+"&"+P2 not in RndListPairs:
                iPairs+=1
                # add negative pairs to the dictionaty
                RndListPairs[P1+"&"+P2]= "0"
                # write negative pair (only m numbers!)
                o.write(P1+s+P2+s+"0\n")
                fLimit=0
            else:
                # increase the limit number
                fLimit+=1

    o.close()
    return


#########################################

##inFile="MI-ProtResults.txt"
##outFile="MI-Prot_Pairs.txt"
##m=2
##
##TIs2PairFile(inFile,outFile,m)
