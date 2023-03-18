from numpy import *

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

# ------------------------------------------------------------------------------
def AverageTIs_PDBHeader(sIn,sOut,HeaderList,iShift):
    # HeaderList = list with user PDB header fields in GUI
    
    # read the input file
    fIn=open(sIn,'r')
    lines=fIn.readlines()
    for i in range(len(lines)):
        OneLine=lines[i]
        if OneLine[-1]=="\n" : lines[i] = OneLine[:-1] # removing strange characters
        if OneLine[-1]=="\r" : lines[i] = OneLine[:-1]
    
    # take the header line
    FirstLine=lines[0]           # get the header

    # numbers of header items = columns
    c=len(FirstLine.split("\t"))
    l=len(lines) # no of lines in the input file FULL header info
    nFields=11 # max no of PDB header fields
    nIndices=c-iShift-nFields   # no of calculated indices
    nUserFields=len(HeaderList) # no of fields to be averaged

    # get the header as a list of [PDB],[fields],[indices labels]
    Headers=[]
    colFirstLine=FirstLine.split("\t")
    # get the PDBs
    temp=[]
    for s in range(iShift):temp.append(colFirstLine[s])
    Headers.append(temp)
    # get the fields
    temp=[]
    FieldPos=[]
    for f in range(nFields):
        for userH in HeaderList: # print only the user header
            if colFirstLine[iShift+f]==userH:
                temp.append(colFirstLine[iShift+f]) # get the value
                FieldPos.append(f) # get the position of the field in full fields
    Headers.append(temp)
    # get the indices
    temp=[]
    for n in range(nIndices):temp.append(colFirstLine[iShift+nFields+n])
    Headers.append(temp)

    # add the future calculations labels: Avgs and Diff
    temp=[]
    for h in range(len(HeaderList)): # for each user field
        IndLabels=Headers[2]
        for n in range(len(IndLabels)): # for each indice label
            temp.append("AvgBy"+HeaderList[h]+"["+IndLabels[n]+"]")
        for n in range(len(IndLabels)): # for each indice label
            temp.append("Diff_"+IndLabels[n]+"-AvgBy"+HeaderList[h])
    Headers.append(temp)

    # get the info in 3 lists: case name, fields, indices
    PDBs=[]     # list with PDBs (col 1)
    Fields=[]   # list with the fields (nFields cols)
    Indices=[]  # list with the indices (nIndices cols)

    # process each line in the file except the header
    for i in range(1,l):
        curRow=(lines[i]).split("\t") # take the values by columns
        # get the PDBs
        temp=[]
        for s in range(iShift):temp.append(curRow[s])
        PDBs.append(temp)
        # get the fields
        temp=[]
        for f in range(nFields):
            for fPos in range(len(FieldPos)): # verify the positions of the user fields in the list
                if f==int(FieldPos[fPos]):
                    temp.append(curRow[iShift+f])
        Fields.append(temp)
        # get the indices
        temp=[]
        for n in range(nIndices):temp.append(curRow[iShift+nFields+n])
        Indices.append(temp) 
        
    # check all the user fields and take the unique values
    UniqueFields=[] # list with list of unique headers in the same order as in HeaderList
    Avgs_Diffs=[] # list with Averages and diffs Indice-Avg to be completed as [avgs1],[diffs1]..etc

    for h in range(len(HeaderList)): # for each user field
        # make unique lists
        HFields=Headers[1] # get only the header fields
        for f in range(len(HFields)): # check with all the list of fields
            if HeaderList[h]==HFields[f]:
                temp=[] # all values of one field
                for k in range(l-1): # l-1 = lines - header=cases
                    kFields=Fields[k]
                    temp.append(kFields[f]) # the fields has the same position in the header with the list of Fields
        uField=unique(temp) # unique list with one field
        #UniqueFields.append(uField) # list with list of unique values for eahc user field

        #print "uField=",uField
        # calculate the averages for each field value and indice
        Avg_Field=  zeros((len(uField),nIndices)) # avgs for a values in the current user field and an indice
        Count_Field= zeros((len(uField),nIndices))

        # take each line of data and compare with unique values of fields
        HFields=Headers[1] # get only the header fields
        for f in range(len(HFields)): # check with all the list of fields
            if HeaderList[h]==HFields[f]: # if user and Full header fields coincide
                for u in range(len(uField)): # for each unique values
                    for k in range(l-1): # l-1 = lines - header=cases
                        kFields=Fields[k]
                        if uField[u]==kFields[f]: # if data fields are the unique values
                            curIndices=Indices[k]
                            for n in range(len(curIndices)): # add all the indices
                                Avg_Field[u][n]+=float(curIndices[n])
                                Count_Field[u][n]+=float(1)
        # calc avgs
        for u in range(len(uField)):
            for n in range(nIndices):
                Avg_Field[u][n]=float(Avg_Field[u][n])/float(Count_Field[u][n])

        # create list Avgs to be added to Avgs_Diffs
        Avgs=[]
        Diffs=[]
        # take each line of data and compare with unique values of fields
        HFields=Headers[1] # get only the header fields
        for f in range(len(HFields)): # check with all the list of fields
            if HeaderList[h]==HFields[f]: # if user and Full header fields coincide
                
                for k in range(l-1): # l-1 = lines - header=cases
                    temp1=[]
                    temp2=[]
                    kFields=Fields[k]
                    for u in range(len(uField)): # for each unique values
                        if uField[u]==kFields[f]: # if data fields are the unique values
                            curIndices=Indices[k]
                            for n in range(len(curIndices)): # add all the indices
                                # print uField[u],curIndices[n],Avg_Field[u][n]
                                temp1.append(Avg_Field[u][n])
                                temp2.append(float(curIndices[n])-Avg_Field[u][n])
                                 
                    Avgs.append(temp1)
                    Diffs.append(temp2)

        Avgs_Diffs.append(Avgs) # add the first list with Avgs
        Avgs_Diffs.append(Diffs) # add the first list with Avgs

    # write the output:
    # Headers
    # PDBs, Fields, Indices, Avgs_Diffs

    fOut=open(sOut,"w")
    sOut=""

    # write Headers
    for item in Headers:
        for details in item:
            sOut+=details+"\t"
    fOut.write(sOut[:-1]+"\n")

    # write PDBs, Fields, Indices, Avgs_Diffs
    sOut=""
    for k in range(l-1): # for each case
        cPDBs=PDBs[k]
        for item in cPDBs:
            sOut+=item+"\t"
        cFields=Fields[k]
        for item in cFields:
            sOut+=item+"\t"
        cIndices=Indices[k]
        for item in cIndices:
            sOut+=item+"\t"
        for items in Avgs_Diffs:         
            cAvgs_Diffs=items[k]
            for item in cAvgs_Diffs:
                sOut+=str(item)+"\t"
            
        sOut=sOut[:-1]+"\n"


    fOut.write(sOut)
    
    fIn.close()
    fOut.close()
    return

# ------------------------------------------------------------------------------
def AverageTIs_PDBHeader_old(sIn,sOut,HeaderList,flag_Header,iShift):
    # uses an output with PDBchain (1), heading from PDBs (17), TIs (90)
    # adds the averages values for each class in each type of header field
    # parameters: input file to modify, output file with the avgs by classes
    # list of the header items, matrix of flags for each header used
    # iShift = number of columns to skip
    
    iHeader=len(HeaderList)        # numbers of possible header items
    
    # read the input file
    fIn=open(sIn,'r')
    lines=fIn.readlines()
    for i in range(len(lines)):
        OneLine=lines[i]
        if OneLine[-1]=="\n" : lines[i] = OneLine[:-1] # removing strange characters
        if OneLine[-1]=="\r" : lines[i] = OneLine[:-1]
    
    l=len(lines)                 # no. of lines in the input file
    FirstLine=lines[0]           # get the header
    c=len(FirstLine.split("\t")) # numbers of header items = columns

    # verify the header choise and create the unique list of header idems = list of lists
    # list of list of headers (non-uniques)
    AllHeaders=[]

    # 11 empty lists for each header type
    HeaderItems=[]
    for Choise in range(iHeader):
        HeaderItems.append([])
    
    # get the user headers excluding the header line
    for i in range(1,l): # process each line in the file
        CurrR=(lines[i]).split("\t") # take the values by columns
        
        # check the choise of the user
        for h in range(iHeader): # verify each header flag            
            if flag_Header[h]==1:
                # get the items of this header
                CurrHeader=HeaderItems[h] # take the list to be completed

                # PDBchain 11*[HeaderFiels] TIs
                CurrHeader.append(CurrR[h+1]) # get the item of a specific header
                HeaderItems[h]=CurrHeader # place back the list for the specific header

    # make unique lists of specific header lists
    for j in range(len(HeaderItems)):
        HeaderItems[j]=unique(HeaderItems[j])
    
    Avg_Header=  zeros((iHeader,len(HeaderItems),c-1-iHeader))  # matrix the avg values of each TI for each Header item
    countHeader= zeros((iHeader,len(HeaderItems),c-1-iHeader))  # matrix with counting for each item in a specific header
    
    # calculate the averages if found the header label
    for i in range(1,l): # process each line in the file
        CurrR=(lines[i]).split("\t") # take the values by columns
        # check the choise of the user
        for h in range(iHeader): # verify each header flag
            if flag_Header[h]==1:
                CurrH=HeaderItems[h] # get the list of unique items of a special header
                # print "CurrH=\n",CurrH
                # process each unique item of header
                for ii in range(len(CurrH)):
                    # print "CurrR[h+1],item=",CurrR[h+1],CurrH[ii]
                    # print CurrR
                    if CurrR[h+1]==CurrH[ii]: # if item found
                        for cc in range(c-iShift-iHeader): #for each TI column
                            Avg_Header[h][ii][cc]+=float(CurrR[cc+iHeader+1])
                            countHeader[h][ii][cc]+=1
    # print HeaderItems

    # create the avg values
    for h in range(iHeader): # verify each header flag
        if flag_Header[h]==1:
            CurrH=HeaderItems[h] # get the list of unique items of a special header
            # process each unique item of header
            for ii in range(len(CurrH)):
                for cc in range(c-iShift-iHeader): #for each TI column
                    if countHeader[h][ii][cc]!=0:
                        Avg_Header[h][ii][cc]=float(Avg_Header[h][ii][cc])/float(countHeader[h][ii][cc])
            
    # write the output file by adding to the header result the columns with the averages
    fOut=open(sOut,"w")
    sOut=""

    for iL in range(l):          # for each old line
        oldLine=lines[iL]        # older header
        lines[iL]=oldLine[:-1]   # remove the end of line
        sOut+=lines[iL]
        CurrR=(lines[iL]).split("\t") # take the values by columns
        
        # if header, we are completing with the new header
        if iL==0: 
            for h in range(iHeader):  # verify each header flag
                if flag_Header[h]==1: # only if header is used
                    for cc in range(c-iShift-iHeader): #for each TI column
                        sOut+="\t"+"AvgBy["+HeaderList[h]+"]for["+CurrR[cc+iHeader+iShift]+"]" # add the new header avgs
            #fOut.write(sOut+"\n")

            # add Diffs TIs - Avgs
            for h in range(iHeader):  # verify each header flag
                if flag_Header[h]==1: # only if header is used
                    for cc in range(c-iShift-iHeader): #for each TI column
                        sOut+="\t"+"Diff["+CurrR[cc+iHeader+iShift]+"]-[AvgBy"+HeaderList[h]+"]" # add the new header avgs
            fOut.write(sOut+"\n")
            
        else: # the TIs lines
            sOut=""
            sOut+=lines[iL]
            for h in range(iHeader): # verify each header flag
                if flag_Header[h]==1:
                    CurrH=HeaderItems[h] # get the list of unique items of a special header
                    # process each unique item of header
                    for ii in range(len(CurrH)):
                        if CurrR[h+1]==CurrH[ii]: # if item found
                            for cc in range(c-iShift-iHeader): #for each TI column
                                sOut+="\t"+str(Avg_Header[h][ii][cc]) # add the new avg of TIs

                    # add differences TIs - Avgs
                    for ii in range(len(CurrH)):
                        if CurrR[h+1]==CurrH[ii]: # if item found
                            for cc in range(c-iShift-iHeader): #for each TI column
                                sOut+="\t"+str(float(CurrR[cc+iHeader+1])-Avg_Header[h][ii][cc]) # add the new avg of TIs
                                
            fOut.write(sOut+"\n")
            
    fOut.close()
    fIn.close()
    return

# --------------------------------------------------------------------------
# averaged by input file
# ---------------------------------------------------------------------------------------------------------------
def GetInputClass(sInFile,Col):
    # list with the unique classes
    ClassList=[]
    
    fIn=open(sInFile,"r") #open file to read
    lines=fIn.readlines()
    
    l=len(lines)                 # no. of lines in the input file
    FirstLine=lines[0]           # get the header
    c=len(FirstLine.split("\t")) # numbers of header items = columns
    # only if there are more columns
    if c>=Col:
        for i in range(l): # process each line in the file
            CurrL=lines[i] # current line
            if CurrL[-1]=="\n" : CurrL = CurrL[:-1] # removing strange characters
            if CurrL[-1]=="\r" : CurrL = CurrL[:-1]
            CurrR=(CurrL).split("\t")  # split each line by columns tab separated
            ClassList.append(CurrR[Col-1]) # add to the non-unique list a specific value for the column Col

    return unique(ClassList) # return the unique classes from the input

# ---------------------------------------------------------------------------------------------------------------
def AverageTIs4Col(sIn,sOut,Classes,Col):
    # transform sIn file in sOut file using the Classes list with unique classes
    # Col = class columns in the input file
    
    iClasses=len(Classes) # numbers of possible class items
    
    # read the input file
    fIn=open(sIn,'r')
    lines=fIn.readlines()
    for line in lines:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
    
    l=len(lines)                 # no. of lines in the input file
    FirstLine=lines[0]           # get the header
    c=len(FirstLine.split("\t")) # numbers of header items = columns

    # make unique lists of specific header lists HeaderItems
    
    Avg_Header=  zeros((iClasses,c-Col))  # matrix the avg values of each TI for each classs
    countHeader= zeros((iClasses,c-Col))  # matrix with counting
    
    # calculate the averages if found the header label
    for i in range(1,l): # process each line in the file
        CurrR=(lines[i]).split("\t") # take the values by columns
        # process each unique class
        for h in range(iClasses): # verify each class
            CurrH=Classes[h] # get the list of unique class

            if CurrR[Col-1]==CurrH: # if item found
                for cc in range(c-Col): #for each TI column
                    Avg_Header[h][cc]+=float(CurrR[cc+Col])
                    countHeader[h][cc]+=1
                    continue

    # create the avg values
    for h in range(iClasses): # verify each header flag
        CurrH=Classes[h] # get the list of unique items of a special header
        # process each unique item of header
        for cc in range(c-Col-iClasses): #for each TI column
            if countHeader[h][cc]!=0:
                Avg_Header[h][cc]=float(Avg_Header[h][cc])/float(countHeader[h][cc])

    # write the output file by adding to the header result the columns with the averages
    fOut=open(sOut,"w")
    sOut=""

    for iL in range(l):          # for each old line
        oldLine=lines[iL]        # older header
        lines[iL]=oldLine[:-1]   # remove the end of line
        sOut+=lines[iL]
        CurrR=(lines[iL]).split("\t") # take the values by columns
        
        # if header, we are completing with the new header
        if iL==0: 
            for cc in range(c-Col): #for each TI column
                sOut+="\t"+"Avg"+CurrR[cc+Col] # add the new header avgs

            # add Diffs: TI-Avg
            for cc in range(c-Col): #for each TI column
                sOut+="\t"+"Diff["+CurrR[cc+Col]+"]-[Avg]" # add the new header avgs
            fOut.write(sOut+"\n")
            
        else: # the TIs lines
            sOut=""
            sOut+=lines[iL]
            for h in range(iClasses): # verify each header flag
                CurrH=Classes[h] # get the list of unique items of a special header
                # process each unique class
                if CurrR[Col-1]==CurrH: # if item found
                    for cc in range(c-Col): #for each TI column
                        sOut+="\t"+str(Avg_Header[h][cc]) # add the new avg of TIs

                    # add Diffs values
                    for cc in range(c-Col): #for each TI column
                        sOut+="\t"+str(float(CurrR[cc+Col])-Avg_Header[h][cc]) # add the new avg of TIs
                        
            fOut.write(sOut+"\n")
            
    fOut.close()
    fIn.close()
    return

# ---------------------------------------------------------------------------------------------------------------
def DrugAverageTIs4Col(sIn,sOut,Classes,Col):
    # transform sIn file in sOut file using the Classes list with unique classes
    # Col = class columns in the input file
    
    iClasses=len(Classes) # numbers of possible class items
    
    # read the input file
    fIn=open(sIn,'r')
    lines=fIn.readlines()
    
    for line in lines:
        if line[-1]=="\n" : line = line[:-1] # removing strange characters
        if line[-1]=="\r" : line = line[:-1]
    
    l=len(lines)                 # no. of lines in the input file
    FirstLine=lines[0]           # get the header
    c=len(FirstLine.split("\t")) # numbers of header items = columns

    # make unique lists of specific header lists HeaderItems
    
    Avg_Header=  zeros((iClasses,c-Col))  # matrix the avg values of each TI for each classs
    countHeader= zeros((iClasses,c-Col))  # matrix with counting
    
    # calculate the averages if found the header label
    for i in range(1,l): # process each line in the file
        CurrLine=lines[i]
        if CurrLine[-1]=="\n" : CurrLine = CurrLine[:-1] # removing strange characters
        if CurrLine[-1]=="\r" : CurrLine = CurrLine[:-1]
        CurrR=CurrLine.split("\t") # take the values by columns
        
        # process each unique class
        for h in range(iClasses): # verify each class
            CurrH=Classes[h] # get the list of unique class

            if CurrR[Col-1]==CurrH: # if item found
                for cc in range(c-Col): #for each TI column
                    Avg_Header[h][cc]+=float(CurrR[cc+Col])
                    countHeader[h][cc]+=1
                    continue

    # create the avg values
    for h in range(iClasses): # verify each header flag
        CurrH=Classes[h] # get the list of unique items of a special header
        # process each unique item of header
        for cc in range(c-Col-iClasses): #for each TI column
            if countHeader[h][cc]!=0:
                Avg_Header[h][cc]=float(Avg_Header[h][cc])/float(countHeader[h][cc])

    # write the output file by adding to the header result the columns with the averages
    fOut=open(sOut,"w")
    sOut=""

    for iL in range(l):          # for each old line
        oldLine=lines[iL]        # older header
        lines[iL]=oldLine[:-1]   # remove the end of line
        sOut+=lines[iL]
        CurrR=(lines[iL]).split("\t") # take the values by columns
        
        # if header, we are completing with the new header
        if iL==0: 
            for cc in range(c-Col): #for each TI column
                sOut+="\t"+"Avg"+CurrR[cc+Col] # add the new header avgs

            # add Diffs: TI-Avg
            for cc in range(c-Col): #for each TI column
                sOut+="\t"+"Diff["+CurrR[cc+Col]+"]-[Avg]" # add the new header avgs
            fOut.write(sOut+"\n")
            
        else: # the TIs lines
            sOut=""
            sOut+=lines[iL]
            for h in range(iClasses): # verify each header flag
                CurrH=Classes[h] # get the list of unique items of a special header
                # process each unique class
                if CurrR[Col-1]==CurrH: # if item found
                    for cc in range(c-Col): #for each TI column
                        sOut+="\t"+str(Avg_Header[h][cc]) # add the new avg of TIs

                    # add Diffs values
                    for cc in range(c-Col): #for each TI column
                        sOut+="\t"+str(float(CurrR[cc+Col])-Avg_Header[h][cc]) # add the new avg of TIs
                        
            fOut.write(sOut+"\n")
            
    fOut.close()
    fIn.close()
    return

############################################################
### main
############################################################
##
##sIn="PROT_SimpleRes.txt"
##sOut="AvgByInput.txt"
##
##sActiv="PDBlist.txt"
##Col=2 # column in the input file of the class (PDBchain\tClass)
##Classes=[]
##Classes=GetInputClass(sActiv,Col)
##
##print Classes
##
##AverageTIs4Col(sIn,sOut,Classes,Col)


##########################################################
# main
##########################################################

##sIn="MI-ProtResults_Header.txt"
##sOut="MI-ProtResults_Header_Avgs.txt"
##
##HeaderList=['head','expression_system','expression_system_taxid','name','chain','organism_scientific','molecule','expression_system_vector_type','ec','organism_common','expression_system_plasmid','engineered','expression_system_strain','cell_line','cellular_location','gene','organism_taxid']
##flag_Header= zeros((len(HeaderList)))  # matrix with flags for each header type, default=0
##
### virtual choise from the interface
##flag_Header[0]=1
##flag_Header[6]=1
##
##AverageTIs_PDBHeader(sIn,sOut,HeaderList,flag_Header)
