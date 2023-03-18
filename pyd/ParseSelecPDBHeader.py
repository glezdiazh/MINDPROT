# parses selective fields from the header of PDB files into a python dictionary

import re

## -----------------------------------------------------------------------------------------
def _get_journal(inl):
    # JRNL        AUTH   L.CHEN,M.DOI,F.S.MATHEWS,A.Y.CHISTOSERDOV,           2BBK   7
    journal=""
    for l in inl:
        if re.search("\AJRNL",l):
            journal+=l[19:72].lower()
    journal=re.sub("\s\s+"," ",journal)
    return journal
## -----------------------------------------------------------------------------------------
def _get_references(inl):
    # REMARK   1 REFERENCE 1                                                  1CSE  11
    # REMARK   1  AUTH   W.BODE,E.PAPAMOKOS,D.MUSIL                           1CSE  12
    references=[]
    actref=""
    for l in inl:        
        if re.search("\AREMARK   1",l):
            if re.search("\AREMARK   1 REFERENCE",l):
                if actref!="":
                    actref=re.sub("\s\s+"," ",actref)
                    if actref!=" ":
                        references.append(actref)
                    actref=""
            else:
                actref+=l[19:72].lower()

    if actref!="":
        actref=re.sub("\s\s+"," ",actref)
        if actref!=" ":
            references.append(actref)
    return references    
## -----------------------------------------------------------------------------------------      
# bring dates to format: 1909-01-08
def _format_date(pdb_date):
    """Converts dates from DD-Mon-YY to YYYY-MM-DD format."""
    date=""
    year=int(pdb_date[7:])
    if year<50:
        century=2000
    else:
        century=1900            
    date=str(century+year)+"-"
    all_months=['xxx','Jan','Feb','Mar','Apr','May','Jun','Jul',\
    'Aug','Sep','Oct','Nov','Dec']        
    month=str(all_months.index(pdb_date[3:6]))
    if len(month)==1:
        month = '0'+month
    date = date+month+'-'+pdb_date[:2]
    return date
## -----------------------------------------------------------------------------------------
def _chop_end_codes(line):
    """Chops lines ending with  '     1CSA  14' and the like."""
    return re.sub("\s\s\s\s+[\w]{4}.\s+\d*\Z","",line)
## -----------------------------------------------------------------------------------------
def _chop_end_misc(line):
    """Chops lines ending with  '     14-JUL-97  1CSA' and the like."""
    return re.sub("\s\s\s\s+.*\Z","",line)
## -----------------------------------------------------------------------------------------
def _nice_case(line):
    """Makes A Lowercase String With Capitals."""
    l=line.lower()
    s=""
    i=0
    nextCap=1
    while i<len(l):
        c=l[i]
        if c>='a' and c<='z' and nextCap:
            c=c.upper()
            nextCap=0
        elif c==' ' or c=='.' or c==',' or c==';' or c==':' or c=='\t' or\
            c=='-' or c=='_':
            nextCap=1            
        s+=c
        i+=1
    return s
## -----------------------------------------------------------------------------------------
def parse_pdb_header(infile):
    """
    Returns the header lines of a pdb file as a dictionary.

    Dictionary keys are: head, deposition_date, release_date, structure_method,
    resolution, structure_reference, journal_reference, author and
    compound.
    """
    header = []
    do_close = False
    if isinstance(infile, basestring):
        f = open(infile,'r')
        do_close = True
    else:
        f = infile
    for l in f:
        record_type=l[0:6]
        if record_type=='ATOM  ' or record_type=='HETATM' or record_type=='MODEL ':
            break
        else:
            header.append(l)    
    if do_close:
        f.close()
    return _parse_pdb_header_list(header)
## -----------------------------------------------------------------------------------------
def _parse_pdb_header_list(header):
    # database fields
    dict={'name':"",
        'head':'',
        'deposition_date' : "1909-01-08",
        'release_date' : "1909-01-08",
        'structure_method' : "unknown",
        'resolution' : 0.0,
        'structure_reference' : "unknown",
        'journal_reference' : "unknown",
        'author' : "",
        'compound':{'1':{'misc':''}},'source':{'1':{'misc':''}}}

    dict['structure_reference'] = _get_references(header)
    dict['journal_reference'] = _get_journal(header)
    comp_molid="1"
    src_molid="1"
    last_comp_key="misc"
    last_src_key="misc"

    for hh in header:
        h=re.sub("[\s\n\r]*\Z","",hh) # chop linebreaks off
        #key=re.sub("\s.+\s*","",h)
        key = h[:6].strip()
        #tail=re.sub("\A\w+\s+\d*\s*","",h)
        tail = h[10:].strip()
        # print key+":"+tail
        
        # From here, all the keys from the header are being parsed
        if key=="TITLE":
            name=_chop_end_codes(tail).lower()
            if 'name' in dict:
                dict['name'] += " "+name
            else:
                dict['name']=name
        elif key=="HEADER":            
            rr=re.search("\d\d-\w\w\w-\d\d",tail)
            if rr!=None:
                dict['deposition_date']=_format_date(_nice_case(rr.group()))
            head=_chop_end_misc(tail).lower()
            dict['head']=head
        elif key=="COMPND":            
            tt=re.sub("\;\s*\Z","",_chop_end_codes(tail)).lower()
            # look for E.C. numbers in COMPND lines
            rec = re.search('\d+\.\d+\.\d+\.\d+',tt)
            if rec:
                dict['compound'][comp_molid]['ec_number']=rec.group()
                tt=re.sub("\((e\.c\.)*\d+\.\d+\.\d+\.\d+\)","",tt)
            tok=tt.split(":")
            if len(tok)>=2:
                ckey=tok[0]
                cval=re.sub("\A\s*","",tok[1])
                if ckey=='mol_id':
                    dict['compound'][cval]={'misc':''}
                    comp_molid=cval
                    last_comp_key="misc"
                else:
                    dict['compound'][comp_molid][ckey]=cval            
                    last_comp_key=ckey
            else:
                dict['compound'][comp_molid][last_comp_key]+=tok[0]+" "
        elif key=="SOURCE":
            tt=re.sub("\;\s*\Z","",_chop_end_codes(tail)).lower()
            tok=tt.split(":")
            # print tok
            if len(tok)>=2:
                ckey=tok[0]
                cval=re.sub("\A\s*","",tok[1])
                if ckey=='mol_id':
                    dict['source'][cval]={'misc':''}
                    comp_molid=cval
                    last_src_key="misc"
                else:
                    dict['source'][comp_molid][ckey]=cval            
                    last_src_key=ckey
            else:
                dict['source'][comp_molid][last_src_key]+=tok[0]+" "
        elif key=="KEYWDS":
            kwd=_chop_end_codes(tail).lower()
            if 'keywords' in dict:
                dict['keywords']+=" "+kwd
            else:
                dict['keywords']=kwd
        elif key=="EXPDTA":
            expd=_chop_end_codes(tail)
            # chop junk at end of lines for some structures
            expd=re.sub('\s\s\s\s\s\s\s.*\Z','',expd)
            # if re.search('\Anmr',expd,re.IGNORECASE): expd='nmr'
            # if re.search('x-ray diffraction',expd,re.IGNORECASE): expd='x-ray diffraction'
            dict['structure_method']=expd.lower()
        elif key=="CAVEAT":
            # make Annotation entries out of these!!!
            pass
        elif key=="REVDAT":
            rr=re.search("\d\d-\w\w\w-\d\d",tail)
            if rr!=None:
                dict['release_date']=_format_date(_nice_case(rr.group()))
        elif key=="JRNL":
            # print key,tail
            if 'journal' in dict:
                dict['journal']+=tail
            else:
                dict['journal']=tail
        elif key=="AUTHOR":
            auth = _nice_case(_chop_end_codes(tail))
            if 'author' in dict:
                dict['author']+=auth
            else:
                dict['author']=auth
        elif key=="REMARK":
            if re.search("REMARK   2 RESOLUTION.",hh):
                r=_chop_end_codes(re.sub("REMARK   2 RESOLUTION.",'',hh))
                r=re.sub("\s+ANGSTROM.*","",r)
                try:
                    dict['resolution']=float(r)
                except:
                    #print 'nonstandard resolution',r
                    dict['resolution']=None
        else:
            # print key
            pass
    if dict['structure_method']=='unknown': 
        if dict['resolution']>0.0: dict['structure_method']='x-ray diffraction'
    return dict
## -----------------------------------------------------------------------------------------
def Get1DictPDBHeader(filename):
    # read the PDB header
    handle = open(filename,'r')
    data_dict = parse_pdb_header(handle)
    handle.close()

    # get selective elements as a dictionary
    SelectHeader={}
    sOut=""
    for field in ['structure_method','head','journal','journal_reference','keywords','name','author','deposition_date','release_date','resolution']:
        if field in data_dict.keys():
            SelectHeader[field]=data_dict[field]
        else:
            SelectHeader[field]="-"
    # compound
    if 'compound' in data_dict.keys():
        comp_list=data_dict['compound']
        comp_dets=comp_list['1']
        for field in ['molecule','engineered','chain','ec']:
            if field in comp_dets.keys():
                SelectHeader[field]=comp_dets[field]
            else:
                SelectHeader[field]="-"      
    # source
    if 'source' in data_dict.keys():
        source_list=data_dict['source']
        source_dets=source_list['1']
        for field in ['expression_system_vector_type','gene','expression_system_taxid','organism_scientific','organ','tissue','expression_system','expression_system_plasmid','expression_system_strain','cell_line','cellular_location','organism_taxid','organism_common']:
            if field in source_dets.keys():
                SelectHeader[field]=source_dets[field]
            else:
                SelectHeader[field]="-"                
    return SelectHeader # return a selected dictionary
## -----------------------------------------------------------------------------------------
def GetSelectivePDBHeader(AllPDBHeader,FieldList=['cellular_location','tissue','organ','organism_scientific','organism_common','ec','expression_system_vector_type','expression_system_taxid','expression_system','engineered','head']):
    # get selective elements as a dictionary
    SelectHeader={}
    for field in FieldList:
        if field in AllPDBHeader.keys():
            SelectHeader[field]=AllPDBHeader[field]
        else:
            SelectHeader[field]="-" 
    return SelectHeader

################################################
### main
################################################

##if __name__=='__main__':
##    # Reads a PDB file passed as argument, parses its header, extracts
##    # some data and returns it as a dictionary.
##    import sys
##    filename = "1A0H.pdb" # sys.argv[1]
##    handle = open(filename,'r')
##    data_dict = parse_pdb_header(handle)
##    handle.close()
##
##    # print the dictionary
##    for k, y in data_dict.iteritems():
##        print "-"*40
##        print k
##        print y
##
##
##
##    print "="*10
##    filename = "1A0H.pdb"
##    print GetSelectivePDBHeader(AllPDBHeader=Get1DictPDBHeader(filename)) 
