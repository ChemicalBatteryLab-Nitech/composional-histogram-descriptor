#! /usr/bin/env python

import re
import os
import sys
import pandas as pd
import numpy as np
import argparse
import math
import itertools
import tqdm
import subprocess as sp

path=os.getcwd()
home = os.environ['HOME']
defvasp ="{}/bin/DefVASP".format(home)     
defelem="{}/bin/DefElem".format(home)
vfsoftware="vfdir/vf.software"
defvasp_am="{}/bin/DefVASP_am".format(home)
defggau="{}/bin/DefGGAU".format(home)

for i in tqdm.tqdm(range(int(1e7))):
    np.pi*np.pi

lsbefore=str(sp.check_output("ls -alh --full-time",shell=True))
ls_bflist =lsbefore.split('\\n')
ls_bflist.pop(0)
 
parser = argparse.ArgumentParser()
parser.add_argument('-aparam1',type=str, default='histdef.dat',help='histdef.dat formatted structure file (def. histdef.dat)')
parser.add_argument('structure',type=str, help='enter a compositon name')
parser.add_argument('-aparam2',help=" -aparam2=sigma\: sigma value for Gaussian broadning\n")
parser.add_argument('-aparam3',help=" -aparam3=T\: output histgram descriptors considering two elements (Def. F)\n")
parser.add_argument('-sparam',default='F',help=" -sparam=T use POSCAR.names file instead of direct input.\n")

args = parser.parse_args()
reference = args.structure
histdef = args.aparam1
aparam2=args.aparam2
aparam3=args.aparam3
sparam=args.sparam
reference = reference.split()

if aparam2 is None:
        aparam2 = ''

if aparam3 is None:
        aparam3 = ''

def defvasptable(): #-----------------------------------------------------------------------       
     
     
#  $defam[$atomnum] gives atomic weight                                                          
#-- ionic energy -------------------------------------------                                     
#   output  $ie[$atomnum][$valence]{type}                                                        
#    type = use, WEL, CRC                                                                        
#   "use" and "WEL" give former in the unit "kJ/mol"  "CRC" gives latter in the unit "eV"        
#----------------------------------------------------------------------------------------- \     
                                                                                                 
     
     
#    print "reading $defvasp\n"                                                                  
     
    INvasp= open(defvasp,'r')     
    vaspdef = INvasp.readlines()   
    for in_vaspdef in vaspdef:     
        in_vaspdef.rstrip('\n')
    INvasp.close()
    dummy=vaspdef[0]
    atomlist=[]
    defpot={}
    defspin={}
    defenergyPBEsol={}
    defenergyPBE={}
    defPQ={}
    defHL={}
    atom2Z={}
    Z2atom=[]
    vaspdef.pop(0)
    for i in vaspdef:     
        field=i.split()
        atomlist.append(field[0])     
        defpot[field[0]]=field[1]     
        defspin[field[0]]=field[2]     
        defenergyPBEsol[field[0]]=field[3]     
        defenergyPBE[field[0]]=field[5]     
        defPQ[field[0]]=field[8]  #principal quantum number
        if len(field) != 9:                  
            defHL[field[0]]=field[9] 
     
        if field[0] != "ZZZ":     
            atom2Z[field[0]]=field[7]
            while len(Z2atom) -1 < int(field[7]):
                Z2atom.append('')    
            Z2atom[int(field[7])]=field[0]     
            
    """
    INie= open("defvasp_ie",'r')     
    iedef = INie.readlines()     
    for in_iedef in iedef:
        in_iedef.restrip('\n')
    INie.close()
     
    for lines in iedef:     
        nele=0     
        field=lines.split(',')     
        atomnum=field[0]     
        atomsym=field[0]     
        atomname=field[0]     
        ietype=field[0]     
     
        for j in field:     
            nele = nele + 1     
            ie[atomnum][nele][ietype]=j    
    """

    defam={}
    INvasp_am=open(defvasp_am,'r')     
    vaspdef = INvasp_am.readlines()     
    for in_vaspdef in vaspdef:     
        in_vaspdef.rstrip('\n')
    INvasp_am.close()
     
     
    for k in vaspdef:     
        field=k.split()     
        defam[field[0]]=field[1]     
         
    
    INggau= open(defggau,'r')     
    ggau = INggau.readlines()     
    for in_defggau in ggau:     
        in_defggau.rstrip('\n')
    INggau.close()
     
    orblist=["s","p","d","f","g"]     
    
    orb={}
    uval={}
    for j in ggau:     
        field=j.split()     
        nob=0     
        for ob in orblist:     
            if ob == field[1]:     
                orb[field[0]]=nob     
                nob = nob + 1     
                uval[field[0]]=field[2]     
    return defvasp,defvasp_am,defggau,atom2Z,Z2atom

def get_posatom():  #------------------------------------------------------- get_posatom    

    vrh = os.system('grep VRH place_of_potcar').split('\n')

    h=0
    for i in vrh:
        re.search(r'(VRHFIN =)(.+)(\:)',i)
        pos=2
        re.sub(r'\s+$','',i)
        posatom[h]=pos
#        print $posatom[$h]."\n";                                                                                                                                                          
        h = h+1
    

    ii=0
    iiiv=1
    for ions in ions:
        iii=1
        while iii<=pos:
            atom_pos[iiiv]=posatom[ii]                                                    
            iiiv = iiiv + 1
            iii = iii + 1
        
        ii = ii + 1
    return ions

def compconv(): #-------------------------------------------------------------------------------                                                                         
#input $incomp or $reference[0]  eg Li0.1CoO2                                                                                                                          
#output  $indexions[$i] -> element                                                                                                                                     
#        $ions[$i]      -> number of ions                                                                                                                              
#        $totalion      -> summation of ions                                                                                                                           
#        $labelcomp     -> composition  Li0.1Co1O2                                                                                                                     
#        $labelcomp2    -> composition  Li 0.1 Co 1 O 2                                                                                                                
#        $labelcomp3    -> normalized version of labelcomp2                                                                                                            
#        $labelcomp4    -> normalized version of labelcomp2 (delite space)                                                                                             
#---------------------------------------------------------------------------------------------                                                                         
    incomp=reference[0]
    
    ions=[]
    totalion=0
    indexions=[]
    factor=[]
    composition=[]
    elemcomp=[]
    nfactor=[]
    nions={}
    labelcomp=""
    labelcomp2=""
    labelcomp3=""
    labelcomp4=""

    if incomp == "" and reference[0]== "":
        print("Usage\: vaspfiles -compconv Arg1 -sparam=ArgS\n")
        print("Usage2\: input \incomp gives \symbolions[$i] \nions[$i]\n")
        print(" Arg1\ or \incomp\: pretty formula, such as Li2MnO3, H2(SO4), LiZr2(PO4)3\n")
        print(" ArgS\: -sparam=detail  outputs molar ratio\n")
        print(" Notice\: if you input composition by argument, do not use \(\), but \[\].\n")
        return
    else:

        if incomp != "":
            compositionprim = incomp
            printoff="T"
        else:
            compositionprim = reference[0]
    
    #z = readdefelem
    re.sub(r'\[','\(',compositionprim)
    re.sub(r'\]','\)',compositionprim)
    if printoff != "T":
        print ("{} input: {}\n".format(composition,compositionprim))
    i=0
    if not  bool(re.search(r'\(/',compositionprim)):
        composition.append(compositionprim)
        factor.append(1)
#       maxflag=0;                                                                                                                                                    
    else:
        if not compositionprim == "":
            if bool(re.search(r'^([A-Z].*?)(\(.*)',compositionprim)):
                search = re.search(r'^([A-Z].*?)(\(.*)',compositionprim)
                ate1 = search.group(1)
                ate2 =search.group(2)
                composition.append(ate1)
                factor.append(1)
                compositionprim = ate2
            elif bool(re.search(r'^\(([A-Z].*?)\)([0-9\.]*)([A-Z].*)',compositionprim)):
                search = re.search(r'^\(([A-Z].*?)\)([0-9\.]*)([A-Z].*)',compositionprim)
                ate1 = search.group(1)
                ate2 = search.group(2)
                ate3 = search.group(3)
                composition.append(ate1)
                factor.append(ate2)
                compositionprim=ate3
            elif bool(re.search(r'^\(([A-Z].*?)\)([0-9\.]*)(\(.*)',compositionprim)):
                search = re.search(r'^\(([A-Z].*?)\)([0-9\.]*)(\(.*)',compositionprim)
                ate1 = search.group(1)
                ate2 =search.group(2)
                ate3  = search.group(3)
                composition.append(ate1)
                if ate2 == "":
                    factor.append(1)
                else:
                    factor.append(ate2)
                compositionprim=ate3
            elif bool(re.search(r'^\(([A-Z].*)\)([0-9\.]*)',compositionprim)):
                search = re.search(r'^\(([A-Z].*)\)([0-9\.]*)',compositionprim)
                ate1 = search.group(1)                
                ate2 = search.group(2)
                composition.append(ate1)
                if ate2 == "":
                    factor[i]=1
                else:
                    factor[i]=ate2
                compositionprim=""
            elif bool(re.search(r'^\(([A-Z].*)\)',compositionprim)):
                search = re.search(r'^\(([A-Z].*)\)',compositionprim)
                ate1 = search.group(1)
                composition.append(ate1)
                factor.append(1)
                compositionprim=""
            else:
                composition.append(compositionprim)
                factor.append(1)
                compositionprim=""
#           print "Index $i    $composition[$i]  $factor[$i]\n";                                                                                                       
            i = i+1
    i=0
    j=0
    
    for compositionflag in composition:                                        
        while compositionflag != "":
            if bool(re.search(r'^([A-Z].*?)([A-Z].*)',compositionflag)):
                search = re.search(r'^([A-Z].*?)([A-Z].*)',compositionflag)
                ate1 = search.group(1)
                ate2 = search.group(2)
                elemcomp.append(ate1)
                nfactor.append(factor[j])
                compositionflag=ate2
            elif bool(re.search(r'^([A-Z].*[0-9\.]*)',compositionflag)):
                search = re.search(r'^([A-Z].*[0-9\.]*)',compositionflag)
                ate1 = search.group(1)
                elemcomp.append(ate1)
                compositionflag=""
                nfactor.append(factor[j])
#           print "$compositionflag $elecomp[$i] $ncomp[$i]\n";                                                                                                        
            i = i + 1
            if i>9999 :
                break
        j =j +1
    i=0
    
    for j in  elemcomp:
#       print "elemcomp $_\n";                                                                                                                                         
        if  bool(re.search(r'^([A-Za-z]*)([0-9\.].*)',str(j))):
            search = re.search(r'^([A-Za-z]*)([0-9\.].*)',str(j))
            ate1 = search.group(1)
            ate2 = search.group(2)
#            $indexions[$i]=$1;                                                                                                                                        
#            $ions[$i]=$2;
            nions[ate1]=0
            nions[ate1]=float(ate2)*nfactor[i]+nions[ate1]                                                                                                                      
        elif bool(re.search('^([A-Za-z]*)',str(j))):
            search = re.search('^([A-Za-z]*)',str(j))
            ate1 = search.group(1)
            nions[ate1]=0
            nions[ate1]=1*nfactor[i]+nions[ate1]
        i = i + 1
    totalion=0
    i=0
    
    for l in  sorted(nions.keys()):
        indexions.append(l)
        ions.append(nions[l])
        if printoff != "T":
            print("{} ".format(l))
            print("{} ".format(nions[l]))
        
        if nions[l] > 0:
            labelcomp=labelcomp+"{} {}".format(l,nions[l])
        
        labelcomp2=labelcomp2+"{} {} ".format(l,nions[l])
        totalion=totalion+nions[l]
        i = i + 1

    i=0
    for m in sorted(nions.keys()):
        normions = []
        for mi in range(len(nions.keys())):
            normions.append(mi)
        normions[i]="%1.6f" % (ions[i]/totalion)
        #normions[i]=print("%1.6f" % (ions[i]/totalion))
        labelcomp3=labelcomp3 + "{} {} ".format(m,normions[i])
        i = i + 1
    
    i=0
    for n in sorted(nions.keys()):
        normions4="%1.6f" % (float(ions[i])/float(totalion))
        if float(normions4) >= 0.0000001:
            labelcomp4=labelcomp4 + "{} {}".format(n,normions4)
        
        i = i + 1


    if printoff != "T":
        print ("\n")
        print ("{} ions   {}\n".format(total,totalion))

   # optional outputs------------------------------------;                                                                                                             

    if bool(re.search('[detail]',sparam,re.I)):

        print("Normalized_composition: ")
        for o in indexions:
            nnions=nions[o]/totalion
            print("{} {} ".format(o,nnions))
        print("\n")

    return ions,compositionprim,indexions,totalion

def readdefelem(atom,defelem): #------------------------------------------------------------------------------
    
    home = os.environ['HOME']
    defelem="{}/bin/DefElem".format(home)
    #print(atom)
    IN =open(defelem,'r')
    defelem=IN.readlines()
    for line in defelem:
        line=line.rstrip('\n')
    IN.close()

    skipline=defelem[0]
    el_dic = {'AN':'AN','EN':'EN','MP':'MP','PN':'PN','PG':'PG','MN':'MN','AW':'AW','AR':'AR','IR':'IR','CoR':'CoR','CrR':'CrR','spdf':'spdf'}
    for e in defelem:
        field=e.split()
        defatom=field[1]
        element = '0'
        def_atom = {defatom : element}
        el_dic['AN'] = dict(def_atom)
        el_dic['EN'] = dict(def_atom)
        el_dic['MP'] = dict(def_atom)
        el_dic['PN'] = dict(def_atom)
        el_dic['PG'] = dict(def_atom)
        el_dic['MN'] = dict(def_atom)
        el_dic['AW'] = dict(def_atom)
        el_dic['AR'] = dict(def_atom)
        el_dic['IR'] = dict(def_atom)
        el_dic['CoR'] = dict(def_atom)
        el_dic['CrR'] = dict(def_atom)
        el_dic['spdf'] = dict(def_atom)

        if atom == defatom:
            el_dic['AN'][atom]=field[0]       #atomic number
            el_dic['EN'][atom]=field[2]
            el_dic['MP'][atom]=field[4]
            el_dic['PN'][atom]=field[5]       #electronegativity;
            el_dic['PG'][atom]=field[6]
            el_dic['MN'][atom]=field[7]
            el_dic['AW'][atom]=field[8]
            el_dic['AR'][atom]=field[9]
            el_dic['IR'][atom]=field[10]
            el_dic['CoR'][atom]=field[11]
            el_dic['CrR'][atom]=field[12]
            el_dic['spdf'][atom]=field[13]
        if atom == defatom:
            return el_dic

def broad(arg0,arg1,arg2,arg3,arg4,arg5,arg6,histdef): #--------------------------------------------------------------- Broad                                                                                                        

    tlimit=0.000001

    if arg0 == "":
        print("Usage  -broad Ref1 Ref2 Ref3 Ref4 -aparam1=RefA \n")
        print("       Ref1: Filename\n")
        print("       Ref2: column index for x (def.1)\n")
        print("       Ref3: column index for y (def.2)\n")
        print("       Ref4: sigma (def.1)\n")
        print("       aparam1: truncation limit to zero (def. tlimit\)\n")
        exit()
    else:
        inpf=arg0
    #print(arg0)
    if arg1 == "":
        rowx=1
    else:
        rowx=arg1

    if arg2 == "":
        rowy=2
    else:
        rowy=arg2

    if arg3 == "":
        sigma=1
    else:
        sigma=arg3

    if not histdef == "":
        tlimit=histdef
#    print "tlimit $tlimit\n";                                                                                                                                                              
    f =open(inpf,'r')
    lines=f.readlines()
    for i in lines:    
        i = i.rstrip('\n')
    f.close()
    n=0
    sumy1=0
    xval=[]
    yval=[]
    #print(rowx)
    for i in lines:
        field=i.split()
        xval.append(field[rowx-1])
        yval.append(field[rowy-1])
        sumy1=sumy1+float(yval[n])
        n = n+1
    #print(xval)

    g= open("out.broad",'w')
    x=0
    try :
        sumy2
    except:
        sumy2=0
    while x<=len(lines)-1:
        sumy=0
        i=1
        while i<=len(lines)-1:
            sumy=sumy+(float(yval[i])/(math.sqrt(2*math.pi*sigma**2)))*math.exp(-((float(xval[i])-float(xval[x]))**2)/(2*sigma**2))*(float(xval[i])-float(xval[i-1]))
            i = i + 1
        sumy2=sumy2+sumy
        ssi=""
        ssi=abs(sumy)
        if ssi < tlimit:
            sumy=0
        
        g.write("{} {}\n".format(xval[x],sumy))
        x = x +1
    g.close()
    
    return "inpf","out.broad"

#    print "check sum of area  inp\:$sumy1 out\:$sumy2\n";                                                                                                                                  
#    print "out\.broad is outputted\n";                                                                                                                                                    
         


def dfmake(arg0,arg1,arg2,arg3,arg4,arg5,arg6,histdef,aparam2): #----------------------------------------------------------------------------------                                                                        

    if arg6 == "":
        print("Usage\: -dfmake Ref1 Ref2 Ref3 Ref4 Ref5 Ref6 Ref7 -aparam1=RefA -aparam=RefB\n")
        print(" Ref1\: data filename\n")
        print(" Ref2\: column number of x-axis\n")
        print(" Ref3\: column number of y-axis (if 0, count frequency instead)\n")
        print(" Ref4\: minimum value of x-axis\n")
        print(" Ref5\: maximum value of x-axis\n")
        print(" Ref6\: bin number\n")
        print(" Ref7\: broadning factor\n")
        print(" -aparam1=RefA\: Normarization factor (def.=1)(Distribution function will be devided by RefA)\n")
        print(" -aparam2=RefB\: Normarization factor (Total area of the histgram will be RefB)\n")
        print("              \:  RefB = Number of lines in Ref1 file, if RefB=line.\n")
        print(" -aparam3=RefC\: if T, Ref2 will be converted to absolute value\n")
#       print(" -aparam2=RefB\: Truncation limit to zero (def. = $tlimit\)\n")                                                                                         
        print(" -o=(filename)\: output filename (def: out.distfunc)\n")
        exit()
    try:
        ofile
    except:
        ofile = ''

    if ofile == "":
        ofile="out.distfunc"
    
    IN= open(arg0)
    datalines=IN.readlines()
    for i in datalines:
        i=i.rstrip('\n')
    #print(datalines)
    
    histdef=1

    clmx=arg1 
    clmy=arg2
    minx=arg3 
    maxx=arg4
    binnum=arg5
    sigma=arg6
    
    if histdef == "":
        norm=1
    else:
        norm=histdef
    histdef = ""

    n=0
    field=[]
    orgx=[]
    orgy=[]
    
    for j in datalines:
        field =j.split()
        if bool(re.search(r"^\#",field[0])):
            continue

        orgx.append(field[int(clmx)-1]) #$orgy[$n]=$field[$clmy-1];
        #print(orgx)                                                                                                         
        if not clmy == 0:
            orgy.append(field[int(clmy)-1])
        else:
            orgy.append(1)
        if bool(re.search(r'[T]',aparam3,re.I)):
            orgy[n]=abs(float(orgy[n]))
        n = n+1
    IN.close()

    totaldata=n-1

    int_=[]
    i=0
    bin0=[]
    bin1=[]
    while i<= int(binnum)-1:
        bin0.append((float(maxx)-float(minx))/float(binnum)*i+float(minx))
        bin1.append((float(maxx)-float(minx))/float(binnum)*(i+1)+float(minx))
        j=0
        while j<=totaldata:
            if float(orgx[j]) >= float(bin0[i]) and float(orgx[j]) < float(bin1[i]):
                while len(int_)-1 < i:
                    int_.append('')
                if int_[i] == '':
                    int_[i]=0
                int_[i]=(float(int_[i])+float(orgy[j]))
            j = j+1
        i = i +1
    
    j=open(ofile,'w')

    if bool(re.search(r'[line]',str(aparam2),re.I)):
        norm=totaldata+1

    elif aparam2 != '':
        if aparam2 > 0:
            sumint=0
            i=0
            while i<=binnum-1:
                sumint=sumint+int_[i]
                i = i + 1
            norm=sumint/aparam2
    
    #print(int_)
                #    print "normalization factor\: $norm\n";                                                                                                                           
    i=0
    while i<= int(binnum)-1:
        while len(int_)-1 < i:
            int_.append('')
        if int_[i] == "":
            int_[i]=0
        int_[i]=int_[i]/norm
        #       $inti=abs($int[$i]);                                                                                                                                           
        #       print "inti $inti\n";                                                                                                                                          
        #       if ($inti < $tlimit){$int[$i] = 0;}                                                                                                               
        j.write("{}  {}\n".format(bin0[i],int_[i]))
        i = i +1
    j.close()
    if sigma !='sigma':
        if float(sigma) > 0:
            arg0 = ofile
            arg1 = 1
            arg2 = 2
            arg3 = sigma
            histdef = ""
            X =broad(arg0,arg1,arg2,arg3,arg4,arg5,arg6,histdef)
            dummy=os.system('mv {} {}'.format(X[1],ofile))
    
    return arg0,arg1,arg2,arg3,arg4,arg5,arg6,ofile
#    print "out.distfunc is outputed.\n";                                                                                                                              



def readposcar(inputf):     #------------------------------------------------------ read poscar (VASP)                                                                                              
     
    POSCAR=[] 
    linesPOSCAR=[] 
    ions=[] 
    totalion=0 
    headPOSCAR5=[]
    
    
    dir_list =os.listdir(os.getcwd())
    for ii in dir_list:
        if ii == inputf:
            with open(inputf,'r')as IN:     
                POSCAR = IN.readlines()     
    orgPOSCAR=POSCAR     
     
     
    k=0                #--------------------------- Remove head-spaces                         
    for l in POSCAR:
        l = re.sub('\s^', '', l)
        POSCAR[k]=l      
        k = k+1      

    slab_fb=[]
    slab_fr=[]
    while len(POSCAR) -1 < 5:
        POSCAR.append('')
    latmul=POSCAR[1]     
    a1 = a2 = a3=POSCAR[2].split()     
    b1 = b2 = b3=POSCAR[3].split()     
    c1 = c2 = c3=POSCAR[4].split()     
    fb_fr = POSCAR[4].split()
    while len(fb_fr)-1 < 4:
        fb_fr.append('')     
    slab_fb = fb_fr[3] 
    slab_fr = fb_fr[4]     
     
    _=POSCAR[5] 
    fieldline5=_.split()
    
    chr5 = fieldline5[0][:2]

    if bool(re.search('[A-Za-z]',chr5)):     
        atomline_ver5=POSCAR[5]     
        ver5name=fieldline5     
        ver5add=1 
        DorC=1     
    else:     
        ver5add=0 
        DorC=0     
        ver5name=[]     
         
    i = 0
    while i <= 6+ver5add:
        headPOSCAR5.append(POSCAR[i])
        i = i + 1
        
    i=5+ver5add     
    ions = POSCAR[i].split()     
     
     
    _=POSCAR[6+ver5add] 
    fieldline6=_.split()     
    chr6=fieldline6[0][:2]     
     
    if bool(re.search('d',chr6, re.I)) or bool(re.search('c',chr6,re.I)):     
        selectiveadd=0     
        DorCchr=chr6     
    elif re.search('s',chr6,re.I):     
        selectiveadd=1     
        _=POSCAR[6+ver5add+selectiveadd]     
        fielddorc=_.split()     
        DorCchr=fielddorc[0][:2]     
         
     
    if DorCchr == "d":     
        DorCchr ="D"     
    if DorCchr == "c":     
        DorCchr ="C"     
     
    i = 0
    indexions=[]
    indexionsend=[]
    try:
        totalion
    except:
        totalion=0
    while i <= len(ions)-1 :     
        indexions.append(totalion+1)
        totalion=totalion+int(ions[i])     
        indexionsend.append(totalion)     
        i = i + 1
     
#    $indexions[$i]=$totalion+1                                                                  
    i=0 
    ii=0 
    headPOSCAR=[]     
    while i<=8:                     # $headPOCAR[n: 0-7] referes standard header of POSCAR v.4                                                                                   
        if not i == 5 and not ver5add == 1 or not i == 6+ver5add and not selectiveadd == 1:              
            ii = ii + 1     
            headPOSCAR.append(POSCAR[i])     
        i = i + 1

     
    counttotalion = 0
    try:
        nototalion
    except:
        nototalion='F'    
    if nototalion == "T":
        i = 7+ver5add+selectiveadd
        while i <= len(POSCAR)-1:
            _= POSCAR[i] 
            field = _.split()
            if field[0] == "":
                break
            counttotalion =counttotalion + 1     
            i = i + 1

        totalion=counttotalion     
        ions=totalion     
         
     
    ii=0
    i = 7+ver5add+selectiveadd
    while i <= totalion+6+ver5add+selectiveadd:      
        ii = ii + 1     
        linesPOSCAR.append(POSCAR[i])     
        i = i + 1
                                           #--- from $linesPOSCAR[1]                      \     
                                                                                                 
                                            #--- $linesPOSCAR[n] refers coordinates of ions #n          
    coodPOSCAR=[]
    coodX=[]
    coodY=[]
    coodZ=[]
    ii=0
    TorF=[]
    coodAtom=[]
    coodAtom2=[]
    coodAtom3=[]
    coodComment=[]
    coodSiteindex=[]
    coodocp=[]
    coodchg=[]
    coodeffchg=[]
    coodCorA=[]
    site=[]
    Idx1=[]
    Idx2=[]
    ionslabel=[]
    i = 7+ver5add+selectiveadd
    while i<=totalion+6+ver5add+selectiveadd:      
        ii = ii + 1
        _=POSCAR[i]     
        coodPOSCAR.append(POSCAR[i])     
        splitline=_.split()      
        if DorCchr == "D":
            k = 0
            while k<=2:     
                if splitline[k] >= 1.0:
                    splitline[k]=splitline[k]-1      
                if splitline[k] <  0.0:
                    splitline[k]=splitline[k]+1      
                k = k +1
        
        coodX.append("%.10f" % float(splitline[0]))      
        coodY.append("%.10f" % float(splitline[1]))      
        coodZ.append("%.10f" % float(splitline[2]))
        #print print print
        if selectiveadd == 1:
            while len(TorF)-1 < 2:
                TorF.append('')
            TorF[0]=splitline[3]      
            TorF[1]=splitline[4]      
            TorF[2]=splitline[5]      
             
        while len(coodAtom)-1 < ii:
            coodAtom.append('')

        coodAtom[ii]=splitline[3+selectiveadd*3]


        if coodAtom[ii] == "\!":
            coodAtom[ii]=splitline[4]
        
        try:
            allowva
        except:
            allowva='F'
        
        if not bool(re.search(r'[T]',allowva,re.I)):
            if bool(re.search(r'([A-Za-z]*)[0-9\+\-\_]',coodAtom[ii])):
                search=re.search(r'([A-Za-z]*)[0-9\+\-\_]',coodAtom[ii])
                ate1 = search.group(1)
                coodAtom[ii]=ate1

        while len(splitline)-1 < 4+selectiveadd*3 :
            splitline.append('')
        
        coodAtom2.append(splitline[4+selectiveadd*3])

        while len(splitline)-1 < 5+selectiveadd:
            splitline.append('')
        
        coodAtom3.append(splitline[5+selectiveadd*3])      
        coodComment.append('')
        try :
            siteindex
        except:
            siteindex=1  
            
        if ii > indexionsend[siteindex]:
            siteindex = siteindex+1
        coodSiteindex.append(siteindex+1)      
        
        sln = 4

        while len(coodComment)-1 < ii:
            coodComment.append('')
        #print(coodComment)
        while sln <= len(splitline) - 1:
            coodComment[ii]="{}".format(coodComment[ii])+" {}".format(splitline[sln])   
            sln=sln+1
        for splitline in splitline:     
            if bool(re.search(r'OC_([0-9\.]*)',_)):
                search=re.search(r'OC_([0-9\.]*)',_)
                ate1 = search.group(1)
                coodocp.append(ate1) 
            if bool(re.search(r'Q_(.*)',_)):
                search=re.search(r'Q_(.*)',_)
                ate1 = search.group(1)
                coodchg.append(ate1) 
                existchg="T" 
                if existeffchg != "T" :
                    ate1 = search.group(1)
                    coodeffchg.append(ate1)
            if bool(re.search(r'QEff_(.*)',_)):
                search=re.search(r'QEff_(.*)',_)
                ate1 = search.group(1)
                coodeffchg.append(ate1) 
                existeffchg="T"
            if bool(re.search(r'CorA_(.*)',_)):
                search=re.search(r'CorA_(.*)',_)
                ate1 = search.group(1)
                coodCorA.append(ate1) 
            if bool(re.search(r'Site_(.*)',_)):
                search=re.search(r'Site_(.*)',_)
                ate1 = search.group(1)
                Site.append(ate1)
            if bool(re.search(r'Idx1_(.*)',_)):
                search=re.search(r'Idx1_(.*)',_)
                ate1 = search.group(1)
                Idx1[ii]=ate1
            if bool(re.search(r'Idx2_(.*)',_)):
                search=re.search(r'Idx2_(.*)',_)
                ate1 = search.group(1)
                Idx2[ii]=ate1      
        i = i + 1
        sln = sln + 1
        
    
    loci = 1
    while loci <= totalion:     
        if coodocp[loci] == "":     
            coodocp[loci] =1
        loci = loci + 1
             
     
    indexsite=1   
    for isite  in ions:     
        ionslabel.append(coodAtom[int(indexsite)])      
        indexsite = float(indexsite)+float(isite)      
         
     
    jelem=0
    outcomp=[]
    for labelelem  in ionslabel:     
        outcomp.append("{}".format(outcomp)+"{}".format(labelelem)+"{}".format(ions[jelem]))      
        jelem = jelem + 1      
     
#    &coodsort                                                                                   
     
    if DorC == 1:     
        for POSCAR in POSCAR:     
            if not n == 5:     
                POSCAR1.append(_)     
            n = n + 1      
        POSCAR=POSCAR1  
        POSCAR1=[]      
         
     
     
    ii=0 
    iii=-1      
    for numa in ions:     
        iii = iii + 1
        i = 1
        while i<= int(numa):     
            ii = ii + 1      
            if not bool(re.search(r'^[A-Z]',coodAtom[ii])):     
                coodAtom[ii] = ver5name[iii]      
            i = i + 1
                     
                 
    if ver5add == 1:
        ionslabel=ver5name      
    
    try:
        chgcaron
    except:
        chgcaron='F'

    if bool(re.search(r'[T]',chgcaron, re.I)):
        #---------------  CHGCAR -----------------------------------------------                                                                                        
     
        #print "ver5add $ver5add\n"                                                              
        #for ($i=1 $i<=10 $i++){                                                                 
        #    print "POSCAR$i\:   $POSCAR[$i]"                                                    
        #}                                                                                       
     
        _=POSCAR[totalion+6+2]      
#       print "$totalion\n"                                                                      
#       print "_ $_\n"                                                                           
        headchgspin=_      
        mesh = mesh.split()      
        totalmesh=mesh[0]*mesh[1]*mesh[2]      
        print ("{}-{} {}  ({}  {}  {})\n".format(total,mesh,totalmesh,mesh[0],mesh[1],mesh[2]))      
     
        count=0  
        chgsum=0      
        maxchg = -9999  
        absmaxchg = -9999      
        minchg = +9999  
        absminchg = +9999      
     
        maxspn = -9999  
        absmaxspn = -9999      
        minspn = +9999  
        absminspn = +9999      
        i = totalion + 6 + 3
        while i <= len(POSCAR) - 1:      
            if POSCAR[i] == "":
                break
            if bool(re.search('[augm]',POSCAR[i])):
                break     
            _=POSCAR[i]      
            field=_.split()      
            for val in field:     
                chgall[count]=val      
                chgsum=chgsum+val      
                count = count + 1      
                if maxchg < val:
                    maxchg =val      
                if absmaxchg < abs(val):
                    absmaxchg=val      
                if minchg > val:
                    minchg =val      
                if absminchg > abs(val):
                    absminchg=val 
            i = i + 1 
            if count >= totalmesh:
                break     
     
        t1mesh=count      
        count=0
        j = i
        while j <= len(POSCAR) - 1:      
            if POSCAR[j] == "headchgspin":     
                if not bool(re.search(r'[augm]',POSCAR[j])):     
                    j = j + 1      
                    _=POSCAR[j]      
                    field=_.split()      
                    for val in field:     
                        spinall[count]=val      
                        spinsum=spinsum+val      
                        count = count + 1      
     
                        if maxspn < val:
                            maxspn =val      
                        if absmaxspn < abs(val):
                            absmaxspn=val      
                        if minspn > val:
                            minspn =val      
                        if absminspn > abs(val):
                            absminspn=val      
     
                         
                    if count >= totalmesh:
                        break     
                     
                break      
            j = j + 1
                 
             
     
        print("Sum-charge chgsum  (total mesh t1mesh)\n")      
        print("  Maximum of charge density maxchg\n")      
        print("  Minimum of charge density minchg\n")      
        print("  Abs(Max-Chg) absmaxchg\n")      
        print("  Abs(Min-Chg) absminchg\n")      
        if count > 10:     
            print("Sum-spin spinsum  (total mesh count)\n")      
            print("  Maximum of spin density maxspn\n")      
            print("  Minimum of spin density minspn\n")      
            print("  Abs(Max-Spin) absmaxspn\n")      
            print("  Abs(Min-Spin) absminspn\n")      
             
        path = 'POTCAR'
        if  os.path.isfile(path):    
            place_of_potcar="./POTCAR"      
            get_posatom()      
            indb=0      
            for numa in ions:
                ind=1
                while ind <= len(numa) - 1:      
                    coodAtom[inda]=posatom[indb]      
                     
                indb = indb + 1    
                 
             
     
        count=0      
        z = 1
        y = 1
        x = 1
        while z<=mesh[2]:      
            while y<=mesh[1]:      
                while x<=mesh[0]:      
                    chgxyz[x][y][z]=chgall[count]      
                    spinxyz[x][y][z]=spinall[count]      
                    count = count + 1
                    x = x + 1
                y = y + 1
            z = z + 1
     
        meshall=count      
        chgsum=chgsum/meshall      
        spinsum=spinsum/meshall      
         
    return ions,coodAtom,totalion,coodocp
     
def compdescript(defelem):

                   
    if reference == "":     
        print("Usage \: vaspfiles -compdescript3 Arg1 Arg2 -sparam=\n")     
        print("  Arg1\: Chemical Formula (e.g. Li1.0Co1.0O2.0) or POSCAR.names formatted file (if -sparam=T)\n")     
#        print "  Arg2\: sigma (def. range/25)\n"                                                                        
        print(" -sparam=T use POSCAR.names file instead of direct input.\n")     
        print(" -aparam1=filename\: histogram format file\n")     
        print(" -aparam2=sigma\: sigma value for Gaussian broadning\n")          
#       print "   Arg3=en  electronegativity table\n"                                                                    
#       print "   Arg3=mp  melting point table\n"                                                                        
#       print "   Arg3=mn  mendeleev number table\n"                                                                     
        print(" Histogram format\: 1:Prop (eg. AN, EN)  2: 3  3: prop_min  4: prop_max  5: bin_num  6: gaussian_sigma\n")     
        return
     
     
     
    filename=["AN","EN","MP","PN","PG","MN","AW","AR","IR","CoR","CrR","spdf"]     
    ini_r=[["AN","EN","MP","PN","PG","MN","AW","AR","IR","CoR","CrR","spdf"]]     

    ini_r = {"AN":"AN","EN":"EN","MP":"MP","PN":"PN","PG":"PG","MN":"MN","AW":"AW","AR":"AR","IR":"IR","CoR":"CoR","CrR":"CrR","spdf":"spdf"}
    ini_r["AN"] = {k:k for k in range(7)}
    ini_r["EN"] = {k:k for k in range(7)}
    ini_r["MP"] = {k:k for k in range(7)}
    ini_r["PN"] = {k:k for k in range(7)}
    ini_r["PG"] = {k:k for k in range(7)}
    ini_r["MN"] = {k:k for k in range(7)}
    ini_r["AW"] = {k:k for k in range(7)}
    ini_r["AR"] = {k:k for k in range(7)}
    ini_r["IR"] = {k:k for k in range(7)}
    ini_r["CoR"] = {k:k for k in range(7)}
    ini_r["CrR"] = {k:k for k in range(7)}
    ini_r["spdf"] = {k:k for k in range(7)}
  
    for i in filename:
        ini_r[i][1]=2

    for k in filename:
        ini_r[k][2]=3

    ini_r['AN'][4]=120
    ini_r['AN'][3]=-0.2*ini_r['AN'][4]
    ini_r['AN'][5]=102
    ini_r['AN'][6]=0

    ini_r["EN"][4]=5.0
    ini_r["EN"][3]=-0.2*ini_r["EN"][4]
    ini_r["EN"][5]=50
    ini_r["EN"][6]=0.2
            
    ini_r["MP"][4]=5000
    ini_r["MP"][3]=-0.2*ini_r["MP"][4]
    ini_r["MP"][5]=50
    ini_r["MP"][6]=ini_r["MP"][4]/25
     
    ini_r["PN"][4]=10
    ini_r["PN"][3]=-0.2*ini_r["PN"][4]
    ini_r["PN"][5]=20
    ini_r["PN"][6]=ini_r["PN"][4]/25
    
    ini_r["PG"][4]=20
    ini_r["PG"][3]=-0.2*ini_r["PG"][4]
    ini_r["PG"][5]=20
    ini_r["PG"][6]=ini_r["PG"][4]/25
                            
    ini_r["MN"][4]=120
    ini_r["MN"][3]=-0.2*ini_r["MN"][4]
    ini_r["MN"][5]=102
    ini_r["MN"][6]=ini_r["MN"][4]/25
                        
    ini_r["AW"][4]=400
    ini_r["AW"][3]=-0.2*ini_r["AW"][4]
    ini_r["AW"][5]=50
    ini_r["AW"][6]=ini_r["AW"][4]/25
                                       
    ini_r["AR"][4]=3.2
    ini_r["AR"][3]=-0.2
    ini_r["AR"][5]=50
    ini_r["AR"][6]=ini_r["AR"][4]/25

    ini_r["IR"][4]=3.2
    ini_r["IR"][3]=-0.2
    ini_r["IR"][5]=50
    ini_r["IR"][6]=ini_r["IR"][4]/25

    ini_r["CoR"][4]=3.2
    ini_r["CoR"][3]=-0.2
    ini_r["CoR"][5]=50
    ini_r["CoR"][6]=ini_r["CoR"][4]/25
    
    ini_r["CrR"][4]=3.2
    ini_r["CrR"][3]=-0.2
    ini_r["CrR"][5]=50
    ini_r["CrR"][6]=ini_r["CrR"][4]/25
                
    ini_r["spdf"][4]=5
    ini_r["spdf"][3]=1
    ini_r["spdf"][5]=4
    ini_r["spdf"][6]=0
    
    if aparam2 != "":     
        for ff in filename:     
            ini_r[ff][6]=aparam2
     
     
    if histdef != "":     
        f= open("histdef.dat")
        linesem = f.readlines()
        for line in linesem:
            line.rstrip('\n')
        
        for k in linesem:
            print ("{}".format(k))     
            field=k.split()  
            if bool(re.search('^\#', field[0])):
                continue     
            label = field[0]
            ie =1
            while ie<=6 :
                ini_r[label][ie]=field[ie]
                ie = ie+1
        f.close()

    aparam=[]
    #aparam3=1
    numatom={}
    if sparam != '':
        if bool(re.search(r'T', sparam, re.I)):     
            YY=defvasptable()     
            inputf=reference[0]     
            XX =readposcar(inputf)     
            i=1     
     
            if XX[1][1] == "":     
                print("POTCAR is read......\n")
                iit=os.system('grep VRH POTCAR')     
                vrh = vrh.split('/\n/,{}'.format(iit))     
                for k in vrh:     
                    k = re.match('(VRHFIN =)(.+)(\:)', s)     
                    atom[i] = 2     
                    numatom[YY[3][atom[i]]]=ions[i-1]     
                    i = i + 1     
                 
            else:     
                ia=1
                while ia<=XX[2]:
                    try:
                        numatom[int(YY[3][XX[1][ia]])]
                    except:
                        numatom[int(YY[3][XX[1][ia]])]=0
                    numatom[int(YY[3][XX[1][ia]])]=float(numatom[int(YY[3][XX[1][ia]])])+1*float(XX[3][ia]) 
                    ia = ia + 1

            f =open("out.composition_label",'w')
            i =1
            while i<=103 :
                try:
                    numatom[i]
                except:
                    numatom[i]=''
                #print(numatom[i])
                if numatom[i] != '':
                    if int(numatom[i]) >0:
                        while len(YY[4])-1 < i:
                            YY[4].append('')
                        try:
                            compname
                        except:
                            compname=''
                        f.write("{}{}".format(YY[4][i],numatom[i]))
                        compname="{}{}{}".format(compname,YY[4][i],numatom[i])      
                    f.write("\n")
                i = i + 1
            f.close()      
            
            print("composition: {}\n".format(compname))      
            reference[0]=compname      
         
     
    incomp=reference[0]     
    x =compconv() 
    y =readdefelem
    reference.append('')
    if reference[1] == "":
        sigma=0.0 
    else:
        sigma=reference[1] 

    #    @filename=("AN","EN","MP","PN","PG","MN","AW","AR","IR","CoR","CrR")                                                                                                                                                                        
    
    PGBox={k: k for k in range(19)}
    def_PGBox={l: l for l in range(19)}
    for ii in range(19):
        PGBox[ii]=dict(def_PGBox)
    
    for FN in filename:     
        if FN == "PG":
            i=0
            while  i<=len(x[0]) - 1:
                ENatom1=x[2][i]      
                inty1=x[0][i]/x[3]
                j = i + 1
                while  j<=len(x[0]) - 1:    
                    if i >= j:
                       continue
                    ENatom2=x[2][j]      
                    inty2=x[0][j]/x[3]
                    if int(y(ENatom1,defelem)[FN][ENatom1]) > int(y(ENatom2,defelem)[FN][ENatom2]):     
                        kk1=y(ENatom2,defelem)[FN][ENatom2]  
                        kk2=y(ENatom1,defelem)[FN][ENatom1]
                    else:     
                        kk2=y(ENatom2,defelem)[FN][ENatom2]  
                        kk1=y(ENatom1,defelem)[FN][ENatom1]
                    
                    PGBox[int(kk1)][int(kk2)]=inty1+inty2+PGBox[int(kk1)][int(kk2)]
                    j = j + 1
                i= i + 1
                          
            g=open("out.PGmatrix",'w')       
            PGmtrx={}
            k=-1
            i=1
            while i<=18:
                j=i
                while j<=18:     
                    k = k + 1      
                    PGmtrx[int(k)]=PGBox[int(i)][int(j)]  
                    g.write("{}  {:.6f}  (pair {} - {})\n".format(k,PGmtrx[k],i,j))
                    j = j + 1
                i = i + 1
            g.close()      
       
        PGBox={k: k for k in range(19)}
        def_PGBox={l: 0 for l in range(19)}
        for ii in range(19):
            PGBox[ii]=dict(def_PGBox)
        if FN == "PN":
            i = 0
            while i<=len(x[0]) -1:     
                ENatom1=x[2][i]      
                inty1=x[0][i]/x[3]
                j=i+1
                while  j<=len(x[0])-1:      
                    if i >= j:
                        continue     
                    ENatom2=x[2][j]      
                    inty2=x[0][j]/x[3]      
     
                    if y(ENatom1,defelem)[FN][ENatom1] > y(ENatom2,defelem)[FN][ENatom2]:     
                        kk1=y(ENatom2,defelem)[FN][ENatom2]  
                        kk2=y(ENatom1,defelem)[FN][ENatom1]      
                    else:     
                        kk2=y(ENatom2,defelem)[FN][ENatom2]  
                        kk1=y(ENatom1,defelem)[FN][ENatom1]      
                    j = j + 1
                    PGBox[int(kk1)][int(kk2)]=inty1+inty2+PGBox[int(kk1)][int(kk2)]      
                i = i + 1
            
     
     
            h = open("out.PNmatrix",'w')      
            k=-1
            i=1
            PGmtrx={m:0 for m in range(8)}
            while i<=7 :
                j=i
                while j<=7:      
                    k = k + 1      
                    PGmtrx[k]=PGBox[i][j]      
                    h.write("{}  {:1.6f}  (pair {} - {})\n".format(k,PGmtrx[k],i,j))
                    j = j + 1
                i = i + 1     
                
            h.close()      
             
     
        #PGBox=[]  
        #PGmtrx=[]      
        if FN == "spdf":     
            i=0
            while i<=len(x[0])-1:      
                ENatom1=x[2][i]      
                inty1=x[0][i]/x[3]
                j=i+1
                while j<=len(x[0])-1:     
                    if i >= j:
                        continue   
                    ENatom2=x[2][j]      
                    inty2=x[0][j]/x[3]      
     
                    if y(ENatom1,defelem)[FN][ENatom1] > y(ENatom2,defelem)[FN][ENatom2]:     
                        kk1=y(ENatom2,defelem)[FN][ENatom2]
                        kk2=y(ENatom1,defelem)[FN][ENatom1]      
                    else:     
                        kk2=y(ENatom2,defelem)[FN][ENatom2]  
                        kk1=y(ENatom1,defelem)[FN][ENatom1]
                        
                    if  kk1 == 'nodata' or kk2 == 'nodata':
                        PGBox[11][11]=inty1+inty2+PGBox[11][11]
                    else:
                        PGBox[int(kk1)][int(kk2)]=inty1+inty2+PGBox[int(kk1)][int(kk2)]
                    j = j+1
                i = i + 1
                 
            l = open("out.SPDFmatrix",'w')      
     
            k=-1
            i=1
            while i<=4:
                j=i
                while j<=4:     
                    k = k + 1      
                    PGmtrx[k]=PGBox[i][j]      
                    l.write("{}  {:1.6f}  (pair {} - {})\n".format(k,PGmtrx[k],i,j))      
                    j = j + 1
                i = i + 1

                 
            l.close()      
             
     
     
        f = open('list.{}-distfunc'.format(FN), 'w')
        f.write("# Atom-Atom characteristics Concentration\n")      
        i=0
        while  i <= len(x[0])-1:
            ENatom=x[2][i]
            inty=x[0][i]/x[3]
            if y(ENatom,defelem)[FN][ENatom] != 'nodata':      
                f.write("{} {:1.6f} {:1.6f}\n".format(ENatom,float(y(ENatom,defelem)[FN][ENatom]),inty))
            else:
                f.write("{} {} {:1.6f}\n".format(ENatom,y(ENatom,defelem)[FN][ENatom],inty))
            i = i + 1
        f.close()     
        if aparam3 == 'T':
            m=open("list.{}{}-distfunc".format(FN,FN),'w')      
            m.write("# Atom*Atom characteristics^2 Concentration^2\n")      
            i=0
            while i<=len(x[0])-1:     
                ENatom1=x[2][i]       
                inty1=x[0][i]/x[3]
                j=i+1
                while j <= len(x[0]) - 1 :     
                    if i >= j:
                        continue     
                    ENatom2=x[2][j]      
                    inty2=x[0][j]/x[3]    
                    m.write("{}-{} {:1.6f} {:1.6f}\n".format(ENatom1,ENatom2,math.log(float(y(ENatom1,defelem)[FN][ENatom1]))+math.log(float(y(ENatom2,defelem)[FN][ENatom2])),inty1+inty2))
                    j = j + 1
                i = i + 1
            m.close()      

            o = open("list.{}-{}-distfunc".format(FN,FN),'w')      
            o.write("# ABS(Atom-Atom) characteristics Concentration\n")      
            i=0
            while i<=len(x[0])-1:      
                ENatom1=x[2][i]      
                inty1=x[0][i]/x[3]
                j=i+1
                while j <= len(x[0])-1:     
                    if i >= j:
                        continue     
                    ENatom2=x[2][j]      
                    inty2=x[0][j]/x[3]    
                    o.write("{}-{} {:1.6f} {:1.6f}\n".format(ENatom1,ENatom2,abs(float(y(ENatom1,defelem)[FN][ENatom1])-float(y(ENatom2,defelem)[FN][ENatom2])),inty1+inty2))      
                    j = j + 1
                i= i + 1
     
            o.close()   




        list_reference = reference[0].split()
        list_reference[0]="list.{}-distfunc".format(FN)
        #reference[0]="list.FN-distfunc"      
        list_reference.append(float(ini_r[FN][1]))      
        list_reference.append(float(ini_r[FN][2]))      
        list_reference.append(float(ini_r[FN][3]))      
        list_reference.append(float(ini_r[FN][4]))      
        list_reference.append(float(ini_r[FN][5]))
        if ini_r[FN][6] =='sigma':
            list_reference.append(ini_r[FN][6])
        else:
            list_reference.append(float(ini_r[FN][6]))    
        #       if ($aparam[3] == 0){$reference[6]=0 } 
        #if aparam3 == 0:
        #    list_reference[6]=0
        outcompfile=[]
        if aparam3 == 'T':
            
            #print(list_reference)
            z = dfmake(*list_reference,histdef,aparam2)
            dummy=os.system('mv {} out.{}_distfunc'.format(z[7],FN))     
            outcompfile.append("out.{}_distfunc".format(FN))

            list_reference[0]="list.{}{}-distfunc".format(FN,FN)     
            list_reference[1]=float(ini_r[FN][1])     
            list_reference[2]=float(ini_r[FN][2])     
            list_reference[3]=-0.2*math.log(float(ini_r[FN][4]))     
            list_reference[4]=2*math.log(float(ini_r[FN][4]))     
            list_reference[5]=float(ini_r[FN][5])
            list_reference[6]=2*math.log(float(ini_r[FN][4]))/float(ini_r[FN][5])      
        
            #if aparam3 == 0:        
            #    list_reference[6]=0
            #print(list_reference)
            z = dfmake(*list_reference,histdef,aparam2)      
            dummy=os.system('mv {} out.{}{}_distfunc'.format(z[7],FN,FN))      
            outcompfile.append("out.{}{}_distfunc".format(FN,FN))   

            list_reference[0]="list.{}-{}-distfunc".format(FN,FN)      
            list_reference[1]=float(ini_r[FN][1])       
            list_reference[2]=float(ini_r[FN][2])      
            list_reference[3]=float(ini_r[FN][3])      
            list_reference[4]=float(ini_r[FN][4])      
            list_reference[5]=float(ini_r[FN][5])
            if ini_r[FN][6] =='sigma':
                list_reference[6]=ini_r[FN][6]
            else:    
                list_reference[6]=float(ini_r[FN][6]) 

            #print('ref6 1st {}\n'.format(list_reference[6]))
            #if aparam3 == 0:
            #list_reference[6]=0
            #print('ref6 2nd {}\n'.format(list_reference[6]))
            #print(list_reference)
            z = dfmake(*list_reference,histdef,aparam2)      
            dummy=os.system('mv {} out.{}-{}_distfunc'.format(z[7],FN,FN))      
            outcompfile.append("out.{}-{}_distfunc".format(FN,FN)) 

        if FN=='AN':
            n= open("out.compdescript",'w')
        
        for fnc in outcompfile:     
            IN=open(fnc)
            lines=IN.readlines()
            for IN_v in lines:
                IN_v.rstrip('\n')
                    
            num=0
            for ik in lines:
                ik=ik.rstrip('\n')
                num = num + 1
                n.write("{}  {}_{}\n".format(ik,fnc,num))
        
     
compdescript(defelem)
     
lsafter=str(sp.check_output("ls -alh --full-time",shell=True))
ls_list =lsafter.split('\\n')
ls_list.pop(0)

for aft in ls_list:
    ruleout=0
    for bfr in ls_bflist:
        if aft == bfr:
            ruleout =ruleout+1
    aft=aft.split()
    if aft[-1] == '.'or aft[-1]=='..' or aft[-1]=="'":
        ruleout=ruleout+1
    if ruleout==0:
        print('Output file: {}'.format(aft[-1]))
     
