#! /usr/bin/env python

import re
import os
import numpy as np
import argparse
import math
import subprocess as sp
from subprocess import check_output
import pandas as pd

path=os.getcwd()

def compconv(reference): #-------------------------------------------------------------------------------                                                                         
#input incomp or reference[0]  eg Li0.1CoO2                                                                                                                          
#output  indexions[i] -> element                                                                                                                                     
#        ions[i]      -> number of ions                                                                                                                              
#        totalion      -> summation of ions                                                                                             
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
            i = i + 1
            if i>9999 :
                break
        j =j +1
    i=0
    
    for j in  elemcomp:                                                                                                                                         
        if  bool(re.search(r'^([A-Za-z]*)([0-9\.].*)',str(j))):
            search = re.search(r'^([A-Za-z]*)([0-9\.].*)',str(j))
            ate1 = search.group(1)
            ate2 = search.group(2)
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
        #if printoff != "T":

            #print("{} ".format(l))
            #print("{} ".format(nions[l]))

        totalion=totalion+nions[l]
        i = i + 1

    if printoff != "T":
        print ("\n")
        print ("total ions   {}\n".format(totalion))

    return ions,compositionprim,indexions,totalion

def readdefelem(atom,elem_csv): #------------------------------------------------------------------------------
    
    """
    IN =open(defelem,'r')
    defelem=IN.readlines()
    for line in defelem:
        line=line.rstrip('\n')
    IN.close()
    """
    #df_elem=df_csv.loc[0:101,'atomic_num.':'spdf-block']
    #print(elem_csv)
    #df_elem=elem_csv.loc[0:101,'atomic_num.':'spdf-block']
    desc_list=elem_csv.columns.tolist()
    elem_list=elem_csv.to_numpy().tolist()
    el_dic={}
    i=0
    i_list=[]
    for desc in desc_list:
        if desc=='element_symbol' or desc=='Endatastat.':
            pass
        else:
            i_list.append(i)
        i+=1
    for field in elem_list:
        defatom=field[1]
        element = '0'
        def_atom = {defatom : element}
        for desc in desc_list:
            if desc=='element_symbol' or desc=='Endatastat.':
                pass
            else:
                el_dic[desc] = dict(def_atom)
                i_list.append(i)
            i+=1

        desc_list=desc_list=[des for des in desc_list if des not in ['element_symbol','Endatastat.']]
        #print(desc_list)
        if atom == defatom:
            for i,desc in zip(i_list,desc_list):
                el_dic[desc][atom]=field[i]       

            return el_dic

def broad(df,arg3,df_hist,j_list_nise): #--------------------------------------------------------------- Broad                                                                                                        

    tlimit=0.000001

    if arg3 == "":
        sigma=1
    else:
        sigma=arg3

    if not df_hist == "":
        tlimit=df_hist                                                                                                                                                              
    #f =open(inpf,'r')
    #lines=f.readlines()
    #print(df)
    col_list2=df.columns.tolist()
    data_list2=[]
    for df_col2 in col_list2: 
        #print(df[df_col2])
        df_data2 = df[df_col2]
        #print(df_data2)
        data_list2.append(df_data2)
    
    n=0
    sumy1=0
    xval=[]
    yval=[]
    for data in data_list2:
        #xval.append(field[0][rowx-1])
        yval.append(data[0])
        sumy1=sumy1+float(yval[n])
        n = n+1

    
    for nise in j_list_nise:
        xval.append(nise)
    #print('xval dayo   :  {}'.format(xval))
    #print('yval dayo   :  {}'.format(yval))

    g_list_index=[]
    g_list_columns=[]
    x=0
    #print(data_list2)
    try :
        sumy2
    except:
        sumy2=0
    while x<=len(data_list2)-1:
        sumy=0
        i=1
        while i<=len(data_list2)-1:
            sumy=sumy+(float(yval[i])/(math.sqrt(2*math.pi*sigma**2)))*math.exp(-((float(xval[i])-float(xval[x]))**2)/(2*sigma**2))*(float(xval[i])-float(xval[i-1]))
            i = i + 1
        sumy2=sumy2+sumy
        ssi=""
        ssi=abs(sumy)
        if ssi < tlimit:
            sumy=0
        
        g_list_index.append("{}".format(sumy))
        #g_list_index.append("{}".format(sumy))
        g_list_columns.append('out_{}'.format(x))
        x = x +1
    g_dict=dict(zip(g_list_columns,g_list_index))
    g_df=pd.DataFrame(g_dict,index=['test'])
    
    return g_df


def dfmake(arg0,arg1,arg2,arg3,arg4,arg5,arg6,aparam2,aparam3,df,df_csv): #----------------------------------------------------------------------------------                                                                        

    """
    try:
        ofile
    except:
        ofile = ''

    if ofile == "":
        ofile="out.distfunc"
    """
    
    col_list=df.columns.tolist()
    data_list=[]
    for df_col in col_list:
        #print(df[df_col])
        df_data=df[df_col]
        data_list.append(df_data)
    #print(data_list)

    """
    IN= open(arg0)
    datalines=IN.readlines()
    for i in datalines:
        i=i.rstrip('\n')
    """

    hist_symbol = '_his'
    column_inc_hist = [column for column in df_csv.columns if hist_symbol in column]    

    df_hist=pd.DataFrame(columns=column_inc_hist)

    df_hist=1

    clmx=arg1 
    clmy=arg2
    minx=arg3 
    maxx=arg4
    binnum=arg5
    sigma=arg6
    
    if df_hist == "":
        norm=1
    else:
        norm=df_hist
    df_hist = ""

    n=0
    field=[]
    orgx=[]
    orgy=[]
    
    for j in data_list:
        field =j.str.split()
        orgx.append(field[0][int(clmx)-1])                                                                                                          
        if not clmy == 0:
            orgy.append(field[0][int(clmy)-1])
        else:
            orgy.append(1)
        if bool(re.search(r'[T]',aparam3,re.I)):
            orgy[n]=abs(float(orgy[n]))
        n = n+1
    

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
            if orgx[j] == 'nodata':
                j = j + 1
                j_df=''
                return arg1,arg2,arg3,arg4,arg5,arg6,j_df
            if float(orgx[j]) >= float(bin0[i]) and float(orgx[j]) < float(bin1[i]):
                while len(int_)-1 < i:
                    int_.append('')
                if int_[i] == '':
                    int_[i]=0
                int_[i]=(float(int_[i])+float(orgy[j]))
            j = j+1
        i = i +1
    j_list_index=[]
    j_list_columns=[]
    j_list_nise=[]
    if bool(re.search(r'[T]',str(aparam2),re.I)):
        norm=totaldata+1
    else:
        pass
    #elif aparam2 != '':
        #if aparam2 > 0:
        #    sumint=0
        #    i=0
        #    while i<=binnum-1:
        #        sumint=sumint+int_[i]
        #        i = i + 1
        #    norm=sumint/aparam2                                                                                                                
    i=0
    while i<= int(binnum)-1:
        while len(int_)-1 < i:
            int_.append('')
        if int_[i] == "":
            int_[i]=0
        if norm != 0:
            int_[i]=int_[i]/norm                                                                                                                     
        j_list_index.append("{}".format(int_[i]))
        j_list_nise.append('{}'.format(bin0[i]))

        #j_list_index.append("{}".format(int_[i]))
        j_list_columns.append('out_{}'.format(i))
        i = i +1
    j_dict=dict(zip(j_list_columns,j_list_index))
    j_df=pd.DataFrame(j_dict,index=['test'])
    
    if sigma !='T':
        if float(sigma) > 0:
            df = j_df
            arg3 = sigma
            df_hist = ""
            X =broad(df,arg3,df_hist,j_list_nise)
            j_df = X
    
    return arg1,arg2,arg3,arg4,arg5,arg6,j_df
     
def compdescript(aparam2,aparam3,reference):
     
    #ini_r=[["AN","EN","MP","PN","PG","MN","AW","AR","IR","CoR","CrR","spdf"]]     
    #csv_list=['atomic_num.','element_symbol','electronegativity','Endatastat.','melting_point','pricipal_quantum_num','group_of_element','mendeleev_num','Atomic_weight','Atomic_Radii','Ionic_Radii','Covalent_Radii','Crystal_Radii','spdf-block']

    df_csv=pd.read_csv('input_data.csv')
    elem_csv=pd.read_csv('DefElem.csv')
    
    filename=elem_csv.columns.tolist()
    filename=[des for des in filename if des not in ['element_symbol','Endatastat.']]

    ans_df_dict={}
    ans_df_hiki_dict={}
    ans_df_kake_dict={}
    df_dict={}
    df_kake_dict={}
    df_hiki_dict={}
      
    if aparam2 == 'T':     
        for ff in filename:
            dict_df_csv = df_csv.to_dict()
            #print('dict_df_csv',dict_df_csv)
            dict_df_csv[ff][5]=aparam2
            #print('dict_df_csv[ff][5]',dict_df_csv[ff][5])
            df_csv = pd.DataFrame(dict_df_csv)
            #df_csv.at[7,ff]=aparam2
     
    hist_symbol = '_his'
    column_inc_hist = [column for column in df_csv.columns if hist_symbol in column]
    
    h_list=[]
    
    if column_inc_hist != '':
        for h_symbol in column_inc_hist:
            h_symbol_list = h_symbol.split('_')
            #print(h_symbol_list)
            h_list.append(h_symbol_list[0])

        for col_h ,hh in zip(column_inc_hist,h_list):
            ie=0
            while ie <= 5:
                df_csv.at[ie,hh]=df_csv.at[ie,col_h]
                ie = ie + 1
    #print(df_csv['AN'])
    """        
    if histdef != "":     
        f= open(histdef)
        linesem = f.readlines()
        for line in linesem:
            line.rstrip('\n')
        
        for k in linesem:   
            field=k.split()  
            if bool(re.search('^\#', field[0])):
                continue     
            label = field[0]
            ie =3
            while ie<=8 :
                df.at[ie,label]
                df.at[label][ie]=field[ie]
                ie = ie+1
    """

    #aparam=[]
    #aparam3=1

    x =compconv(reference) 
    reference.append('')
    #if reference[1] == "":
    #    sigma=0.0 
    #else:
    #    sigma=reference[1] 

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
                    
                    if int(readdefelem(ENatom1,elem_csv)[FN][ENatom1]) > int(readdefelem(ENatom2,elem_csv)[FN][ENatom2]):     
                        kk1=readdefelem(ENatom2,elem_csv)[FN][ENatom2]  
                        kk2=readdefelem(ENatom1,elem_csv)[FN][ENatom1]
                    else:     
                        kk2=readdefelem(ENatom2,elem_csv)[FN][ENatom2]  
                        kk1=readdefelem(ENatom1,elem_csv)[FN][ENatom1]
                    
                    PGBox[int(kk1)][int(kk2)]=inty1+inty2+PGBox[int(kk1)][int(kk2)]
                    j = j + 1
                i= i + 1
                          
            #g=open("out.PGmatrix",'w')
            k_list_columns=[]
            k_list_index=[]       
            PGmtrx={}
            k=-1
            i=1
            while i<=18:
                j=i
                while j<=18:     
                    k = k + 1      
                    PGmtrx[int(k)]=PGBox[int(i)][int(j)]  
                    k_list_index.append("{:.6f}".format(PGmtrx[k]))
                    k_list_columns.append('PNmatrix_{}'.format(k))
                    j = j + 1
                i = i + 1
            k_dict=dict(zip(k_list_columns,k_list_index))
            k_df=pd.DataFrame(k_dict,index=['test'])    
       
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
                    if readdefelem(ENatom1,elem_csv)[FN][ENatom1] > readdefelem(ENatom2,elem_csv)[FN][ENatom2]:     
                        kk1=readdefelem(ENatom2,elem_csv)[FN][ENatom2]  
                        kk2=readdefelem(ENatom1,elem_csv)[FN][ENatom1]      
                    else:     
                        kk2=readdefelem(ENatom2,elem_csv)[FN][ENatom2]  
                        kk1=readdefelem(ENatom1,elem_csv)[FN][ENatom1]      
                    j = j + 1
                    PGBox[int(kk1)][int(kk2)]=inty1+inty2+PGBox[int(kk1)][int(kk2)]      
                i = i + 1
            h_list_index=[]
            h_list_columns=[]     
            k=-1
            i=1
            PGmtrx={m:0 for m in range(8)}
            while i<=7 :
                j=i
                while j<=7:      
                    k = k + 1      
                    PGmtrx[k]=PGBox[i][j]      
                    h_list_index.append("{:1.6f}".format(PGmtrx[k]))
                    h_list_columns.append("PGmatrix_{}".format(k))
                    j = j + 1
                i = i + 1     
            h_dict=dict(zip(h_list_columns,h_list_index))
            h_df=pd.DataFrame(h_dict,index=['test'])    
                     
     
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
     
                    if readdefelem(ENatom1,elem_csv)[FN][ENatom1] > readdefelem(ENatom2,elem_csv)[FN][ENatom2]:     
                        kk1=readdefelem(ENatom2,elem_csv)[FN][ENatom2]
                        kk2=readdefelem(ENatom1,elem_csv)[FN][ENatom1]      
                    else:     
                        kk2=readdefelem(ENatom2,elem_csv)[FN][ENatom2]  
                        kk1=readdefelem(ENatom1,elem_csv)[FN][ENatom1]
                        
                    if  kk1 == 'nodata' or kk2 == 'nodata':
                        j=j+1
                        return 'owari'
                        #PGBox[11][11]=inty1+inty2+PGBox[11][11]
                        
                    else:
                        PGBox[int(kk1)][int(kk2)]=inty1+inty2+PGBox[int(kk1)][int(kk2)]
                    j = j+1
                i = i + 1
                 
            #l = open("out.SPDFmatrix",'w')      
     
            k=-1
            i=1
            l_list_index=[]
            l_list_columns=[]
            while i<=4:
                j=i
                while j<=4:     
                    k = k + 1      
                    PGmtrx[k]=PGBox[i][j]      
                    l_list_index.append("{:1.6f}".format(PGmtrx[k]))
                    l_list_columns.append('SPDFmatrix_{}'.format(k))      
                    j = j + 1
                i = i + 1
            l_dict=dict(zip(l_list_columns,l_list_index))
            l_df=pd.DataFrame(l_dict,index=['test']) 
             
        #'list.{}-distfunc'.format(FN), )
        #f_list.write("# Atom-Atom characteristics Concentration\n")      
        i=0
        f_list_index=[]
        f_list_columns=[]
        while  i <= len(x[0])-1:
            ENatom=x[2][i]
            inty=x[0][i]/x[3]
            if readdefelem(ENatom,elem_csv)[FN][ENatom] != 'nodata':      
                f_list_index.append("{} {:1.6f} {:1.6f}".format(ENatom,float(readdefelem(ENatom,elem_csv)[FN][ENatom]),inty))
                f_list_columns.append('out.{}_{}'.format(FN,i))
            else:
                i=i+1
                return 'owari'
                #f_list_index.append("{} {} {:1.6f}".format(ENatom,readdefelem(ENatom,elem_csv)[FN][ENatom],inty))
                #f_list_columns.append('out.{}_{}'.format(FN,i))
            i = i + 1
        f_dict=dict(zip(f_list_columns,f_list_index))
        f_df = pd.DataFrame(f_dict,index=['test'])
        ans_df_dict[FN]=f_df
        if aparam3 == 'T':
            #m=open("list.{}{}-distfunc".format(FN,FN),'w')      
            #m.write("# Atom*Atom characteristics^2 Concentration^2\n")
            m_list_index=[]
            m_list_columns=[]      
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
                    if readdefelem(ENatom1,elem_csv)[FN][ENatom1] != 'nodata' and readdefelem(ENatom2,elem_csv)[FN][ENatom2] != 'nodata':
                        m_list_index.append("{}-{} {:1.6f} {:1.6f}".format(ENatom1,ENatom2,math.log(float(readdefelem(ENatom1,elem_csv)[FN][ENatom1]))+math.log(float(readdefelem(ENatom2,elem_csv)[FN][ENatom2])),inty1+inty2))
                        m_list_columns.append('out.{}{}_{}'.format(FN,FN,i))
                        #m_list_index.append("{}-{} {} {:1.6f}".format(ENatom1,ENatom2,readdefelem(ENatom1,df_csv)[FN][ENatom1],inty1+inty2))
                        #m_list_columns.append('out.{}{}_{}'.format(FN,FN,i))
                        #m_list_index.append("{}-{} {} {:1.6f}".format(ENatom1,ENatom2,readdefelem(ENatom2,df_csv)[FN][ENatom2],inty1+inty2))
                        #m_list_columns.append('out.{}{}_{}'.format(FN,FN,i))
                    else:
                        return 'owari'
                        #m_list_index.append("{}-{} {} {:1.6f}".format(ENatom1,ENatom2,readdefelem(ENatom2,df_csv)[FN][ENatom2],inty1+inty2))
                        #m_list_columns.append('out.{}{}_{}'.format(FN,FN,i))
                    j = j + 1
                i = i + 1
            #m.close()
            m_dict=dict(zip(m_list_columns,m_list_index))
            m_df=pd.DataFrame(m_dict,index=['test']) 
            ans_df_kake_dict[FN]=m_df   
            
            #o = open("list.{}-{}-distfunc".format(FN,FN),'w')      
            #o.write("# ABS(Atom-Atom) characteristics Concentration\n")      
            i=0
            o_list_index=[]
            o_list_columns=[]
            while i<=len(x[0])-1:      
                ENatom1=x[2][i]      
                inty1=x[0][i]/x[3]
                j=i+1
                while j <= len(x[0])-1:     
                    if i >= j:
                        continue     
                    ENatom2=x[2][j]      
                    inty2=x[0][j]/x[3]  
                    if readdefelem(ENatom1,elem_csv)[FN][ENatom1] != 'nodata' and readdefelem(ENatom2,elem_csv)[FN][ENatom2] != 'nodata':
                        o_list_index.append("{}-{} {} {:1.6f}".format(ENatom1,ENatom2,abs(float(readdefelem(ENatom1,elem_csv)[FN][ENatom1])-float(readdefelem(ENatom2,elem_csv)[FN][ENatom2])),inty1+inty2))      
                        o_list_columns.append('out.{}-{}_{}'.format(FN,FN,i))
                    elif readdefelem(ENatom1,elem_csv)[FN][ENatom1] == 'nodata' and readdefelem(ENatom2,elem_csv)[FN][ENatom2] != 'nodata':
                        j=j+1
                        return 'owari'
                        #o_list_index.append("{}-{} {} {:1.6f}".format(ENatom1,ENatom2,readdefelem(ENatom1,elem_csv)[FN][ENatom1],inty1+inty2))      
                        #o_list_columns.append('out.{}-{}_{}'.format(FN,FN,i))
                    elif readdefelem(ENatom1,elem_csv)[FN][ENatom1] != 'nodata' and readdefelem(ENatom2,elem_csv)[FN][ENatom2] == 'nodata':
                        j=j+1
                        return 'owari'
                        #o_list_index.append("{}-{} {} {:1.6f}".format(ENatom1,ENatom2,readdefelem(ENatom2,elem_csv)[FN][ENatom2],inty1+inty2))      
                        #o_list_columns.append('out.{}-{}_{}'.format(FN,FN,i))
                    else:
                        j=j+1
                        return 'owari'
                        #o_list_index.append("{}-{} {} {:1.6f}".format(ENatom1,ENatom2,readdefelem(ENatom2,elem_csv)[FN][ENatom2],inty1+inty2))      
                        #o_list_columns.append('out.{}-{}_{}'.format(FN,FN,i))
                    j = j + 1
                i= i + 1
            o_dict=dict(zip(o_list_columns,o_list_index))
            o_df=pd.DataFrame(o_dict,index=['test'])
            ans_df_hiki_dict[FN]=o_df
            #o.close()   
        
        list_reference = reference[0].split()
        #list_reference[0]="list.{}-distfunc".format(FN) 
        list_reference.append(float(df_csv.at[0,FN]))      
        list_reference.append(float(df_csv.at[1,FN]))      
        list_reference.append(float(df_csv.at[2,FN]))      
        list_reference.append(float(df_csv.at[3,FN]))      
        list_reference.append(float(df_csv.at[4,FN]))
        #print(list_reference)
        if df_csv.at[5,FN] =='T':
            list_reference.append(df_csv.at[5,FN])
        else:
            list_reference.append(float(df_csv.at[5,FN]))    
        #       if ($aparam[3] == 0){$reference[6]=0 } 
        #if aparam3 == 0:
        #    list_reference[6]=0
        #print(df_csv[FN])
        #print(list_reference)
        z = dfmake(*list_reference,aparam2,aparam3,ans_df_dict[FN],df_csv)
        if z[6].empty :
            return 'owari'
        old_col=list(z[6].columns)
        new_col=[]
        for col in range(len(old_col)):
            new_col.append('{}_{}'.format(FN,col))
        z[6].set_axis(new_col, axis='columns',inplace=True)
        df_dict[FN]=z[6]
        #os.system('mv {} out.{}'.format(z[7],FN))     
        #outcompfile.append("out.{}".format(FN))
        #outcompfile.append(z[7])

        if aparam3 == 'T':        
            #list_reference[0]="list.{}{}-distfunc".format(FN,FN)     
            list_reference[1]=float(df_csv.at[0,FN])     
            list_reference[2]=float(df_csv.at[1,FN])     
            list_reference[3]=-0.2*math.log(float(df_csv.at[3,FN]))     
            list_reference[4]=2*math.log(float(df_csv.at[3,FN]))     
            list_reference[5]=float(df_csv.at[4,FN])
            list_reference[6]=2*math.log(float(df_csv.at[3,FN]))/float(df_csv.at[4,FN])      
            #if aparam3 == 0:        
            #    list_reference[6]=0
            z = dfmake(*list_reference,aparam2,aparam3,ans_df_kake_dict[FN],df_csv)
            if z[6].empty :
                return 'owari'
            old_col=list(z[6].columns)
            new_col=[]
            for col in range(len(old_col)):
                new_col.append('{}{}_{}'.format(FN,FN,col))
            z[6].set_axis(new_col, axis='columns',inplace=True)
            df_kake_dict[FN] =z[6]


            #os.system('mv {} out.{}{}'.format(z[7],FN,FN))      
            #outcompfile.append("out.{}{}".format(FN,FN))   
            #outcompfile.append(z[7]) 

            #list_reference[0]="list.{}-{}-distfunc".format(FN,FN)      
            list_reference[1]=float(df_csv.at[0,FN])       
            list_reference[2]=float(df_csv.at[1,FN])      
            list_reference[3]=float(df_csv.at[2,FN])      
            list_reference[4]=float(df_csv.at[3,FN])      
            list_reference[5]=float(df_csv.at[4,FN])
            if df_csv.at[5,FN] =='T':
                list_reference[6]=df_csv.at[5,FN]
            else:    
                list_reference[6]=float(df_csv.at[5,FN]) 

            #if aparam3 == 0:
            #list_reference[6]=0
            z = dfmake(*list_reference,aparam2,aparam3,ans_df_hiki_dict[FN],df_csv)
            if z[6].empty :
                return 'owari'
            old_col=list(z[6].columns)
            new_col=[]
            for col in range(len(old_col)):
                new_col.append('{}-{}_{}'.format(FN,FN,col))
            z[6].set_axis(new_col, axis='columns',inplace=True)     
            df_hiki_dict[FN]=z[6] 
            #os.system('mv {} out.{}-{}'.format(z[7],FN,FN))      
            #outcompfile.append("out.{}-{}".format(FN,FN)) 
            #outcompfile.append(z[7]) 

    return k_df,h_df,l_df,df_dict,df_kake_dict,df_hiki_dict,filename


     
#compdescript(defelem)

"""  
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
"""     
