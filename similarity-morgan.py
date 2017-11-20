#!/usr/bin/env python
# -*- coding: utf-8 -*-

##    Description    Tool for generating a SDFile from diverse sources
##                   
##    Authors:       Inés Martínez (mmartinez4@imim.es)
##                   Manuel Pastor (manuel.pastor@upf.edu)
##
##    Copyright 2015 Manuel Pastor
##
##    This file is part of PhiTools
##
##    PhiTools is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation version 3.
##
##    PhiTools is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with PhiTools.  If not, see <http://www.gnu.org/licenses/>

import urllib
import os
import sys
import getopt
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from rdkit import DataStructs
from standardise import standardise

def getSimilarity(data1,data2,numMFPS):

    suppl=Chem.SDMolSupplier(data1)
    mF_model = []
    mF_DB = [] ## morganFP from Drug Bank
    name_model = []
    id_model = []
    for m in suppl:
        try:
            mf = AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=1024)
            mF_model.append(mf)
            name = m.GetProp('_Name')
            if len(name) < 1 :
                name = m.GetProp('name')
            name_model.append(name)
            code = m.GetProp('molID')
            id_model.append(code)
            
        except:
            return ('error_fingerprinting')
    
    suppl2=Chem.SDMolSupplier(data2)
    counter = 0
    simStruct = []

    print 'Number of Compounds\tid_Tox21\tname_Tox21\tmolID-CS1\tname_CS1\tsimilarity'
    for m in suppl2:
        print str(m.GetProp('Substance_Name')
        max_similarity = 0
        max_similarity_id = ''
        max_similarity_name = ''
        try:
            mf = AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=1024)
            for j in range(len(mF_model)):
                sim = DataStructs.FingerprintSimilarity(mF_model[j],mf)
                if sim>max_similarity:
                   max_similarity = sim
                   max_similarity_id = id_model[j]
                   max_similarity_name = name_model[j]
        except:
            max_similarity = 0
            max_similarity_id = ''  

##        print max_similarity, m.GetProp('Substance_Name'),
##        print max_similarity, m.GetProp('Substance_Name')

        if max_similarity > numMFPS:
            simStruct.append(m)
            counter =counter+1
            print str(counter)+'\t'+str(m.GetProp('ToxCast_chid'))+'\t',
            print str(m.GetProp('Substance_Name'))+'\t',
            print str(max_similarity_id)+'\t'+ str(max_similarity_name)+'\t',
            print str(max_similarity)
        else:
            #print
            pass

    #print 'total molecules found ', counter
    
    return (simStruct)


def writeSDF(simStruct,sdfname):

    nameSDF = Chem.SDWriter(sdfname)
    
    for m in simStruct:
        nameSDF.write(m)    

def usage ():
    """Prints in the screen the command syntax and argument"""
    print 'getSimilarity -m model.sdf -d database.sdf [-o output.sdf] -s'
    print '\n\t -m model.sdf model structures to compare with'
    print '\t -d database.sdf (molecules from DrugBank Database)'
    print '\t -o output.sdf (output SDFile)'
    print '\n\t -s fingerprint score from 0 to 1, used to filter molecules'
    sys.exit(1)

def main ():
    sdfname = 'output.sdf'
    data1 = None
    data2 = None
    
    try:
       opts, args = getopt.getopt(sys.argv[1:],'o:m:d:s:')
    except getopt.GetoptError:
        print "Arguments not recognized"
        usage()
    
    if not len(opts):
        print "Arguments not recognized"
        usage()

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-o':
                sdfname = arg
            elif opt in '-m':
                data1 = arg
            elif opt in '-d':
                data2 = arg
            elif opt in '-s':
                numMFPS = arg ### range number to filter Morganfingerprints scores
                numMFPS = float(numMFPS)
                
    if data1 ==None or data2 == None: usage()

    simStruct = getSimilarity(data1, data2,numMFPS)
    writeSDF(simStruct, sdfname)

if __name__ == '__main__':    
    main()
