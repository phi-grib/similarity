##    Description    Similarity Tool for comparing two datasets 
##                   and generating a SDFile with similar compounds 
##                   
##    Authors:       Kevin Pinto Gil (kevin.pinto@upf.edu)
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


# # Libraries to import to run similarity analysis

import urllib
import os
import sys
import getopt
import numpy as np
from cmath import *

import warnings


from rdkit import Chem
from rdkit.Chem import AllChem
with warnings.catch_warnings():
     warnings.simplefilter("ignore")
from rdkit import DataStructs


def getSimilarity(md,db,cutoff,outSim,dbIDname,mdIDname):
    ''' Protocol to be followed before similarity analysis:
        - Standardize field names in model molecules to be distinguised with Database field names ( e.g. name for name-model) 
        - Standardize molecules to eliminate contra ions, duplicates, molecules containing metal ions.... 
        - Double check manually if the molecules excluded are well excluded or they need to be included. 
        - Obtain 2D or 3D coordinates 
        - Protonate structures at the same pH (e.g. 7.4 pH) to be comparable
        - Calculate final smile, inchi and inchikey

        This functions reads two SD files where:
        - one contains a Reference Database
        - Model database
        Then a similiraty analysis using Morgan Fingerprints is performed at the cutoff provided

        The output will be the molecules from the reference database similar at the Model database compounds
        at the cutoff provided.'''
    
    ### Loading Model SDfile
    
    suppl=Chem.SDMolSupplier(md)
    
    mF_model = [] ## will contain model morgan fingerprints 
    id_model = [] ## molID of the model created by me ( e.g. mol00001)

    mdictnames = {}
    mdictvalues = {}
    
    for m in suppl:
        
        try:
            mf = AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=1024)
            mF_model.append(mf)
            id_model.append(m.GetProp(mdIDname))      
            
            mnames = [] # model names
            mnewnames = [] # model names  plus MODEL_ to be diferentiated with database ones
            mvalues = [] # model values
            for i in m.GetPropNames():
                mnames.append(i)
                mnewnames.append('MODEL_'+i) ### modify 'MODEL_' for the name that you want to 
            for i in mnames:   
                mvalues.append(m.GetProp(i))
            mdictnames[m.GetProp(mdIDname)] = mnewnames    
            mdictvalues[m.GetProp(mdIDname)] = mvalues  
        except:
            return ('error_fingerprinting')
   
    ### Loading Reference Database SD file and running the similiraty analysis. 
        
    suppl2=Chem.SDMolSupplier(db)
    counter = 0
    simStruct = []
    
    fo= open (outSim,'w')

    print 'numOFcompounds\trefID\tmodelID\tsimilarity'    
    for m in suppl2:       
#         if m is None: continue
        max_similarity = 0
        max_similarity_id = ''

        try:
            mf = AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=1024)
            for j in xrange(len(mF_model)):
                sim = DataStructs.FingerprintSimilarity(mF_model[j],mf)
                if sim>max_similarity:
                    max_similarity = sim
                    max_similarity_id = id_model[j]
        except:
            max_similarity = 0
            max_similarity_id= ''

        if max_similarity > cutoff and max_similarity < 1:
            simStruct.append(m)
            parent = Chem.MolToMolBlock(m)
            fo.write(parent)
            pnames = []
            for i in m.GetPropNames():
                pnames.append(i)

            pvalues = []

            for i in pnames:
                pvalues.append(m.GetProp(i))
                        
            for pn,pv in zip(pnames,pvalues):        
                fo.write('>  <'+str(pn)+'>\n'+str(pv)+'\n\n')

            for sn, sv in zip(mdictnames[max_similarity_id], mdictvalues[max_similarity_id]):
                fo.write('>  <'+str(sn)+'>\n'+str(sv)+'\n\n')
            fo.write('>  <similarity>\n'+str(max_similarity)+'\n\n')
            fo.write('$$$$\n')
            
            counter =counter+1
            print str(counter)+'\t'+str(m.GetProp(dbIDname))+'\t',
            print str(max_similarity_id)+'\t',
            print str(max_similarity)

        else:
            pass

    print 'total molecules found = ', counter
    fo.close()
    return (simStruct)



def usage ():
    
    """Prints in the screen the command syntax and argument"""
    print 'getSimilarity -a model.sdf -b database.sdf -c cutoff -d refIDname -e modelIDname -o output.sdf'
    print '\n\t -a ReferenceDatabase.sdf (molecules from Reference database e.g. DrugBank Database)'
    print '\n\t -b Modeldatabase.sdf (model structures to compare with)'
    print '\n\t -c cutoff similarity distance by tanimoto Morgan fingerprints (score from 0 to 1, used to filter molecules)'
    print '\n\t -d drugbankID (RefDatabase field compound name) One must provide an Identifier from the Reference database'
    print '\n\t -e molID (Model field compound name) One must provide an Identifier from the model database' 
    print '\n\t -o output.sdf (output SDFile)'

    ## db =  '/media/sf_users/eutoxrisk/biosimilarity_task/5-CASE1-EU-ToxRisk/7-Tox21/6-moka3D/tox21-3D-moka.sdf' ## database
    ## md = '/media/sf_users/eutoxrisk/biosimilarity_task/5-CASE1-EU-ToxRisk/4-moka3D/steatosis-3D-moka-plus-activity.sdf' ## model
    ## cutoff = 0.6 ### cutoff similarity distance by tanimoto Morgan fingerprints
    ## outSim = 'similars-Kevin.sdf' ### output name similars
    ## idName = 'molID'
    ## simStruct = getSimilarity(md,db,cutoff,outSim, idName)


    sys.exit(1)

def main ():
    
    db =  None
    md = None
    cutoff = None
    dbIDname = None
    mdIDname = None 
    outSim = None
   

    try:
       opts, args = getopt.getopt(sys.argv[1:],'a:b:c:d:e:o:')
    except getopt.GetoptError:
        print "Arguments not recognized"
        usage()
    
    if not len(opts):
        print "Arguments not recognized"
        usage()

    if len(opts)>0:
        for opt, arg in opts:
            if opt in '-a':
                db = arg
            elif opt in '-b':
                md = arg
            elif opt in '-c':
                cutoff = arg ### range number to filter Morganfingerprints scores
                cutoff = float(cutoff)
            elif opt in '-d':
                dbIDname = arg
            elif opt in '-e':
                mdIDname = arg
            elif opt in '-o':
                outSim = arg
                
    if db ==None or md == None or cutoff == None or dbIDname ==None or mdIDname == None or outSim == None: usage()

    simStruct = getSimilarity(md,db,cutoff,outSim,dbIDname,mdIDname)


if __name__ == '__main__':    
    main()


