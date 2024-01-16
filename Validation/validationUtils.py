#!/usr/bin/env python

import numpy as np
import uproot
from ROOT import TFile, TTree, TChain
import array
#INCLUDE 0 for eventIndex when counting all vars!
inputVars = [
    'rawEnergy',
    'r9HLT',
    'sigmaIEtaIEta',
    'etaWidth',
    'phiWidth',
    's4',
    'eta',
    'hOvrE',
    'ecalPFIso'
]

extraVars = ['trkIsoPho','et','energy','phi','pdgID']
singleVars = ['mass','nEgs','passFailStd','passFailL1Single','passFailL1Double','triggerBits','nEvent','nEgsPassing','runNr']

#----------------------------------------------------------------------

#def load_file(input_file, geoSelection = None, ptCuts = None, dptCuts = None):
def load_file(inputFile, geoSelection = None, ptCuts = None):


    """input_file should be a uproot object corresponding to a ROOT file

    :return: input_values, target_values, orig_weights, train_weights, pt, scEta, input_var_names
    """
    inputValuesSig = []
    targetValuesSig = []


    inputValuesBkg = []
    targetValuesBkg = []


    # names of variables used as BDT input
    inputVarNames = []

    # original weights without pt/eta reweighting
    # we can use these also for evaluation
    origWeightsSig = []
    origWeightsBkg = []

    varValuesSig = []
    varValuesBkg = []
    varNames = []

    isFirstVar = True
        
    fileIn = uproot.open(inputFile)
    
    nEventsSig = 0
    nEventsBkg = 0
    
    for treeName, label in [
        ('sigTree', 1),
        ('bkgTree', 0)
        ]:

        tree = fileIn[treeName]
        
        mask = None

        if not mask is None:
            indices = mask
        else:
            indices = np.ones(len(tree.arrays()['eta']), dtype = 'bool')
                
        
        countingArray = tree.arrays()['eta'][indices]
        
        if label == 0:
            nEventsBkg = len(countingArray)
        if label == 1:
            nEventsSig = len(countingArray)

        for varname in ['index'] + inputVars + extraVars + singleVars:
            isFirstProc = True
        
            varArray = []
            if varname != 'index':
                varArray = tree.arrays()[varname][indices]
                
            allPhotonsVarArray = []
         
            if varname not in singleVars and varname != 'index':
                for i in range(0,len(varArray)):
                    for j in range(0,len(varArray[i])):
                        allPhotonsVarArray.append(varArray[i][j])
                
                if label == 0:
                    varValuesBkg.append(allPhotonsVarArray)
                if label == 1:
                    varValuesSig.append(allPhotonsVarArray)

                if isFirstProc:
                    varNames.append(varname)

                if varname in inputVars:
                    # BDT input variable
                    if label == 0:
                        inputValuesBkg.append(allPhotonsVarArray)
                    if label == 1:
                        inputValuesSig.append(allPhotonsVarArray)

                    if isFirstProc:
                        inputVarNames.append(varname)
            
            elif varname in singleVars:
                for i  in range(0,len(countingArray)):
                    for j in range(0,len(countingArray[i])):
                        allPhotonsVarArray.append(varArray[i])
                if label == 0:
                    varValuesBkg.append(allPhotonsVarArray)
                if label == 1:
                    varValuesSig.append(allPhotonsVarArray)

            elif varname == 'index':
                indexTmp = 0
                indexArray = []
                for i  in range(0,len(countingArray)):
                    for j in range(0,len(countingArray[i])):
                        indexArray.append(indexTmp)
                        indexTmp += 1
                if label == 0:
                    varValuesBkg.append(indexArray)
                if label == 1:
                    varValuesSig.append(indexArray)
                    

            isFirstProc = False

        isFirstVar = False

    inputValuesSig = np.vstack(inputValuesSig)
    varValuesSig = np.vstack(varValuesSig)

    inputValuesBkg = np.vstack(inputValuesBkg)
    varValuesBkg = np.vstack(varValuesBkg)

    varNames = np.hstack(varNames)

    return inputValuesSig, inputValuesBkg, varValuesSig, varValuesBkg, inputVarNames, varNames, nEventsSig, nEventsBkg
