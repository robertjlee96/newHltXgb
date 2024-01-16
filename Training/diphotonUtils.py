#!/usr/bin/env python
import re
import numpy as np
import uproot
from ROOT import TFile, TTree, TChain
from array import array
import xml.etree.cElementTree as ET
regex_float_pattern = r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'


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

extraVars = ['trkIsoPho','et','energy']
singleVars = ['mass','nEgs','passFailStd','passFailL1Single','passFailL1Double','triggerBits']

#----------------------------------------------------------------------

#def load_file(input_file, geoSelection = None, ptCuts = None, dptCuts = None):
def loadFile(inputFile, geoSelection = None, ptCuts = None, massCuts = None):


    """input_file should be a uproot object corresponding to a ROOT file

    :return: input_values, target_values, orig_weights, train_weights, pt, scEta, input_var_names
    """
    inputValuesSig = []
    inputValuesSigLead = []
    inputValuesSigSub = []
    
    targetValuesSig = []

    inputValuesBkg = []
    inputValuesBkgLead = []
    inputValuesBkgSub = []
    
    targetValuesBkg = []


    # names of variables used as BDT input
    inputVarNames = []

    # original weights without pt/eta reweighting
    # we can use these also for evaluation
    origWeightsSig = []
    origWeightsBkg = []

    varValuesSig = []
    varValuesSigLead = []
    varValuesSigSub = []
    
    varValuesBkg = []
    varValuesBkgLead = []
    varValuesBkgSub = []
    
    varNames = []

    isFirstVar = True

    for varname in inputVars + extraVars + singleVars + ['index']:

        thisValuesSig = []
        thisValuesBkg = []

        isFirstProc = True

        for treeName, label in [
            ('sigTree', 1),
            ('bkgTree', 0)
        ]:

            tree = inputFile[treeName]
           
            #mask = []
            
            mask = massCuts(tree)
            #print mask
            
#            for i in range(0,len(tree.array('eta'))):
#                if tree.array('mass')[i][0] > 60:
#                    mask.append(True)
#                else:
#                    mask.append(False)

            if not mask is None:
                indices = mask
            else:
                indices = np.ones(len(tree.arrays()[varname]), dtype = 'bool')
                
           
            if varname != 'index' and varname not in singleVars:
                if label == 0:
                    #varValuesBkg.append([(item[0],item[1]) for item in tree.arrays()[varname][indices]])
                    varValuesBkgLead.append([item[0] for item in tree.arrays()[varname][indices]])
                    varValuesBkgSub.append([item[1] for item in tree.arrays()[varname][indices]])
                if label == 1:
                    #varValuesSig.append([(item[0],item[1]) for item in tree.arrays()[varname][indices]])
                    varValuesSigLead.append([item[0] for item in tree.arrays()[varname][indices]])
                    varValuesSigSub.append([item[1] for item in tree.arrays()[varname][indices]])


                if isFirstProc:
                    varNames.append(varname)
#
                if varname in inputVars:
                    # BDT input variable
                    if label == 0:
                        #inputValuesBkg.append([(item[0],item[1]) for item in tree.arrays()[varname][indices]])
                        inputValuesBkgLead.append([item[0] for item in tree.arrays()[varname][indices]])
                        inputValuesBkgSub.append([item[1] for item in tree.arrays()[varname][indices]])
                    if label == 1:
                        #inputValuesSig.append([(item[0],item[1]) for item in tree.arrays()[varname][indices]])
                        inputValuesSigLead.append([item[0] for item in tree.arrays()[varname][indices]])
                        inputValuesSigSub.append([item[1] for item in tree.arrays()[varname][indices]])

                    if isFirstProc:
                        inputVarNames.append(varname)
                        
            elif varname in singleVars:
                if label == 0:
                    #varValuesBkg.append([(item[0],item[1]) for item in tree.arrays()[varname][indices]])
                    varValuesBkgLead.append([item for item in tree.arrays()[varname][indices]])
                    varValuesBkgSub.append([item for item in tree.arrays()[varname][indices]])
                if label == 1:
                    #varValuesSig.append([(item[0],item[1]) for item in tree.arrays()[varname][indices]])
                    varValuesSigLead.append([item for item in tree.arrays()[varname][indices]])
                    varValuesSigSub.append([item for item in tree.arrays()[varname][indices]])
                        
            elif varname == 'index':
                if label == 0:
                    indicesBkg = range(0,len(tree.arrays()['eta'][indices]))
                    #indicesArray = [val for pair in zip(indicesBkg, indicesBkg) for val in pair]
                    varValuesBkgLead.append(indicesBkg)
                    varValuesBkgSub.append(indicesBkg)
                if label == 1:
                    indicesSig = range(0,len(tree.arrays()['eta'][indices]))
                    indicesArray = [val for pair in zip(indicesSig, indicesSig) for val in pair]
                    varValuesSigLead.append(indicesSig)
                    varValuesSigSub.append(indicesSig)
                

            # append target values and weights
            if isFirstVar:
#                #thisWeights =  tree.array('weight')[indices]
                thisWeights =  np.ones(len(mask))[indices]
                if label == 0:
                    targetValuesBkg.append(np.ones(len(inputValuesBkgLead[-1])) * 1)
                    origWeightsBkg.append(thisWeights)
                    
                if label == 1:
                    targetValuesSig.append(np.ones(len(inputValuesSigLead[-1])) * 1)
                    origWeightsSig.append(thisWeights)
    
            isFirstProc = False

        if isFirstVar:
            if label == 0:
                targetValuesBkg = np.hstack(targetValuesBkg)
                origWeightsBkg = np.hstack(origWeightsBkg)
           
            if label == 1:
                targetValuesSig = np.hstack(targetValuesSig)
                origWeightsSig = np.hstack(origWeightsSig)

        isFirstVar = False




#    inputValuesSig = np.vstack(inputValuesSig)
#    varValuesSig = np.vstack(varValuesSig)

#    inputValuesBkg = np.vstack(inputValuesBkg)
#    varValuesBkg = np.vstack(varValuesBkg)
#
#
    inputValuesSig = np.concatenate((inputValuesSigLead,inputValuesSigSub), axis=0)
    varValuesSig = np.concatenate((varValuesSigLead,varValuesSigSub), axis=0)

    inputValuesBkg = np.concatenate((inputValuesBkgLead,inputValuesBkgSub), axis=0)
    varValuesBkg = np.concatenate((varValuesBkgLead,varValuesBkgSub), axis=0)

    varNames = np.hstack(varNames)

    return inputValuesSig, inputValuesBkg, targetValuesSig, targetValuesBkg, origWeightsSig, origWeightsBkg, varValuesSig, varValuesBkg, inputVarNames, varNames

def buildBranch(treeName,vars,scores,outFileRoot):
    xgbScore = array('f',[0.0,0.0])
    rawEnergy = array('f',[0.0,0.0])
    r9HLT = array('f',[0.0,0.0])
    sigmaIEtaIEta = array('f',[0.0,0.0])
    etaWidth = array('f',[0.0,0.0])
    phiWidth = array('f',[0.0,0.0])
    s4 = array('f',[0.0,0.0])
    trkIsoPho = array('f',[0.0,0.0])
    eta = array('f',[0.0,0.0])
    hOvrE = array('f',[0.0,0.0])
    ecalPFIso = array('f',[0.0,0.0])
    et = array('f',[0.0,0.0])
    energy = array('f',[0.0,0.0])
    mass = array('f',[0.0,0.0])
    nEgs = array('f',[0.0,0.0])
    triggerBits = array('f',[0.0,0.0])
    passFailStd = array('i',[0])
    passFailL1Single = array('i',[0])
    passFailL1Double = array('i',[0])

    tree = TTree(treeName,treeName)
    tree.Branch( 'xgbScore', xgbScore, 'xgbScore[2]/F' )
    tree.Branch( 'rawEnergy', rawEnergy, 'rawEnergy[2]/F' )
    tree.Branch( 'r9HLT', r9HLT, 'r9HLT[2]/F' )
    tree.Branch( 'sigmaIEtaIEta', sigmaIEtaIEta, 'sigmaIEtaIEta[2]/F' )
    tree.Branch( 'etaWidth', etaWidth, 'etaWidth[2]/F' )
    tree.Branch( 'phiWidth', phiWidth, 'phiWidth[2]/F' )
    tree.Branch( 's4', s4, 's4[2]/F' )
    tree.Branch( 'trkIsoPho', trkIsoPho, 'trkIsoPho[2]/F' )
    tree.Branch( 'eta', eta, 'eta[2]/F' )
    tree.Branch( 'hOvrE', hOvrE, 'hOvrE[2]/F' )
    tree.Branch( 'ecalPFIso', ecalPFIso, 'ecalPFIso[2]/F' )
    tree.Branch( 'et', et, 'et[2]/F' )
    tree.Branch( 'energy', energy, 'energy[2]/F' )
    tree.Branch( 'mass', mass, 'mass[2]/F' )
    tree.Branch( 'nEgs', nEgs, 'nEgs[2]/F' )
    tree.Branch( 'triggerBits', triggerBits, 'triggerBits[2]/F' )
    tree.Branch( 'passFailStd', passFailStd, 'passFailStd/I' )
    tree.Branch( 'passFailL1Single', passFailL1Single, 'passFailL1Single/I' )
    tree.Branch( 'passFailL1Double', passFailL1Double, 'passFailL1Double/I' )

    nEvents = int(len(vars)/2)

    for n in range(0,nEvents):
        #Order indices by Et, selecting Lead vs Sub for tree here
        nLead = -100000
        nSub = -100000
        nFirst = 2*n
        nSecond = 2*n + 1

        if vars[nFirst][10] > vars[nSecond][10]:
            nLead = nFirst
            nSub = nSecond
        else:
            nLead = nSecond
            nSub = nFirst

        xgbScore[0] = scores[nLead][0]
        xgbScore[1] = scores[nSub][0]
        rawEnergy[0] = vars[nLead][0]
        rawEnergy[1] = vars[nSub][0]
        r9HLT[0] = vars[nLead][1]
        r9HLT[1] = vars[nSub][1]
        sigmaIEtaIEta[0] = vars[nLead][2]
        sigmaIEtaIEta[1] = vars[nSub][2]
        etaWidth[0] = vars[nLead][3]
        etaWidth[1] = vars[nSub][3]
        phiWidth[0] = vars[nLead][4]
        phiWidth[1] = vars[nSub][4]
        s4[0] = vars[nLead][5]
        s4[1] = vars[nSub][5]
        trkIsoPho[0] = vars[nLead][6]
        trkIsoPho[1] = vars[nSub][6]
        eta[0] = vars[nLead][7]
        eta[1] = vars[nSub][7]
        hOvrE[0] = vars[nLead][8]
        hOvrE[1] = vars[nSub][8]
        ecalPFIso[0] = vars[nLead][9]
        ecalPFIso[1] = vars[nSub][9]
        et[0] = vars[nLead][10]
        et[1] = vars[nSub][10]
        energy[0] = vars[nLead][11]
        energy[1] = vars[nSub][11]
        mass[0] = vars[nLead][12]
        mass[1] = vars[nSub][12]
        nEgs[0] = vars[nLead][13]
        nEgs[1] = vars[nSub][13]
        passFailStd[0] = int(vars[nLead][14])
        passFailL1Single[0] = int(vars[nLead][15])
        passFailL1Double[0] = int(vars[nLead][16])
        triggerBits[0] = int(vars[nLead][17])
        triggerBits[1] = int(vars[nLead][17])

        tree.Fill()
    outFileRoot.Write()
    return
    
def build_tree(xgtree, base_xml_element, var_indices,var_list):
    parent_element_dict = {'0':base_xml_element}
    pos_dict = {'0':'s'}
    
    varListGenNames = ['f0','f1','f2','f3','f4','f5','f6','f7','f8']
    
#    varListBothNames = {}
#
#    for i in range(0,len(varListGenNames)):
#        var = varListGenNames[i]
#        varListBothNames[var] = var_list
        
    for line in xgtree.split('\n'):
        if not line: continue
        if ':leaf=' in line:
            #leaf node
            result = re.match(r'(\t*)(\d+):leaf=({0})$'.format(regex_float_pattern), line)
            if not result:
                print(line)
            depth = result.group(1).count('\t')
            inode = result.group(2)
            res = result.group(3)
            node_elementTree = ET.SubElement(parent_element_dict[inode], "Node", pos=str(pos_dict[inode]),
                                             depth=str(depth), NCoef="0", IVar="-1", Cut="0.0e+00", cType="1", res=str(res), rms="0.0e+00", purity="0.0e+00", nType="-99")
        else:
            #\t\t3:[var_topcand_mass<138.19] yes=7,no=8,missing=7
            result = re.match(r'(\t*)([0-9]+):\[(?P<var>.+)<(?P<cut>{0})\]\syes=(?P<yes>\d+),no=(?P<no>\d+)'.format(regex_float_pattern),line)
            if not result:
                print(line)
            depth = result.group(1).count('\t')
            inode = result.group(2)
            var = result.group('var')
            cut = result.group('cut')
            lnode = result.group('yes')
            rnode = result.group('no')
            pos_dict[lnode] = 'l'
            pos_dict[rnode] = 'r'
#            indexTmp = var_indicies.index(var)
#            newName = var_indices[
#            print "var_indices = {}".format(var_indices)
#            print "var_list = {}".format(var_list)
#            print "var = {}".format(var)
            varIndex = varListGenNames.index(var)
#            print "varIndex = {}".format(varIndex)
            varCorr = var_list[varIndex][0]
#            print "varCorr = {}".format(varCorr)
            #print "indexTmp = {}".format(var_indices[var])
        
            node_elementTree = ET.SubElement(parent_element_dict[inode], "Node",
                pos=str(pos_dict[inode]), depth=str(depth),
                NCoef="0", IVar=str(var_indices[varCorr]),
                Cut=str(cut), cType="1", res="0.0e+00", rms="0.0e+00", purity="0.0e+00", nType="0")
            parent_element_dict[lnode] = node_elementTree
            parent_element_dict[rnode] = node_elementTree
        


def convert_model(model, input_variables, output_xml):
    NTrees = len(model)
    var_list = input_variables
#    print "var_list Within Convert Model = {}".format(var_list)
    var_indices = {}
    
    # <MethodSetup>
    MethodSetup = ET.Element("MethodSetup", Method="BDT::BDT")

    # <Variables>
    Variables = ET.SubElement(MethodSetup, "Variables", NVar=str(len(var_list)))
    for ind, val in enumerate(var_list):
#        print "ind, val = {}, {}".format(ind,val)
        name = val[0]
        var_type = val[1]
        var_indices[name] = ind
        Variable = ET.SubElement(Variables, "Variable", VarIndex=str(ind), Type=val[1],
            Expression=name, Label=name, Title=name, Unit="", Internal=name,
            Min="0.0e+00", Max="0.0e+00")
#    print "var_indices = {}".format(var_indices)
    GeneralInfo = ET.SubElement(MethodSetup, "GeneralInfo")
    Info_Creator = ET.SubElement(GeneralInfo, "Info", name="Creator", value="xgboost2TMVA")
    Info_AnalysisType = ET.SubElement(GeneralInfo, "Info", name="AnalysisType", value="Classification")

    # <Options>
    Options = ET.SubElement(MethodSetup, "Options")
    Option_NodePurityLimit = ET.SubElement(Options, "Option", name="NodePurityLimit", modified="No").text = "5.00e-01"
    Option_BoostType = ET.SubElement(Options, "Option", name="BoostType", modified="Yes").text = "Grad"
    
    # <Weights>
    Weights = ET.SubElement(MethodSetup, "Weights", NTrees=str(NTrees), AnalysisType="1")
    
    for itree in range(NTrees):
        BinaryTree = ET.SubElement(Weights, "BinaryTree", type="DecisionTree", boostWeight="1.0e+00", itree=str(itree))
        build_tree(model[itree], BinaryTree, var_indices,var_list)
        
    tree = ET.ElementTree(MethodSetup)
    tree.write(output_xml)
    # format it with 'xmlli
