#!/usr/bin/env python -W ignore::DeprecationWarning
import numpy as np
import uproot
from sklearn.model_selection import train_test_split
#from sklearn.utils.fixes import _Iterable as Iterable, _Sized as Sized
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import accuracy_score
from ROOT import TFile, TTree, TGraph
from ROOT import gROOT, AddressOf
import xgboost as xgb
from xgboost import XGBClassifier
import matplotlib
matplotlib.use('agg') #Sketchy way of using alternate backend for problematic Tkinter package (FIX?)
import matplotlib.pyplot as plt
import diphotonUtils
from diphotonUtils import buildBranch
#import oldXGBoost2tmva
import time,pickle
from tqdm import tqdm
import sys
from tempfile import TemporaryFile
import time
from array import array
#import joblib
import gc
from operator import itemgetter

np.random.seed(1337)

inFileName = "../../CMSSW_12_4_8/src/allPhotonsNTuples/GGH13andEphemD_GenInfoCUT_DiphotonTrain_1115.root"
#inFileName = "../../CMSSW_12_4_8/src/allPhotonsNTuples/GGH13andDataD_AllCombined_DiphotonTrain_0808.root"
outNamesGen = []
outNamesGen.append("0111/M7L25_GGH13andDataD_NoTrkIso_M60_PdgIDCut_0111")

modelNamesB = []
modelNamesB.append("barrelOut/" + outNamesGen[0] + "_Barrel")
#modelNamesB.append("../Training/barrelOut/Model_MD13LR07_M90_Barrel_0418.Model")
#modelNamesB.append("../Training/barrelOut/Model_MD17LR05_M90_Barrel_0418.Model")
#
modelNamesE = []
modelNamesE.append("endcapOut/" + outNamesGen[0] + "_Endcap")
#modelNamesE.append("../Training/endcapOut/Model_MD13LR07_M90_Endcap_0418.model")
#modelNamesE.append("../Training/endcapOut/Model_MD17LR05_M90_Endcap_0418.model")
#
fNames = []
fNames.append("trainingNTuples/" + outNamesGen[0] +".root")
#fNames.append("validationNTuples/0425/GGHandData_MD13LR07")
#fNames.append("validationNTuples/0425/GGHandData_MD17LR05")

fin = uproot.open(inFileName)
prompt = fin['sigTree']
fake = fin['bkgTree']

## for barrel
geometry_selection = lambda tree: np.abs(tree.array('eta')) < 2.5
ptCut = lambda tree: tree.array('et') > 14.25
massCut = lambda tree: tree.arrays()['mass'] > 60.0

#diphotonUtils.load_file(fin,geometry_selection,ptCut,massCut)
#Initial script gives event-based arrays. One element will have 38 vars, 19 for each photon (all lead, then all sub)
inputValuesSig, inputValuesBkg, targetValuesSigX, targetValuesBkgX, origWeightsSig, origWeightsBkg, varValuesSig, varValuesBkg, inputVarNames, varNames = diphotonUtils.loadFile(fin,geometry_selection,ptCut,massCut)

weightsSig = np.ones(len(inputValuesSig[0]))
weightsBkg = np.ones(len(inputValuesBkg[0]))

targetValuesSig = np.ones(len(inputValuesSig[0]))
targetValuesBkg = np.zeros(len(inputValuesBkg[0]))

#Add together sig and bkg so we can do Test/Train split EVENT-BY-EVENT
inputTotal = np.concatenate((inputValuesSig,inputValuesBkg),axis = 1)
weightsTotal = np.concatenate((weightsSig,weightsBkg),axis = 0)
targetsTotal = np.concatenate((targetValuesSig,targetValuesBkg),axis = 0)
varsTotal = np.concatenate((varValuesSig,varValuesBkg),axis = 1)

#Do Test/Train split
#THIS RETURNS A TRANSPOSED VERSION OF INPUT AND VARIABLE ARRRAYS
xTrain, xTest, wTrain, dummy, yTrain, yTest, dummy, wTest, varValuesTrain, varValuesTest = train_test_split(inputTotal.T,weightsTotal,targetsTotal,weightsTotal,varsTotal.T,test_size=0.25)

#Split into lead/sub arrays, so we can make single-photon eta cuts (to split into B and E)
#xLeadTrain = xTrain[:,0:10]
#xSubTrain = xTrain[:,10:20]
#xLeadTest = xTest[:,0:10]
#xSubTest = xTest[:,10:20]
print(xTrain[0,0:20])
xLeadTrain = xTrain[:,0:9]
xSubTrain = xTrain[:,9:18]
xLeadTest = xTest[:,0:9]
xSubTest = xTest[:,9:18]

varsLeadTrain = varValuesTrain[:,0:19]
varsSubTrain = varValuesTrain[:,19:38]
varsLeadTest = varValuesTest[:,0:19]
varsSubTest = varValuesTest[:,19:38]

yLeadTrain = yTrain
ySubTrain = yTrain
yLeadTest = yTest
ySubTest = yTest

wLeadTrain = wTrain
wSubTrain = wTrain
wLeadTest = wTest
wSubTest = wTest

#Make Eta cuts on single photon arrays to put each photon into B or E
xLeadTrainB = xLeadTrain[np.abs(xLeadTrain[:,6]) < 1.5]
xSubTrainB = xSubTrain[np.abs(xSubTrain[:,6]) < 1.5]
xLeadTrainE = xLeadTrain[np.abs(xLeadTrain[:,6]) > 1.5]
xSubTrainE = xSubTrain[np.abs(xSubTrain[:,6]) > 1.5]

xLeadTestB = xLeadTest[np.abs(xLeadTest[:,6]) < 1.5]
xSubTestB = xSubTest[np.abs(xSubTest[:,6]) < 1.5]
xLeadTestE = xLeadTest[np.abs(xLeadTest[:,6]) > 1.5]
xSubTestE = xSubTest[np.abs(xSubTest[:,6]) > 1.5]

varsLeadTrainB = varsLeadTrain[np.abs(xLeadTrain[:,6]) < 1.5]
varsSubTrainB = varsSubTrain[np.abs(xSubTrain[:,6]) < 1.5]
varsLeadTrainE = varsLeadTrain[np.abs(xLeadTrain[:,6]) > 1.5]
varsSubTrainE = varsSubTrain[np.abs(xSubTrain[:,6]) > 1.5]

varsLeadTestB = varsLeadTest[np.abs(xLeadTest[:,6]) < 1.5]
varsSubTestB = varsSubTest[np.abs(xSubTest[:,6]) < 1.5]
varsLeadTestE = varsLeadTest[np.abs(xLeadTest[:,6]) > 1.5]
varsSubTestE = varsSubTest[np.abs(xSubTest[:,6]) > 1.5]

yLeadTrainB = yLeadTrain[np.abs(xLeadTrain[:,6]) < 1.5]
ySubTrainB = ySubTrain[np.abs(xSubTrain[:,6]) < 1.5]
yLeadTrainE = yLeadTrain[np.abs(xLeadTrain[:,6]) > 1.5]
ySubTrainE = ySubTrain[np.abs(xSubTrain[:,6]) > 1.5]

yLeadTestB = yLeadTest[np.abs(xLeadTest[:,6]) < 1.5]
ySubTestB = ySubTest[np.abs(xSubTest[:,6]) < 1.5]
yLeadTestE = yLeadTest[np.abs(xLeadTest[:,6]) > 1.5]
ySubTestE = ySubTest[np.abs(xSubTest[:,6]) > 1.5]

wLeadTrainB = wLeadTrain[np.abs(xLeadTrain[:,6]) < 1.5]
wSubTrainB = wSubTrain[np.abs(xSubTrain[:,6]) < 1.5]
wLeadTrainE = wLeadTrain[np.abs(xLeadTrain[:,6]) > 1.5]
wSubTrainE = wSubTrain[np.abs(xSubTrain[:,6]) > 1.5]

wLeadTestB = wLeadTest[np.abs(xLeadTest[:,6]) < 1.5]
wSubTestB = wSubTest[np.abs(xSubTest[:,6]) < 1.5]
wLeadTestE = wLeadTest[np.abs(xLeadTest[:,6]) > 1.5]
wSubTestE = wSubTest[np.abs(xSubTest[:,6]) > 1.5]

#Now combine lead and sublead for each eta region, making final training sample
xTrainB = np.concatenate((xLeadTrainB,xSubTrainB),axis = 0)
xTrainE = np.concatenate((xLeadTrainE,xSubTrainE),axis = 0)
xTestB = np.concatenate((xLeadTestB,xSubTestB),axis = 0)
xTestE = np.concatenate((xLeadTestE,xSubTestE),axis = 0)

varsTrainB = np.concatenate((varsLeadTrainB,varsSubTrainB),axis = 0)
varsTrainE = np.concatenate((varsLeadTrainE,varsSubTrainE),axis = 0)
varsTestB = np.concatenate((varsLeadTestB,varsSubTestB),axis = 0)
varsTestE = np.concatenate((varsLeadTestE,varsSubTestE),axis = 0)

yTrainB = np.concatenate((yLeadTrainB,ySubTrainB),axis = 0)
yTrainE = np.concatenate((yLeadTrainE,ySubTrainE),axis = 0)
yTestB = np.concatenate((yLeadTestB,ySubTestB),axis = 0)
yTestE = np.concatenate((yLeadTestE,ySubTestE),axis = 0)

wTrainB = np.concatenate((wLeadTrainB,wSubTrainB),axis = 0)
wTrainE = np.concatenate((wLeadTrainE,wSubTrainE),axis = 0)
wTestB = np.concatenate((wLeadTestB,wSubTestB),axis = 0)
wTestE = np.concatenate((wLeadTestE,wSubTestE),axis = 0)

mdB = 7
lrB = 0.25
nEstB = 1500

mdE = mdB
lrE = lrB
nEstE = nEstB

#bkgToSig = (158104 + 116852)/(27776 + 10900)
bkgToSig = 1
MDS = 0

startTime = time.time()

modelB = XGBClassifier(max_depth = mdB, learning_rate = lrB, early_stopping_rounds=50, scale_pos_weight = bkgToSig, subsample = 0.9, max_delta_step = MDS, min_child_weight = 0.0, reg_alpha = 0.0, reg_lambda = 4.0, colsample_bytree = 0.65, verbosity = 0, n_estimators=nEstB, nthread = 16, feature_names = inputVarNames, eval_metric=["auc", "logloss"],)
#modelB = XGBClassifier(max_depth = mdB, learning_rate = lrB, scale_pos_weight = bkgToSig, subsample = 0.9, max_delta_step = MDS, min_child_weight = 0.0, reg_alpha = 0.0, reg_lambda = 4.0, colsample_bytree = 0.65, verbosity = 0, n_estimators=nEstB, nthread = 16, feature_names = inputVarNames, tree_method = 'hist')
evalSetB = [(xTrainB, yTrainB), (xTestB, yTestB)]
modelB.fit(xTrainB, yTrainB, sample_weight = wTrainB, eval_set=evalSetB, verbose=1000)

modelE = XGBClassifier(max_depth = mdE, learning_rate = lrE, early_stopping_rounds=50, scale_pos_weight = bkgToSig, subsample = 0.9, max_delta_step = MDS, min_child_weight = 0.0, reg_alpha = 0.0, reg_lambda = 4.0, colsample_bytree = 0.65, verbosity = 0, n_estimators=nEstB, nthread = 16, feature_names = inputVarNames, eval_metric=["auc", "logloss"],)
#modelE = XGBClassifier(max_depth = mdE, learning_rate = lrE, scale_pos_weight = bkgToSig, subsample = 0.9, max_delta_step = MDS, min_child_weight = 0.0, reg_alpha = 0.0, reg_lambda = 4.0, colsample_bytree = 0.65, verbosity = 0, n_estimators=nEstB, nthread = 16, feature_names = inputVarNames, tree_method = 'hist')
evalSetE = [(xTrainE, yTrainE), (xTestE, yTestE)]
modelE.fit(xTrainE, yTrainE, sample_weight = wTrainE, eval_set=evalSetE, verbose=1000)

modelB.save_model(modelNamesB[0]+".bin")
modelE.save_model(modelNamesE[0]+".bin")

gc.collect()
pickle.dump(modelB, open((modelNamesB[0]+".Model"), "wb"))
pickle.dump(modelE, open((modelNamesE[0]+".Model"), "wb"))

diphotonUtils.convert_model(modelB.get_booster().get_dump(), input_variables = inputVarNames, output_xml=(modelNamesB[0]+".xml"))
diphotonUtils.convert_model(modelE.get_booster().get_dump(), input_variables = inputVarNames, output_xml=(modelNamesE[0]+".xml"))

##18 is the event index determined before splitting. Allows Lead and Sub to be recombined after evaluation
##So XGB results are an array with the score result and the event index to be recombined into diphoton event
yPredTrainB = np.stack((modelB.predict_proba(xTrainB)[:,1],varsTrainB[:,18]),axis=1)
yPredTestB = np.stack((modelB.predict_proba(xTestB)[:,1],varsTestB[:,18]),axis=1)

yPredTrainE = np.stack((modelE.predict_proba(xTrainE)[:,1],varsTrainE[:,18]),axis=1)
yPredTestE = np.stack((modelE.predict_proba(xTestE)[:,1],varsTestE[:,18]),axis=1)

#Split into sig/bkg arrays by label
yPredTrainSigB = yPredTrainB[yTrainB[:] == 1]
yPredTrainBkgB = yPredTrainB[yTrainB[:] != 1]
yPredTrainSigE = yPredTrainE[yTrainE[:] == 1]
yPredTrainBkgE = yPredTrainE[yTrainE[:] != 1]

yPredTestSigB = yPredTestB[yTestB[:] == 1]
yPredTestBkgB = yPredTestB[yTestB[:] != 1]
yPredTestSigE = yPredTestE[yTestE[:] == 1]
yPredTestBkgE = yPredTestE[yTestE[:] != 1]

varsTrainSigB = varsTrainB[yTrainB[:] == 1]
varsTrainBkgB = varsTrainB[yTrainB[:] != 1]
varsTrainSigE = varsTrainE[yTrainE[:] == 1]
varsTrainBkgE = varsTrainE[yTrainE[:] != 1]

varsTestSigB = varsTestB[yTestB[:] == 1]
varsTestBkgB = varsTestB[yTestB[:] != 1]
varsTestSigE = varsTestE[yTestE[:] == 1]
varsTestBkgE = varsTestE[yTestE[:] != 1]

#Add barrel and endcap for sig and bkg seperately
yPredTrainSigTotal = np.concatenate((yPredTrainSigB,yPredTrainSigE),axis = 0)
yPredTrainBkgTotal = np.concatenate((yPredTrainBkgB,yPredTrainBkgE),axis = 0)

yPredTestSigTotal = np.concatenate((yPredTestSigB,yPredTestSigE),axis = 0)
yPredTestBkgTotal = np.concatenate((yPredTestBkgB,yPredTestBkgE),axis = 0)

varsTrainSigTotal = np.concatenate((varsTrainSigB,varsTrainSigE), axis = 0)
varsTrainBkgTotal = np.concatenate((varsTrainBkgB,varsTrainBkgE), axis = 0)

varsTestSigTotal = np.concatenate((varsTestSigB,varsTestSigE), axis = 0)
varsTestBkgTotal = np.concatenate((varsTestBkgB,varsTestBkgE), axis = 0)

#Sort the sig and bkg by event index, so events with same index are next to each other
yPredTrainSigTotal = sorted(yPredTrainSigTotal, key=itemgetter(1))
yPredTrainBkgTotal = sorted(yPredTrainBkgTotal, key=itemgetter(1))

yPredTestSigTotal = sorted(yPredTestSigTotal, key=itemgetter(1))
yPredTestBkgTotal = sorted(yPredTestBkgTotal, key=itemgetter(1))

varsTrainSigTotal = sorted(varsTrainSigTotal, key=itemgetter(18))
varsTrainBkgTotal = sorted(varsTrainBkgTotal, key=itemgetter(18))

varsTestSigTotal = sorted(varsTestSigTotal, key=itemgetter(18))
varsTestBkgTotal = sorted(varsTestBkgTotal, key=itemgetter(18))

print("Sample Sizes Train:")
print('Train Signal Barrel =  {0:.0f}'.format(len(yPredTrainSigB)))
print('Train Signal Endcap =  {0:.0f}'.format(len(yPredTrainSigE)))
print('Train Bkg Barrel =  {0:.0f}'.format(len(yPredTrainBkgB)))
print('Train Bkg Endcap =  {0:.0f}'.format(len(yPredTrainBkgE)))
print("Sample Sizes Test:")
print('Test Signal Barrel =  {0:.0f}'.format(len(yPredTestSigB)))
print('Test Signal Endcap =  {0:.0f}'.format(len(yPredTestSigE)))
print('Test Bkg Barrel =  {0:.0f}'.format(len(yPredTestBkgB)))
print('Test Bkg Endcap =  {0:.0f}'.format(len(yPredTestBkgE)))


endTime = time.time()
timeSpent = endTime - startTime
timeSpentMin = timeSpent/60.0
print('Run Time(min) =  {0:.4f}'.format(timeSpentMin))
        
outFileRoot = TFile(fNames[0],'recreate')

buildBranch("sigTrain",varsTrainSigTotal,yPredTrainSigTotal,outFileRoot)
buildBranch("bkgTrain",varsTrainBkgTotal,yPredTrainBkgTotal,outFileRoot)
buildBranch("sigTest",varsTestSigTotal,yPredTestSigTotal,outFileRoot)
buildBranch("bkgTest",varsTestBkgTotal,yPredTestBkgTotal,outFileRoot)

outFileRoot.Close()
