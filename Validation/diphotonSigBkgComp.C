#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCut.h"
#include "TPaletteAxis.h"
#include "TMVA/Reader.h"
#include "./pairSelection.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
double newGetWidth(TH1F *histo);
void diphotonSigBkgComp(){
    gROOT->Reset();
    gStyle->SetPalette(1);
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitle(0);
    //gErrorIgnoreLevel = kFatal;
    
    string fileName = "validationNTuples/0115/M7L25_GGH13andDataD_DiphotonV6_NoTrkIso_M60.root";
    //string fileName = "validationNTuples/0725/M9LR15_GGH13andData_NoTrkIso_M60_0725Validated.root";
    string genTitleStringSignal = "GGH Signal ";
    string genTitleStringBkg = "Data ";
    string dirStr ="sigBkgVarPlots/0116/MD7L25_GGH13andDataD_NoTrkIso_M60_Added/";
    string plotType = "Added";
    //string plotType = "Separate";
    //string dirStr ="varPlots/0130/DataRelaxedUnseededChoose2/";
        
    int binFactor = 1;
    
    int nEta;
    
    nEta = 4;
    string etaLabels[4] = {"Barrel-Barrel","Barrel-Endcap","Endcap-Endcap","All #eta"};
    string etaFLabels[4] = {"_BB","_BE","_EE","_All"};

    //New cuts (11/09) chosen by hand (just Barrel and Endcap)
    double leadCuts1[2] = {0.85,0.90};
    double subCuts1[2] = {0.0,0.0};//If First score above leadCuts1, second score must be above subCuts1
    double leadCuts2[2] = {0.75,0.8};
    double subCuts2[2] = {0.015,0.02};//If First score between leadCuts1 and leadCuts2, second score must be above subCuts2
    double subCuts3[2] = {0.075,0.075};//If First score below leadCuts2, second score must be above subCuts3
    
    string outNameGen = dirStr;
    
    int nVars = 21;
    string varNames[21] = {"nRun","nEgs","nEgsPassing","triggerBits","mass","dR","et","energy","rawEnergy","eta","etaWidth","phi","phiWidth","sigmaIEtaIEta","ecalPFIso","r9HLT","s4","hOvrE","trkIsoPho","xgbScore","etOvrM"};
    double limsLow[21] = {369800,1,1,0.0,0.0,0.0,0.0,0.0,0.0,-3,0.0,-3.15,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    double limsHigh[21] = {371000,10,10,4,250,6,300,600,600,3,0.1,3.15,0.2,0.1,40.0,1.2,1.2,1.0,40.0,1.0,1.5};
    
    int nBinsSig[21] = {1200,9,9,4,50000,250,250,250,250,75,250,200,250,250,250,250,250,250,250,250,250};
    int nBinsBkg[21] = {1200,9,9,4,125,250,125,250,250,75,250,250,250,250,250,125,250,250,250,250,125};
    
    for (int i = 3; i < nVars; i++)nBinsSig[i] = nBinsSig[i]/binFactor;
    for (int i = 3; i < nVars; i++)nBinsBkg[i] = nBinsBkg[i]/binFactor;
    
    TCanvas *can = new TCanvas ("can","can",10,10,1600,900);
    for(int i = 0; i < nVars; i++){
    //for(int i = 4; i < 5; i++){//Mass
    //for(int i = 9; i < 10; i++){//Eta

        TH1F *hSigAllLead[4];
        TH1F *hSigAllSub[4];
        TH1F *hBkgAllLead[4];
        TH1F *hBkgAllSub[4];
        
        TH1F *hSigFailXGBLead[4];
        TH1F *hSigFailXGBSub[4];
        TH1F *hBkgFailXGBLead[4];
        TH1F *hBkgFailXGBSub[4];
        
        TH1F *hSigPassTightXGBLead[4];
        TH1F *hSigPassTightXGBSub[4];
        TH1F *hBkgPassTightXGBLead[4];
        TH1F *hBkgPassTightXGBSub[4];
        
        TH1F *hSigPassTightXGBTightMassLead[4];
        TH1F *hSigPassTightXGBTightMassSub[4];
        TH1F *hBkgPassTightXGBTightMassLead[4];
        TH1F *hBkgPassTightXGBTightMassSub[4];

        TH1F *hSigPassStdLead[4];
        TH1F *hSigPassStdSub[4];
        TH1F *hBkgPassStdLead[4];
        TH1F *hBkgPassStdSub[4];
        
        for (int e = 0; e < nEta; e++){
            hSigAllLead[e] = new TH1F ("hSigAllLead","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hSigAllSub[e] = new TH1F ("hSigAllSub","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hBkgAllLead[e] = new TH1F ("hBkgAllLead","",nBinsBkg[i],limsLow[i],limsHigh[i]);
            hBkgAllSub[e] = new TH1F ("hBkgAllSub","",nBinsBkg[i],limsLow[i],limsHigh[i]);

            hSigFailXGBLead[e] = new TH1F ("hSigFailXGBLead","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hSigFailXGBSub[e] = new TH1F ("hSigFailXGBSub","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hBkgFailXGBLead[e] = new TH1F ("hBkgFailXGBLead","",nBinsBkg[i],limsLow[i],limsHigh[i]);
            hBkgFailXGBSub[e] = new TH1F ("hBkgFailXGBSub","",nBinsBkg[i],limsLow[i],limsHigh[i]);

            hSigPassTightXGBLead[e] = new TH1F ("hSigPassTightXGBLead","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hSigPassTightXGBSub[e] = new TH1F ("hSigPassTightXGBSub","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hBkgPassTightXGBLead[e] = new TH1F ("hBkgPassTightXGBLead","",nBinsBkg[i],limsLow[i],limsHigh[i]);
            hBkgPassTightXGBSub[e] = new TH1F ("hBkgPassTightXGBSub","",nBinsBkg[i],limsLow[i],limsHigh[i]);
            
            hSigPassTightXGBTightMassLead[e] = new TH1F ("hSigPassTightXGBTightMassLead","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hSigPassTightXGBTightMassSub[e] = new TH1F ("hSigPassTightXGBTightMassSub","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hBkgPassTightXGBTightMassLead[e] = new TH1F ("hBkgPassTightXGBTightMassLead","",nBinsBkg[i],limsLow[i],limsHigh[i]);
            hBkgPassTightXGBTightMassSub[e] = new TH1F ("hBkgPassTightXGBTightMassSub","",nBinsBkg[i],limsLow[i],limsHigh[i]);
           
            hSigPassStdLead[e] = new TH1F ("hSigPassStdLead","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hSigPassStdSub[e] = new TH1F ("hSigPassStdSub","",nBinsSig[i],limsLow[i],limsHigh[i]);
            hBkgPassStdLead[e] = new TH1F ("hBkgPassStdLead","",nBinsBkg[i],limsLow[i],limsHigh[i]);
            hBkgPassStdSub[e] = new TH1F ("hBkgPassStdSub","",nBinsBkg[i],limsLow[i],limsHigh[i]);
        }
        
        int nPassEta = 0;
        int nFailEta = 0;
        int nPassEtaSig = 0;
        int nFailEtaSig = 0;
        
        TFile *f = new TFile(fileName.c_str());
        
        string varNameTmp;
        varNameTmp = varNames[i];
        string varNameTmpHere = varNameTmp;
        int isSinglePhoton = 0;
        if (varNameTmp == "nRun" || varNameTmp == "nEgs" || varNameTmp == "nEgsPassing" || varNameTmp == "triggerBits" || varNameTmp == "mass" || varNameTmp == "dR")isSinglePhoton = 1;
        if (isSinglePhoton ==  1 || varNameTmp == "etOvrM")varNameTmpHere = "eta";
        
        //FIRST handle Signal
        TTreeReader inReaderSig("sigTree", f);
        inReaderSig.Restart();
        
        TTreeReaderArray<Float_t> inVarSig(inReaderSig, varNameTmpHere.c_str());
        TTreeReaderArray<Float_t> inESig(inReaderSig,"energy");
        TTreeReaderArray<Float_t> inEtSig(inReaderSig,"et");
        TTreeReaderArray<Float_t> inEtaSig(inReaderSig, "eta");
        TTreeReaderArray<Float_t> inPhiSig(inReaderSig, "phi");
        TTreeReaderArray<Float_t> inR9Sig(inReaderSig,"r9HLT");
        TTreeReaderArray<Float_t> inHovrESig(inReaderSig,"hOvrE");
        TTreeReaderArray<Float_t> inSigIEtaIEtaSig(inReaderSig,"sigmaIEtaIEta");
        TTreeReaderArray<Float_t> inPhoIsoSig(inReaderSig,"trkIsoPho");
        TTreeReaderArray<Float_t> inEcalIsoSig(inReaderSig,"ecalPFIso");
        TTreeReaderArray<Float_t> inXGBSig(inReaderSig, "xgbScore");
        TTreeReaderArray<int> inPdgIDSig(inReaderSig,"pdgID");

        TTreeReaderValue<int> inNEgsSig(inReaderSig,"nEgs");
        TTreeReaderValue<int> inNRunSig(inReaderSig,"nRun");
        TTreeReaderValue<int> inNEgsPassingSig(inReaderSig,"nEgsPassing");
        TTreeReaderValue<int> inTriggerBitsSig(inReaderSig,"triggerBits");
        TTreeReaderValue<int> inPassFailStdSig(inReaderSig,"passFailStd");
        TTreeReaderValue<int> inPassFailDoubleL1Sig(inReaderSig,"passFailL1Double");
        TTreeReaderValue<int> inPassFailSingleL1Sig(inReaderSig,"passFailL1Single");
        
        //NORMAL ETA THRESHOLDS ARE 1.444 and 1.556
        int nTotalSig = 0;
        int nDisagreeSig = 0;
        int nNoneChosenSig = 0;
        int nFailingDRSig = 0;
        int nFailingM50Sig = 0;

        while (inReaderSig.Next()) {
            int iOutXGBSig,jOutXGBSig;
            double massCalcXGBSig = calcMassXGB(&inXGBSig,0.25,0.25,&inEtSig,&inEtaSig,&inPhiSig,*inNEgsSig,&iOutXGBSig,&jOutXGBSig);
            double drXGBSig = calcDR(inEtaSig[iOutXGBSig], inEtaSig[jOutXGBSig], inPhiSig[iOutXGBSig], inPhiSig[jOutXGBSig]);
            if (iOutXGBSig == -1 || jOutXGBSig == -1)nNoneChosenSig+=1;
            //if(inEtSig[iOutXGBSig] > 14.25 && inEtSig[jOutXGBSig] > 14.25 && abs(inEtaSig[iOutXGBSig]) < 2.55 && abs(inEtaSig[jOutXGBSig]) < 2.55 && iOutXGBSig != -1 && jOutXGBSig != -1){
            if(inPdgIDSig[iOutXGBSig] == 22 && inPdgIDSig[jOutXGBSig] == 22 && inEtSig[iOutXGBSig] > 14.25 && inEtSig[jOutXGBSig] > 14.25 && abs(inEtaSig[iOutXGBSig]) < 2.55 && abs(inEtaSig[jOutXGBSig]) < 2.55 && iOutXGBSig != -1 && jOutXGBSig != -1){//Added pdgID requirement for signal
                int indexHigh,indexLow;
                if (inXGBSig[iOutXGBSig] > inXGBSig[jOutXGBSig]){
                    indexHigh = iOutXGBSig;
                    indexLow = jOutXGBSig;
                }
                else if (inXGBSig[iOutXGBSig] < inXGBSig[jOutXGBSig]){
                    indexHigh = jOutXGBSig;
                    indexLow = iOutXGBSig;
                }
        
                int indexLead,indexSub;
                if (inEtSig[iOutXGBSig] > inEtSig[jOutXGBSig]){
                    indexLead = iOutXGBSig;
                    indexSub = jOutXGBSig;
                }
                else if (inEtSig[iOutXGBSig] < inEtSig[jOutXGBSig]){
                    indexLead = jOutXGBSig;
                    indexSub = iOutXGBSig;
                }
                double outOne,outTwo;
                //if (varNameTmp != "mass" && varNameTmp != "nEgs" && varNameTmp != "triggerBits" && varNameTmp != "dR" && varNameTmp != "etOvrM"){
                if (isSinglePhoton == 0 && varNameTmp != "etOvrM"){
                    outOne = inVarSig[indexLead];
                    outTwo = inVarSig[indexSub];
                }
                else if (varNameTmp == "mass"){
                    outOne = massCalcXGBSig;
                    outTwo = massCalcXGBSig;
                }
                else if (varNameTmp == "nRun"){
                    outOne = (double)(*inNRunSig);
                    outTwo = (double)(*inNRunSig);
                    
                }
                else if (varNameTmp == "nEgs"){
                    outOne = (double)(*inNEgsSig);
                    outTwo = (double)(*inNEgsSig);
                    
                }
                else if (varNameTmp == "nEgsPassing"){
                    outOne = (double)(*inNEgsPassingSig);
                    outTwo = (double)(*inNEgsPassingSig);
                }
                else if (varNameTmp == "triggerBits"){
                    outOne = (double)(*inTriggerBitsSig);
                    outTwo = (double)(*inTriggerBitsSig);
                    
                }
                else if (varNameTmp == "dR"){
                    outOne = drXGBSig;
                    outTwo = drXGBSig;
                    
                }
                else if (varNameTmp == "etOvrM"){
                    outOne = inEtSig[indexLead]/massCalcXGBSig;
                    outTwo = inEtSig[indexSub]/massCalcXGBSig;
                }
                else if (varNameTmp == "xgbScore"){
                    outOne = inXGBSig[indexHigh];
                    outTwo = inXGBSig[indexLow];
                }
                
                int e;//Combined eta index -- 0=BB, 1=BE, 2=EE
                int e1,e2;//Individual eta index -- 0=Barrel, 1=Endcap
                if (abs(inEtaSig[indexHigh]) < 1.5 && abs(inEtaSig[indexLow]) < 1.5){
                    e = 0;
                    e1 = 0;//Barrel
                    e2 = 0;//Barrel
                }
                if (abs(inEtaSig[indexHigh]) < 1.5 && abs(inEtaSig[indexLow]) > 1.5){
                    e = 1;
                    e1 = 0;//Barrel
                    e2 = 1;//Endcap
                }
                if(abs(inEtaSig[indexHigh]) > 1.5 && abs(inEtaSig[indexLow]) < 1.5){
                    e = 1;
                    e1 = 1;//Endcap
                    e2 = 0;//Barrel
                }
                if (abs(inEtaSig[indexHigh]) > 1.5 && abs(inEtaSig[indexLow]) > 1.5){
                    e = 2;
                    e1 = 1;//Endcap
                    e2 = 1;//Endcap
                }
                hSigAllLead[e]->Fill(outOne);
                hSigAllSub[e]->Fill(outTwo);
                
                if (drXGBSig < 0.5)nFailingDRSig+=1;
                if (drXGBSig > 0.5 && massCalcXGBSig < 50)nFailingM50Sig+=1;
                
                if((inXGBSig[indexHigh] > leadCuts1[e1] && inXGBSig[indexLow] > subCuts1[e2])
                   || (inXGBSig[indexHigh] > leadCuts2[e1] && inXGBSig[indexLow] > subCuts2[e2])
                   || (inXGBSig[indexHigh] < leadCuts2[e1] && inXGBSig[indexLow] > subCuts3[e2])
                   ){
                    if(massCalcXGBSig > 50.0){
                        hSigPassTightXGBLead[e]->Fill(outOne);
                        hSigPassTightXGBSub[e]->Fill(outTwo);
                    }
                    if(massCalcXGBSig > 95.0){
                        hSigPassTightXGBTightMassLead[e]->Fill(outOne);
                        hSigPassTightXGBTightMassSub[e]->Fill(outTwo);
                    }
                }
                else{
                    hSigFailXGBLead[e]->Fill(outOne);
                    hSigFailXGBSub[e]->Fill(outTwo);
                }

            }//XGB Eta selection
            
            int iOutHLTSig,jOutHLTSig;
            double massCalcHLTSig = calcMassHLT(&inEtSig, &inEtaSig, &inPhiSig, &inR9Sig, &inHovrESig, &inSigIEtaIEtaSig, &inPhoIsoSig, &inEcalIsoSig, *inNEgsSig, &iOutHLTSig, &jOutHLTSig);
            double drHLTSig = calcDR(inEtaSig[iOutHLTSig], inEtaSig[jOutHLTSig], inPhiSig[iOutHLTSig], inPhiSig[jOutHLTSig]);
            //if(iOutHLTSig != -1 && jOutHLTSig != -1){
            if(inPdgIDSig[iOutHLTSig] == 22 && inPdgIDSig[jOutHLTSig] == 22 && iOutHLTSig != -1 && jOutHLTSig != -1){//Added pdgID requirement
                double outOne,outTwo;
                if (isSinglePhoton == 0 && varNameTmp != "etOvrM"){
                    outOne = inVarSig[iOutHLTSig];
                    outTwo = inVarSig[jOutHLTSig];
                }
                else if (varNameTmp == "mass"){
                    outOne = massCalcHLTSig;
                    outTwo = massCalcHLTSig;
                }
                else if (varNameTmp == "nRun"){
                    outOne = (double)(*inNRunSig);
                    outTwo = (double)(*inNRunSig);
                }
                else if (varNameTmp == "nEgs"){
                    outOne = (double)(*inNEgsSig);
                    outTwo = (double)(*inNEgsSig);
                    
                }
                else if (varNameTmp == "nEgsPassing"){
                    outOne = (double)(*inNEgsPassingSig);
                    outTwo = (double)(*inNEgsPassingSig);
                }
                else if (varNameTmp == "triggerBits"){
                    outOne = (double)(*inTriggerBitsSig);
                    outTwo = (double)(*inTriggerBitsSig);
                    
                }
                else if (varNameTmp == "dR"){
                    outOne = drHLTSig;
                    outTwo = drHLTSig;
                }
                else if (varNameTmp == "etOvrM"){
                    outOne = inEtSig[iOutHLTSig]/massCalcHLTSig;
                    outTwo = inEtSig[jOutHLTSig]/massCalcHLTSig;
                }
                int e;
                if (abs(inEtaSig[iOutHLTSig]) < 1.5 && abs(inEtaSig[jOutHLTSig]) < 1.5)e = 0;
                if ((abs(inEtaSig[iOutHLTSig]) > 1.5 && abs(inEtaSig[jOutHLTSig]) < 1.5)
                    || (abs(inEtaSig[iOutHLTSig]) < 1.5 && abs(inEtaSig[jOutHLTSig]) > 1.5))e = 1;
                if (abs(inEtaSig[iOutHLTSig]) > 1.5 && abs(inEtaSig[jOutHLTSig]) > 1.5)e = 2;
                
                if(*inPassFailStdSig == 1 && (*inPassFailSingleL1Sig == 1 || *inPassFailDoubleL1Sig == 1)){
                    hSigPassStdLead[e]->Fill(outOne);
                    hSigPassStdSub[e]->Fill(outTwo);
                }//If HLT == 1 && L1 == 1
            }//HLT Eta Selection if
            nTotalSig+=1;
            if (iOutHLTSig != iOutXGBSig || jOutHLTSig != jOutXGBSig)nDisagreeSig+=1;
        }//While inReaderSig.Next()
        
        
        //THEN handle Bkg
        TTreeReader inReaderBkg("bkgTree", f);
        inReaderBkg.Restart();
        
        TTreeReaderArray<Float_t> inVarBkg(inReaderBkg, varNameTmpHere.c_str());
        TTreeReaderArray<Float_t> inEBkg(inReaderBkg,"energy");
        TTreeReaderArray<Float_t> inEtBkg(inReaderBkg,"et");
        TTreeReaderArray<Float_t> inEtaBkg(inReaderBkg, "eta");
        TTreeReaderArray<Float_t> inPhiBkg(inReaderBkg, "phi");
        TTreeReaderArray<Float_t> inR9Bkg(inReaderBkg,"r9HLT");
        TTreeReaderArray<Float_t> inHovrEBkg(inReaderBkg,"hOvrE");
        TTreeReaderArray<Float_t> inSigIEtaIEtaBkg(inReaderBkg,"sigmaIEtaIEta");
        TTreeReaderArray<Float_t> inPhoIsoBkg(inReaderBkg,"trkIsoPho");
        TTreeReaderArray<Float_t> inEcalIsoBkg(inReaderBkg,"ecalPFIso");
        TTreeReaderArray<Float_t> inXGBBkg(inReaderBkg, "xgbScore");
        
        TTreeReaderValue<int> inNEgsBkg(inReaderBkg,"nEgs");
        TTreeReaderValue<int> inNRunBkg(inReaderBkg,"nRun");
        TTreeReaderValue<int> inNEgsPassingBkg(inReaderBkg,"nEgsPassing");
        TTreeReaderValue<int> inTriggerBitsBkg(inReaderBkg,"triggerBits");
        TTreeReaderValue<int> inPassFailStdBkg(inReaderBkg,"passFailStd");
        TTreeReaderValue<int> inPassFailDoubleL1Bkg(inReaderBkg,"passFailL1Double");
        TTreeReaderValue<int> inPassFailSingleL1Bkg(inReaderBkg,"passFailL1Single");
        //NORMAL ETA THRESHOLDS ARE 1.444 and 1.556
        
        int nTotalBkg = 0;
        int nDisagreeBkg = 0;
        int nNoneChosenBkg = 0;
        int nFailingDRBkg = 0;
        int nFailingM50Bkg = 0;
        while (inReaderBkg.Next()) {
            int iOutXGBBkg,jOutXGBBkg;
            double massCalcXGBBkg = calcMassXGB(&inXGBBkg,0.25,0.25,&inEtBkg,&inEtaBkg,&inPhiBkg,*inNEgsBkg,&iOutXGBBkg,&jOutXGBBkg);
            double drXGBBkg = calcDR(inEtaBkg[iOutXGBBkg], inEtaBkg[jOutXGBBkg], inPhiBkg[iOutXGBBkg], inPhiBkg[jOutXGBBkg]);
            if (iOutXGBBkg == -1 || jOutXGBBkg == -1)nNoneChosenBkg+=1;
            if(inEtBkg[iOutXGBBkg] > 14.25 && inEtBkg[jOutXGBBkg] > 14.25 && abs(inEtaBkg[iOutXGBBkg]) < 2.55 && abs(inEtaBkg[jOutXGBBkg]) < 2.55 && iOutXGBBkg != -1 && jOutXGBBkg != -1){
                int indexHigh,indexLow;
                if (inXGBBkg[iOutXGBBkg] > inXGBBkg[jOutXGBBkg]){
                    indexHigh = iOutXGBBkg;
                    indexLow = jOutXGBBkg;
                }
                else if (inXGBBkg[iOutXGBBkg] < inXGBBkg[jOutXGBBkg]){
                    indexHigh = jOutXGBBkg;
                    indexLow = iOutXGBBkg;
                }

                int indexLead,indexSub;
                if (inEtBkg[iOutXGBBkg] > inEtBkg[jOutXGBBkg]){
                    indexLead = iOutXGBBkg;
                    indexSub = jOutXGBBkg;
                }
                else if (inEtBkg[iOutXGBBkg] < inEtBkg[jOutXGBBkg]){
                    indexLead = jOutXGBBkg;
                    indexSub = iOutXGBBkg;
                }
                
                double outOne,outTwo;
                if (isSinglePhoton == 0 && varNameTmp != "etOvrM"){
                    outOne = inVarBkg[indexLead];
                    outTwo = inVarBkg[indexSub];
                }
                else if (varNameTmp == "mass"){
                    outOne = massCalcXGBBkg;
                    outTwo = massCalcXGBBkg;
                }
                else if (varNameTmp == "nRun"){
                    outOne = (double)(*inNRunBkg);
                    outTwo = (double)(*inNRunBkg);
                    
                }
                else if (varNameTmp == "nEgs"){
                    outOne = (double)(*inNEgsBkg);
                    outTwo = (double)(*inNEgsBkg);
                    
                }
                else if (varNameTmp == "nEgsPassing"){
                    outOne = (double)(*inNEgsPassingBkg);
                    outTwo = (double)(*inNEgsPassingBkg);
                }
                else if (varNameTmp == "triggerBits"){
                    outOne = (double)(*inTriggerBitsBkg);
                    outTwo = (double)(*inTriggerBitsBkg);
                }
                else if (varNameTmp == "dR"){
                    outOne = drXGBBkg;
                    outTwo = drXGBBkg;
                }
                else if (varNameTmp == "etOvrM"){
                    outOne = inEtBkg[indexLead]/massCalcXGBBkg;
                    outTwo = inEtBkg[indexSub]/massCalcXGBBkg;
                }
                else if (varNameTmp == "xgbScore"){
                    outOne = inXGBSig[indexHigh];
                    outTwo = inXGBSig[indexLow];
                }
                int e;//Combined eta index -- 0=BB, 1=BE, 2=EE
                int e1,e2;//Individual eta index -- 0=Barrel, 1=Endcap
                if (abs(inEtaBkg[indexHigh]) < 1.5 && abs(inEtaBkg[indexLow]) < 1.5){
                    e = 0;
                    e1 = 0;//Barrel
                    e2 = 0;//Barrel
                }
                if (abs(inEtaBkg[indexHigh]) < 1.5 && abs(inEtaBkg[indexLow]) > 1.5){
                    e = 1;
                    e1 = 0;//Barrel
                    e2 = 1;//Endcap
                }
                if(abs(inEtaBkg[indexHigh]) > 1.5 && abs(inEtaBkg[indexLow]) < 1.5){
                    e = 1;
                    e1 = 1;//Endcap
                    e2 = 0;//Barrel
                }
                if (abs(inEtaBkg[indexHigh]) > 1.5 && abs(inEtaBkg[indexLow]) > 1.5){
                    e = 2;
                    e1 = 1;//Endcap
                    e2 = 1;//Endcap
                }
                hBkgAllLead[e]->Fill(outOne);
                hBkgAllSub[e]->Fill(outTwo);
                
                if (drXGBBkg < 0.5)nFailingDRBkg+=1;
                if (drXGBBkg > 0.5 && massCalcXGBBkg < 50)nFailingM50Bkg+=1;
                
                if((inXGBBkg[indexHigh] > leadCuts1[e1] && inXGBBkg[indexLow] > subCuts1[e2])
                   || (inXGBBkg[indexHigh] > leadCuts2[e1] && inXGBBkg[indexLow] > subCuts2[e2])
                   || (inXGBBkg[indexHigh] < leadCuts2[e1] && inXGBBkg[indexLow] > subCuts3[e2])
                   ){
                    if(massCalcXGBBkg > 50.0){
                        hBkgPassTightXGBLead[e]->Fill(outOne);
                        hBkgPassTightXGBSub[e]->Fill(outTwo);
                    }
                    if(massCalcXGBBkg > 95.0){
                        hBkgPassTightXGBTightMassLead[e]->Fill(outOne);
                        hBkgPassTightXGBTightMassSub[e]->Fill(outTwo);
                    }
                }
                else{
                    hBkgFailXGBLead[e]->Fill(outOne);
                    hBkgFailXGBSub[e]->Fill(outTwo);
                }
            }//XGB Eta selection
            int iOutHLTBkg,jOutHLTBkg;
            double massCalcHLTBkg = calcMassHLT(&inEtBkg, &inEtaBkg, &inPhiBkg, &inR9Bkg, &inHovrEBkg, &inSigIEtaIEtaBkg, &inPhoIsoBkg, &inEcalIsoBkg, *inNEgsBkg, &iOutHLTBkg, &jOutHLTBkg);
            double drHLTBkg = calcDR(inEtaBkg[iOutHLTBkg], inEtaBkg[jOutHLTBkg], inPhiBkg[iOutHLTBkg], inPhiBkg[jOutHLTBkg]);
            if(iOutHLTBkg != -1 && jOutHLTBkg != -1){
                double outOne,outTwo;
                if (isSinglePhoton == 0 && varNameTmp != "etOvrM"){
                    outOne = inVarBkg[iOutHLTBkg];
                    outTwo = inVarBkg[jOutHLTBkg];
                }
                else if (varNameTmp == "mass"){
                    outOne = massCalcHLTBkg;
                    outTwo = massCalcHLTBkg;
                }
                else if (varNameTmp == "nRun"){
                    outOne = (double)(*inNRunBkg);
                    outTwo = (double)(*inNRunBkg);
                    
                }
                else if (varNameTmp == "nEgs"){
                    outOne = (double)(*inNEgsBkg);
                    outTwo = (double)(*inNEgsBkg);
                }
                else if (varNameTmp == "nEgsPassing"){
                    outOne = (double)(*inNEgsPassingBkg);
                    outTwo = (double)(*inNEgsPassingBkg);
                }
                else if (varNameTmp == "triggerBits"){
                    outOne = (double)(*inTriggerBitsBkg);
                    outTwo = (double)(*inTriggerBitsBkg);
                }
                else if (varNameTmp == "dR"){
                    outOne = drHLTBkg;
                    outTwo = drHLTBkg;
                    
                }
                else if (varNameTmp == "etOvrM"){
                    outOne = inEtBkg[iOutHLTBkg]/massCalcHLTBkg;
                    outTwo = inEtBkg[jOutHLTBkg]/massCalcHLTBkg;
                }
                int e;
                if (abs(inEtaBkg[iOutHLTBkg]) < 1.5 && abs(inEtaBkg[jOutHLTBkg]) < 1.5)e = 0;
                if ((abs(inEtaBkg[iOutHLTBkg]) > 1.5 && abs(inEtaBkg[jOutHLTBkg]) < 1.5)
                    || (abs(inEtaBkg[iOutHLTBkg]) < 1.5 && abs(inEtaBkg[jOutHLTBkg]) > 1.5))e = 1;
                if (abs(inEtaBkg[iOutHLTBkg]) > 1.5 && abs(inEtaBkg[jOutHLTBkg]) > 1.5)e = 2;
                if(*inPassFailStdBkg == 1 && (*inPassFailSingleL1Bkg == 1 || *inPassFailDoubleL1Bkg == 1)){
                    hBkgPassStdLead[e]->Fill(outOne);
                    hBkgPassStdSub[e]->Fill(outTwo);
                }
            }//HLT Eta Selection if
            nTotalBkg+=1;
            if (iOutHLTBkg != iOutXGBBkg || jOutHLTBkg != jOutXGBBkg)nDisagreeBkg+=1;
        }//While inReaderBkg.Next()
        cout<<"Sig N Failing DR Cut = "<<nFailingDRSig<<endl;
        cout<<"Sig N Failing M>50 Cut = "<<nFailingM50Sig<<endl;
        cout<<"Bkg N Failing DR Cut = "<<nFailingDRBkg<<endl;
        cout<<"Bkg N Failing M>50 Cut = "<<nFailingM50Bkg<<endl;

        
        
        for (int p = 0; p < 3; p++){
               
            hSigAllLead[3]->Add(hSigAllLead[p]);
            hSigAllSub[3]->Add(hSigAllSub[p]);
            hBkgAllLead[3]->Add(hBkgAllLead[p]);
            hBkgAllSub[3]->Add(hBkgAllSub[p]);
            
            hSigFailXGBLead[3]->Add(hSigFailXGBLead[p]);
            hSigFailXGBSub[3]->Add(hSigFailXGBSub[p]);
            hBkgFailXGBLead[3]->Add(hBkgFailXGBLead[p]);
            hBkgFailXGBSub[3]->Add(hBkgFailXGBSub[p]);
                
            hSigPassTightXGBLead[3]->Add(hSigPassTightXGBLead[p]);
            hSigPassTightXGBSub[3]->Add(hSigPassTightXGBSub[p]);
            hBkgPassTightXGBLead[3]->Add(hBkgPassTightXGBLead[p]);
            hBkgPassTightXGBSub[3]->Add(hBkgPassTightXGBSub[p]);
            
            hSigPassTightXGBTightMassLead[3]->Add(hSigPassTightXGBTightMassLead[p]);
            hSigPassTightXGBTightMassSub[3]->Add(hSigPassTightXGBTightMassSub[p]);
            hBkgPassTightXGBTightMassLead[3]->Add(hBkgPassTightXGBTightMassLead[p]);
            hBkgPassTightXGBTightMassSub[3]->Add(hBkgPassTightXGBTightMassSub[p]);

            hSigPassStdLead[3]->Add(hSigPassStdLead[p]);
            hSigPassStdSub[3]->Add(hSigPassStdSub[p]);
            hBkgPassStdLead[3]->Add(hBkgPassStdLead[p]);
            hBkgPassStdSub[3]->Add(hBkgPassStdSub[p]);
        }
        
        string varNameTmp2;
        if (varNames[i] == "mass")varNameTmp2 = "M_{#gamma#gamma}";
        else if (varNames[i] == "r9HLT")varNameTmp2 = "R_{9}";
        else if (varNames[i] == "sigmaIEtaIEta")varNameTmp2 = "#sigma_{i#eta i#eta}";
        else if (varNames[i] == "etaWidth")varNameTmp2 = "#eta Width";
        else if (varNames[i] == "phiWidth")varNameTmp2 = "#phi Width";
        else if (varNames[i] == "s4")varNameTmp2 = "S_{4}";
        else if (varNames[i] == "eta")varNameTmp2 = "#eta_{SC}";
        else if (varNames[i] == "et")varNameTmp2 = "E_{T}";
        else if (varNames[i] == "dR")varNameTmp2 = "#Delta r";
        else if (varNames[i] == "etOvrM")varNameTmp2 = "E_{T}/M";
        else varNameTmp2 = varNames[i];

        //Now loop over all filled histograms, plotting and saving each seperately
        for (int e = 0; e < nEta; e++){
            string outName = outNameGen + varNames[i] + etaFLabels[e];

            can->Clear();
            
            can->Divide(2,1);
            
            string plotTitle1 = genTitleStringSignal + etaLabels[e] + ";" + varNameTmp2;
            THStack *hStack1 = new THStack("hStack1",plotTitle1.c_str());
            
            string plotTitle2 = genTitleStringBkg + etaLabels[e] + ";" + varNameTmp2;
            THStack *hStack2 = new THStack("hStack2",plotTitle2.c_str());
            
            TLegend *legend1;
            if (varNames[i] == "r9HLT"|| varNames[i] == "s4") legend1 = new TLegend(0.1,0.65,0.5,0.9,"","brNDC");
            else legend1 = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
            
            TLegend *legend2;
            if (varNames[i] == "r9HLT"|| varNames[i] == "s4") legend2 = new TLegend(0.1,0.65,0.5,0.9,"","brNDC");
            else legend2 = new TLegend(0.50,0.65,0.90,0.9,"","brNDC");
            
            string label1A = "All";

            string label1PLoose = "#splitline{Failing XGBScore Cuts or M < 50";
            string label1PTight = "#splitline{Passing XGBScore Cuts (M > 50)";
            string label1PTightTightMass = "#splitline{Passing XGBScore Cuts (M > 95) ";
            string label1PStd = "#splitline{Passing Std. HLT (M > 95) ";
            
            string label2A = "All";
            string label2PLoose = "#splitline{Failing XGBScore Cuts or M < 50";
            string label2PTight = "#splitline{Passing XGBScore Cuts (M > 50)";
            string label2PTightTightMass = "#splitline{Passing XGBScore Cuts (M > 95) ";
            string label2PStd = "#splitline{Passing Std. HLT (M > 95) ";
            
            if (plotType == "Added"){
                
                TH1F *hSigAll = (TH1F*)hSigAllLead[e]->Clone();
                TH1F *hBkgAll = (TH1F*)hBkgAllLead[e]->Clone();
                
                TH1F *hSigFailXGB = (TH1F*)hSigFailXGBLead[e]->Clone();
                TH1F *hBkgFailXGB = (TH1F*)hBkgFailXGBLead[e]->Clone();
                
                TH1F *hSigPassTightXGB = (TH1F*)hSigPassTightXGBLead[e]->Clone();
                TH1F *hBkgPassTightXGB = (TH1F*)hBkgPassTightXGBLead[e]->Clone();
                
                TH1F *hSigPassTightXGBTightMass = (TH1F*)hSigPassTightXGBTightMassLead[e]->Clone();
                TH1F *hBkgPassTightXGBTightMass = (TH1F*)hBkgPassTightXGBTightMassLead[e]->Clone();
                
                TH1F *hSigPassStd = (TH1F*)hSigPassStdLead[e]->Clone();
                TH1F *hBkgPassStd = (TH1F*)hBkgPassStdLead[e]->Clone();
              
                if (isSinglePhoton == 0){
                    hSigAll->Add(hSigAllSub[e]);
                    hBkgAll->Add(hBkgAllSub[e]);
                    hSigFailXGB->Add(hSigFailXGBSub[e]);
                    hBkgFailXGB->Add(hBkgFailXGBSub[e]);
                    hSigPassTightXGB->Add(hSigPassTightXGBSub[e]);
                    hBkgPassTightXGB->Add(hBkgPassTightXGBSub[e]);
                    hSigPassTightXGBTightMass->Add(hSigPassTightXGBTightMassSub[e]);
                    hBkgPassTightXGBTightMass->Add(hBkgPassTightXGBTightMassSub[e]);
                    hSigPassStd->Add(hSigPassStdSub[e]);
                    hBkgPassStd->Add(hBkgPassStdSub[e]);
                }
                
                char buff1[100], buff2[100], buff3[100], buff4[100], buff5[100], buff6[100], buff7[100], buff8[100], buff9[100], buff10[100], buff11[100], buff12[100];
                
                snprintf(buff1, sizeof(buff1), "(%0.1f)", hSigAll->GetSumOfWeights());
                string nEvents1 = buff1;
                label1A += nEvents1;
                
                snprintf(buff2, sizeof(buff2), "}{(%0.1f)(%0.2f)", hSigFailXGB->GetSumOfWeights(),100*(hSigFailXGB->GetSumOfWeights()/hSigAll->GetSumOfWeights()));
                string nEvents1PLoose = buff2;
                label1PLoose += nEvents1PLoose;
                
                snprintf(buff3, sizeof(buff3), "}{(%0.1f)(%0.2f)", hSigPassTightXGB->GetSumOfWeights(),100*(hSigPassTightXGB->GetSumOfWeights()/hSigAll->GetSumOfWeights()));
                string nEvents1PTight = buff3;
                label1PTight += nEvents1PTight;
                
                snprintf(buff9, sizeof(buff9), "}{(%0.1f)(%0.2f)", hSigPassTightXGBTightMass->GetSumOfWeights(),100*(hSigPassTightXGBTightMass->GetSumOfWeights()/hSigAll->GetSumOfWeights()));
                string nEvents1PTightTightMass = buff9;
                label1PTightTightMass += nEvents1PTightTightMass;
                
                snprintf(buff4, sizeof(buff4), "}{(%0.1f)(%0.2f)", hSigPassStd->GetSumOfWeights(),100*(hSigPassStd->GetSumOfWeights()/hSigAll->GetSumOfWeights()));
                string nEvents1PStd = buff4;
                label1PStd += nEvents1PStd;
                
                //Ripped from dnaPlotSideBySide
                if(varNames[i] == "mass"){
                    char buff13[100], buff14[100], buff15[100], buff16[100], buff17[100];
                    
                    //snprintf(buff11, sizeof(buff11), "(Width = %0.2f)",getExactWidth(hSigAll[e]));
                    snprintf(buff13, sizeof(buff13), "(2#sigma_{eff} = %0.2f)",newGetWidth(hSigAll));
                    string widthAll = buff13;
                    label1A += widthAll;
                    
                    //snprintf(buff12, sizeof(buff12), "(Width = %0.2f)",getExactWidth(hSigPresel[e]));
                    //snprintf(buff14, sizeof(buff14), "(2#sigma_{eff} = %0.2f)",newGetWidth(hSigFailXGB));
                    //string widthFail = buff14;
                    //label1PLoose += widthFail;
                    
                    //snprintf(buff13, sizeof(buff13), "(Width = %0.2f)",getExactWidth(hSigIDMVA[e]));
                    snprintf(buff15, sizeof(buff15), "(2#sigma_{eff} = %0.2f)",newGetWidth(hSigPassTightXGB));
                    string widthIDMVATight = buff15;
                    label1PTight += widthIDMVATight;
                    
                    //snprintf(buff14, sizeof(buff14), "(Width = %0.2f)",getExactWidth(hSigIDMVA2[e]));
                    snprintf(buff16, sizeof(buff16), "(2#sigma_{eff} = %0.2f)",newGetWidth(hSigPassTightXGBTightMass));
                    string widthIDMVA2 = buff16;
                    label1PTightTightMass += widthIDMVA2;
                    
                    //snprintf(buff15, sizeof(buff15), "(Width = %0.2f)",getExactWidth(hSigIDMVA3[e]));
                    snprintf(buff17, sizeof(buff17), "(2#sigma_{eff} = %0.2f)",newGetWidth(hSigPassStd));
                    string widthPresel = buff17;
                    label1PStd += widthPresel;
                }
                
                if(e == 3 && varNames[i] == "mass"){
                    hSigAll->Rebin(100);
                    hSigFailXGB->Rebin(100);
                    hSigPassTightXGB->Rebin(100);
                    hSigPassTightXGBTightMass->Rebin(100);
                    hSigPassStd->Rebin(100);

                }
                
                label1PLoose += "}";
                label1PTight += "}";
                label1PTightTightMass += "}";
                label1PStd += "}";
                
                
                TH1F *hSigAllDraw = DrawOverflow(hSigAll);
                TH1F *hStackHistoSig = (TH1F*)hSigAllDraw->Clone();
                hStackHistoSig->Reset();
                hStackHistoSig->GetXaxis()->SetRange(0,hStackHistoSig->GetNbinsX());
                hStack1->SetHistogram(hStackHistoSig);
                
                hStack1->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
                hStack1->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
                hStack1->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
                hStack1->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
                
                hStack1->Add(hSigAllDraw);
                legend1->AddEntry(hSigAllDraw,label1A.c_str(),"pl");
                
                TH1F *hSigFailDraw = DrawOverflow(hSigFailXGB);
                hSigFailDraw->SetLineColor(2);
                hStack1->Add(hSigFailDraw);
                legend1->AddEntry(hSigFailDraw,label1PLoose.c_str(),"pl");
                
                TH1F *hSigPassStdDraw = DrawOverflow(hSigPassStd);
                hSigPassStdDraw->SetLineColor(3);
                hStack1->Add(hSigPassStdDraw);
                legend1->AddEntry(hSigPassStdDraw,label1PStd.c_str(),"pl");

                
                TH1F *hSigPassTightTightMassDraw = DrawOverflow(hSigPassTightXGBTightMass);
                hSigPassTightTightMassDraw->SetLineColor(4);
                hStack1->Add(hSigPassTightTightMassDraw);
                legend1->AddEntry(hSigPassTightTightMassDraw,label1PTightTightMass.c_str(),"pl");
                
                TH1F *hSigPassTightDraw = DrawOverflow(hSigPassTightXGB);
                hSigPassTightDraw->SetLineColor(4);
                hSigPassTightDraw->SetLineStyle(3);
                hStack1->Add(hSigPassTightDraw);
                legend1->AddEntry(hSigPassTightDraw,label1PTight.c_str(),"pl");
                
                can->cd(1);
                gPad->SetGrid();
                if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho" || varNames[i] == "hOvrE")gPad->SetLogy();
                else gPad->SetLogy(0);
                
                hStack1->Draw("nostackhist");
                legend1->Draw("same");
                
                snprintf(buff5, sizeof(buff5), "(%0.1f)", hBkgAll->GetSumOfWeights());
                string nEvents2 = buff5;
                label2A += nEvents2;
                
                snprintf(buff6, sizeof(buff6), "}{(%0.1f)(%0.2f)}", hBkgFailXGB->GetSumOfWeights(),100*(hBkgFailXGB->GetSumOfWeights()/hBkgAll->GetSumOfWeights()));
                string nEvents2PLoose = buff6;
                label2PLoose += nEvents2PLoose;
                
                snprintf(buff7, sizeof(buff7), "}{(%0.1f)(%0.2f)}", hBkgPassTightXGB->GetSumOfWeights(),100*(hBkgPassTightXGB->GetSumOfWeights()/hBkgAll->GetSumOfWeights()));
                string nEvents2PTight = buff7;
                label2PTight += nEvents2PTight;
                
                snprintf(buff10, sizeof(buff10), "}{(%0.1f)(%0.2f)}", hBkgPassTightXGBTightMass->GetSumOfWeights(),100*(hBkgPassTightXGBTightMass->GetSumOfWeights()/hBkgAll->GetSumOfWeights()));
                string nEvents2PTightTightMass = buff10;
                label2PTightTightMass += nEvents2PTightTightMass;
                
                snprintf(buff8, sizeof(buff8), "}{(%0.1f)(%0.2f)}", hBkgPassStd->GetSumOfWeights(),100*(hBkgPassStd->GetSumOfWeights()/hBkgAll->GetSumOfWeights()));
                string nEvents2PStd = buff8;
                label2PStd += nEvents2PStd;

                TH1F *hBkgAllDraw = DrawOverflow(hBkgAll);
                TH1F *hStackHistoBkg = (TH1F*)hBkgAllDraw->Clone();
                hStackHistoBkg->Reset();
                hStackHistoBkg->GetXaxis()->SetRange(0,hStackHistoBkg->GetNbinsX());
                hStack2->SetHistogram(hStackHistoBkg);
                
                hStack2->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
                hStack2->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
                hStack2->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
                hStack2->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
                
                hStack2->Add(hBkgAllDraw);
                legend2->AddEntry(hBkgAllDraw,label2A.c_str(),"pl");
                
                TH1F *hBkgFailXGBDraw = DrawOverflow(hBkgFailXGB);
                hBkgFailXGBDraw->SetLineColor(2);
                legend2->AddEntry(hBkgFailXGBDraw,label2PLoose.c_str(),"pl");
                hStack2->Add(hBkgFailXGBDraw);
                
                TH1F *hBkgPassStdDraw = DrawOverflow(hBkgPassStd);
                hBkgPassStdDraw->SetLineColor(3);
                hStack2->Add(hBkgPassStdDraw);
                legend2->AddEntry(hBkgPassStdDraw,label2PStd.c_str(),"pl");
                
                TH1F *hBkgPassTightXGBTightMassDraw = DrawOverflow(hBkgPassTightXGBTightMass);
                hBkgPassTightXGBTightMassDraw->SetLineColor(4);
                hStack2->Add(hBkgPassTightXGBTightMassDraw);
                legend2->AddEntry(hBkgPassTightXGBTightMassDraw,label2PTightTightMass.c_str(),"pl");
                
                TH1F *hBkgPassTightXGBDraw = DrawOverflow(hBkgPassTightXGB);
                hBkgPassTightXGBDraw->SetLineColor(4);
                hBkgPassTightXGBDraw->SetLineStyle(3);
                hStack2->Add(hBkgPassTightXGBDraw);
                legend2->AddEntry(hBkgPassTightXGBDraw,label2PTight.c_str(),"pl");
                
                can->cd(2);
                gPad->SetGrid();
                //if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho")gPad->SetLogy();
                //else gPad->SetLogy(0);
                gPad->SetLogy();
                gPad->SetGrid();
                
                hStack2->Draw("nostackhist");
                legend2->Draw("same");
                
                
                can->SaveAs((outName+".png").c_str());
                can->SaveAs((outName+".root").c_str());
                
                can->Clear();
                
                f->Close();
                f->Delete();
            }//For "Added" Plotting
            
            
            if (plotType == "Separate"){
                
                char buff1[100], buff2[100], buff3[100], buff4[100], buff5[100], buff6[100], buff7[100], buff8[100], buff9[100], buff10[100], buff11[100], buff12[100];
                
                snprintf(buff1, sizeof(buff1), "(%0.1f)", hSigAllLead[e]->GetSumOfWeights());
                string nEvents1 = buff1;
                label1A += nEvents1;
                
                snprintf(buff2, sizeof(buff2), "}{(%0.1f)(%0.2f)}", hSigFailXGBLead[e]->GetSumOfWeights(),100*(hSigFailXGBLead[e]->GetSumOfWeights()/hSigAllLead[e]->GetSumOfWeights()));
                string nEvents1PLoose = buff2;
                label1PLoose += nEvents1PLoose;
                
                snprintf(buff3, sizeof(buff3), "}{(%0.1f)(%0.2f)}", hSigPassTightXGBLead[e]->GetSumOfWeights(),100*(hSigPassTightXGBLead[e]->GetSumOfWeights()/hSigAllLead[e]->GetSumOfWeights()));
                string nEvents1PTight = buff3;
                label1PTight += nEvents1PTight;
                
                snprintf(buff4, sizeof(buff4), "}{(%0.1f)(%0.2f)}", hSigPassTightXGBTightMassLead[e]->GetSumOfWeights(),100*(hSigPassTightXGBTightMassLead[e]->GetSumOfWeights()/hSigAllLead[e]->GetSumOfWeights()));
                string nEvents1PTightTightMass = buff4;
                label1PTightTightMass += nEvents1PTightTightMass;
                
                snprintf(buff5, sizeof(buff5), "}{(%0.1f)(%0.2f)}", hSigPassStdLead[e]->GetSumOfWeights(),100*(hSigPassStdLead[e]->GetSumOfWeights()/hSigAllLead[e]->GetSumOfWeights()));
                string nEvents1PStd = buff5;
                label1PStd += nEvents1PStd;
                
                TH1F *hSigAllLeadDraw = DrawOverflow(hSigAllLead[e]);
                TH1F *hStackHistoSig = (TH1F*)hSigAllLeadDraw->Clone();
                hStackHistoSig->Reset();
                hStackHistoSig->GetXaxis()->SetRange(0,hStackHistoSig->GetNbinsX());
                hStack1->SetHistogram(hStackHistoSig);
                
                hStack1->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
                hStack1->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
                hStack1->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
                hStack1->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
                
                hStack1->Add(hSigAllLeadDraw);
                legend1->AddEntry(hSigAllLeadDraw,label1A.c_str(),"pl");
                TH1F *hSigAllSubDraw = DrawOverflow(hSigAllSub[e]);
                hSigAllSubDraw->SetLineStyle(3);
                hStack1->Add(hSigAllSubDraw);
                
                TH1F *hSigFailLeadDraw = DrawOverflow(hSigFailXGBLead[e]);
                hSigFailLeadDraw->SetLineColor(2);
                hStack1->Add(hSigFailLeadDraw);
                legend1->AddEntry(hSigFailLeadDraw,label1PLoose.c_str(),"pl");
                TH1F *hSigFailSubDraw = DrawOverflow(hSigFailXGBSub[e]);
                hSigFailSubDraw->SetLineColor(2);
                hSigFailSubDraw->SetLineStyle(3);
                hStack1->Add(hSigFailSubDraw);
                
                TH1F *hSigPassStdLeadDraw = DrawOverflow(hSigPassStdLead[e]);
                hSigPassStdLeadDraw->SetLineColor(3);
                hStack1->Add(hSigPassStdLeadDraw);
                legend1->AddEntry(hSigPassStdLeadDraw,label1PStd.c_str(),"pl");
                TH1F *hSigPassStdSubDraw = DrawOverflow(hSigPassStdSub[e]);
                hSigPassStdSubDraw->SetLineColor(3);
                hSigPassStdSubDraw->SetLineStyle(3);
                hStack1->Add(hSigPassStdSubDraw);

                TH1F *hSigPassTightXGBTightMassLeadDraw = DrawOverflow(hSigPassTightXGBTightMassLead[e]);
                hSigPassTightXGBTightMassLeadDraw->SetLineColor(4);
                hStack1->Add(hSigPassTightXGBTightMassLeadDraw);
                legend1->AddEntry(hSigPassTightXGBTightMassLeadDraw,label1PTightTightMass.c_str(),"pl");
                TH1F *hSigPassTightXGBTightMassSubDraw = DrawOverflow(hSigPassTightXGBTightMassSub[e]);
                hSigPassTightXGBTightMassSubDraw->SetLineColor(4);
                hSigPassTightXGBTightMassSubDraw->SetLineStyle(3);
                hStack1->Add(hSigPassTightXGBTightMassSubDraw);
                
                
                can->cd(1);
                gPad->SetGrid();
                if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho")gPad->SetLogy();
                else gPad->SetLogy(0);
                
                hStack1->Draw("nostackhist");
                legend1->Draw("same");
                
                snprintf(buff7, sizeof(buff7), "(%0.1f)", hBkgAllLead[e]->GetSumOfWeights());
                string nEvents2 = buff7;
                label2A += nEvents2;
                
                snprintf(buff8, sizeof(buff8), "}{(%0.1f)(%0.2f)}", hBkgFailXGBLead[e]->GetSumOfWeights(),100*(hBkgFailXGBLead[e]->GetSumOfWeights()/hSigAllLead[e]->GetSumOfWeights()));
                string nEvents2PLoose = buff8;
                label2PLoose += nEvents2PLoose;
                
                snprintf(buff9, sizeof(buff9), "}{(%0.1f)(%0.2f)}", hBkgPassTightXGBLead[e]->GetSumOfWeights(),100*(hBkgPassTightXGBLead[e]->GetSumOfWeights()/hSigAllLead[e]->GetSumOfWeights()));
                string nEvents2PTight = buff9;
                label2PTight += nEvents2PTight;
                
                snprintf(buff10, sizeof(buff10), "}{(%0.1f)(%0.2f)}", hBkgPassTightXGBTightMassLead[e]->GetSumOfWeights(),100*(hBkgPassTightXGBTightMassLead[e]->GetSumOfWeights()/hSigAllLead[e]->GetSumOfWeights()));
                string nEvents2PTightTightMass = buff10;
                label2PTightTightMass += nEvents2PTightTightMass;
                
                snprintf(buff11, sizeof(buff11), "}{(%0.1f)(%0.2f)}", hBkgPassStdLead[e]->GetSumOfWeights(),100*(hBkgPassStdLead[e]->GetSumOfWeights()/hSigAllLead[e]->GetSumOfWeights()));
                string nEvents2PStd = buff11;
                label2PStd += nEvents2PStd;
                
                TH1F *hBkgAllLeadDraw = DrawOverflow(hBkgAllLead[e]);
                TH1F *hStackHistoBkg = (TH1F*)hBkgAllLeadDraw->Clone();
                hStackHistoBkg->Reset();
                hStackHistoBkg->GetXaxis()->SetRange(0,hStackHistoBkg->GetNbinsX());
                hStack2->SetHistogram(hStackHistoBkg);
                
                hStack2->GetHistogram()->GetXaxis()->SetTitle(varNames[i].c_str());
                hStack2->GetHistogram()->GetXaxis()->SetLabelSize(0.025);
                hStack2->GetHistogram()->GetXaxis()->SetTitleSize(0.03);
                hStack2->GetHistogram()->GetYaxis()->SetLabelSize(0.025);
                
                hStack2->Add(hBkgAllLeadDraw);
                legend2->AddEntry(hBkgAllLeadDraw,label2A.c_str(),"pl");
                TH1F *hBkgAllSubDraw = DrawOverflow(hBkgAllSub[e]);
                hBkgAllSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgAllSubDraw);
                
                TH1F *hBkgFailLeadDraw = DrawOverflow(hBkgFailXGBLead[e]);
                hBkgFailLeadDraw->SetLineColor(2);
                hStack2->Add(hBkgFailLeadDraw);
                legend2->AddEntry(hBkgFailLeadDraw,label2PLoose.c_str(),"pl");
                TH1F *hBkgFailSubDraw = DrawOverflow(hBkgFailXGBSub[e]);
                hBkgFailSubDraw->SetLineColor(2);
                hBkgFailSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgFailSubDraw);
                
                TH1F *hBkgPassStdLeadDraw = DrawOverflow(hBkgPassStdLead[e]);
                hBkgPassStdLeadDraw->SetLineColor(3);
                hStack2->Add(hBkgPassStdLeadDraw);
                legend2->AddEntry(hBkgPassStdLeadDraw,label2PStd.c_str(),"pl");
                TH1F *hBkgPassStdSubDraw = DrawOverflow(hBkgPassStdSub[e]);
                hBkgPassStdSubDraw->SetLineColor(3);
                hBkgPassStdSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgPassStdSubDraw);

                TH1F *hBkgPassTightXGBTightMassLeadDraw = DrawOverflow(hBkgPassTightXGBTightMassLead[e]);
                hBkgPassTightXGBTightMassLeadDraw->SetLineColor(4);
                hStack2->Add(hBkgPassTightXGBTightMassLeadDraw);
                legend2->AddEntry(hBkgPassTightXGBTightMassLeadDraw,label2PTightTightMass.c_str(),"pl");
                TH1F *hBkgPassTightXGBTightMassSubDraw = DrawOverflow(hBkgPassTightXGBTightMassSub[e]);
                hBkgPassTightXGBTightMassSubDraw->SetLineColor(4);
                hBkgPassTightXGBTightMassSubDraw->SetLineStyle(3);
                hStack2->Add(hBkgPassTightXGBTightMassSubDraw);
                
                can->cd(2);
                gPad->SetGrid();
                //if (varNames[i] == "ecalPFIso" || varNames[i] == "trkIsoPho")gPad->SetLogy();
                //else gPad->SetLogy(0);
                gPad->SetGrid();
                gPad->SetLogy();
                hStack2->Draw("nostackhist");
                legend2->Draw("same");
                
                
                can->SaveAs((outName+".png").c_str());
                can->SaveAs((outName+".root").c_str());
                
                can->Clear();
                
                f->Close();
                f->Delete();
            }//For "Separate" Plotting
        }
        for (int p = 0; p < 4; p++){
            hSigAllLead[p]->Delete();
            hSigAllSub[p]->Delete();
            hSigFailXGBLead[p]->Delete();
            hSigFailXGBSub[p]->Delete();
            hSigPassTightXGBLead[p]->Delete();
            hSigPassTightXGBSub[p]->Delete();
            hSigPassTightXGBTightMassLead[p]->Delete();
            hSigPassTightXGBTightMassSub[p]->Delete();
            hSigPassStdLead[p]->Delete();
            hSigPassStdSub[p]->Delete();
            
            hBkgAllLead[p]->Delete();
            hBkgAllSub[p]->Delete();
            hBkgFailXGBLead[p]->Delete();
            hBkgFailXGBSub[p]->Delete();
            hBkgPassTightXGBLead[p]->Delete();
            hBkgPassTightXGBSub[p]->Delete();
            hBkgPassTightXGBTightMassLead[p]->Delete();
            hBkgPassTightXGBTightMassSub[p]->Delete();
            hBkgPassStdLead[p]->Delete();
            hBkgPassStdSub[p]->Delete();
        }
    }// variable loop (i)
}//Function definition

double newGetWidth(TH1F *histo){
    double desiredProp = 0.68;
    double totalInt = histo->Integral();
    int peakBin = histo->GetMaximumBin();
    cout<<"peakBin = "<<peakBin<<" with value of "<<histo->GetBinCenter(peakBin)<<endl;
    int nBinsEachWay;
    double finalProp;
    for (int i = 0; i < 5000; i++){
        double tmpProp = histo->Integral(peakBin-i,peakBin+i)/totalInt;
        if (tmpProp >= desiredProp){
            nBinsEachWay = i;
            finalProp = tmpProp;
            break;
        }
    }
    double binSize = histo->GetBinWidth(peakBin);
    double widthInGeV = (nBinsEachWay*2 + 1)*binSize;
    cout<<"Final Width = "<<widthInGeV<<" With proportion of events = "<<finalProp<<endl;
    return widthInGeV;
}
