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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
string ToStringWithPrecision(double inVal, const int n = 2);
TH1F *DrawOverflow(TH1F* h);
double calcDR(float eta1, float eta2, float phi1, float phi2);
double calcMassXGB(TTreeReaderArray<float>* xgbScore,double xgbScoreB, double xgbScoreE,TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, int nEgsIn, int* pho1, int* pho2);
double calcMassHLT(TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, TTreeReaderArray<float>* r9, TTreeReaderArray<float>* hOvrE, TTreeReaderArray<float>* sIeIe, TTreeReaderArray<float>* phoIso, TTreeReaderArray<float>* ecalIso, int nEgsIn, int* pho1, int* pho2);

TH1F *DrawOverflow(TH1F* h){
    int nBins    = h->GetNbinsX()+1;
    Double_t *xbins= new Double_t[nBins+1];
    for (UInt_t i=0;i<nBins;i++)
        xbins[i]=h->GetBinLowEdge(i+1);
    xbins[nBins]=xbins[nBins-1]+h->GetBinWidth(nBins);
    //book a temporary histogram having extra bins for overflows
    TH1F *htmp = new TH1F(h->GetName(), h->GetTitle(), nBins, xbins);
    htmp->Sumw2();
    //fill the new histogram including the overflows
    for (UInt_t i=1; i<=nBins; i++) {
        htmp->SetBinContent(htmp->FindBin(htmp->GetBinCenter(i)),h->GetBinContent(i));
    }
    //htmp->SetBinContent(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinContent(0));
    htmp->SetBinContent(0, h->GetBinContent(0));
    htmp->GetXaxis()->SetRange(0,nBins);
    return htmp;
}

string ToStringWithPrecision(double inVal, const int n = 2){
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << inVal;
    return std::move(out).str();
}
double calcMassXGB(TTreeReaderArray<float>* xgbScore,double xgbCutB, double xgbCutE, TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, int nEgsIn, int* pho1, int* pho2){
    double mass = -999.0;
    int iOut = -1;
    int jOut = -1;
    double scoreHighOut = -1;
    double scoreLowOut = -1;
    int bothPass = 0;
    int someSelectionMade = 0;
    if(nEgsIn > 1){
        for(int i = 0; i < nEgsIn;i++){
            //if ((*et)[i] > 14.25 && fabs((*eta)[i]) < 2.55){
            if (fabs((*eta)[i]) < 2.55 && (*et)[i] > 14.25 ){
                ROOT::Math::PtEtaPhiMVector v1((*et)[i],(*eta)[i],(*phi)[i],0.0);
                for(int j = 0; j < nEgsIn; j ++){
                    if (j != i){
                        //if(fabs((*eta)[j]) < 2.55 && (*et)[j] > 14.25){
                        if (fabs((*eta)[j]) < 2.55 && (*et)[j] > 14.25){
                            ROOT::Math::PtEtaPhiMVector v2((*et)[j],(*eta)[j],(*phi)[j],0.0);
                            auto v3 = v1+v2;
                            double massTmp = v3.M();
                            double dR = calcDR((*eta)[i],(*eta)[j],(*phi)[i],(*phi)[j]);
                            //double dRTmp = calcDR((*eta)[i], (*eta)[j], (*phi)[i], (*phi)[j])
                            if(massTmp > 95 && dR > 0.5){
                                double scoreHighTmp,scoreLowTmp,indexHighTmp,indexLowTmp;
                                if ((*xgbScore)[i] > (*xgbScore)[j]){
                                    scoreHighTmp = (*xgbScore)[i];
                                    indexHighTmp = i;
                                    scoreLowTmp = (*xgbScore)[j];
                                    indexLowTmp = j;
                                }
                                else if ((*xgbScore)[i] < (*xgbScore)[j]){
                                    scoreHighTmp = (*xgbScore)[j];
                                    indexHighTmp = j;
                                    scoreLowTmp = (*xgbScore)[i];
                                    indexLowTmp = i;
                                }
                                
                                if (scoreHighTmp >= scoreHighOut && scoreLowTmp >= scoreLowOut){
                                    bothPass = 1;
                                    scoreHighOut = scoreHighTmp;
                                    scoreLowOut = scoreLowTmp;
                                    iOut = indexHighTmp;
                                    jOut = indexLowTmp;
                                    mass = massTmp;
                                    someSelectionMade = 1;
                                }//Save pair if they pass mass cuts and have highest scores
                            }//Passing highest mass and dR cut
                        }//If et and eta good [j]
                    }//j != i
                }//j Loop
            }//If et and eta good [i]
        }//i Loop for passing both
        if(bothPass != 1){
            for(int i = 0; i < nEgsIn;i++){
                //if ((*et)[i] > 14.25 && fabs((*eta)[i]) < 2.55){
                if (fabs((*eta)[i]) < 2.55 && (*et)[i] > 14.25 ){
                    ROOT::Math::PtEtaPhiMVector v1((*et)[i],(*eta)[i],(*phi)[i],0.0);
                    for(int j = i+1; j < nEgsIn; j ++){
                        //if(fabs((*eta)[j]) < 2.55 && (*et)[j] > 14.25){
                        if (fabs((*eta)[j]) < 2.55 && (*et)[j] > 14.25){
                            ROOT::Math::PtEtaPhiMVector v2((*et)[j],(*eta)[j],(*phi)[j],0.0);
                            auto v3 = v1+v2;
                            double massTmp = v3.M();
                            double dR = calcDR((*eta)[i],(*eta)[j],(*phi)[i],(*phi)[j]);
                            //double dRTmp = calcDR((*eta)[i], (*eta)[j], (*phi)[i], (*phi)[j])
                            if(massTmp > 50 && dR > 0.5){
                                double scoreHighTmp,scoreLowTmp,indexHighTmp,indexLowTmp;
                                if ((*xgbScore)[i] > (*xgbScore)[j]){
                                    scoreHighTmp = (*xgbScore)[i];
                                    indexHighTmp = i;
                                    scoreLowTmp = (*xgbScore)[j];
                                    indexLowTmp = j;
                                }
                                else if ((*xgbScore)[i] < (*xgbScore)[j]){
                                    scoreHighTmp = (*xgbScore)[j];
                                    indexHighTmp = j;
                                    scoreLowTmp = (*xgbScore)[i];
                                    indexLowTmp = i;
                                }
                                
                                if (scoreHighTmp >= scoreHighOut && scoreLowTmp >= scoreLowOut){
                                    bothPass = 1;
                                    scoreHighOut = scoreHighTmp;
                                    scoreLowOut = scoreLowTmp;
                                    iOut = indexHighTmp;
                                    jOut = indexLowTmp;
                                    mass = massTmp;
                                    someSelectionMade = 1;
                                }//Save pair if they pass mass cuts and have highest scores
                            }//Passing Medium mass and dR cut
                        }//If et and eta good [j]
                    }//j Loop
                }//If et and eta good [i]
            }//i Loop for passing both
        }//If first iteration fails, try again with lowered mass cut
        if(someSelectionMade != 1){
            double scoreLead = -1;
            double scoreSub = -1;
            for(int i = 0; i < nEgsIn;i++){
                if(fabs((*eta)[i]) < 2.55 && (*et)[i] > 14.25){
                    double scoreLeadTmp = (*xgbScore)[i];
                    if(scoreLeadTmp > scoreLead){
                        scoreLead = scoreLeadTmp;
                        for(int j = i+1; j < nEgsIn; j ++){
                            if(fabs((*eta)[j]) < 2.55 && (*et)[j] > 14.25){
                                double scoreSubTmp = (*xgbScore)[j];
                                if(scoreSubTmp > scoreSub){
                                    scoreSub = scoreSubTmp;
                                    iOut = i;
                                    jOut = j;
                                    someSelectionMade = 1;
                                }//If score of j is larger than highest sub score
                            }//Et and Eta cuts on second photon
                        }//j Loop
                    }//If score of i is larger than highest lead score
                }//Et and Eta cuts on first photon
            }//i Loop for not passing both
            ROOT::Math::PtEtaPhiMVector v1((*et)[iOut],(*eta)[iOut],(*phi)[iOut],0.0);
            ROOT::Math::PtEtaPhiMVector v2((*et)[jOut],(*eta)[jOut],(*phi)[jOut],0.0);
            auto v3 = v1+v2;
            mass = v3.M();
        }//If bothPass != 1
    }//nEgs
    if(someSelectionMade != 1 || nEgsIn < 2){
        iOut = -1;
        jOut = -1;
        mass = -999;
    }
    *pho1 = iOut;
    *pho2 = jOut;
    return mass;
}

//double calcMassXGB(TTreeReaderArray<float>* xgbScore,double xgbCutB, double xgbCutE, TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, int nEgsIn, int* pho1, int* pho2){
//    double mass = -999.0;
//    int iOut = -1;
//    int jOut = -1;
//    int bothPass = 0;
//    int someSelectionMade = 0;
//    if(nEgsIn > 1){
//        for(int i = 0; i < nEgsIn;i++){
//            //if ((*et)[i] > 14.25 && fabs((*eta)[i]) < 2.55){
//            if (fabs((*eta)[i]) < 2.55 && (*et)[i] > 14.25 && ((abs((*eta)[i]) < 1.5 && (*xgbScore)[i] > xgbCutB) ||(abs((*eta)[i]) > 1.5 && (*xgbScore)[i] > xgbCutE))){
//                ROOT::Math::PtEtaPhiMVector v1((*et)[i],(*eta)[i],(*phi)[i],0.0);
//                for(int j = i+1; j < nEgsIn; j ++){
//                    //if(fabs((*eta)[j]) < 2.55 && (*et)[j] > 14.25){
//                    if (fabs((*eta)[j]) < 2.55 && (*et)[j] > 14.25 && ((abs((*eta)[j]) < 1.5 && (*xgbScore)[j] > xgbCutB) ||(abs((*eta)[j]) > 1.5 && (*xgbScore)[j] > xgbCutE))){
//                        ROOT::Math::PtEtaPhiMVector v2((*et)[j],(*eta)[j],(*phi)[j],0.0);
//                        auto v3 = v1+v2;
//                        double massTmp = v3.M();
//                        double dR = calcDR((*eta)[i],(*eta)[j],(*phi)[i],(*phi)[j]);
//                        //double dRTmp = calcDR((*eta)[i], (*eta)[j], (*phi)[i], (*phi)[j])
//                        if(massTmp > mass){
//                            mass = massTmp;
//                            if (massTmp > 60 && dR > 0.5){
//                                bothPass = 1;
//                                iOut = i;
//                                jOut = j;
//                                someSelectionMade = 1;
//                            }//Mass cut for photon selection. If mass is less than 60, don't use this pair.
//                        }//massTmp > mass
//                    }//If et and eta good [j]
//                }//j Loop
//            }//If et and eta good [i]
//        }//i Loop for passing both
//        if(bothPass != 1){
//            double scoreLead = -1;
//            double scoreSub = -1;
//            for(int i = 0; i < nEgsIn;i++){
//                if(fabs((*eta)[i]) < 2.55 && (*et)[i] > 14.25){
//                    double scoreLeadTmp = (*xgbScore)[i];
//                    if(scoreLeadTmp > scoreLead){
//                        scoreLead = scoreLeadTmp;
//                        for(int j = i+1; j < nEgsIn; j ++){
//                            if(fabs((*eta)[j]) < 2.55 && (*et)[j] > 14.25){
//                                double scoreSubTmp = (*xgbScore)[j];
//                                if(scoreSubTmp > scoreSub){
//                                    scoreLead = scoreSubTmp;
//                                    double dR = calcDR((*eta)[i],(*eta)[j],(*phi)[i],(*phi)[j]);
//                                    if (dR > 0.5){
//                                        iOut = i;
//                                        jOut = j;
//                                        someSelectionMade = 1;
//                                    }
//                                }//If score oof j is larger than highest sub score
//                            }//Et and Eta cuts on second photon
//                        }//j Loop
//                    }//If score of i is larger than highest lead score
//                }//Et and Eta cuts on first photon
//            }//i Loop for not passing both
//            ROOT::Math::PtEtaPhiMVector v1((*et)[iOut],(*eta)[iOut],(*phi)[iOut],0.0);
//            ROOT::Math::PtEtaPhiMVector v2((*et)[jOut],(*eta)[jOut],(*phi)[jOut],0.0);
//            auto v3 = v1+v2;
//            mass = v3.M();
//        }//If bothPass != 1
//    }//nEgs
//    if(someSelectionMade != 1 || nEgsIn < 2){
//        iOut = -1;
//        jOut = -1;
//        mass = -999;
//    }
//    *pho1 = iOut;
//    *pho2 = jOut;
//    return mass;
//}
double calcMassHLT(TTreeReaderArray<float>* et, TTreeReaderArray<float>* eta, TTreeReaderArray<float>* phi, TTreeReaderArray<float>* r9, TTreeReaderArray<float>* hOvrE, TTreeReaderArray<float>* sIeIe, TTreeReaderArray<float>* phoIso, TTreeReaderArray<float>* ecalIso, int nEgsIn, int* pho1, int* pho2){
    double mass = -999.0;
    int iOut = -1;
    int jOut = -1;
    int bothPass = 0;
    int hltPassFail = 0;
    if(nEgsIn > 1){
        for(int i = 0; i < nEgsIn;i++){
            //if ((*et)[i] > 14.25 && fabs((*eta)[i]) < 2.55){
            if ((*et)[i] > 30.0 && fabs((*eta)[i]) < 2.55 && ((fabs((*eta)[i]) < 1.479 && (*r9)[i] > 0.5 && (*hOvrE)[i] < 0.12 && !((*r9)[i] < 0.85 && ((*sIeIe)[i] > 0.015 || (*phoIso)[i] > 6))) || (fabs((*eta)[i]) > 1.479 && (*r9)[i] > 0.8 && (*hOvrE)[i] < 0.10 && !((*r9)[i] < 0.90 && ((*sIeIe)[i] > 0.035 || (*phoIso)[i] > 6))))){
                ROOT::Math::PtEtaPhiMVector v1((*et)[i],(*eta)[i],(*phi)[i],0.0);
                for(int j = 0; j < nEgsIn; j ++){
                    if(j != i){
                        //if(fabs((*eta)[j]) < 2.5 && (*et)[j] > 14.25){
                        if ((*et)[j] > 22.0 && fabs((*eta)[j]) < 2.55 && ((fabs((*eta)[j]) < 1.479 && (*r9)[j] > 0.5 && (*hOvrE)[j] < 0.12 && !((*r9)[j] < 0.85 && ((*sIeIe)[j] > 0.015 || (*phoIso)[j] > 6))) || (fabs((*eta)[j]) > 1.479 && (*r9)[j] > 0.8 && (*hOvrE)[j] < 0.10 && !((*r9)[j] < 0.90 && ((*sIeIe)[j] > 0.035 || (*phoIso)[j] > 6))))){
                            ROOT::Math::PtEtaPhiMVector v2((*et)[j],(*eta)[j],(*phi)[j],0.0);
                            auto v3 = v1+v2;
                            double massTmp = v3.M();
                            //double dRTmp = calcDR((*eta)[i], (*eta)[j], (*phi)[i], (*phi)[j])
                            if(massTmp > mass){
                                //cout<<"(i,j) = ("<<i<<","<<j<<"), mass = "<<ma
                                mass = massTmp;
                                bothPass = 1;
                                iOut = i;
                                jOut = j;
                            }//massTmp > mass
                        }//If et and eta good [j]
                    }
                }//j Loop
            }//If et and eta good [i]
        }//i Loop for passing both
        if (bothPass == 1){
            for(int i = 0; i < nEgsIn; i ++){
                if ((*et)[i] > 22.0 && fabs((*eta)[i]) < 2.55 && ((fabs((*eta)[i]) < 1.5 && (*r9)[i] > 0.5 && (*hOvrE)[i] < 0.12) || (fabs((*eta)[i]) > 1.5 && (*r9)[i] > 0.8 && (*hOvrE)[i] < 0.10))){
                    ROOT::Math::PtEtaPhiMVector v1((*et)[i],(*eta)[i],(*phi)[i],0.0);
                    for(int j = 0; j < nEgsIn; j ++){
                        if (j != i){
                            if ((*et)[j] > 22.0 && fabs((*eta)[j]) < 2.55 && ((fabs((*eta)[j]) < 1.5 && (*r9)[j] > 0.5 && (*hOvrE)[j] < 0.12 && !((*r9)[j] < 0.85 && ((*sIeIe)[j] > 0.015 || (*phoIso)[j] > 6 || (*ecalIso)[j] > 6))) || (fabs((*eta)[j]) > 1.5 && (*r9)[j] > 0.8 && (*hOvrE)[j] < 0.10 && !((*r9)[j] < 0.90 && ((*sIeIe)[j] > 0.035 || (*phoIso)[j] > 6 || (*ecalIso)[j] > 6))))){
                                ROOT::Math::PtEtaPhiMVector v2((*et)[j],(*eta)[j],(*phi)[j],0.0);
                                auto v3 = v1+v2;
                                double massTmp = v3.M();
                                if(massTmp > mass){
                                    mass = massTmp;
                                }
                            }
                        }
                    }
                }
            }
        }
    }//nEgs
    if(bothPass != 1){
        mass = -999;
        iOut = -1;
        jOut = -1;
    }
    *pho1 = iOut;
    *pho2 = jOut;
    return mass;
}
double calcDR(float eta1, float eta2, float phi1, float phi2){
    double DR = 0.0;
    double dE = eta1-eta2;
    double dP = phi1-phi2;
    if (dP > M_PI) dP -= 2*M_PI;
    if (dP < -M_PI) dP += 2*M_PI;
    
    DR = sqrt(dE*dE + dP*dP);
    
    return DR;
}
