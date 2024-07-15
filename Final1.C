#define Final1_cxx
#include "Final1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

const int npTbins = 8;
const int ndphibins = 2;

const int massBins = 20;
const double massLow = 2;
const double massHigh = 5;
//double pTbins[npTbins+1] = {0,0.5,1,1.5,2,2.5,3,4,5};
const int file = 0;
const int inout = 0;
int ptBins;

TH1D *h_[2][2][npTbins];

void Final1::Loop()
{

for(int k = 0; k<ndphibins; k++){
for(int j = 0; j<2; j++){
for(int i = 0; i<npTbins; i++){
        h_[k][j][i] = new TH1D(Form("h_%d%d%d",k,j,i), "", massBins, massLow,massHigh);
        h_[k][j][i]->Sumw2();
}}}

const Float_t pi = TMath::Pi();

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<nentries<<endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

if (jentry%150000==0)
{
cout<<jentry<<" of "<<nentries<<endl;
}

      Float_t psi_s = CompactEPArray_psi[17];   //TRUE
      Float_t psi_n = CompactEPArray_psi[18];   //TRUE
//      Float_t psi_n = CompactEPArray_psi[17];   //SAME_ARM (Systematic)
//      Float_t psi_s = CompactEPArray_psi[18];   //SAME_ARM (Systematic)
//      Float_t psi_s = CompactEPArray_psi[16];   //CNT (Systematic)
//      Float_t psi_n = CompactEPArray_psi[16];   //CNT (Systematic)
//      Float_t psi_s = CompactEPArray_psi[10];   //BBC REACTION PLANE (Systematic)
//      Float_t psi_n = CompactEPArray_psi[11];   //BBC REACTION PLANE (Systematic)



      Float_t Evt_Cent = Evt_Cent[nentries];

      for (Int_t idmu = 0; idmu<nDiMuons; idmu++){

        Float_t mass = DiMuons_mass[idmu];
        Bool_t se = DiMuons_same_event[idmu];
        Short_t charge = DiMuons_charge[idmu];
        Float_t dcar = DiMuons_dca_r[idmu];
        Float_t px = DiMuons_Px[idmu];
        Float_t py = DiMuons_Py[idmu];
        Float_t rapidity = DiMuons_rapidity[idmu];
        Short_t lg0 = DiMuons_Tr0_lastgap[idmu];
        Short_t lg1 = DiMuons_Tr1_lastgap[idmu];
        Float_t chi = DiMuons_Evt_vtxchi2[idmu];
        Float_t px0 = DiMuons_Tr0_px[idmu];
        Float_t py0 = DiMuons_Tr0_py[idmu];
        Float_t rap0 = DiMuons_Tr0_rapidity[idmu];
        Float_t rap1 = DiMuons_Tr1_rapidity[idmu];
        Float_t px1 = DiMuons_Tr1_px[idmu];
        Float_t py1 = DiMuons_Tr1_py[idmu];
        Float_t pz0 = DiMuons_Tr0_pz[idmu];
        Float_t pz1 = DiMuons_Tr1_pz[idmu];
        Float_t DG0_0 = DiMuons_Tr0_DG0[idmu];
        Float_t DG0_1 = DiMuons_Tr1_DG0[idmu];
        Float_t DDG0_0 = DiMuons_Tr0_DDG0[idmu];
        Float_t DDG0_1 = DiMuons_Tr1_DDG0[idmu];
        Float_t trchi2_0 = DiMuons_Tr0_trchi2[idmu];
        Float_t trchi2_1 = DiMuons_Tr1_trchi2[idmu];
        Float_t idchi2_0 = DiMuons_Tr0_idchi2[idmu];
        Float_t idchi2_1 = DiMuons_Tr1_idchi2[idmu];
        Float_t fvtxchi2_0 = DiMuons_Tr0_chi2_fvtxmutr[idmu];
        Float_t fvtxchi2_1 = DiMuons_Tr1_chi2_fvtxmutr[idmu];
        Float_t fvtxnhits_0 = DiMuons_Tr0_nhits_fvtx[idmu];
        Float_t fvtxnhits_1 = DiMuons_Tr1_nhits_fvtx[idmu];
        Float_t dcar0 = DiMuons_Tr0_dca_r[idmu];
        Float_t dcar1 = DiMuons_Tr1_dca_r[idmu];
        Float_t xst01 = DiMuons_Tr0_xst1[idmu];
        Float_t xst11 = DiMuons_Tr1_xst1[idmu];
        Float_t yst01 = DiMuons_Tr0_yst1[idmu];
        Float_t yst11 = DiMuons_Tr1_yst1[idmu];

        Short_t ntrhits0 = DiMuons_Tr0_ntrhits[idmu];
        Short_t ntrhits1 = DiMuons_Tr1_ntrhits[idmu];
        Short_t nidhits0 = DiMuons_Tr0_nidhits[idmu];
        Short_t nidhits1 = DiMuons_Tr1_nidhits[idmu];

        Float_t pt = TMath::Sqrt(px*px + py*py);
        Float_t pt0 = TMath::Sqrt(px0*px0 + py0*py0);
        Float_t pt1 = TMath::Sqrt(px1*px1 + py1*py1);
        Float_t p0 = TMath::Sqrt(px0*px0 + py0*py0 + pz0*pz0);
        Float_t p1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);

        Int_t octant0 = ((int((atan2(xst01,yst01)+TMath::Pi())/(TMath::Pi()/8))+1)/2)%8;
        Int_t octant1 = ((int((atan2(xst11,yst11)+TMath::Pi())/(TMath::Pi()/8))+1)/2)%8;

        Float_t DG0cut0 = p0*DG0_0;
        Float_t DG0cut1 = p1*DG0_1;

        Float_t DDG0cut0 = p0*DDG0_0;
        Float_t DDG0cut1 = p1*DDG0_1;

        Float_t phi = atan2(py,px);

        Float_t weight = 1./(0.562425*exp(-0.5*((Evt_Cent-24.9)/(31.99))**2));

        if (rapidity<0){
        Float_t dphi = TMath::Abs(phi-psi_n);
        } else{
        Float_t dphi = TMath::Abs(phi-psi_s);
        }


if (dphi<=7 && charge==0 && octant0!=octant1 && chi<3 && TMath::Abs(rapidity)<=2.2 && TMath::Abs(rapidity)>=1.2 && mass>=2 && mass<=5 && TMath::Abs(rap0)<2.2 && TMath::Abs(rap0)>1.2 && TMath::Abs(rap1)<2.2 && TMath::Abs(rap1)>1.2 && Evt_Cent<=60 && Evt_Cent>=10 && px!=-9999.900 && py!=-9999.900 && (lg0==3||lg0==4) && (lg1==3||lg1==4) && pt0>=1. && pt1>=1. && DG0cut0<60 && DG0cut1<60 && DDG0cut0<40 && DDG0cut1<40 && trchi2_0<10 && trchi2_1<10 && idchi2_0<3 && idchi2_1<3 && ntrhits0>=6 && ntrhits1>=6 && nidhits0>3 && nidhits1>3 && TMath::Abs(pz0)>=3. && TMath::Abs(pz1)>=3.) {

//(dphi<=pi/4||(dphi>=3*pi/4 && dphi<=5*pi/4) ) in plane
//((dphi>=pi/4 && dphi<=3*pi/4) || dphi>=5*pi/4) out of plane

        if(pt<=1.0){ptBins=0;}
        else if(pt>1.0 && pt<=2.0){ptBins=1;}
        else if(pt>2.0 && pt<=3.0){ptBins=2;}
        else if(pt>3.0 && pt<=5.0){ptBins=3;}


        if(se){int evtType=0;}
        else{evtType=1;}

        if(dphi<=pi/4||(dphi>=3*pi/4 && dphi<=5*pi/4)){int dphiType=0;}
        else{dphiType=1;}

        switch(ptBins){
        case ptBins:
        switch(evtType){
        case evtType:
        switch(dphiType){
        h_[dphiType][evtType][ptBins]->Fill(mass,weight);
        break;
        }
        break;
        }
        break;
        }

}
}
}

for(int k = 0; k<ndphibins; k++){
for(int j = 0; j<2; j++){
for(int i = 0; i<npTbins; i++){
        h_[k][j][i]->SaveAs(Form("../WEIGHT/TRUE/%d%d%d_%d.root",k,j,i,file));
}
}
}
