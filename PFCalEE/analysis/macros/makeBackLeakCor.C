#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TMath.h"

int makeBackLeakCor(const unsigned nLayers,
		    const unsigned nBack,
		    const unsigned iSR,
		    const unsigned pT,
		    const unsigned eta,
		    const unsigned pu,
		    const double offset,
		    const double calib,
		    double & backLeakCor,
		    TCanvas * mycE2D,
		    TTree *ltree,
		    TFile *outfile,
		    TGraphErrors * corrBackLeakFit,
		    const TString & plotDir
		    ){

  double Eval = E(pT,eta);
  
  mycE2D->cd();
  gPad->SetRightMargin(0.15);
  gPad->SetLogz();

  std::ostringstream hName;
  hName.str("");
  hName << "energy" << pT << "_vsBackFraction";

  const int enRegions = 4;
  const int nxbins[enRegions] = {5,7,6,10};
  const int nybins[enRegions] = {50,50,50,50};
  double en0xbins[nxbins[0]+1] = {0.0001,0.003,0.005,0.007,0.02,0.05};
  double en1xbins[nxbins[1]+1] = {0.0001,0.003,0.005,0.007,0.01,0.02,0.03,0.05};
  double en2xbins[nxbins[2]+1] = {0.0001,0.004,0.007,0.01,0.015,0.03,0.08};
  double en3xbins[nxbins[3]+1] = {0.0001,0.004,0.005,0.006,0.008,0.010,0.013,0.016,0.02,0.04,0.11};
  double *xbins[enRegions];
  xbins[0] = en0xbins;
  xbins[1] = en1xbins;
  xbins[2] = en2xbins;
  xbins[3] = en3xbins;

  double minEnFrac[enRegions];
  double maxEnFrac[enRegions];
  if(nBack==2) {
    minEnFrac[0] = 0.85;
    minEnFrac[1] = 0.90;
    minEnFrac[2] = 0.93;
    minEnFrac[3] = 0.95;

    maxEnFrac[0] = 1.2;
    maxEnFrac[1] = 1.1;
    maxEnFrac[2] = 1.4;
    maxEnFrac[3] = 1.04;
  }
  else if(nBack==3) {
    minEnFrac[0] = 0.80;
    minEnFrac[1] = 0.90;
    minEnFrac[2] = 0.88;
    minEnFrac[3] = 0.9;

    maxEnFrac[0] = 1.2;
    maxEnFrac[1] = 1.1;
    maxEnFrac[2] = 1.1;
    maxEnFrac[3] = 1.1;
  }
  else if(nBack==4) {
    minEnFrac[0] = 0.75;
    minEnFrac[1] = 0.85;
    minEnFrac[2] = 0.85;
    minEnFrac[3] = 0.85;

    maxEnFrac[0] = 1.25;
    maxEnFrac[1] = 1.15;
    maxEnFrac[2] = 1.15;
    maxEnFrac[3] = 1.15;
  }
  else {
  std::cout << "fix!" << std::endl;
  std::exit(0);
}

    
  TH2F *p_ErecovsEback;
  unsigned enIdx = 0;
  if(pT > 101.) enIdx = 3;
  else if(pT > 31.) enIdx = 2;
  else if(pT > 11.) enIdx = 1;

  p_ErecovsEback = new TH2F(hName.str().c_str(), "", nxbins[enIdx], xbins[enIdx], nybins[enIdx],
			    minEnFrac[enIdx]*Eval, maxEnFrac[enIdx]*Eval);
  
  if (!p_ErecovsEback){
    std::cout << " -- ERROR, pointer for histogram " << hName.str() << " is null." << std::endl;
    return 1;
  }

  std::ostringstream lName;
  std::string lNameTot,lNameBack;
  getTotalEnergyString(nLayers,nBack,lNameTot,lNameBack,iSR);

  lName.str("");
  lName << std::setprecision(6);
  lName << "(";
  lName << lNameTot;
  lName << " - " << offset << ")/" << calib;
  lName << ":(" << lNameBack << ")/(" << lNameTot << ")";
  lName << ">>"+hName.str();  
  ltree->Draw(lName.str().c_str(),"","colz");

  std::cout << lName.str()  << std::endl;

  p_ErecovsEback->SetTitle(";f_{back} #equiv E_{back}/E_{tot};E_{tot} [GeV]");
  p_ErecovsEback->Draw("colz");

  std::unordered_map<unsigned,float> radius_map;
  radius_map[4] = 26.0;

  char buf1[500];
  sprintf(buf1,"#gamma, PU 0");
  TLatex lat1;
  lat1.SetTextSize(0.03);
  lat1.DrawLatexNDC(0.15,0.9,buf1);
  sprintf(buf1,"r = %3.0f mm", radius_map[iSR]);
  lat1.DrawLatexNDC(0.15,0.87,buf1);
  sprintf(buf1,(("|#eta| = " + std::to_string(eta)).c_str()));
  lat1.DrawLatexNDC(0.15,0.84,buf1);
  sprintf(buf1,(("pT = " + std::to_string(pT) + "GeV").c_str()));
  lat1.DrawLatexNDC(0.15,0.81,buf1);
  lat1.SetTextSize(0.05);
  lat1.DrawLatexNDC(0.13,0.95,"HGCAL G4 standalone");
  lat1.DrawLatexNDC(0.71,0.95,("nBack="+std::to_string(nBack)).c_str());

  //normalize histogram by x bin width before the fit
  for(int i(1); i<=nxbins[enIdx]; ++i) {
    for(int j(1); j<=nybins[enIdx]; ++j) {
      double cont = p_ErecovsEback->GetBinContent(i, j);
      double width = p_ErecovsEback->GetXaxis()->GetBinWidth(i);
      p_ErecovsEback->SetBinContent(i, j, cont/width);
    }
  }
      
  lName << "_pfx";
  TProfile *tmpProf = p_ErecovsEback->ProfileX(lName.str().c_str());
  tmpProf->SetMarkerStyle(20);
  tmpProf->SetMarkerColor(1);
  //tmpProf->SetMarkerSize(5);
  tmpProf->SetLineColor(1);
  tmpProf->Draw("PEsame");
  tmpProf->Fit("pol1","","same");
  
  //tmpProf->SetStats(1);
  //gStyle->SetOptFit(1111);
  TF1 *fitcor = (TF1*)tmpProf->GetFunction("pol1");
  if (!fitcor) {
    std::cout << " Fit failed for back leakage correction" << std::endl;
    mycE2D->Update();
    backLeakCor = 0;
    if (Eval < 1000) {
      //for low energies: ignore the fit and return success....
      return 0;
    }
    return 1;
  }
  else {
    backLeakCor = fitcor->GetParameter(1);
    std::cout << " ---- back leakage correction factor: " << backLeakCor << " +/- " << fitcor->GetParError(1) << std::endl;
    char buf[500];
    TLatex lat;
    lat.SetTextSize(0.03);
    sprintf(buf,"E=%3.3f #times f_{back} + %3.3f",backLeakCor,fitcor->GetParameter(0));
    lat.DrawLatexNDC(0.45,0.90,buf);
		    
    Int_t np=corrBackLeakFit->GetN();
    //if (!dovsE) corrBackLeakFit[iSR]->SetPoint(np,genEn[iE],backLeakCor[oldIdx[iE]][iSR]);
    corrBackLeakFit->SetPoint(np,Eval,backLeakCor);
    corrBackLeakFit->SetPointError(np,0.0,fitcor->GetParError(1));
    
  }

  if (system(TString("mkdir -p ")+plotDir+TString("/BackCor2D"))) return 1;
  
  std::ostringstream saveName;
  saveName.str("");
  saveName << plotDir << "BackCor2D/ErecovsbackFraction_eta" << eta << "_pu" << pu;
  saveName << "_E" << pT << "_SR" << iSR;
  mycE2D->Update();
  mycE2D->Print((saveName.str()+".png").c_str());
  mycE2D->Print((saveName.str()+".C").c_str());
  
  
  outfile->cd();
  p_ErecovsEback->Write();
  tmpProf->Write();

  
  return 0;
};

int plotBackLeakFit(const TString & plotDir,
		    TGraphErrors *corrBackLeakFit,
		    const unsigned eta,
		    const unsigned pu){
  
  TCanvas *myc = new TCanvas("mycBF","mycBF",1);
  myc->cd();
  corrBackLeakFit->SetTitle(";E [GeV];back cor");
  corrBackLeakFit->Draw("APE");
  myc->Update();
  std::ostringstream lsave;
  lsave.str("");

  if (system(TString("mkdir -p ")+plotDir+TString("/BackCor"))) return 1;

  lsave << plotDir << "/BackCor/";
  lsave << "BackLeakCor";
  lsave << "_eta" << eta << "_pu" << pu ;
  lsave << "_vsE";
  myc->Print((lsave.str()+".png").c_str());
  myc->Print((lsave.str()+".C").c_str());

  return 0;
};
