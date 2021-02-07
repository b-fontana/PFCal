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

bool plotResolution(TGraphErrors *gr,TPad *pad,
		    const unsigned pu,
		    const double & stoch0,
		    const double & const0,
		    const double & noise0,
		    double & stoch,double & stochErr, 
		    double & constant, double & constErr,
		    double & noise,double & noiseErr,
		    const bool dovsE,
		    const bool doRaw){

  pad->cd();
  gr->SetTitle("");
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(1);
  gr->SetLineColor(1);
  gr->GetXaxis()->SetLabelSize(0.06);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(1.);
  gr->GetYaxis()->SetTitleOffset(1.);
  gr->SetMinimum(0);
  gr->SetMaximum(0.3);
  gr->Draw("ap");
  if (!dovsE) gr->GetXaxis()->SetTitle("p_{T} (GeV)");
  else gr->GetXaxis()->SetTitle("E (GeV)");
  gr->GetYaxis()->SetTitle("#sigma/E");

  TF1 *fitFunc = new TF1(doRaw?"resoRaw":"reso","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());

  std::cout << " -- Setting inital function parameters to: " << stoch0 << " " << const0 << " " << noise0/2. << std::endl;

  fitFunc->SetParameter(0,stoch0);
  fitFunc->SetParLimits(0,0,1);
  fitFunc->SetParameter(1,const0);
  fitFunc->SetParLimits(1,0,1);
  fitFunc->SetParameter(2,pu==0?noise0:noise0/2.);
  //fitFunc->SetParLimits(2,0,noise0);

  if (pu==0) 
    fitFunc->FixParameter(2,noise0);
  if (pu>0) fitFunc->FixParameter(1,const0);
  
  //std::cout << " Initial params: "  << fitFunc->GetParameter(0) << " "<< fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2)
  //<< std::endl;

  int status = gr->Fit(fitFunc,"BIME");

  if (fitFunc->GetChisquare()/fitFunc->GetNDF()>50){
    fitFunc->ReleaseParameter(1);
    fitFunc->SetParameters(stoch0,const0,noise0);
    if (!dovsE) status = gr->Fit(fitFunc,"BIME");//,"",5,110);
    else status = gr->Fit(fitFunc,"BIME");//,"",10,400);
  }

  fitFunc->SetLineColor(6);
  fitFunc->Draw("same");

  stoch = fitFunc->GetParameter(0);
  stochErr = fitFunc->GetParError(0);
  constant = fitFunc->GetParameter(1);
  constErr = fitFunc->GetParError(1);
  noise = fitFunc->GetParameter(2);
  noiseErr = fitFunc->GetParError(2);

  char buf[500];
  TLatex lat;
  lat.SetTextSize(0.04);
  lat.SetTextColor(6);
  if (!dovsE) sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{p_{T}}} #oplus c #oplus #frac{n}{p_{T}}");
  else sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c #oplus #frac{n}{E}");

  lat.DrawLatexNDC(0.5,0.85,buf);
  sprintf(buf,"s=%3.3f #pm %3.3f",stoch,stochErr);
  lat.DrawLatexNDC(0.5,0.75,buf);
  sprintf(buf,"c=%3.3f #pm %3.3f",constant,constErr);
  lat.DrawLatexNDC(0.5,0.65,buf);
  sprintf(buf,"n=%3.3f #pm %3.3f",noise,noiseErr);
  lat.DrawLatexNDC(0.5,0.55,buf);
  //sprintf(buf,"status = %d, #chi^{2}/N = %3.1f/%d = %3.1f",status,fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  sprintf(buf,"#chi^{2}/N = %3.1f/%d = %3.1f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatexNDC(0.5,0.45,buf);
  
  if (status != 0 && status != 4000) {
    std::cout << " -- Fit failed with status " << status << std::endl;
    return false;
  }

  return true;
};


int makeResolution(const bool dovsE,
		   const bool doRaw,
		   const bool doBackLeakCor,
		   const unsigned eta,
		   const unsigned pu,
		   const double radius,
		   TGraphErrors *resoRecoFit,
		   const double sigmaStochRef,
		   const double sigmaConstRef,
		   const double noiseRef,
		   double & sigmaStoch,
		   double & sigmaStochErr,
		   double & sigmaConst,
		   double & sigmaConstErr,
		   double & sigmaNoise,
		   double & sigmaNoiseErr,
		   TString plotDir,
		   TFile *foutEfit)
{
  double etaval = eta/10.;
  
  TCanvas *mycR = new TCanvas("mycR","mycReso",1);
  
  TPad *lpad = (TPad*)(mycR->cd());

  double stoch0 = pu==0? (dovsE?0.25 : 0.14) : sigmaStochRef;
  double const0 = pu==0? 0.01 : sigmaConstRef;

  double noise0 = noiseRef;
  
  bool success = plotResolution(resoRecoFit,lpad,
				pu,
				stoch0,const0,noise0,
				sigmaStoch,
				sigmaStochErr,
				sigmaConst,
				sigmaConstErr,
				sigmaNoise,
				sigmaNoiseErr,
				dovsE,doRaw);
  lpad->cd();
  char buf[500];
  sprintf(buf,"#gamma #eta=%3.1f + PU %d",etaval,pu);
  TLatex lat;
  lat.SetTextSize(0.04);
  lat.DrawLatexNDC(0.16,0.85,buf);
  sprintf(buf,"r = %3.0f mm",radius);
  lat.DrawLatexNDC(0.16,0.75,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
  
  if (!success) {
    return 1;	    
  }
  
  
  mycR->Update();
  if (system(TString("mkdir -p ")+plotDir+TString("/Reso"))) return 1;

  std::ostringstream lsave;
  lsave.str("");
  lsave << plotDir << "/Reso/";
  if (doRaw) lsave << "ResoRaw";
  else lsave << "Reso";
  if (doBackLeakCor) lsave << "BackLeakCor";
  lsave << "_eta" << eta << "_pu" << pu ;
  if (dovsE) lsave << "_vsE";
  mycR->Print((lsave.str()+".png").c_str());
  mycR->Print((lsave.str()+".C").c_str());
  
  foutEfit->cd();
  resoRecoFit->Write();

  return 0;

};//main
