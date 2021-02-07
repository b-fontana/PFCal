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


std::vector<double> GetCalib(/*const std::string file*/){
  //TBD from file
  std::vector<double> calib;
  calib.resize(6,1);
  return calib;
};


TPad* plotCalibration(TGraphErrors *& gr,
		      TPad *pad,
		      bool doRatio, 
		      TGraphErrors *& grDelta,
		      std::string unit, 
		      double & calib,double & calibErr, 
		      double & offset, double & offsetErr,
		      const bool dovsE){

  TPad *upper = 0;
  TPad *lower = 0;

  if (!doRatio) pad->cd();
  else {
    pad->Clear();
    upper = plot_ratio(pad, true);
    lower = plot_ratio(pad, false);
    upper->cd();
  }
  gr->SetTitle("");
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(1);
  gr->SetLineColor(1);
  gr->GetXaxis()->SetLabelSize(0.0);
  gr->GetXaxis()->SetTitleSize(0.0);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(7);
  gr->GetYaxis()->SetTitleOffset(0.9);
  
  gr->Draw("ap");

  if (!doRatio){
    if (!dovsE) gr->GetXaxis()->SetTitle("p_{T} (GeV)");
    else gr->GetXaxis()->SetTitle("E (GeV)");
  }
  else gr->GetXaxis()->SetTitle("");
  
  gr->GetYaxis()->SetTitle(("Average energy deposited ("+unit+")").c_str()); 
  char buf[500];
  TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
  fitFunc->SetLineColor(6);

  if (dovsE) gr->Fit(fitFunc,"IME","same");
  else gr->Fit(fitFunc,"IME","same");
  TLatex lat;
  lat.SetTextColor(6);
  lat.SetTextSize(0.05);
  if (!dovsE) sprintf(buf,"<E> #propto a + b #times p_{T} ");
  else sprintf(buf,"<E> #propto a + b #times E ");
  lat.DrawLatexNDC(0.16,0.85,buf);
  sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),unit.c_str());
  lat.DrawLatexNDC(0.16,0.75,buf);
  sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),unit.c_str());
  lat.DrawLatexNDC(0.16,0.65,buf);
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatexNDC(0.16,0.55,buf);

  calib = fitFunc->GetParameter(1);
  offset = fitFunc->GetParameter(0);
  calibErr = fitFunc->GetParError(1);
  offsetErr = fitFunc->GetParError(0);

  if (doRatio){
    //draw deltaE/E vs E
    lower->cd();
    gPad->SetLogx(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    double loffset = fitFunc->GetParameter(0);
    double lslope = fitFunc->GetParameter(1);
    double range = 0.05;
    if (unit=="GeV") {
      loffset=0;
      lslope=1;
      range = 0.1;
    }

    //fill delta
    for (int ip(0);ip<gr->GetN();++ip){
      double x=0;
      double y=0;
      gr->GetPoint(ip,x,y);
      grDelta->SetPoint(ip,x,((y-loffset)/lslope-x)/x);
      double err = gr->GetErrorY(ip)/lslope*1./x;
      grDelta->SetPointError(ip,0,err);
      std::cout << "Calib " << ip << " Egen=" << x << " Erec=" << y << " delta=" << ((y-loffset)/lslope-x)/x << std::endl;
    }
    grDelta->SetTitle("");
    grDelta->SetMinimum(-1.*range);
    grDelta->SetMaximum(range);
    grDelta->GetXaxis()->SetLabelSize(0.15);
    grDelta->GetXaxis()->SetTitleSize(0.15);
    grDelta->GetYaxis()->SetLabelSize(0.12);
    grDelta->GetYaxis()->SetTitleSize(0.15);
    grDelta->GetXaxis()->SetTitleOffset(0.8);
    grDelta->GetYaxis()->SetTitleOffset(0.4);
    grDelta->SetMarkerStyle(20);
    grDelta->SetMarkerColor(1);
    grDelta->SetLineColor(1);
    
    grDelta->Draw("ap");
    //grDelta->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
    if (!dovsE) grDelta->GetXaxis()->SetTitle("p_{T} (GeV)");
    else grDelta->GetXaxis()->SetTitle("E (GeV)");
    grDelta->GetYaxis()->SetTitle("(#Delta E)/E");

    TLine *line = new TLine(grDelta->GetXaxis()->GetXmin(),0,grDelta->GetXaxis()->GetXmax(),0);
    line->SetLineColor(2);//kYellow+4);
    line->Draw();
    

    lower->Update();

  }

  return upper;
};


int makeCalibration(const bool doRaw,
		    const bool doBackLeakCor,
		    const unsigned eta,
		    const unsigned pu,
		    const double radius,
		    TGraphErrors *calibRecoFit,
		    TGraphErrors *calibRecoDelta,
		    double & calib,
		    double & calibErr,
		    double & offset,
		    double & offsetErr,
		    TString plotDir,
		    TFile *foutEfit){
  
  double etaval = eta/10.;

  TCanvas *mycC = new TCanvas("mycC","mycCalib",1);
  
  std::string unit = "MIPS";
  if (!doRaw) unit = "GeV";

  TPad *lpad = (TPad*)(mycC->cd());
  TPad *upper = plotCalibration(calibRecoFit,lpad,
				true,calibRecoDelta,
				unit,
				calib,calibErr,
				offset,offsetErr,
				true);
  upper->cd();
  char buf[500];
  sprintf(buf,"#gamma #eta=%3.1f + PU %d",etaval,pu);
  TLatex lat;
  lat.SetTextSize(0.05);
  lat.DrawLatexNDC(0.5,0.17,buf);
  sprintf(buf,"r = %3.0f mm",radius);
  lat.DrawLatexNDC(0.5,0.3,buf);
  lat.SetTextSize(0.04);
  lat.DrawLatexNDC(0.6,0.96,"HGCAL G4 standalone");
  

  mycC->Update(); 
  if (system(TString("mkdir -p ")+plotDir+TString("/Calib"))) return 1;

  std::ostringstream lsave;
  lsave.str("");
  lsave << plotDir << "/Calib/";
  if (doRaw) lsave << "CalibMipToGeV";
  else if (doBackLeakCor) lsave << "CalibBackLeakCor";
  else lsave << "Calib";
  lsave << "_eta" << eta << "_pu" << pu ;
  //if (dovsE) 
  lsave << "_vsE";
  mycC->Print((lsave.str()+".png").c_str());
  mycC->Print((lsave.str()+".C").c_str());

  foutEfit->cd();
  calibRecoFit->Write();
  calibRecoDelta->Write();

  return 0;

};


