#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"

#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

int main(int argc, char** argv){//main  

  if (argc < 4) {
    std::cout << " Usage: " 
	      << argv[0] << " <nEvts to process (0=all)>"
	      << " <full path to input file>"
	      << " <full path to output file>" 
	      << " <optional: debug (default=0)>"
	      << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////
  //// Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////
  bool concept = true;

  double minX=-1700,maxX=1700;
  double minY=-1700,maxY=1700;
  double minZ=3170,maxZ=5070;
  //double minX=-510,maxX=510;
  //double minY=-510,maxY=510;
  //double minZ=-1000,maxZ=1000;

  unsigned nX=(maxX-minX)/10,nY=(maxY-minY)/10;
  //if (nLayers == nEcalLayers) {
  //maxZ = 3500;
  //}
  unsigned nZ=maxZ-minZ;

  //////////////////////////////////////////////////////////
  //// End Hardcoded config ////////////////////////////////////
  //////////////////////////////////////////////////////////

  const unsigned pNevts = atoi(argv[1]);
  std::string filePath = argv[2];
  std::string outPath = argv[3];

  if (argc >4) debug = atoi(argv[4]);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  //initialise calibration class

  HGCSSCalibration mycalib(filePath,concept,calibrate);
  const unsigned nLayers = mycalib.nLayers();//64;//Calice 54;//Scint 9;//HCAL 33;//All 64
  const unsigned nEcalLayers = mycalib.nEcalLayers();//31;
  const unsigned nFHcalLayers = mycalib.nFHcalLayers();//concept 24;//calice 47

  HGCSSDigitisation myDigitiser;
  //myDigitiser.setRandomSeed(1);

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- nEcalLayers = " << nEcalLayers 
	    << ", mip weights = " << mycalib.mipWeight(0) << " " << mycalib.mipWeight(1) << " " << mycalib.mipWeight(2)
	    << ", GeV weights = " << mycalib.gevWeight(0) << " offset " << mycalib.gevOffset(0)
	    << std::endl
	    << " -- nFHcalLayers  = " << nFHcalLayers 
	    << ", mip weights = " << mycalib.mipWeight(3) << " " << mycalib.mipWeight(4) << " " << mycalib.mipWeight(5) 
	    << ", GeV weights = " << mycalib.gevWeight(1) << " " << mycalib.gevWeight(2) << " offsets: " << mycalib.gevOffset(1) << " " <<   mycalib.gevOffset(2)
	    << std::endl
	    << " -- conversions: HcalToEcalConv = " <<mycalib.HcalToEcalConv() << " BHcalToFHcalConv = " << mycalib.BHcalToFHcalConv() << std::endl
	    << " -- Elim for Cglobal = " << pElim[0] << "-" << pElim[nLimits-1] << " Mips" << std::endl
	    << " -- HCAL time threshold: " <<  mycalib.hcalTimeThreshold() << " ns " << std::endl
	    << " ----------------------------------- " << std::endl;

  if (applyDigi){
    myDigitiser.Print(std::cout);
  }
  else {
    std::cout << " -- No digi applied -- " << std::endl;
  }

  std::cout << " ----------------------------------- " << std::endl
	    << " -----------------------------------" << std::endl;

  std::cout << " -- Output file directory is : " << plotDir << std::endl;

  TFile *outputFile;
  if (!lightOutput) outputFile = TFile::Open(plotDir+"/CalibHistosNew_"+pSuffix+".root","RECREATE");
  else outputFile = TFile::Open(plotDir+"/CalibHistosNew_"+pSuffix+"_light.root","RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << plotDir << "/CalibHistos.root cannot be opened. Please create output directory : " << plotDir << "  Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }


  std::cout << " -- 2-D histograms: " << std::endl
	    << " -- X: " << nX << " " << minX << " " << maxX << std::endl
	    << " -- Y: " << nY << " " << minY << " " << maxY << std::endl
	    << " -- Z: " << nZ << " " << minZ << " " << maxZ << std::endl
    ;

  TH2F *p_xy[nGenEn][nLayers];
  TH3F *p_xyz[nGenEn];
  TH2F *p_xy_evt[nLayers];
  TH3F *p_xyz_evt = 0;
  TH2F *p_recoxy[nGenEn][nLayers];
  TH3F *p_recoxyz[nGenEn];
  TH2F *p_recoxy_evt[nLayers];

  //std::cout << " -- Checking layer-eta conversion:" << std::endl;

  for (unsigned iL(0); iL<nLayers; ++iL){
    p_recoxy_evt[iL] = 0;
    p_xy_evt[iL] = 0;
    //std::cout << iL << " " << zlay[iL] << " " 
    //<< llrBinWidth[iL] << " y=0@eta3 " 
    //<< getEta(0,"eta30",zlay[iL]) << std::endl;
  }

  TH3F *p_recoxyz_evt = 0;
  TH2F *p_EvsLayer[nGenEn];
  TH1F *p_time[nGenEn];
  TH1F *p_Etot[nGenEn][nLayers];
  TH1F *p_Efrac[nGenEn][nLayers];
  TH1F *p_Etotal[nGenEn][nSections];
  TH1F *p_EfracLateral[nGenEn][nSections][nLat];
  TH2F *p_HCALvsECAL[nGenEn];
  TH2F *p_BHCALvsFHCAL[nGenEn];
  TH1F *p_Cglobal[nGenEn][nLimits];
  TH1F *p_Eshower[nGenEn][nLimits];
  TH2F *p_EvsCglobal[nGenEn][nLimits];


  TH1F *p_Ereco[nGenEn][nSmear];
  TH1F *p_nSimHits[nGenEn];
  TH1F *p_nRecHits[nGenEn];

  TH1F *p_Edensity[nGenEn][nLayers];
  TH2F *p_occupancy[nOcc][nLayers];

  std::ostringstream lName;

  for (unsigned iO(0); iO<nOcc; ++iO){
    for (unsigned iL(0); iL<nLayers; ++iL){
      lName.str("");
      lName << "p_occupancy_" << occThreshold[iO] << "_" << iL;
      p_occupancy[iO][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",totalWidth/10.,-totalWidth/2.,totalWidth/2.,totalWidth/10.,-totalWidth/2.,totalWidth/2.);
      //p_occupancy[iO][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",100,getEta(totalWidth/2.,pEta),getEta(-totalWidth/2.,pEta),100,);
    }
  }


  bool isG4Tree = true;
  
  //for basic control plots
  TCanvas *mycE = new TCanvas("mycE","mycE",1500,1000);
  mycE->Divide(5,2);
  gStyle->SetOptStat("eMRuoi");


  double EmipMean[nGenEn];


  for (unsigned iE(0); iE<nGenEn; ++iE){

    if (isMB) genEn[iE] = iE;

    EmipMean[iE] = 1.;


    std::cout << "- Processing energy : " << genEn[iE] << std::endl;
    TString genEnStr;
    genEnStr += genEn[iE];
    std::ofstream hitsOut;
    if (saveTxtOutput) hitsOut.open(plotDir+"/hits_"+genEnStr+"GeV.dat");
    std::ofstream optim;
    if (saveTxtOutput) optim.open(plotDir+"/optimisation_"+genEnStr+"GeV.dat");
    std::ofstream edensity;
    genEnStr = "";
    genEnStr += iE;
    if (saveLLR) edensity.open(plotDir+"/edensity_"+genEnStr+".dat");

    double E1=0;
    double E2=0;
    double E3=0;


    std::ostringstream input;
    input << filePath ;
    //if (filePath.find("pi+")!=filePath.npos) input << "/pt_" ;
    //HGcal_version_8_e50.root
    if (pEOS) {
      input << "_e" ;
      input << genEn[iE] << fileName;
    }
    else {
      input << "/e_" ;
      input << genEn[iE];
      input << "/" << fileName;
    }
    TFile *inputFile = TFile::Open(input.str().c_str());

    if (!inputFile) {
      std::cout << " -- Error, input file " << input.str() << " cannot be opened. Trying next one..." << std::endl;
      continue;
    }

    TTree *lTree = (TTree*)inputFile->Get("HGCSSTree");
    if (!lTree){
      std::cout << " -- Error, tree HGCSSTree cannot be opened. Trying RecoTree instead." << std::endl;
      isG4Tree = false;
      lTree = (TTree*)inputFile->Get("RecoTree");
      if (!lTree){
	std::cout << " -- Error, tree RecoTree cannot be opened either. Exiting..." << std::endl;
	return 1;
      }
    }

    HGCSSEvent * event = 0;
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;

    if (isG4Tree){
      lTree->SetBranchAddress("HGCSSEvent",&event);
      lTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
      lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    }
    else {
      //lTree->SetBranchAddress("event",&event);
      //lTree->SetBranchAddress("volX0",&volX0);
      //lTree->SetBranchAddress("volLambda",&volLambda);
      lTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
      lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
    }

    const unsigned nEvts = ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;

    std::cout << "- Processing = " << nEvts  << " events out of " << lTree->GetEntries() << std::endl;
    
    TFile *puFile = 0;
    TTree *lPuTree = 0;
    std::vector<HGCSSSimHit> * pusimhitvec = 0;
    std::vector<HGCSSRecoHit> * purechitvec = 0;

    unsigned nPuEvts= 0;

    if (overlayPU){
      std::string inStr;// = "/afs/cern.ch/work/a/amagnan/public/";
      
      //if (isG4Tree) inStr += "HGCalEEGeant4/"+pVersion+"/"+pScenario+"/PedroPU/"+pEta+"/PFcal.root";
      //else  inStr += "HGCalEEDigi/"+pVersion+"/"+pScenario+"/PedroPU/"+pEta+"/DigiPFcal.root";
      if (isG4Tree) inStr += "140PU/PFcal_140PU_EEHE_full.root";
      puFile = TFile::Open(inStr.c_str());
      if (!puFile) {
	std::cout << " -- Error, input file " << inStr << " for PU cannot be opened. Exiting..." << std::endl;
	return 1;
      }
      else std::cout << " -- using file " << inStr << " for PU." << std::endl;
      lPuTree = (TTree*)puFile->Get("PUTree");
      if (!lPuTree){
	std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
	return 1;
      }

      if (isG4Tree){
	lPuTree->SetBranchAddress("HGCSSSimHitVec",&pusimhitvec);
      }
      else {
	lPuTree->SetBranchAddress("HGCSSSimHitVec",&pusimhitvec);
	lPuTree->SetBranchAddress("HGCSSRecoHitVec",&purechitvec);
      }
      nPuEvts = static_cast<unsigned>(lPuTree->GetEntries());
      //nPuEvts =  static_cast<unsigned>(lPuTree->GetEntries());

      std::cout << "- For PU: " << nPuEvts  << " single-interaction events are available." << std::endl;
      
    }


    //Initialise histos
    //necessary to have overflows ?
    gStyle->SetOptStat(1111111);
    double Etot[nLayers];
    double eLLR[nLayers];
    double zlay[nLayers];

    for (unsigned iL(0);iL<nLayers;++iL){
      Etot[iL] = 0;
      eLLR[iL] = 0;
      zlay[iL] = 0;
    }
    double Etotal[nSections];
    TH2F * EmipHits = new TH2F("EmipHits",";x(mm);y(mm)",34,-510,510,34,-510,510);

    double Elateral[nSections][nLat];
    for (unsigned iD(0); iD<nSections; ++iD){
      Etotal[iD] = 0;
      p_Etotal[iE][iD] = 0;
      for (unsigned iR(0); iR<nLat; ++iR){
	Elateral[iD][iR] = 0;
      }
    }

    unsigned nTotal = 0;
    unsigned nTotalSignal = 0;
    double Ereco[nSmear];
    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
      Ereco[iSmear]= 0;
    }
 

    lName.str("");
    lName << "p_xyz_" << genEn[iE];
    p_xyz[iE] = new TH3F(lName.str().c_str(),";z(mm);x(mm);y(mm)",
			 nZ,minZ,maxZ,
			 nX,minX,maxX,
			 nY,minY,maxY);
    if (!isG4Tree){
      lName.str("");
      lName << "p_recoxyz_" << genEn[iE];
      p_recoxyz[iE] = new TH3F(lName.str().c_str(),";z(mm);x(mm);y(mm)",
			       nZ,minZ,maxZ,
			       nX,minX,maxX,
			       nY,minY,maxY);
    }
    else p_recoxyz[iE] = 0;


    for (unsigned iLim(0); iLim<nLimits;++iLim){
      lName.str("");
      lName << "p_Cglobal_" << genEn[iE] << "_" << iLim;
      p_Cglobal[iE][iLim] =  new TH1F(lName.str().c_str(),";C_{global}",1000,0,2);
      
      lName.str("");
      lName << "p_Eshower_" << genEn[iE] << "_" << iLim;
      p_Eshower[iE][iLim] = new TH1F(lName.str().c_str(),";E_{shower} (GeV)",1000,0,2*genEn[iE]);

      lName.str("");
      lName << "p_EvsCglobal_" << genEn[iE] << "_" << iLim;
      p_EvsCglobal[iE][iLim] =  new TH2F(lName.str().c_str(),";C_{global};Etot (GeV)",
					 1000,0,2,
					 1000,0,2*genEn[iE]);

    }
    lName.str("");
    lName << "p_time_" << genEn[iE];
    p_time[iE] = new TH1F(lName.str().c_str(),";G4 time (ns)",1000,0,1000);

    p_nSimHits[iE] = 0;
    p_EvsLayer[iE] = 0;
    p_nRecHits[iE] = 0;
    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
      p_Ereco[iE][iSmear] = 0;
    }
    for (unsigned iR(0); iR<nLat; ++iR){
      lName.str("");
      lName << "p_EfracLateral_" << genEn[iE] << "_ECAL_" << iR;
      p_EfracLateral[iE][0][iR] =  new TH1F(lName.str().c_str(),";E_{total}^{x #times x}/E_{total}",101,0,1.01);
      lName.str("");
      lName << "p_EfracLateral_" << genEn[iE] << "_FHCAL_" << iR;
      p_EfracLateral[iE][1][iR] =  new TH1F(lName.str().c_str(),";E_{total}^{x #times x}/E_{total}",101,0,1.01);
      lName.str("");
      lName << "p_EfracLateral_" << genEn[iE] << "_BHCAL_" << iR;
      p_EfracLateral[iE][2][iR] =  new TH1F(lName.str().c_str(),";E_{total}^{x #times x}/E_{total}",101,0,1.01);
      lName.str("");
      lName << "p_EfracLateral_" << genEn[iE] << "_ALL_" << iR;
      p_EfracLateral[iE][3][iR] =  new TH1F(lName.str().c_str(),";E_{total}^{x #times x}/E_{total}",101,0,1.01);
    }
    for (unsigned iL(0); iL<nLayers; ++iL){
      lName.str("");
      lName << "p_xy_" << genEn[iE] << "_" << iL;
      p_xy[iE][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
			      nX,minX,maxX,
			      nY,minY,maxY);
      if (!isG4Tree){
	lName.str("");
	lName << "p_recoxy_" << genEn[iE] << "_" << iL;
	p_recoxy[iE][iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
				    nX,minX,maxX,
				    nY,minY,maxY);

      }
      else p_recoxy[iE][iL] = 0;
      Etot[iL] = 0;
      lName.str("");
      lName << "p_Efrac_" << genEn[iE] << "_" << iL;
      p_Efrac[iE][iL] = new TH1F(lName.str().c_str(),";integrated E_{layer}/E_{total}",101,0,1.01);
      lName.str("");
      lName << "p_Edensity_" << genEn[iE] << "_" << iL;
      p_Edensity[iE][iL] = new TH1F(lName.str().c_str(),";#eta",130,1.8,4.4);
    }



    bool firstEvent = true;
    unsigned ipuevt = 0;


    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
      if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
      else if (ievt%1000 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
      
      unsigned lEvt = ievt;
      bool proceed = evtList.size()==0?true:false;
      for (unsigned ie(0); ie<evtList.size();++ie){
	if (lEvt==evtList[ie]) {
	  proceed=true;
	  break;
	}
      }

      if (!proceed) continue;

      lTree->GetEntry(ievt);

      if (debug){
	if (isG4Tree) {
	  std::cout << "... Size of hit vector " << (*simhitvec).size() << std::endl;
	}
	else {
	  std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< std::endl;
	}
      }
      
      double Emax = 0;

      if (saveEventByEvent){
	
	if (p_xyz_evt) p_xyz_evt->Delete();
	if (!doOverview) p_xyz_evt = new TH3F("p_xyz",";z(mm);x(mm);y(mm)",nZ/2,minZ,maxZ,nX/2,minX,maxX,nY/2,minY,maxY);
	else p_xyz_evt = new TH3F("p_xyz",";z(mm);x(mm);y(mm)",nZ/20,minZ,maxZ,nX/8,minX,maxX,nY/8,minY,maxY);
	if (!isG4Tree) {
	  if (p_recoxyz_evt) p_recoxyz_evt->Delete();
	  p_recoxyz_evt = new TH3F("p_recoxyz",";z(mm);x(mm);y(mm)",nZ/2,minZ,maxZ,nX/2,minX,maxX,nY/2,minY,maxY);
	}
	
	for (unsigned iL(0); iL<nLayers; ++iL){
	  lName.str("");
	  lName << "p_xy_" << iL;
	  if (p_xy_evt[iL]) p_xy_evt[iL]->Delete();
	  p_xy_evt[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",nX,minX,maxX,nY,minY,maxY);
	  if (!isG4Tree){
	    lName.str("");
	    lName << "p_recoxy_" << iL;
	    if (p_recoxy_evt[iL]) p_recoxy_evt[iL]->Delete();
	    p_recoxy_evt[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",nX/2,minX,maxX,nY/2,minY,maxY);

	  }
	}
      }

      //add PU hits
      if (overlayPU){//add PU hits
	ipuevt = lRndm.Integer(nPuEvts);
	//get random PU events among available;
	std::cout << " -- PU Random number for entry " << ievt << " is " << ipuevt << std::endl;
	
	lPuTree->GetEntry(ipuevt);
	if (debug>1) std::cout << " -- Number of pu hits: " << (*pusimhitvec).size() << std::endl;
	for (unsigned iH(0); iH<(*pusimhitvec).size(); ++iH){//loop on hits
	  HGCSSSimHit lHit = (*pusimhitvec)[iH];
	  unsigned layer = lHit.layer();
	  //select "flavour"
	  if (doProtons){
	    //if ((lHit.nHadrons() == 0 && lHit.nMuons() == 0)) continue;
	    //HACK !!
	    if (lHit.nProtons() == 0) continue;
	  }
	  else if (doNeutrons){
	    //if (lHit.nGammas()==0 && lHit.nElectrons()==0 ) continue;
	    //HACK !!
	    if (lHit.nNeutrons()==0 ) continue;
	  }
	    
	  double posx = useOldEncoding? lHit.get_x_old():lHit.get_x();
	  double posy = useOldEncoding? lHit.get_y_old():lHit.get_y();
	  double posz = lHit.get_z();
	  
	  double radius = sqrt(posx*posx+posy*posy);
	  if (radius < minRadius) continue;
	  //double posx = lHit.get_y();
	  //double posy = lHit.get_x();
	  double energy = lHit.energy();
	  if (debug>1) {
	    std::cout << " --  SimHit " << iH << "/" << (*pusimhitvec).size() << " --" << std::endl
		      << " --  position x,y " << posx << "," << posy << std::endl;
	    lHit.Print(std::cout);
	  }
	  double weightedE = energy;
	  if (applyWeight) weightedE *= mycalib.MeVToMip(layer,true);
	  else weightedE *= mycalib.MeVToMip(layer);
	  //reject lower region for boundary effects
	  if (doFiducialCuts){
	    if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
		 (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	    //region left and right also for boundary effects
	    if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
		 (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	  }
	  p_xy[iE][layer]->Fill(posx,posy,energy);
	  p_xyz[iE]->Fill(posz,posx,posy,energy);
	    
	  if (saveEventByEvent){
	    bool pass = (layer < nEcalLayers) ||
	      (layer>= nEcalLayers && 
	       layer<nEcalLayers+nFHcalLayers &&
	       (!concept || (concept && layer%2==1))
	       ) ||
	      layer>=nEcalLayers+nFHcalLayers;
	    if (pass){
	      double tmpE = weightedE;
	      p_xy_evt[layer]->Fill(posx,posy,tmpE);
	      p_xyz_evt->Fill(posz,posx,posy,tmpE);
	    }
	  }
	  //correct for time of flight
	  double lRealTime = mycalib.correctTime(lHit.time(),posx,posy,posz);
	  p_time[iE]->Fill(lRealTime);
	  //restrict in y and x to have just the influence of the eta bin wanted
	  //(at high eta, surface of the cone in x is smaller than detector size)
	  if (fabs(posx)<signalRegionInX/2.){
	    double lEtaTmp = getEta(posy,pEta,posz);
	    p_Edensity[iE][layer]->Fill(lEtaTmp,energy);
	  }
	  Etot[layer] += energy;
	  if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << Etot[layer];
	    
	  mycalib.incrementEnergy(layer,weightedE,lHit.time(),posx,posy);
	  
	  for (unsigned iR(0); iR<nLat; ++iR){
	    if (radius<latSize[iR]) mycalib.incrementBlockEnergy(layer,weightedE,Elateral[0][iR],Elateral[1][iR],Elateral[2][iR]);
	  }

	  //calculate mean Emip
	  EmipHits->Fill(posx,posy,weightedE);
	}//loop on hits
	  
	  
	nTotal += (*pusimhitvec).size();
	if (debug)  std::cout << std::endl;
	if (!isG4Tree){
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){	  
	    Ereco[iSmear] = 0;
	  }
	  for (unsigned iH(0); iH<(*purechitvec).size(); ++iH){//loop on rechits
	    HGCSSRecoHit lHit = (*purechitvec)[iH];
	    if (debug>1) {
	      std::cout << " --  RecoHit " << iH << "/" << (*purechitvec).size() << " --" << std::endl
			<< " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
	      lHit.Print(std::cout);
	    }
	    double posx = lHit.get_x();
	    double posy = lHit.get_y();
	    double posz = lHit.get_z();
	    if (doFiducialCuts){
	      //reject lower region for boundary effects
	      if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
		   (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	      //region left and right also for boundary effects
	      if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
		   (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	    }
	      
	    double energy = lHit.energy();//in MIPs already
	    unsigned layer = lHit.layer();
	    p_recoxy[iE][layer]->Fill(posx,posy,energy);
	    p_recoxyz[iE]->Fill(posz,posx,posy,energy);
	    for (unsigned iO(0); iO<nOcc; ++iO){
	      if (energy > occThreshold[iO]) p_occupancy[iO][layer]->Fill(posx,posy);
	    }
	    double weightedE = energy;
	    //if (applyWeight)  weightedE *= mycalib.MeVToMip(layer);
	    for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	      double smearFactor = lRndm.Gaus(0,smearFrac[iSmear]);
	      Ereco[iSmear] += weightedE*(1+smearFactor);
	    }
	      
	    if (saveEventByEvent){
	      p_recoxy_evt[layer]->Fill(posx,posy,energy);
	      p_recoxyz_evt->Fill(posz,posx,posy,energy);
	    }
	      
	      
	  }
	}
      }//add PU hits

      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	HGCSSSimHit lHit = (*simhitvec)[iH];
	unsigned layer = lHit.layer();

	//select "flavour"
	if (doProtons){
	  //if ((lHit.nHadrons() == 0 && lHit.nMuons() == 0)) continue;
	  //HACK !!
	  if (lHit.nProtons()==0 ) continue;
	}
	else if (doNeutrons){
	  //HACK !!
	  if (lHit.nNeutrons()==0 ) continue;
	  //if (lHit.nGammas()==0 && lHit.nElectrons()==0 ) continue;
	}

	double posx = useOldEncoding? lHit.get_x_old():lHit.get_x();
	double posy = useOldEncoding? lHit.get_y_old():lHit.get_y();
	double posz = lHit.get_z();
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!! WARNING, quick and dirty fix !!!!!!!!!!!!!!!!!!
	//double posx = lHit.get_y();
	//double posy = lHit.get_x();
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	double radius = sqrt(posx*posx+posy*posy);
	if (radius < minRadius) continue;
	
	double energy = lHit.energy();
	if (debug>1) {
	  std::cout << " --  SimHit " << iH << "/" << (*simhitvec).size() << " --" << std::endl
		    << " --  position x,y " << posx << "," << posy << std::endl;
	  lHit.Print(std::cout);
	}
	if (layer<10) E1 += energy;
	else if (layer<20) E2 += energy;
	else E3 += energy;

	if (energy>0 && saveTxtOutput) {
	  if (isG4Tree) hitsOut << static_cast<unsigned>(ievt/nLayers);
	  else hitsOut << ievt;
	  hitsOut << " " << layer << " " << posx << " " << posy << " " << energy << std::endl;
	}

	double weightedE = energy;
	if (applyWeight) weightedE *= mycalib.MeVToMip(layer,true);
	else weightedE *= mycalib.MeVToMip(layer);

	if (doFiducialCuts){
	  //reject lower region for boundary effects
	  if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
	       (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	  //region left and right also for boundary effects
	  if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
	       (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	}

	p_xy[iE][layer]->Fill(posx,posy,energy);
	p_xyz[iE]->Fill(posz,posx,posy,energy);
	if (saveEventByEvent){
	  bool pass = (layer < nEcalLayers) ||
	    (layer>= nEcalLayers && 
	     layer<nEcalLayers+nFHcalLayers &&
	     (!concept || (concept && layer%2==1))
	     ) ||
	    layer>=nEcalLayers+nFHcalLayers;
	  if (pass){
	    double tmpE = weightedE;
	    p_xy_evt[layer]->Fill(posx,posy,tmpE);
	    p_xyz_evt->Fill(posz,posx,posy,tmpE);
	  }
	}
	//correct for time of flight
	double lRealTime = mycalib.correctTime(lHit.time(),posx,posy,posz);
	p_time[iE]->Fill(lRealTime);

	if (fabs(posx)<signalRegionInX/2.){
	  double lEtaTmp = getEta(posy,pEta,posz);
	  p_Edensity[iE][layer]->Fill(lEtaTmp,energy);
	  if (lEtaTmp>=etaMinLLR && lEtaTmp<etaMaxLLR) eLLR[layer]+=energy;
	}
	Etot[layer] += energy;
	if (debug>1) std::cout << "-hit" << iH << "-" << layer << " " << energy << " " << Etot[layer];
	    
	mycalib.incrementEnergy(layer,weightedE,lHit.time(),posx,posy);
	
	for (unsigned iR(0); iR<nLat; ++iR){
	  if (radius<latSize[iR]) mycalib.incrementBlockEnergy(layer,weightedE,Elateral[0][iR],Elateral[1][iR],Elateral[2][iR]);
	}
	
	//calculate mean Emip
	EmipHits->Fill(posx,posy,weightedE);
 	zlay[layer] = posz;

      }//loop on hits

      nTotal += (*simhitvec).size();
      nTotalSignal += (*simhitvec).size();

      if (debug)  std::cout << std::endl;
      if (!isG4Tree){
	if (!overlayPU) {
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){	  
	    Ereco[iSmear] = 0;
	  }
	}
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
	  HGCSSRecoHit lHit = (*rechitvec)[iH];
	  if (debug>1) {
 	    std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		      << " --  position x,y " << lHit.get_x() << "," << lHit.get_y() << std::endl;
	    lHit.Print(std::cout);
	  }
	    
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  double posz = lHit.get_z();

	  if (doFiducialCuts){
	    //reject lower region for boundary effects
	    if ( (pVersion.find("23") != pVersion.npos && posy < -200) ||
		 (pVersion.find("23") == pVersion.npos && posy < -70)) continue;
	    //region left and right also for boundary effects
	    if ( (pVersion.find("23") != pVersion.npos && fabs(posx) > 200) ||
		 (pVersion.find("23") == pVersion.npos && fabs(posx) > 60)) continue;
	  }

	  double energy = lHit.energy();//in MIP already...
	  unsigned layer = lHit.layer();
	  for (unsigned iO(0); iO<nOcc; ++iO){
	    if (energy > occThreshold[iO]) p_occupancy[iO][layer]->Fill(posx,posy);
	  }
	  double weightedE = energy;
	  //if (applyWeight)  weightedE *= mycalib.MeVToMip(layer);
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	    double smearFactor = lRndm.Gaus(0,smearFrac[iSmear]);
	    Ereco[iSmear] += weightedE*(1+smearFactor);
	  }

	  p_recoxy[iE][layer]->Fill(posx,posy,energy);
	  p_recoxyz[iE]->Fill(posz,posx,posy,energy);
	  if (saveEventByEvent){
	    p_recoxy_evt[layer]->Fill(posx,posy,energy);
	    p_recoxyz_evt->Fill(posz,posx,posy,energy);
	  }
	}
	if (firstEvent){
	  for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	    lName.str("");
	    lName << "p_Ereco_" << genEn[iE] << "_smear" << iSmear;
	    p_Ereco[iE][iSmear] = new TH1F(lName.str().c_str(),";Reco Etotal (MIPs)",1000,0,2*genEn[iE]);
	  }
	  lName.str("");
	  lName << "p_nRecHits_" << genEn[iE];
	  unsigned nRecHits = (*rechitvec).size();
	  if (overlayPU) nRecHits += (*purechitvec).size();
	  p_nRecHits[iE] = new TH1F(lName.str().c_str(),"; nRecHits",1000,static_cast<unsigned>(nRecHits/10.),nRecHits*6);
	}
	for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	  p_Ereco[iE][iSmear]->Fill(Ereco[iSmear]);
	  Ereco[iSmear] = 0;
	}
	p_nRecHits[iE]->Fill((*rechitvec).size());
	if (debug) std::cout << "... recoE = " << Ereco << std::endl;
      }


      //fill histos
      //if (debug) std::cout << " -- Filling histograms..." << std::endl;
      if (applyDigi){
	myDigitiser.digitiseECAL(mycalib.getECAL2DHistoVec());
	myDigitiser.digitiseFHCAL(mycalib.getFHCAL2DHistoVec());
	myDigitiser.digitiseBHCAL(mycalib.getBHCAL2DHistoVec());
      }
      
      Etotal[0] = mycalib.getECALenergy(mipThresh);
      Etotal[1] = mycalib.getFHCALenergy(mipThresh);
      Etotal[2] = mycalib.getBHCALenergy(mipThresh);
      Etotal[3] = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2]);
      
      if (addECALenergyLoss) {
	//remove ECAL energy
	mimicECALenergyLoss(Etotal[3],lRndm,Etotal[0]);
	Etotal[3] = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2]);
      }
      
      double EReco_FHCAL = mycalib.recoEnergyUncor(Etotal[0],Etotal[1],Etotal[2],true);
      //bool doFill = fillHistoConditions(selectEarlyDecay,EReco_FHCAL,Etotal[3]);
      double Etmp = 0;
      for (unsigned iL(0);iL<6;++iL){//loop on layers
	Etmp += Etot[iL]*mycalib.MeVToMip(iL,true);
      }
      bool doFill = true;
      if (selectEarlyDecay && Etotal[0]+Etotal[1]+Etotal[2] > 0 && Etmp/(Etotal[0]+Etotal[1]+Etotal[2])<0.08) doFill = false;
      if (selectPiHCAL && Etotal[0]>10) doFill = false;
      
      Etmp = 0;
      
      if (saveTxtOutput) optim << E1 << " " << E2 << " " << E3 << std::endl;
      E1=0;
      E2=0;
      E3=0;
      
      if (firstEvent){
	for (unsigned iL(0);iL<nLayers;++iL){
	  if (Etot[iL]>Emax) Emax = Etot[iL];
	}
	for (unsigned iL(0);iL<nLayers;++iL){
	  lName.str("");
	  lName << "p_Etot_" << genEn[iE] << "_" << iL;
	  p_Etot[iE][iL] = new TH1F(lName.str().c_str(),";G4 Etot (MeV)",1000,0,genEn[iE]);
	}
	lName.str("");
	lName << "p_EvsLayer_" << genEn[iE];
	p_EvsLayer[iE] = new TH2F(lName.str().c_str(),";layer;E (MeV)",nLayers,0,nLayers,1000,0,genEn[iE]);
	
	lName.str("");
	lName << "p_Etotal_" << genEn[iE] << "_ECAL";
	p_Etotal[iE][0] = new TH1F(lName.str().c_str(),";G4 E ECAL (MeV)",
				   1000,0,2*genEn[iE]);
	p_Etotal[iE][0]->StatOverflows();
	
	lName.str("");
	lName << "p_Etotal_" << genEn[iE] << "_FHCAL";
	p_Etotal[iE][1] = new TH1F(lName.str().c_str(),";G4 E FHCAL (MeV)",
				   1000,0,2*genEn[iE]);
	
	p_Etotal[iE][1]->StatOverflows();
	
	lName.str("");
	lName << "p_Etotal_" << genEn[iE] << "_BHCAL";
	p_Etotal[iE][2] = new TH1F(lName.str().c_str(),";G4 E BHCAL (MeV)",
				   1000,0,2*genEn[iE]);
	
	p_Etotal[iE][2]->StatOverflows();
	
	lName.str("");
	lName << "p_Etotal_" << genEn[iE] << "_ECALHCAL";
	p_Etotal[iE][3] = new TH1F(lName.str().c_str(),";G4 Etotal ECAL+HCAL (MeV)",1000,0,2*genEn[iE]);
	p_Etotal[iE][3]->StatOverflows();
	
	lName.str("");
	lName << "p_HCALvsECAL_" << genEn[iE];
	p_HCALvsECAL[iE] = new TH2F(lName.str().c_str(),";G4 E ECAL; G4 E HCAL",
				    1000,0,2*genEn[iE],
				    1000,0,2*genEn[iE]);
	p_HCALvsECAL[iE]->StatOverflows();
	
	lName.str("");
	lName << "p_BHCALvsFHCAL_" << genEn[iE];
	p_BHCALvsFHCAL[iE] = new TH2F(lName.str().c_str(),";G4 E FHCAL; G4 E BHCAL",
				      1000,0,2*genEn[iE],
				      1000,0,2*genEn[iE]);
	p_HCALvsECAL[iE]->StatOverflows();
	
	lName.str("");
	lName << "p_nSimHits_" << genEn[iE];
	p_nSimHits[iE] = new TH1F(lName.str().c_str(),"; nSimHits",1000,static_cast<unsigned>(nTotal/10.),nTotal*10);
      }
      
      for (unsigned iL(0);iL<nLayers;++iL){//loop on layers
	if (debug) std::cout << " -- Layer " << iL << " total E = " << Etot[iL] << std::endl;
	if (doFill) {
	  p_Etot[iE][iL]->Fill(Etot[iL]*mycalib.MeVToMip(iL));
	  p_EvsLayer[iE]->Fill(iL,Etot[iL]*mycalib.MeVToMip(iL));
	  Etmp += Etot[iL]*mycalib.MeVToMip(iL,true);
	  if (Etotal[0]+Etotal[1]+Etotal[2] > 0) p_Efrac[iE][iL]->Fill(Etmp/(Etotal[0]+Etotal[1]+Etotal[2]));
	  else p_Efrac[iE][iL]->Fill(0);
	}
	Etot[iL] = 0;
	
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	//For LLR study
	//get total density in one eta bin
	if (saveLLR){// && iL>0 && iL<30){
	  double llrBinWidth = getBinWidth(etaMinLLR,etaMaxLLR,"eta35",zlay[iL]);
	  double density = eLLR[iL]/signalRegionInX*1/llrBinWidth;
	  if (isG4Tree) edensity << static_cast<unsigned>(ievt/nLayers);
	  else edensity << ievt;
	  edensity << " " << iL << " " << density << std::endl;
	  eLLR[iL] = 0;
	}
	//end-For LLR study
	//////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
	
      }//loop on layers
      if (doFill) {
	
	//fill Cglobal histos
	unsigned nHitsCountNum[nLimits];
	unsigned nHitsCountDen = 0;
	
	double Cglobal[nLimits];
	for (unsigned iLim(0); iLim<nLimits;++iLim){
	  nHitsCountNum[iLim] = 0;
	  Cglobal[iLim] = 1.;
	  for (int ix(1); ix<EmipHits->GetNbinsX()+1; ++ix){
	    for (int iy(1); iy<EmipHits->GetNbinsY()+1; ++iy){
	      double eTmp = EmipHits->GetBinContent(ix,iy);
	      if (eTmp<pElim[iLim]) nHitsCountNum[iLim]++;
	      if (iLim==0 && eTmp<EmipMean[iE]) nHitsCountDen++;
	    }
	  }
	  if (nHitsCountDen>0) Cglobal[iLim] = nHitsCountNum[iLim]*1.0/nHitsCountDen;
	  p_Cglobal[iE][iLim]->Fill(Cglobal[iLim]);
	  p_EvsCglobal[iE][iLim]->Fill(Cglobal[iLim],EReco_FHCAL);
	  p_Eshower[iE][iLim]->Fill(mycalib.hcalShowerEnergy(Cglobal[0],Etotal[1],Etotal[2]));
	}
	
	p_HCALvsECAL[iE]->Fill(Etotal[0]*(calibrate?mycalib.gevWeight(0):1),
			       mycalib.hcalShowerEnergy(Cglobal[0],Etotal[1],Etotal[2]));
	
	
	p_BHCALvsFHCAL[iE]->Fill(Cglobal[0]*Etotal[1]*(calibrate?mycalib.gevWeight(1):1),Etotal[2]*(calibrate?mycalib.gevWeight(2):1));
	
	// std::cout << lEvt << " " 
	// 	    << nTotHits << " "
	// 	    << EmipMean << " " 
	// 	    << pElim << " "
	// 	    << nHitsCount[0] << " "
	// 	    << nHitsCount[1] << " "
	// 	    << Cglobal << std::endl;
	
	for (unsigned iD(0); iD<nSections; ++iD){
	  p_Etotal[iE][iD]->Fill(Etotal[iD]*(calibrate?mycalib.gevWeight(iD):1.));
	  for (unsigned iR(0); iR<nLat; ++iR){
	    if (iD==nSections-1){
	      double totTmp = mycalib.showerEnergy(Cglobal[0],Etotal[0],Etotal[1],Etotal[2]);
	      double latTmp = mycalib.showerEnergy(Cglobal[0],Elateral[0][iR],Elateral[1][iR],Elateral[2][iR]);
	      if (totTmp > 0) p_EfracLateral[iE][iD][iR]->Fill(latTmp/totTmp);
	    }
	    else if (Etotal[iD]>0) p_EfracLateral[iE][iD][iR]->Fill(Elateral[iD][iR]/Etotal[iD]);
	  }
	}
	
	p_nSimHits[iE]->Fill(nTotal);
	
	if (saveEventByEvent && nTotalSignal >= nMaxHits){
	  std::ostringstream evtName;
	  evtName << plotDir << "/CalibHistos_E" << genEn[iE] << "_evt" << lEvt << ".root";
	  TFile *outputEvt = TFile::Open(evtName.str().c_str(),"RECREATE");
	  
	  if (!outputEvt) {
	    std::cout << " -- Error, output file for evt " << lEvt << " cannot be opened. Please create output directory : " << plotDir << ". Exiting..." << std::endl;
	    return 1;
	  }
	  else {
	    std::cout << " -- Opening event file for evt " << lEvt << std::endl;
	  }
	  outputEvt->cd();
	  for (unsigned iL(0); iL<nLayers; ++iL){
	    p_xy_evt[iL]->Write();
	    if (!isG4Tree) p_recoxy_evt[iL]->Write();
	  }
	  if (p_xyz_evt) p_xyz_evt->Write();
	  if (!isG4Tree) p_recoxyz_evt->Write();
	  outputEvt->Close();
	}
      }//if dofill

      for (unsigned iD(0); iD<nSections; ++iD){
	Etotal[iD] = 0;
	for (unsigned iR(0); iR<nLat; ++iR){
	  Elateral[iD][iR] = 0;
	}
      }
      
      EmipHits->Reset();
      mycalib.reset2DHistos();
      nTotal = 0;
      nTotalSignal = 0;
      firstEvent = false;
      
    }//loop on entries
  
    if (saveTxtOutput) {
      hitsOut.close();
      optim.close();
    }
    if (saveLLR) edensity.close();
    
    for (unsigned iL(0); iL<nLayers; ++iL){
      outputFile->cd();
      p_xy[iE][iL]->Write();
      if (!isG4Tree) p_recoxy[iE][iL]->Write();
      if (!lightOutput) p_Edensity[iE][iL]->Write();
      if (!lightOutput) p_Etot[iE][iL]->Write();
      p_Efrac[iE][iL]->Write();
    }
    p_time[iE]->Write();
    if (!lightOutput) p_xyz[iE]->Write();
    p_EvsLayer[iE]->Write();
    if (!lightOutput) {
      for (unsigned iLim(0); iLim<nLimits;++iLim){
	p_Cglobal[iE][iLim]->Write();
	p_EvsCglobal[iE][iLim]->Write();
	p_Eshower[iE][iLim]->Write();
      }
    }
    for (unsigned iD(0); iD<nSections; ++iD){
      p_Etotal[iE][iD]->Write();
      if (!lightOutput) {
	for (unsigned iR(0); iR<nLat; ++iR){
	  p_EfracLateral[iE][iD][iR]->Write();
	}
      }
    }
    p_HCALvsECAL[iE]->Write();
    p_BHCALvsFHCAL[iE]->Write();
    p_nSimHits[iE]->Write();
    if (!isG4Tree) {
      if (!lightOutput){
	for (unsigned iSmear(0); iSmear<nSmear;++iSmear){
	  p_Ereco[iE][iSmear]->Write();
	}
      }
      p_nRecHits[iE]->Write();
      if (!lightOutput) p_recoxyz[iE]->Write();
    }

    std::cout << " -- Summary of energies ECAL: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE][0]->GetEntries() 
	      << " mean " << p_Etotal[iE][0]->GetMean() 
	      << " rms " << p_Etotal[iE][0]->GetRMS() 
	      << " underflows " << p_Etotal[iE][0]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE][0]->GetBinContent(p_Etotal[iE][0]->GetNbinsX()+1)
	      << std::endl;

    std::cout << " -- Summary of energies HCAL section 1: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE][1]->GetEntries() 
	      << " mean " << p_Etotal[iE][1]->GetMean() 
	      << " rms " << p_Etotal[iE][1]->GetRMS() 
	      << " underflows " << p_Etotal[iE][1]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE][1]->GetBinContent(p_Etotal[iE][1]->GetNbinsX()+1)
	      << std::endl;

    std::cout << " -- Summary of energies HCAL section 2: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE][2]->GetEntries() 
	      << " mean " << p_Etotal[iE][2]->GetMean() 
	      << " rms " << p_Etotal[iE][2]->GetRMS() 
	      << " underflows " << p_Etotal[iE][2]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE][2]->GetBinContent(p_Etotal[iE][2]->GetNbinsX()+1)
	      << std::endl;
 
    std::cout << " -- Summary of energies ECAL+HCAL: " << std::endl
	      << " ---- SimHits: entries " << p_Etotal[iE][3]->GetEntries() 
	      << " mean " << p_Etotal[iE][3]->GetMean() 
	      << " rms " << p_Etotal[iE][3]->GetRMS() 
	      << " underflows " << p_Etotal[iE][3]->GetBinContent(0)
	      << " overflows " << p_Etotal[iE][3]->GetBinContent(p_Etotal[iE][3]->GetNbinsX()+1)
	      << std::endl;
    if (!isG4Tree) {
      std::cout << " ---- RecHits: entries " << p_Ereco[iE][0]->GetEntries() 
		<< " mean " << p_Ereco[iE][0]->GetMean() 
		<< " rms " << p_Ereco[iE][0]->GetRMS() 
		<< " underflows " << p_Ereco[iE][0]->GetBinContent(0)
		<< " overflows " << p_Ereco[iE][0]->GetBinContent(p_Ereco[iE][0]->GetNbinsX()+1)
		<< std::endl;
    }
    std::cout << " -- Summary of hits: " << std::endl
	      << " ---- SimHits: entries " << p_nSimHits[iE]->GetEntries() 
	      << " mean " << p_nSimHits[iE]->GetMean() 
	      << " rms " << p_nSimHits[iE]->GetRMS() 
	      << " underflows " << p_nSimHits[iE]->GetBinContent(0)
	      << " overflows " << p_nSimHits[iE]->GetBinContent(p_nSimHits[iE]->GetNbinsX()+1)
	      << std::endl;
    if (!isG4Tree) {
      std::cout << " ---- RecHits: entries " << p_nRecHits[iE]->GetEntries() 
		<< " mean " << p_nRecHits[iE]->GetMean() 
		<< " rms " << p_nRecHits[iE]->GetRMS() 
		<< " underflows " << p_nRecHits[iE]->GetBinContent(0)
		<< " overflows " << p_nRecHits[iE]->GetBinContent(p_nRecHits[iE]->GetNbinsX()+1)
		<< std::endl;
    }

    // //create summary control plots
    // mycE->cd(iE+1);
    // p_Etotal[iE]->Rebin(10);
    // p_Etotal[iE]->Draw();
    // TLatex lat;
    // char buf[500];
    // sprintf(buf,"Egen = %d GeV",genEn[iE]);
    // lat.DrawLatex(p_Etotal[iE]->GetMean(),p_Etotal[iE]->GetMaximum()*1.1,buf);

  }//loop on energies

 if (!isG4Tree && !lightOutput) {
   outputFile->cd();
   for (unsigned iO(0); iO<nOcc; ++iO){
     for (unsigned iL(0); iL<nLayers; ++iL){
       p_occupancy[iO][iL]->Write();
     }
   }
 }

  std::cout << "\\hline \n";
  std::cout << " Energy (GeV) & Eecal (RMS) & Ehcal1 (RMS) & Ehcal2 (RMS) & Etot (RMS) "
    // << " & nSimHits (RMS) "
	    << "\\\\ \n"; 
  std::cout << "\\hline \n";
  for (unsigned iE(0); iE<nGenEn; ++iE){//loop on energies
    if (!p_Etotal[iE][0] || !p_Etotal[iE][1] || !p_Etotal[iE][2]) continue;
    std::cout << genEn[iE] << " & "
	      << p_Etotal[iE][0]->GetMean() << " (" << p_Etotal[iE][0]->GetRMS() 
	      << ") & "
	      << p_Etotal[iE][1]->GetMean() << " (" << p_Etotal[iE][1]->GetRMS() 
	      << ") & "
	      << p_Etotal[iE][2]->GetMean() << " (" << p_Etotal[iE][2]->GetRMS() 
	      << ") & "
	      << p_Etotal[iE][3]->GetMean() << " (" << p_Etotal[iE][3]->GetRMS() 
        //<< ") & "
      //<< p_nSimHits[iE]->GetMean() << " (" << p_nSimHits[iE]->GetRMS() 
	      << ") \\\\ \n ";
  }
  std::cout << "\\hline \n";



  // std::ostringstream saveName;
  // saveName.str("");
  // saveName << plotDir << "/ControlG4Etotal";
  // mycE->Update();
  // mycE->Print((saveName.str()+".png").c_str());
  // mycE->Print((saveName.str()+".pdf").c_str());

  // for (unsigned iE(0); iE<nGenEn; ++iE){//loop on energies
  //   p_nSimHits[iE]->Draw();
  //   p_nSimHits[iE]->Rebin(10);
  //   TLatex lat;
  //   char buf[500];
  //   sprintf(buf,"Egen = %d GeV",genEn[iE]);
  //   lat.DrawLatex(p_nSimHits[iE]->GetMean(),p_nSimHits[iE]->GetMaximum()*1.1,buf);
    
  // }//loop on energies

  // saveName.str("");
  // saveName << plotDir << "/ControlG4nSimHits";
  // mycE->Update();
  // mycE->Print((saveName.str()+".png").c_str());
  // mycE->Print((saveName.str()+".pdf").c_str());

  
  outputFile->Write();

  return 0;


}//main
