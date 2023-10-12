#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<numeric>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

//ORDER MATTERS...
#include "TDRStyle.h"
#include "plotUtilities.C"
#include "fitEnergy.C"
#include "makeBackLeakCor.C"
#include "makePuSigma.C"
#include "makeCalibration.C"
#include "makeResolution.C"

#include "../include/InputParserBase.h"

class InputParserPlotEGReso: public InputParserBase {
public:
  explicit InputParserPlotEGReso(int argc, char** argv): InputParserBase(argc, argv) { run(); }

  std::string tag() const { return tag_; }
  std::vector<float> etas() const { return etas_; }
  std::vector<std::string> versions() const { return versions_; }
  std::vector<unsigned> thicknesses() const { return thicknesses_; }
  std::vector<unsigned> signalRegions() const { return signalRegions_; }
  unsigned nBack() const { return nBack_; }
  bool useSigmaEff() const { return useSigmaEff_; }

private:
  std::string tag_;
  std::vector<float> etas_;
  std::vector<unsigned> thicknesses_;
  std::vector<unsigned> signalRegions_;
  std::vector<std::string> versions_;
  unsigned nBack_;
  bool useSigmaEff_;

  void set_args_final_()
  {
    tag_ = chosen_args_["--tag"];
    nBack_ = static_cast<unsigned>( std::stoi(chosen_args_["--nBack"]) );
    useSigmaEff_ = chosen_args_["--useSigmaEff"]=="true" ? true : false;
    for(auto&& x: chosen_args_v_["--etas"])
      etas_.push_back( std::stof(x) );
    for(auto&& x: chosen_args_v_["--thicknesses"])
      thicknesses_.push_back( static_cast<unsigned>(std::stoi(x)) );
    for(auto&& x: chosen_args_v_["--signalRegions"])
      signalRegions_.push_back( static_cast<unsigned>(std::stoi(x)) );
    versions_ = chosen_args_v_["--versions"];
  }
  
  void set_args_options_()
  {
    required_args_ = { "--nBack", "--thicknesses", "--signalRegions",
		       "--versions", "--tag", "--etas" };
    optional_args_ = {"--useSigmaEff"};

    valid_args_["--nBack"] = {"0", "1", "2", "3", "4"};
    valid_args_v_["--thicknesses"] = {"100", "200", "300"};
    valid_args_v_["--signalRegions"] = {"0", "1", "2", "3", "4", "5"};
    valid_args_v_["--versions"] = {"60", "70", "80"};
    free_args_ = {"--tag"};
    free_args_v_ = {"--etas"};
  }
};

int makeEfit(const bool useSigmaEff,
	     const bool dovsE,
	     const bool doRaw,
	     const bool doBackLeakCor,
	     const unsigned nBack,
	     const unsigned pu,
	     const unsigned ICval,
	     const std::string scenario,
	     const TString version,
	     const unsigned nLayers,
	     const unsigned eta,
	     const unsigned pT,
	     const unsigned iSR,
	     const double radius,
	     const double & offset,
	     const double & calib,
	     const double & backLeakCor,
	     TCanvas * mycE,
	     TTree *ltree,
	     TFile *foutEfit,
	     TString plotDir,
	     TGraphErrors *calibRecoFit,
	     TGraphErrors *resoRecoFit
	     ){//main


  double Eval = E(pT,eta);
  double etaval = eta/10.;
  double etaerr = 0;

  const unsigned nEvtMin = 150;

  TString pSuffix = doBackLeakCor?"backLeakCor":"";


  double fitQual = iSR==0?50:30;

  
  TH1F *p_chi2ndf;
  std::ostringstream label;
  label << "p_chi2ndf_E_" << pT << "_SR" << iSR;
  p_chi2ndf = new TH1F(label.str().c_str(),";#chi^{2}/N;entries",500,0,50);
  p_chi2ndf->StatOverflows();
 
  
  std::string unit = "MIPS";
  if (!doRaw) unit = "GeV";
  
  //identify valid energy values
  
  TH1F *p_Ereco;
  TH2F *p_ErecovsEback;
	      
  mycE->cd();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
		
  std::ostringstream lName;
  std::string lNameTot,lNameBack;

  getTotalEnergyString(nLayers,nBack,lNameTot,lNameBack,iSR);

  //std::cout << " Check totE = " << lNameTot << std::endl;
  //std::cout << " Check backE = " << lNameBack << std::endl;

  lName.str("");
  lName << std::setprecision(6);
  
  if (!doRaw) lName << "(";
  lName << lNameTot;
  if (!doRaw) lName << " - " << offset << ")/" << calib;

  if (doBackLeakCor){// && E>300) {
    lName << " - " << backLeakCor << "*(" << lNameBack << ")/(" << lNameTot << ")";
  }
  
  //low energy cut
  std::ostringstream ecut;
  ecut << lName.str() << ">0.01";
  
  //cut low outliers....
  if (iSR>0 && !doRaw){
    ecut << " && " << lName.str() << ">0.7*" << Eval;
  }

  ltree->Draw(lName.str().c_str(),ecut.str().c_str(),"");
  TH1F  *hTmp =  (TH1F*)gPad->GetPrimitive("htemp");
  if (!hTmp) {
    std::cout << " -- Problem with tree Draw for "
	      <<  lName.str() << " cut " << ecut.str()
	      << " -> not drawn. Exiting..." << std::endl;
    return 1;
  }
  
  //get mean and rms, and choose binning accordingly then redraw.
  unsigned nrms = 5;
  TH1F *hE = new TH1F("hE","",50,hTmp->GetMean()-nrms*hTmp->GetRMS(),hTmp->GetMean()+nrms*hTmp->GetRMS());
  ltree->Draw((lName.str()+">>hE").c_str(),ecut.str().c_str(),"");

  lName.str("");
  lName << "energy" << pT << "_SR" << iSR ;
  p_Ereco = (TH1F*)hE->Clone(lName.str().c_str()); // 1D
  
  if (!p_Ereco){
    std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
    return 1;
  }
  std::cout << " --- Reco E = entries " << p_Ereco->GetEntries() 
	    << " mean " << p_Ereco->GetMean() 
	    << " rms " << p_Ereco->GetRMS()
    	    << " nbinsx " << p_Ereco->GetNbinsX()
	    << " overflows " << p_Ereco->GetBinContent(p_Ereco->GetNbinsX()+1)
	    << std::endl;
  
  
  p_Ereco->SetTitle((";E ("+unit+");events").c_str());
  
  //take min 20 bins
  //if(p_Ereco->GetNbinsX()>40) p_Ereco->Rebin(2);
  
  //skip data with too little stat
  if (p_Ereco->GetEntries()<nEvtMin) {
    gPad->Clear();
    return 1;
    
  }
  
  
  TPad *lpad = (TPad*)(mycE->cd());
  FitResult lres;
  if (fitEnergy(p_Ereco,lpad,unit,lres,iSR,useSigmaEff)!=0) return 1;
  lpad->cd();
  char buf[500];
  sprintf(buf,"#gamma p_{T}=%d GeV + PU %d",pT,pu);
  TLatex lat;
  lat.SetTextSize(0.04);
  lat.DrawLatexNDC(0.2,0.96,buf);
  sprintf(buf,"#eta=%3.1f, r = %3.0f mm",etaval,radius);
  lat.SetTextSize(0.04);
  lat.DrawLatexNDC(0.2,0.85,buf);
  lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
  p_chi2ndf->Fill(lres.chi2/lres.ndf);
  
  //filter out bad points
  if (lres.chi2/lres.ndf > fitQual) {
    std::cout << " --- INFO! Point Egen=" 
	      << pT 
	      << " eta=" << etaval
	      << " pu=" << pu
	      << " skipped, chi2/ndf = "
	      << lres.chi2/lres.ndf
	      << std::endl;
    return 1;
  }

  Int_t np=calibRecoFit->GetN();
  //if (!dovsE) calibRecoFit[iSR]->SetPoint(np,genEn[iE],lres.mean);
  //else 
  calibRecoFit->SetPoint(np,Eval,lres.mean);
  calibRecoFit->SetPointError(np,0.0,lres.meanerr);
  
  //use truth: residual calib won't affect resolution....
  double reso = fabs(lres.sigma/lres.mean);
  //if (pu>0) reso = fabs(lres.sigma/(dovsE?E(genEn[iE],eta[ieta]):genEn[iE]));
  if (pu>0) reso = fabs(lres.sigma/Eval);

  np=resoRecoFit->GetN();
  if (!dovsE) resoRecoFit->SetPoint(np,pT,reso);
  else resoRecoFit->SetPoint(np,Eval,reso);
  double errFit = reso*sqrt(pow(lres.sigmaerr/lres.sigma,2)+pow(lres.meanerr/lres.mean,2));
  resoRecoFit->SetPointError(np,0,errFit);
  
  std::cout << " -- calib/reso value for point " << np << " " << Eval << " " << lres.mean << " +/- " << lres.meanerr 
	    << " sigma " << lres.sigma << " +/- " << lres.sigmaerr
	    << " reso " << reso << " +/- " << errFit
	    << std::endl;


  if (system(TString("mkdir -p ")+plotDir+TString("/Ereco"))) return 1;
  std::ostringstream saveName;
  saveName.str("");
  saveName << plotDir << "/Ereco/Ereco_eta" << eta << "_pu" << pu;
  if (doRaw) saveName << "raw";
  saveName << "_E" << pT << "_SR" << iSR << pSuffix;
  mycE->Update();
  mycE->Print((saveName.str()+".png").c_str());
  mycE->Print((saveName.str()+".C").c_str());

  //    drawChi2(p_chi2ndf);
	    
  foutEfit->cd();
  p_chi2ndf->Write();
  p_Ereco->Write();
  
  return 0;
  
};//makeEfit

int plotEGReso(const InputParserPlotEGReso& ip) {
 
  SetTdrStyle();
  
  const std::array<unsigned,1> ICvals = {3}; //0,1,2,3,4,5,10,15,20,50};

  const bool dovsE = true;
  const bool doBackLeakCor = true;
  const bool redoCalib = true;
  const bool redoLeakCor = true;

  const std::array<unsigned,1> pu = {0};//,140,200};
  const unsigned nPu = pu.size();

  const std::vector<unsigned> thicknesses = ip.thicknesses();
  const unsigned nS = thicknesses.size();
  std::vector<std::string> scenarios(nS);
  for(unsigned inS(0); inS<nS; ++inS)
    scenarios[inS] = std::string("model2/gamma/") + std::to_string(thicknesses[inS]) + "u/";
  
  const std::array<unsigned,10> genEnAll = {5,10,20,30,40,60,80,100,150,200};
  //const std::array<unsigned,10> genEnAll = {150};
  const unsigned nGenEnAll = genEnAll.size();

  const unsigned nEvtMin = 150;

  const std::array<double,6> radius = {13,15,20,23,26,53};
  const unsigned nSR = radius.size(); 
  std::array<double,nSR> noise100, noise200, noise300;

  const std::array<unsigned,nSR> ncellsSmall = {7,13,19,31,37,151};
  const std::array<unsigned,nSR> ncellsLarge = {7,7,13,19,19,85};

  for (auto&& iSR: ip.signalRegions()) {
    noise100[iSR] = sqrt(pow(sqrt(10*ncellsSmall[iSR])*0.00192,2)+pow(sqrt(10*ncellsSmall[iSR])*0.00241,2)+pow(sqrt(8*ncellsSmall[iSR])*0.00325,2));
    noise200[iSR] = sqrt(pow(sqrt(10*ncellsLarge[iSR])*0.00097,2)+pow(sqrt(10*ncellsLarge[iSR])*0.00121,2)+pow(sqrt(8*ncellsLarge[iSR])*0.00164,2));
    noise300[iSR] = sqrt(pow(sqrt(10*ncellsLarge[iSR])*0.00049,2)+pow(sqrt(10*ncellsLarge[iSR])*0.00062,2)+pow(sqrt(8*ncellsLarge[iSR])*0.00083,2));
    std::cout << " Noise vals for iSR " << iSR << " = " << noise100[iSR] << " " << noise200[iSR] << " " << noise300[iSR] << " GeV." << std::endl;
  }
  //canvas so they are created only once
  TCanvas *mycE[nGenEnAll];
  for (unsigned iE(0); iE<nGenEnAll;++iE){
    std::ostringstream lName;
    lName << "mycE" << genEnAll[iE];
    mycE[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
  }

  TCanvas *mycE2D[nGenEnAll];
  if (redoLeakCor){
    for (unsigned iE(0); iE<nGenEnAll;++iE){
      std::ostringstream lName;
      lName << "mycE2D" << genEnAll[iE];
      mycE2D[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1);
    }
  }

  //strip the 'noSigmaEff' substring from the tag name
  std::string tag_modif;
  std::string ext("-noSigmaEff");
  if(ip.tag().size() > ext.size() and ip.tag().substr(ip.tag().size() - ext.size()) == ext)
    tag_modif = ip.tag().substr(0, ip.tag().size() - ext.size());
  else
    tag_modif = ip.tag();
  TString saveDir = "/eos/user/b/bfontana/www/RemoveLayers/" + ip.tag() + "/";

  for (auto&& icval: ICvals){//loop on intercalib

    for (auto&& version: ip.versions()) {
      //TString baseDir = "/afs/cern.ch/work/b/bfontana/PFCalEEAna/HGCalTDR/gitv" + version + "/";
      TString baseDir = "/afs/cern.ch/work/b/bfontana/PFCalEEAna/HGCalTDR/git" + tag_modif + "/";
      const unsigned nLayers = version=="60" ? 28 : 26;

      for (auto&& scenario: scenarios) {

	TString inputDir = baseDir+"version"+version+"/"+scenario;

	for (auto&& eta: ip.etas()) {
	  unsigned etax10 = static_cast<unsigned>(eta*10.);
	  
	  TTree *ltree[nPu][nGenEnAll];

	  double sigmaStoch[nPu];
	  double sigmaConst[nPu];
	  double sigmaNoise[nPu];
	  double sigmaStochErr[nPu];
	  double sigmaConstErr[nPu];
	  double sigmaNoiseErr[nPu];

	  for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	    sigmaStoch[ipu] = 0;
	    sigmaConst[ipu] = 0;
	    sigmaNoise[ipu] = 0;
	    sigmaStochErr[ipu] = 0;
	    sigmaConstErr[ipu] = 0;
	    sigmaNoiseErr[ipu] = 0;
	    
	    //loop over energies
	    //check energy point is valid
	    bool skip[nGenEnAll];
	    unsigned nValid = 0;
	    TFile *inputFile[nGenEnAll];
	    for (unsigned iE(0); iE<nGenEnAll; ++iE){
	      skip[iE] = false;
	      inputFile[iE] = 0;
	      ltree[ipu][iE] = 0;
	    }
	    for (unsigned iE(0); iE<nGenEnAll; ++iE){
	      std::ostringstream linputStr;
	      linputStr << inputDir ;
	      /*      
	      linputStr << "eta" << etax10 << "_et" << genEnAll[iE];// << "_pu" << pu[ipu];
	      if (pu[ipu]>0) linputStr << "_Pu" << pu[ipu];
	      linputStr << "_IC" << icval;// << "_Si2";
	      */
	      linputStr << "ana_npuvtx" << pu[ipu] << "_ic" << icval << "_en" << genEnAll[iE] << "_eta" << etax10;

	      linputStr << ".root";
	      inputFile[iE] = TFile::Open(linputStr.str().c_str());
	      if (!inputFile[iE]) {
		std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
		//	    return 1;
		skip[iE] = true;
	      }
	      else {
		inputFile[iE]->cd("Energies");
		ltree[ipu][iE] = (TTree*)gDirectory->Get("Ereso");
		
		if (!ltree[ipu][iE]){
		  std::cout << " -- File " << inputFile[iE]->GetName() << " sucessfully opened but tree Ereso not found! Skipping." << std::endl;
		  skip[iE] = true;
		} else { 
		  std::cout << " -- File " << inputFile[iE]->GetName() << " sucessfully opened." << std::endl
			    << " ---- Tree found with " << ltree[ipu][iE]->GetEntries() << " entries." << std::endl;
		  if (ltree[ipu][iE]->GetEntries()<nEvtMin) {
		    std::cout << " -- Tree has only " << ltree[ipu][iE]->GetEntries() << " entries, skipping..." << std::endl;
		    skip[iE] = true;
		  } 
		  else nValid++;
		}
	      }
	    }//loop on energies

	    
	    if (nValid < 1) return 1;

	    unsigned newidx = 0;
	    const unsigned nGenEn = nValid;
	    unsigned genEn[nGenEn];
	    unsigned oldIdx[nGenEn];
	    for (unsigned iE(0); iE<nGenEnAll; ++iE){	  
	      if (!skip[iE]) {
		genEn[newidx]=genEnAll[iE];
		oldIdx[newidx]=iE;
		newidx++;
	      }
	    }

	    std::vector<double> noisePu;

	    for (auto&& iSR: ip.signalRegions()) { //loop on signal regions
	      std::cout << " - Processing signal region: " << iSR << " with size " << radius[iSR] << std::endl;

	      TString plotDir = saveDir + "version" + version + "/model2/gamma/SR" + std::to_string(iSR) + "/";
	      if (system(TString("mkdir -p ")+plotDir)) return 1;

	      double calib = 0, offset = 0, calibErr = 0, offsetErr = 0;
	      
	      for (unsigned iT(redoCalib && pu[ipu]==0 ? 0 : 1); iT<2; ++iT){//loop on calib type
		bool doRaw = iT==0?true:false;
		std::ostringstream loutputStr,linputStr;
		loutputStr << plotDir << "IC" << icval << "_pu" << pu[ipu] << "_SR" << iSR << "_Eta" << etax10; 
		if (dovsE) loutputStr  << "_vsE";
		if (doBackLeakCor) loutputStr  << "_backLeakCor";
		linputStr << loutputStr.str();
		if (doRaw) loutputStr << "_raw";
		loutputStr  << ".root";
		TFile *foutEfit = TFile::Open(loutputStr.str().c_str(),"UPDATE");
		if (!foutEfit){
		  std::cout << " -- Could not create output file " << loutputStr.str() << ". Exit..." << std::endl;
		  return 1;
		}

		TFile *finEfit = 0;
		if (!redoCalib){
		  linputStr << "_raw";
		  linputStr  << ".root";
		  finEfit = TFile::Open(linputStr.str().c_str());
		  if (!finEfit){
		    std::cout << " -- Could not open input file " << linputStr.str() << ". Exit..." << std::endl;
		    return 1;
		  }
		}


		std::cout << " -- Output histograms will be saved in: " << foutEfit->GetName() << std::endl;

		////
		///////// add flag to get graphs from file instead of new
		////
		foutEfit->cd();
		TGraphErrors *corrBackLeakFit = 0;
		if (!doRaw && redoLeakCor) {
		  corrBackLeakFit = new TGraphErrors();
		  corrBackLeakFit->SetName("corrBackLeakFit");
		}
		else if (!doRaw && doBackLeakCor){
		  //get from file
		  foutEfit->cd();
		  corrBackLeakFit = (TGraphErrors*)gDirectory->Get("corrBackLeakFit");
		  if (!corrBackLeakFit) {
		    std::cout << " -- back leak corr TGraph not found! " << std::endl;
		    return 1;
		  }
		}
		/////
		TGraphErrors *calibRecoFit = 0;
		TGraphErrors *calibRecoDelta = 0;
		TGraphErrors *resoRecoFit = 0;
		if (!redoCalib){
		  //get calib and offset from file....
		  std::cout << " -- Reading input calibration histograms from: " << finEfit->GetName() << std::endl;
		  finEfit->cd();
		  calibRecoFit = (TGraphErrors*)gDirectory->Get("calibRecoFitRaw");
		  calibRecoDelta = (TGraphErrors*)gDirectory->Get("calibRecoDeltaRaw");
		  if (!calibRecoFit || !calibRecoDelta){
		    std::cout << " -- Calib TGraphs not found !" << std::endl;
		    return 1;
		  }
		  TF1 *fitFunc = calibRecoFit->GetFunction("calib");
		  calib = fitFunc->GetParameter(1);
		  offset = fitFunc->GetParameter(0);
		  calibErr = fitFunc->GetParError(1);
		  offsetErr = fitFunc->GetParError(0);
		  
		  std::cout << " -- Getting calibration factors from file: " << std::endl
			    << " slope " << calib << " +/- " << calibErr
			    << " offset " << offset << " +/- " << offsetErr
			    << std::endl;
		  calibRecoFit->Delete();
		  calibRecoFit = new TGraphErrors();
		  calibRecoFit->SetName("calibRecoFit");
		  calibRecoDelta->Delete();
		  calibRecoDelta =  new TGraphErrors();
		  calibRecoDelta->SetName("calibRecoDelta");

		}
		else {
		  foutEfit->cd();
		  if (calibRecoFit) calibRecoFit->Delete();
		  calibRecoFit = new TGraphErrors();
		  calibRecoFit->SetName("calibRecoFitRaw");
		  if (calibRecoDelta) calibRecoDelta->Delete();
		  calibRecoDelta =  new TGraphErrors();
		  calibRecoDelta->SetName("calibRecoDeltaRaw");
		}
		
		if (resoRecoFit) resoRecoFit->Delete();
		resoRecoFit = new TGraphErrors();
		if (doRaw) resoRecoFit->SetName("resoRecoFitRaw");
		else resoRecoFit->SetName("resoRecoFit");



		if (doRaw) std::cout << " -- Doing RAW calibration fit " << std::endl;
		else std::cout << " -- Check linearity fit " << std::endl;
		//fit of rawE
		
		
		///////////////////////////////////////////////////////
		////////////// TO BE UPDATED IF PU NEEDED /////////////
		///////////////////////////////////////////////////////
		//fill sigma PU for all SR for 40 GeV E
		if (iT>0 && pu[ipu]>0){// && pT>5 && pT<40) {
		  unsigned iE = 6;
		  std::vector<double> calib = GetCalib(foutEfit->GetName());
		  bool success = retrievePuSigma(pu[ipu],version,nLayers,
						 ltree[0][oldIdx[iE]], ltree[ipu][oldIdx[iE]],
						 calib,
						 eta,
						 noisePu);
		  if (!success) return 1;
		}
		///////////////////////////////////////////////////////
		///////////////////////////////////////////////////////
		///////////////////////////////////////////////////////


		for (unsigned iE(0); iE<nGenEn; ++iE){
		  
		  std::cout << " --- Processing energy : " << genEn[iE] 
			    << std::endl;
		  
		  inputFile[iE]->cd("Energies");
		  std::cout << " -- Tree entries for eta=" << etax10 << " pu=" << pu[ipu] << " : " << ltree[ipu][oldIdx[iE]]->GetEntries() << std::endl;

		  double backLeakCor = 0;

		  if (redoLeakCor && !doRaw) {
		    int fail = makeBackLeakCor(nLayers, ip.nBack(),
					       iSR, genEn[iE], etax10, pu[ipu],
					       offset, calib,
					       backLeakCor,
					       mycE2D[oldIdx[iE]],
					       ltree[ipu][oldIdx[iE]],
					       foutEfit,
					       corrBackLeakFit,
					       plotDir) ;
		    if (doBackLeakCor && fail){
		      std::cout << " -- Fit failed for back leakage correction, for Egen " << genEn[iE] << std::endl;
		      return 1;
		    }
		  }
		  else {
		    //get it from some saved input file !


		  }
		  
		  int success = makeEfit(ip.useSigmaEff(), dovsE, doRaw,
					 doBackLeakCor, ip.nBack(),
					 pu[ipu], icval,
					 scenario, version, nLayers,
					 etax10,genEn[iE], iSR, radius[iSR],
					 offset,calib,backLeakCor,
					 mycE[oldIdx[iE]],
					 ltree[ipu][oldIdx[iE]],
					 foutEfit,
					 plotDir,
					 calibRecoFit,
					 resoRecoFit
					 );
		}//loop on energies

		//get calib
		makeCalibration(doRaw,doBackLeakCor,
				etax10,pu[ipu], iSR, radius[iSR],
				calibRecoFit,calibRecoDelta,
				calib,calibErr,
				offset,offsetErr,
				plotDir,foutEfit
				);

		foutEfit->cd();
		if (redoLeakCor && !doRaw) {
		  corrBackLeakFit->Write();
		  plotBackLeakFit(plotDir,corrBackLeakFit,etax10,pu[ipu]);
		}
		double noise = eta==2.4f ? noise100[iSR] : eta==2.0f ? noise200[iSR] : eta==1.7f ? noise300[iSR] : 0 ;
		double noiseRef = pu[ipu]==0? noise : noisePu[iSR];

		makeResolution(dovsE,doRaw,doBackLeakCor,
			       etax10,pu[ipu],iSR,radius[iSR],
			       version,
			       resoRecoFit,
			       sigmaStoch[0],
			       sigmaConst[0],
			       noiseRef,
			       sigmaStoch[ipu],sigmaStochErr[ipu],
			       sigmaConst[ipu],sigmaConstErr[ipu],
			       sigmaNoise[ipu],sigmaNoiseErr[ipu],
			       plotDir,foutEfit
			       );


		foutEfit->Write();
	      }//loop on calib type

	    }//loop on SR
	  }//loop on pu
	}//loop on eta

      }//loop on scenario
	 
    }//loop on version
  }//loop on intercalib


  return 0;

};//main

//compile with `g++ plotEGReso_standalone.cc $(root-config --glibs --cflags --libs) -o plotEGReso`
//run with `./plotEGReso --thicknesses 1 200 --tag V08-08-00 --nBack 2 --etas 3 1.7 2.0 2.4 --versions 1 60 --signalRegions 1 4` (example)
int main(int argc, char** argv)
{
  InputParserPlotEGReso ip(argc, argv);

  plotEGReso(ip);
  return 0;
}
