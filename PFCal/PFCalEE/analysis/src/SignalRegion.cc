#include "SignalRegion.hh"
#include "HGCSSEvent.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "utilities.h"

SignalRegion::SignalRegion(const std::string inputFolder,
			   const unsigned nLayers,
			   const std::vector<double> & zpos,
			   const unsigned nevt,
			   HGCSSGeometryConversion & geomConv,
			   const unsigned versionNumber,
			   const bool doHexa,
			   const unsigned g4trackID){

  g4trackID_ = g4trackID;
  doHexa_ = doHexa;
  nSR_ = 6;
  radius_[0] = 9;
  radius_[1] = 14;
  radius_[2] = 21;
  radius_[3] = 25;
  radius_[4] = 37;
  radius_[5] = 53;
  nevt_ = nevt;
  inputFolder_ = inputFolder;
  nLayers_ = nLayers;
  //geomConv_ = geomConv;
  geomConv_.setCellSize(geomConv.cellSize());
  geomConv_.hexagonMap(*geomConv.hexagonMap());
  geomConv_.diamondMap(*geomConv.diamondMap());
  geomConv_.triangleMap(*geomConv.triangleMap());
  geomConv_.squareMap(*geomConv.squareMap());
  
  geomConv_.copyhexaGeom(geomConv.hexaGeom);
  geomConv_.copydiamGeom(geomConv.diamGeom);
  geomConv_.copytriangleGeom(geomConv.triangleGeom);
  geomConv_.copysquareGeom(geomConv.squareGeom);

  firstEvent_ = true;
  nSkipped_ = 0;

  zPos_ = zpos;
  puDensity_ = HGCSSPUenergy();
  fixForPuMixBug_ =  false;

}

SignalRegion::SignalRegion(const std::string inputFolder,
			   const unsigned nLayers,
			   const std::vector<double> & zpos,
			   const unsigned nevt,
			   HGCSSGeometryConversion & geomConv,
			   const HGCSSPUenergy & puDensity,
			   const bool applyPuMixFix,
			   const unsigned versionNumber,
			   const bool doHexa,
			   const unsigned g4trackID){


  SignalRegion(inputFolder,nLayers,zpos,nevt,geomConv,versionNumber,doHexa,g4trackID);
  puDensity_ = puDensity;
  fixForPuMixBug_ = applyPuMixFix;

}

SignalRegion::~SignalRegion(){
}

bool SignalRegion::initialiseFitPositions(){

  std::ifstream fxypos;
  std::ostringstream finname;
  finname << inputFolder_ << "/accuratePos.dat";
  fxypos.open(finname.str());
  if (!fxypos.is_open()){
    std::cout << " Cannot open input file " << finname.str() << "! Exiting..." << std::endl;
    return false;
  }

  std::cout << " -- Accurate positions found in file " << finname.str() << std::endl;
  
  //all events did not pass the chi2 fit: fill only those found.
  //keep failed ones to emptyvec so they are ignored afterwards
  accurateFit_.clear();
  FitResult res;
  accurateFit_.resize(nevt_,res);

  unsigned nfound = 0;

  while (!fxypos.eof()){
    unsigned eventIndex = nevt_;
    double xpos(0),ypos(0),xangle(0),yangle(0);
    double fitMatrix[4] = {0,0,0,0};
    double xpost(0),ypost(0),xanglet(0),yanglet(0);
    fxypos >> eventIndex >> xpos >> fitMatrix[0] >> xangle >> fitMatrix[1] >> ypos >> fitMatrix[2] >> yangle >> fitMatrix[3] >> xpost >> xanglet >> ypost >> yanglet;
    //testing for nan
    if ( eventIndex != eventIndex || xpos != xpos || fitMatrix[0]!=fitMatrix[0] || xangle!=xangle || fitMatrix[1]!=fitMatrix[1] || ypos!=ypos || fitMatrix[2]!=fitMatrix[2] || yangle!=yangle || fitMatrix[3]!=fitMatrix[3]){
      std::cout << " Found nan ! Fix code !" << std::endl;
      std::cout << eventIndex << " " << xpos << " " << fitMatrix[0] << " " << xangle << " " << fitMatrix[1] << " " << ypos << " " << fitMatrix[2] << " " << yangle << " " << fitMatrix[3]<< std::endl;
      exit(1);
    }
    if (eventIndex<nevt_) {
      accurateFit_[eventIndex].pos_x = xpos;
      accurateFit_[eventIndex].pos_y = ypos;
      accurateFit_[eventIndex].pos_z = 0;
      accurateFit_[eventIndex].tanangle_x = xangle;
      accurateFit_[eventIndex].tanangle_y = yangle;
      accurateFit_[eventIndex].found=true;
      nfound++;
    }
    else break;
  }
  //if not all events found
  if (nfound < nevt_) {
    std::cout << " Warning, file " << finname.str() << " contains only " << nfound 
	      << " events, program running on " << nevt_ 
	      << "." << std::endl;
    if (nfound*1./nevt_ < 0.5) return false;
  }
  std::cout << " -- Now filling signal region histograms..." << std::endl;
  return true;

}

void SignalRegion::initialise(TFile *outputFile,
			      const std::string outputDir){

  //mycalib_ = mycalib;
  outputDir_ = outputDir;
  setOutputFile(outputFile);
  initialiseHistograms();
}

const FitResult & SignalRegion::getAccurateFit(const unsigned ievt) const{
  return accurateFit_[ievt];
}

ROOT::Math::XYZPoint SignalRegion::getAccuratePos(const unsigned ievt, const unsigned iL) const{
  return getAccuratePos(accurateFit_[ievt],iL);
}

ROOT::Math::XYZPoint SignalRegion::getAccuratePos(const FitResult& fit, const unsigned iL) const{
  return ROOT::Math::XYZPoint(fit.pos_x + fit.tanangle_x*(zPos_[iL]-fit.pos_z), fit.pos_y + fit.tanangle_y*(zPos_[iL]-fit.pos_z), zPos_[iL]);
}

Direction SignalRegion::getAccurateDirection(const unsigned ievt) const{
  return Direction(accurateFit_[ievt].tanangle_x,accurateFit_[ievt].tanangle_y);
}


bool SignalRegion::setTruthInfo(const std::vector<HGCSSGenParticle> & genvec,
				const HGCSSEvent & event,
				const int G4TrackID){

  bool found = false;
  for (unsigned iP(0); iP<genvec.size(); ++iP){//loop on gen particles    
    //if (debug_>1 && genvec.size()!= 1) 
    if (genvec[iP].pdgid()==22 && genvec[iP].trackID()<500) {
      std::cout << " iP " <<  iP << " ";
      genvec[iP].Print(std::cout);
    }
    else continue;
    if (genvec[iP].trackID()==G4TrackID){
      found = true;
      //set direction
      truthDir_ = Direction(genvec[iP].px()/genvec[iP].pz(),genvec[iP].py()/genvec[iP].pz());
      double xvtx=event.vtx_x(),yvtx=event.vtx_y();
      double zvtx = event.vtx_z();//genvec[iP].z()-( (genvec[iP].x()+genvec[iP].y()-xvtx-yvtx)/(truthDir_.tanangle_x+truthDir_.tanangle_y) );
      //double xvtx = genvec[iP].x()-((genvec[iP].z()-zvtx)*truthDir_.tanangle_x);
      truthVtx_ = ROOT::Math::XYZPoint(xvtx,yvtx,zvtx);
      //in GeV
      trueE_ = genvec[iP].E()/1000.;
      //if (debug_) 
      std::cout << " Found truth vertex for particle pdgid " << genvec[iP].pdgid() << " pos at: x=" 
		<< xvtx << ", y" << yvtx << ", z=" << zvtx 
		<< " x,ypos was: " 
		<< genvec[iP].x() 
		<< " " << genvec[iP].y() 
		<< std::endl;
      
      break;
    }
    
  }
  if (!found){
    std::cout << " - Info: no photon G4trackID=" << G4TrackID << " found, already converted or outside of acceptance ..." << std::endl;
    //std::cout << " -- Number of genparticles: " << genvec.size() << std::endl;
    //for (unsigned iP(0); iP<genvec.size(); ++iP){//loop on gen particles 
    //std::cout << " --- particle " << iP << std::endl;
    //genvec[iP].Print(std::cout);
    //}
  }
  return found;
}


bool SignalRegion::fillEnergies(const unsigned ievt,
				const HGCSSEvent & event,
				const std::vector<HGCSSGenParticle> & genvec,
				const std::vector<HGCSSSamplingSection> & ssvec,
				const std::vector<HGCSSSimHit> & simhitvec,
				const std::vector<HGCSSRecoHit> & rechitvec,
				const unsigned nPuVtx){
   bool found = false;
   FitResult fit;
   for (unsigned iP(0); iP<genvec.size(); ++iP){//loop on gen particles    
     if (genvec[iP].trackID()==g4trackID_){
       found = true;
       //set direction
       trueE_ = genvec[iP].E()/1000.;
       fit.pos_x = genvec[iP].x();
       fit.pos_y = genvec[iP].y();
       fit.pos_z = genvec[iP].z();
       trueEta_ = genvec[iP].eta();
       truePhi_ = genvec[iP].phi();
       
       fit.tanangle_x = genvec[iP].px()/genvec[iP].pz();
       fit.tanangle_y = genvec[iP].py()/genvec[iP].pz();
       //to do in opposite region to get pile-up only...
       //fit.tanangle_x = -1.*genvec[iP].px()/genvec[iP].pz();
       //fit.tanangle_y = -1.*genvec[iP].py()/genvec[iP].pz();
       fit.found = true;
       //std::cout << " True particle pdgid=" << genvec[iP].pdgid() << " E = " << trueE_ << " z eta phi = " << genvec[iP].z() << " " << genvec[iP].eta() << " " << genvec[iP].phi() << std::endl;
     }
   }
   if (!found){
     std::cout << " - Info: no photon G4trackID=" << g4trackID_ << " found, already converted or outside of acceptance ..." << std::endl;
     fit.found = false;
   }

   vtxX_ = event.vtx_x();
   vtxY_ = event.vtx_y();
   vtxZ_ = event.vtx_z();
 
   return fillEnergies(ievt,ssvec,simhitvec,rechitvec,nPuVtx,fit);
}

bool SignalRegion::fillEnergies(const unsigned ievt,
				const std::vector<HGCSSSamplingSection> & ssvec,
				const std::vector<HGCSSSimHit> & simhitvec,
				const std::vector<HGCSSRecoHit> & rechitvec,
				const unsigned nPuVtx){
  const FitResult & fit = accurateFit_[ievt];
  return fillEnergies(ievt,ssvec,simhitvec,rechitvec,nPuVtx,fit);
}

bool SignalRegion::fillEnergies(const unsigned ievt,
				const std::vector<HGCSSSamplingSection> & ssvec,
				const std::vector<HGCSSSimHit> & simhitvec,
				const std::vector<HGCSSRecoHit> & rechitvec,
				const unsigned nPuVtx,
				const FitResult & fit){

  if(!fit.found) {
    std::cout << " -- Event " << ievt << " skipped, accurate position not found." << std::endl;
    nSkipped_++;
    //fill tree to find correspondance between noPu and PU...
    outtree_->Fill();
    return false;
  }
  //initialise accuratepos per layer
  std::vector<ROOT::Math::XYZPoint> eventPos;
  eventPos.resize(nLayers_,ROOT::Math::XYZPoint(0,0,0));
  for (unsigned iL(0); iL<nLayers_;++iL){
    eventPos[iL] = getAccuratePos(fit,iL);

    //std::cout << " Layer " << iL << " best pos = " << eventPos[iL].X() << " " << eventPos[iL].Y() << " " << eventPos[iL].Z() << std::endl;
  }

  return fillEnergies(ievt,ssvec,simhitvec,rechitvec,nPuVtx,eventPos);

}

bool SignalRegion::fillEnergies(const unsigned ievt,
				const std::vector<HGCSSSamplingSection> & ssvec,
				const std::vector<HGCSSSimHit> & simhitvec,
				const std::vector<HGCSSRecoHit> & rechitvec,
				const unsigned nPuVtx,
				const std::vector<ROOT::Math::XYZPoint> & eventPos){
  
 
  //fill weights for first event only: same in all events
  if (firstEvent_){
    std::cout << " -- Absorber weights used for total energy:" << std::endl;
    for(unsigned iL(0); iL<nLayers_; iL++){
      //double w = absWeight(iL);
      //double w = ssvec[iL].voldEdx()/ssvec[1].voldEdx();
      //double w = ssvec[iL].volX0trans()/ssvec[1].volX0trans();
      //absweight_[iL] = 10.0166;
      absweight_[iL]=iL<(nLayers_-2) ? (ssvec[iL].voldEdx()+ssvec[iL+1].voldEdx())/2. : absweight_[nLayers_-3];
      std::cout << " - Layer " << iL << " wdEdx=" << ssvec[iL].voldEdx()  << " Wx0=" << ssvec[iL].volX0trans() << " W valeri scheme " << absweight_[iL] << std::endl;
      //if (fabs(zPos_[iL]-ssvec[iL].sensitiveZ()) > 0.5){
      //std::cout << " -- error with zPos: not aligned with value in tree." << std::endl;
	//zPos_[iL]=ssvec[iL].sensitiveZ();
      //}
    }
    //absweight_[0] = 20.3628;
    //absweight_[nLayers_-1] = 13.0629;
    for(unsigned iL(0); iL<nLayers_; iL++){
      std::cout << zPos_[iL] << " " << ssvec[iL].sensitiveZ() << std::endl;
    }
    
    firstEvent_=false;
  }

  if (absweight_.size()!=nLayers_) {
    std::cout << " -- Error! Not all layers found! Only: " << absweight_.size() << ". Fix code." << std::endl;
    exit(1);
  }
  
    //initialise values for current event
  evtIdx_ = ievt;
  totalE_ = 0;
  wgttotalE_ = 0;
  mR68_ = 0;
  mR90_ = 0;

  for (unsigned iL(0); iL<nLayers_;++iL){
    dR68_[iL] = 0;
    dR90_[iL] = 0;
    E68_[iL] = 0;
    E90_[iL] = 0;
    E100_[iL] = 0;
    for (unsigned iSR(0);iSR<nSR_;++iSR){
      energySR_[iL][iSR] = 0;
      subtractedenergySR_[iL][iSR] = 0;
      maxhitEoutside_[iL][iSR] = 0;
    }
    for (unsigned idx(0);idx<9;++idx){
      Exy_[iL][idx] = 0;
    }
  }
  double refx[nLayers_],refy[nLayers_];

  for (unsigned iL(0); iL<nLayers_;++iL){
    int refid = 0;
    if (doHexa_) refid = geomConv_.hexagonMap()->FindBin(eventPos[iL].X(),eventPos[iL].Y());
    else refid = geomConv_.squareMap()->FindBin(eventPos[iL].X(),eventPos[iL].Y());
    refx[iL] = doHexa_ ? geomConv_.hexaGeom[refid].first : geomConv_.squareGeom[refid].first;
    refy[iL] = doHexa_ ? geomConv_.hexaGeom[refid].second : geomConv_.squareGeom[refid].second;
  }
  
  //std::cout << " -- Accurate direction for evt " << ievt << ": " << std::endl;
  //getAccurateDirection(ievt).Print();
  
  //double etacor = 1./tanh(getAccurateDirection(ievt).eta());

  //get event-by-event PU
  //get PU contrib from elsewhere in the event
  //loop over phi with same etamax
  //take average per layer: not all 9 cells of 3*3 area have hits...
  /*	std::vector<double> puE;
	puE.resize(nLayers_,0);
	if (nPuVtx>0){
	unsigned nRandomCones = 50;
	double phistep = TMath::Pi()/nRandomCones;
	for (unsigned ipm(0);ipm<nRandomCones;++ipm){
	std::vector<double> xmaxrc;
	xmaxrc.resize(nLayers_,0);
	std::vector<double> ymaxrc;
	ymaxrc.resize(nLayers_,0);
	double phirc = phimax-TMath::Pi();
	if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
	if (ipm%2==0) phirc += ipm/2*phistep+phistep/2.;
	else  phirc = phirc - ipm/2*phistep-phistep/2.;
	if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
	//take from geom to not be biased by hit having PU, because
	//not from geom means find cell with a hit closest to maxpos...
	getMaximumCellFromGeom(phirc,etamax,xmaxrc,ymaxrc);
	getPuContribution(rechitvec,xmaxrc,ymaxrc,puE);
	}
	
	//normalise to one cell: must count cells with 0 hit !
	//use cell size at etamax...
	for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
	unsigned nCells = nRandomCones*nSR_*geomConv_.cellSize()/geomConv_.cellSize(iL,etamax);
	puE[iL] = puE[iL]/nCells;
	}
	
	}//if PU
  */
  
  
  // Define different signal region and sum over energy
  //double maxE[nLayers_];
  //double maxH[nLayers_];
  //for (unsigned iL(0);iL<nLayers_;++iL){
  //maxE[iL] = 0;
  //maxH[iL] = 0;
  //}

  std::vector<MyRecoHit> lhitvec[nLayers_];
  std::vector<MyRecoHit> lhitvectotal;

  for (unsigned iH(0); iH<rechitvec.size(); ++iH){//loop on hits
    const HGCSSRecoHit & lHit = rechitvec[iH];
    
    unsigned layer = lHit.layer();
    if (layer >= nLayers_) {
      continue;
    }
    double leta = lHit.eta();
    //std::cout << " hit eta=" << leta << std::endl;
    //not interested in hits outside of acceptance...
    //but need physics eta, line below works only for vtx 0 0 0...
    //if (leta<1.4 || leta>3.0) continue;

    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;
    double energy = lHit.energy();
    //bool notNoiseOnly = fabs(lHit.noiseFraction()-1)>0.001;
    //double etacor = fabs(1./tanh(leta));

    //if (energy>maxE[layer]) {
    //maxE[layer] = energy;
    //maxH[layer] = iH;
    //}

    totalE_ += energy;
    wgttotalE_ += energy*absweight_[layer];    
    
    double lradius = sqrt(pow(posx,2)+pow(posy,2));
    double puE = 0;
    if (nPuVtx>0) puE = puDensity_.getDensity(leta,layer,geomConv_.cellSizeInCm(layer,lradius),nPuVtx);
    double subtractedenergy = std::max(0.,energy - puE);
    //double halfCell = 0.5*geomConv_.cellSize(layer,lradius);
    //hexagons are side up....
    //std::cout << geomConv_.cellSize(layer,lradius) << std::endl;
    double distance = sqrt(3.)*geomConv_.cellSize(layer,lradius);
    double halfCelly = 0.5*distance;
    //double halfCellx = doHexa_?geomConv_.cellSize(layer,lradius):geomConv_.cellSize(layer,lradius)/2.;
    double halfCellx = doHexa_?geomConv_.cellSize(layer,lradius):5;
    //std::cout << " halfcell = " << halfCell << std::endl;

    //SR0-4

    double dx = posx-refx[layer];
    double dy = posy-refy[layer];
    double dr = sqrt(dx*dx+dy*dy);

    MyRecoHit ltmpHit;
    ltmpHit.dR = dr;
    ltmpHit.E = energy*absweight_[layer];
    lhitvec[layer].push_back(ltmpHit);
    lhitvectotal.push_back(ltmpHit);
    E100_[layer] += energy*absweight_[layer];
    

    for (unsigned isr(0); isr<nSR_;++isr){
      //double dx = isr%2==0? posx-refx[layer] : posx-eventPos[layer].x();
      //double dy = isr%2==0? posy-refy[layer] : posy-eventPos[layer].y();
      /*if (isr==nSR_-1 && dr<50. ){
	energySR_[layer][isr] += energy;
	subtractedenergySR_[layer][isr] += subtractedenergy;
	}
	else if ( (doHexa_ && ((isr%2==0 && dr<(isr*halfCellx+0.1)) || (isr%2==1 && dr <= ((isr+1)*halfCelly) && (fabs(dx) <= ((isr+1)*halfCellx)) && (fabs(dy) <= ((isr+1)*halfCelly))))) || 
	(!doHexa_ && (fabs(dx) <= ((isr+1)*halfCelly)) && (fabs(dy) <= ((isr+1)*halfCelly))) )*/
      if ((doHexa_ && nSR_<=6 && dr<radius_[isr]) || 
	(!doHexa_ && (fabs(dx) <= ((isr+1)*halfCelly)) && (fabs(dy) <= ((isr+1)*halfCelly))) )
	{
	  energySR_[layer][isr] += energy;//*absweight_[layer]*etacor;
	  subtractedenergySR_[layer][isr] += subtractedenergy;//*absweight_[layer];//*etacor;
	//save energy in 1+6 cells
	  if (isr== (doHexa_?0:2) ){
	  int ix = doHexa_ ? dx/halfCellx : dx/(2*halfCelly);
	  int iy = doHexa_ ? dy/(1.5*halfCelly) : dy/(2*halfCelly);
	  unsigned idx = 0;
	  if (halfCellx>1) {
	    if (((doHexa_ && (ix > 2 || ix < -2)) || (!doHexa_ && (ix > 1 || ix < -1))) || (iy>1 || iy<-1)) {
	      std::cout << " error, isr = " << isr << " check ix=" << ix << " iy=" << iy << " posx,y-max=" << dx << " " << dy << " step " ;
	      if (!doHexa_) std::cout << (isr+1)*halfCelly;
	      else std::cout << isr*halfCelly+0.1;
	      std::cout << std::endl;
	      continue;
	    }
	    else {
	      if (!doHexa_) idx = 3*(iy+1)+(ix+1);	    
	      else {
		if (ix==-1 && iy==-1) idx=0;
		else if (ix==1 && iy==-1) idx = 1;
		else if (ix==-2 && iy==0) idx = 2;
		else if (ix==0 && iy==0) idx = 3;
		else if (ix==2 && iy==0) idx = 4;
		else if (ix==-1 && iy==1) idx = 5;
		else if (ix==1 && iy==1) idx = 6;
	      }
	    }
	    Exy_[layer][idx] = energy;	  
	  }
	}

      }
      else {
	if (energy>maxhitEoutside_[layer][isr]) maxhitEoutside_[layer][isr] = energy;
      }
    }//loop on SR
  }//loop on hits

  //fill 68 and 90% containment radius
  for (unsigned iL(0);iL<nLayers_;++iL){
    std::sort(lhitvec[iL].begin(), lhitvec[iL].end(), customdRsort);
    double Esum = 0;
    //std::cout << " -- Layer: " << iL << " nH = " << lhitvec[iL].size() << std::endl;
    for (unsigned iH(0); iH<lhitvec[iL].size(); ++iH){
      //std::cout << " ---- " << iH << " dR " << lhitvec[iL][iH].dR 
      //	<< " Esum " << Esum
      //	<< " E100 " << E100_[iL] 
      //	<< " (68%=" << 0.68*E100_[iL] << ")"
      //	<< " E68 " << E68_[iL]
      //	<< " dR68 " << dR68_[iL]
      //	<< std::endl;
      //if (lhitvec[iL][iH].dR>0) 
      Esum += lhitvec[iL][iH].E;
      p_EsumfracvsdR[iL]->Fill(lhitvec[iL][iH].dR,Esum/E100_[iL]);
      p_EvsdR[iL]->Fill(lhitvec[iL][iH].dR,lhitvec[iL][iH].E);
      p_dR[iL]->Fill(lhitvec[iL][iH].dR,lhitvec[iL][iH].E);
      if (Esum<0.68*E100_[iL]) {
	dR68_[iL] = lhitvec[iL][iH].dR;
	E68_[iL] = Esum;
      }
      if (Esum<0.9*E100_[iL]){
	dR90_[iL] = lhitvec[iL][iH].dR;
	E90_[iL] = Esum;
      }
    }
  }

  //fill moliere radius total
  std::sort(lhitvectotal.begin(), lhitvectotal.end(), customdRsort);
  double Esum = 0;
  for (unsigned iH(0); iH<lhitvectotal.size(); ++iH){
    //std::cout << iL << " " << iH << " " << lhitvec[iL][iH].dR << std::endl;
    Esum += lhitvectotal[iH].E;
    if (Esum<0.68*wgttotalE_) {
      mR68_ = lhitvectotal[iH].dR;
    }
    if (Esum<0.9*wgttotalE_){
      mR90_ = lhitvectotal[iH].dR;
    }
  }

  /*for (unsigned iL(0);iL<nLayers_;++iL){
    const HGCSSRecoHit & lmaxHit = rechitvec[maxH[iL]];
    std::cout << "max hit layer " << iL << " x=" 
	      << lmaxHit.get_x() << " y="
	      << lmaxHit.get_y() << " eta="
	      << lmaxHit.eta() << " E="
	      << lmaxHit.energy()
	      << std::endl;
	      }*/

  outputFile_->cd(outputDir_.c_str());
  fillHistograms();
  outtree_->Fill();
  return true;
}

void SignalRegion::finalise(){
  outputFile_->Flush();
  std::cout << " -- Histograms for signal regions have been filled !" << std::endl;
  std::cout << " -- Number of skipped events: " << nSkipped_ << std::endl;
}


void SignalRegion::initialiseHistograms(){

    outputFile_->cd(outputDir_.c_str());

    outtree_ = new TTree("Ereso","Tree to save energies in signal regions");


    outtree_->Branch("eventIndex",&evtIdx_);
    outtree_->Branch("rawEtotal",&totalE_);
    outtree_->Branch("wgtEtotal",&wgttotalE_);
    outtree_->Branch("trueE",&trueE_);
    outtree_->Branch("vtxX",&vtxX_);
    outtree_->Branch("vtxY",&vtxY_);
    outtree_->Branch("vtxZ",&vtxZ_);
    outtree_->Branch("trueEta",&trueEta_);
    outtree_->Branch("truePhi",&truePhi_);
    outtree_->Branch("mR68",&mR68_);
    outtree_->Branch("mR90",&mR90_);

    std::vector<double> emptyvec;
    emptyvec.resize(nSR_,0);
    energySR_.resize(nLayers_,emptyvec);
    subtractedenergySR_.resize(nLayers_,emptyvec);
    maxhitEoutside_.resize(nLayers_,emptyvec);
    absweight_.resize(nLayers_,1);
    dR68_.resize(nLayers_,1);
    dR90_.resize(nLayers_,1);
    E68_.resize(nLayers_,1);
    E90_.resize(nLayers_,1);
    E100_.resize(nLayers_,1);

    if (zPos_.size()<nLayers_) {
      std::cout << " -- ERROR! z positions not filled properly. Size of vector " << zPos_.size() << " nLayers_= " << nLayers_ << ". Exiting..." << std::endl;
      exit(1);
    }

    std::vector<double> init;
    init.resize(9,0);
    Exy_.resize(nLayers_,init);

    std::ostringstream label;
    for (unsigned iL(0); iL<nLayers_;++iL){
      label.str("");
      label << "absweight_" << iL;
      outtree_->Branch(label.str().c_str(),&absweight_[iL]);

      label.str("");
      label << "zpos_" << iL;
      outtree_->Branch(label.str().c_str(),&zPos_[iL]);

      label.str("");
      label << "rho68_" << iL;
      outtree_->Branch(label.str().c_str(),&dR68_[iL]);
      label.str("");
      label << "rho90_" << iL;
      outtree_->Branch(label.str().c_str(),&dR90_[iL]);

      label.str("");
      label << "E68_" << iL;
      outtree_->Branch(label.str().c_str(),&E68_[iL]);
      label.str("");
      label << "E90_" << iL;
      outtree_->Branch(label.str().c_str(),&E90_[iL]);
      label.str("");
      label << "E100_" << iL;
      outtree_->Branch(label.str().c_str(),&E100_[iL]);

      for (unsigned iSR(0);iSR<nSR_;++iSR){
	label.str("");
	label << "energy_" << iL << "_SR" << iSR;
	outtree_->Branch(label.str().c_str(),&energySR_[iL][iSR]);
	label.str("");
	label << "subtractedenergy_" << iL << "_SR" << iSR;
	outtree_->Branch(label.str().c_str(),&subtractedenergySR_[iL][iSR]);
	label.str("");
	label << "maxhitEoutside_" << iL << "_SR" << iSR;
	outtree_->Branch(label.str().c_str(),&maxhitEoutside_[iL][iSR]);
      }
      for (unsigned iy(0);iy<3;++iy){
	for (unsigned ix(0);ix<3;++ix){
	  unsigned idx = 3*iy+ix;
	  label.str("");     
	  label << "E_" << iL << "_" << idx;
	  outtree_->Branch(label.str().c_str(),&Exy_[iL][idx]);
	}
      }
    }

    //outtree_->SetCacheSize(100000000);
    //outtree_->SetCacheLearnEntries(100);
    //tree_ptr->LoadTree(evt);

    p_rawEtotal = new TH1F("p_rawEtotal", "Total E (MIP)", 5000,0,200000);
    p_wgtEtotal = new TH1F("p_wgtEtotal", "Total weighted E (MIP)",5000, 0, 200000);

    //p_rawESR.resize(nSR_,0);
    p_wgtESR.resize(nSR_,0);
    //p_rawSubtractESR.resize(nSR_,0);
    p_wgtSubtractESR.resize(nSR_,0);

    p_EsumfracvsdR.resize(nLayers_,0);
    p_EvsdR.resize(nLayers_,0);
    p_dR.resize(nLayers_,0);
    for (unsigned iL(0); iL<nLayers_;++iL){
      label.str("");
      label << "EsumfracvsdR_" << iL;
      p_EsumfracvsdR[iL] = new TH2F(("p_"+label.str()).c_str(),";dR (mm);Esum/Etot;events", 2000,0,200,1000,0,1);
      p_EsumfracvsdR[iL]->StatOverflows();
      label.str("");
      label << "EvsdR_" << iL;
      p_EvsdR[iL] = new TH2F(("p_"+label.str()).c_str(),";dR (mm);Ehit (mips);events", 2000,0,200,5000,0,5000);
      p_EvsdR[iL]->StatOverflows();   
      label.str("");
      label << "dR_" << iL;
      p_dR[iL] = new TH1F(("p_"+label.str()).c_str(),";dR (mm);E-weighted events", 2000,0,200);
      p_dR[iL]->StatOverflows();   
    }


    for (unsigned iSR(0);iSR<nSR_;++iSR){
      //label.str("");
      //label << "rawESR" << iSR;
      //p_rawESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR} (MIPs);events", 5000,0,200000);
      //p_rawESR[iSR]->StatOverflows();

      label.str("");
      label << "wgtESR" << iSR;
      p_wgtESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR} (MIPs);events", 5000,0,200000);
      p_wgtESR[iSR]->StatOverflows();
      
      //label.str("");
      //label << "rawSubtractESR" << iSR;
      //p_rawSubtractESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR}^{PUsubtr} (MIPs);events", 5000,0,200000);
      //p_rawSubtractESR[iSR]->StatOverflows();
      label.str("");
      label << "wgtSubtractESR" << iSR;
      p_wgtSubtractESR[iSR] = new TH1F(("p_"+label.str()).c_str(),";E_{SR}^{PUsubtr} (MIPs);events", 5000,0,200000);
      p_wgtSubtractESR[iSR]->StatOverflows();
    }//loop on sr

}

void SignalRegion::fillHistograms(){
        
  outputFile_->cd(outputDir_.c_str());

  p_rawEtotal->Fill(totalE_);
  p_wgtEtotal->Fill(wgttotalE_);

  for (unsigned iSR(0);iSR<nSR_;++iSR){
    //Fill energy without PU subtraction
    bool subtractPU = false;
    //p_rawESR[iSR]->Fill( getEtotalSR(iSR, subtractPU));

    double wgtESR = 0;
    for(unsigned iL(0); iL < nLayers_;iL++){
      wgtESR += getSR(iSR, iL, subtractPU);
    }
    p_wgtESR[iSR]->Fill(wgtESR);

    //Fill energy after PU subtraction
    subtractPU = true;
    //p_rawSubtractESR[iSR]->Fill( getEtotalSR(iSR, subtractPU));
    double wgtSubtractESR = 0;
    for(unsigned iL(0); iL < nLayers_;iL++){
      wgtSubtractESR += getSR(iSR, iL, subtractPU);
    }
    p_wgtSubtractESR[iSR]->Fill(wgtSubtractESR);

  }//loop on SR

}




