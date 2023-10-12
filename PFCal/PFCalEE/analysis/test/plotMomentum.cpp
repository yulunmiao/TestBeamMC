#include <iostream>
#include <sstream>
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

using namespace HepMC;

int main(int argc, char ** argv) {

  std::cout << " -- Starting program..." << std::endl;

  const unsigned nZ = 100;
  const unsigned z0 = 100;
  const unsigned step = 2;

  //TFile *fout = TFile::Open("PLOTS/output_hepmc_vtxorig_MB.root","RECREATE");
  //TFile *fout = TFile::Open("PLOTS/output_hepmc_vtxmodif.root","RECREATE");
  TFile *fout = TFile::Open("PLOTS/output_hepmc_hggmodif.root","RECREATE");
  fout->cd();
  TH1F *hvtx_x = new TH1F("hvtx_x",";x (mm);vertices",1000,-10,10);
  TH1F *hvtx_y = new TH1F("hvtx_y",";y (mm);vertices",1000,-10,10);
  TH1F *hvtx_z = new TH1F("hvtx_z",";z (mm);vertices",1000,-500,500);
  TH1F *hvtx_t = new TH1F("hvtx_t",";t (ns);vertices",20000,-100,100);
  TH2F *hvtx_tvsz = new TH2F("hvtx_tvsz",";z (mm); t (ps);vertices",1000,-200,200,1200,-600,600);
  TH1F *hvtx_t_z[nZ];
  for (unsigned i(0);i<nZ;++i){
    std::ostringstream label;
    int zmin = -1*z0+i*step;
    int zmax = -1*z0+(i+1)*step;
    label << "hvtx_t_z_" << zmin << "_" << zmax;
    hvtx_t_z[i] = new TH1F(label.str().c_str(),";t (ps);vertices",120,-600,600);
  }
  TH1F *hProton = new TH1F("hProton",";p (GeV);particles",1000,0,100);
  TH1F *hNeutron = new TH1F("hNeutron",";p (GeV);particles",1000,0,100);
  TH1F *hPipm = new TH1F("hPipm",";p (GeV);particles",1000,0,100);
  TH1F *hProtonLog = new TH1F("hProtonLog",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *hNeutronLog = new TH1F("hNeutronLog",";log(p) (log(GeV));particles",100,-2,2);
  TH1F *hPipmLog = new TH1F("hPipmLog",";log(p) (log(GeV));particles",100,-2,2);

  const unsigned nFiles = 1;//6;

  for (unsigned i(0);i<nFiles;++i){

    // specify an input file
    std::ostringstream lname;
    //lname << "/afs/cern.ch/work/p/pdauncey/public/Pythia140305_";
    //if (i==0) lname << "000000";
    //else if (i<10) lname << "00000" << i;
    //else if (i<100) lname << "0000" << i;
    //else if (i<1000) lname << "000" << i;
    //else if (i<10000) lname << "00" << i;
    //else if (i<100000) lname << "0" << i;
    //else lname << i;
    //lname << ".dat";
    //lname << "/afs/cern.ch/work/a/amagnan/public/HepMCFiles/ggHgg_origVtx.dat";
    //lname << "/afs/cern.ch/work/a/amagnan/public/HepMCFiles/ggHgg_modifyVtx.dat";
    //lname << "/afs/cern.ch/work/a/amagnan/public/HepMCFiles/vertexHLLHC.dat";
    //lname << "/afs/cern.ch/work/p/pdauncey/public/Pythia140305_000000.dat";
    lname << "/afs/cern.ch/work/a/amagnan/public/HepMCFiles/ggHgg_1428658356.dat";
    HepMC::IO_GenEvent ascii_in(lname.str().c_str(),std::ios::in);
    
    // get the first event
    HepMC::GenEvent* evt = ascii_in.read_next_event();
    // loop until we run out of events
    while ( evt ) {
      unsigned ievt =  evt->event_number();
      if (ievt%10000==0) std::cout << "Processing Event Number "
				 << ievt
				 << std::endl;


      GenVertex * parent = 0;
      HepMC::GenEvent::vertex_const_iterator q = evt->vertices_begin();

      //for (; q != evt->vertices_end(); ++q ){

	//if ((*q)->position().x()!=0 || (*q)->position().y()!=0) continue;

      double z = (*q)->position().z();
      double t = (*q)->position().t();
      hvtx_x->Fill((*q)->position().x());
      hvtx_y->Fill((*q)->position().y());
      hvtx_z->Fill(z);
      hvtx_t->Fill(t);
      //hvtx_tvsz->Fill(z,t*1000);
      hvtx_tvsz->Fill(z,t);

      if (fabs(z)<z0){
	unsigned idx = static_cast<unsigned>((z+z0)*1./step);
	if (idx>(nZ-1)) continue;
	hvtx_t_z[idx]->Fill(t*1000);
      }

	/*std::cout << " -- vtx pos: " << (*q)->position().x() << " " << (*q)->position().y() << " " << (*q)->position().z() << " nParticles: in=" << (*q)->particles_in_size() << " " << (*q)->particles_out_size()
		  << std::endl;
	for ( HepMC::GenVertex::particles_in_const_iterator p
		= (*q)->particles_in_const_begin(); p != (*q)->particles_in_const_end(); ++p ){

	  std::cout << " ---- in particle " << (*p)->pdg_id() << " status " << (*p)->status()
		    << std::endl;
	  
	}
 	for ( HepMC::GenVertex::particles_out_const_iterator p
		= (*q)->particles_out_const_begin(); p != (*q)->particles_out_const_end(); ++p ){

	  std::cout << " ---- out particle " << (*p)->pdg_id() << " status " << (*p)->status()
		    << std::endl;
	  
	}*/
	//}

      /*
      // analyze the event
      HepMC::GenEvent::particle_const_iterator lPart = evt->particles_begin();
      //unsigned counter = 0;
      
      //std::cout << " -- Number of particles: " << evt->particles_size()
      //<< std::endl;
      
      for (; lPart!=evt->particles_end();++lPart){
	//std::cout << counter << " " 
	//<< (*lPart)->pdg_id() <<  " " 
	//<< (*lPart)->status() 
	//<< std::endl;
	//counter++;
	if ((*lPart)->status()!=1) continue;
	double p = sqrt(pow((*lPart)->momentum().px(),2)+pow((*lPart)->momentum().py(),2)+pow((*lPart)->momentum().pz(),2));
	//if (fabs((*lPart)->momentum().eta())>2.8 && fabs((*lPart)->momentum().eta())<3.0){
	if (fabs((*lPart)->momentum().eta())<2.5){
	  if (fabs((*lPart)->pdg_id())==2212) {hProton->Fill(p);hProtonLog->Fill(log10(p));}
	  else if (fabs((*lPart)->pdg_id())==2112) {hNeutron->Fill(p);hNeutronLog->Fill(log10(p));}
	  else if (fabs((*lPart)->pdg_id())==211) {hPipm->Fill(p);hPipmLog->Fill(log10(p));}
	}
      }
      */

      // delete the created event from memory
      delete evt;
      // read the next event
      ascii_in >> evt;
      ievt++;
    }
  }

  fout->Write();
  
  return 0;

}//main




