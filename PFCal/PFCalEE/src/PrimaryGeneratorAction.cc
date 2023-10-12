//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "HepMCG4AsciiReader.hh"
#include "HepMCG4PythiaInterface.hh"

#define PI 3.1415926535

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(G4int mod, double eta)
{
  model_ = mod;
  eta_ = eta;
  G4int n_particle = 1;

  // default generator is particle gun.
  currentGenerator= particleGun= new G4ParticleGun(n_particle);
  currentGeneratorName= "particleGun";
  hepmcAscii= new HepMCG4AsciiReader();
#ifdef G4LIB_USE_PYTHIA
  pythiaGen= new HepMCG4PythiaInterface();
#else
  pythiaGen= 0;
#endif
  gentypeMap["particleGun"]= particleGun;
  gentypeMap["hepmcAscii"]= hepmcAscii;
  gentypeMap["pythia"]= pythiaGen;

  Detector = (DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();
 
  //create a messenger for this class
  gunMessenger = new PrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="e-");
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticleEnergy(10.*GeV);
  G4double position = -0.5*(Detector->GetWorldSizeZ());
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,position));
  
  //G4cout << " -- Gun position set to: 0,0," << position << G4endl;

  rndmFlag = "off";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete hepmcAscii;
  delete pythiaGen;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  // 
  G4double x0 = 0.*cm, y0 = 0.*cm;
  G4double z0 = -0.5*(Detector->GetWorldSizeZ());

  switch(model_) {
  case 2:
    //smear within 1cm...
    z0 = (G4RandGauss::shoot(0.,5.))*cm; break;
  case 4:
    x0 = (G4RandFlat::shoot(0.,460)-170)*mm; break;
    //y0 = (G4RandFlat::shoot(0.,10.)-5)*mm;
  case 3:
    x0 = (G4RandFlat::shoot(0.,250.)-125)*mm;
    //start from near the bottom in y to have enough space for all layers when shooting with an angle...
    y0 = (G4RandFlat::shoot(0.,200.)-250)*mm;
    break;
  case 5:
    //x0 = (G4RandFlat::shoot(0.,30)-15)*mm;
    //y0 = (G4RandFlat::shoot(0.,30.)-15)*mm;
    //update to cover full hexagon
    //size=6.4mm
    x0 = (G4RandFlat::shoot(0.,1.)-0.5)*mm;
    y0 = (G4RandFlat::shoot(0.,1.)-0.5)*mm;
    break;
  default: break;
  }

  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  G4cout << " -- Gun position set to: " << x0 << "," << y0 << "," << z0 << G4endl;

  if (model_ == 2) {
    G4double theta0 = 2*atan(exp(-1*eta_));
    G4double phi0 = (G4RandFlat::shoot(0.,2*PI));
    particleGun->SetParticleMomentumDirection(G4ThreeVector(cos(phi0)*sin(theta0), sin(phi0)*sin(theta0), cos(theta0))); 
  } else {
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));
  }
  
  if(currentGenerator){
    currentGenerator->GeneratePrimaryVertex(anEvent);
  }
  else
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries",
                "PrimaryGeneratorAction001", FatalException,
                "generator is not instanciated." );


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

