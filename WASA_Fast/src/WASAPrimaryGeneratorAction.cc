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
/// \file WASAPrimaryGeneratorAction.cc
/// \brief Implementation of the WASAPrimaryGeneratorAction class

#include "WASAPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "WASAPrimaryParticleInformation.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#define pi 3.14159265

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAPrimaryGeneratorAction::WASAPrimaryGeneratorAction() {
  G4int n_particle = 1;
  fParticleGPS = new G4GeneralParticleSource();
  //fParticleGPS = new G4ParticleGun();
 // fParticleGPS->SetNumberOfParticles(n_particle);

  //G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //G4String particleName;
  //G4ParticleDefinition* particle = 
    //particleTable->FindParticle( particleName = "gamma" );
  //fParticleGPS->SetParticleDefinition( particle );

  //fParticleGPS->SetParticleMomentumDirection( G4ThreeVector( 0.0, 1.0, 1.0 ) );
  //fParticleGPS->SetParticleEnergy( 0.5*GeV );
  //fParticleGPS->SetParticlePosition( G4ThreeVector( 0.0, 0.0, 0.0 ) );
  //fParticleGPS->SetParticleTime(0.0*ns);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAPrimaryGeneratorAction::~WASAPrimaryGeneratorAction() {
  delete fParticleGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WASAPrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent ) {

//Added for random direction

  //Primiary verices along the a disc of radius 50cm
  //G4double x; G4double y; G4double z;
  //G4double radius = 0.25*G4UniformRand();
  //G4double  phi = CLHEP::twopi*G4UniformRand();
  //x = radius*std::cos(phi)*m; 
  //y = radius*std::sin(phi)*m; 
  //z = 0.*m; 

  //G4double cosTheta = 2*G4UniformRand() - 1.;
  //G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  //G4double ux = sinTheta*std::cos(phi),
    //       uy = sinTheta*std::sin(phi),
      //     uz = cosTheta;
  //fParticleGPS->SetParticlePosition(G4ThreeVector(0, 0, 0));
  //fParticleGPS->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
  fParticleGPS->GeneratePrimaryVertex(anEvent);

  // Loop over the vertices, and then over primary particles,
  // and for each primary particle create an info object, in
  // which to store "Monte Carlo true" information.
  // This approach could appear unnecessarily heavy in the present case
  // of a trivial particle gun generator, but it is useful in the more
  // realistic case of a Monte Carlo event generator like Pythia8.
  G4int count_particles = 0;
  for ( G4int ivtx = 0; ivtx < anEvent->GetNumberOfPrimaryVertex(); ivtx++ ) {
    for ( G4int ipp = 0; ipp < anEvent->GetPrimaryVertex( ivtx )->GetNumberOfParticle();
          ipp++ ) {
      G4PrimaryParticle* primary_particle = 
        anEvent->GetPrimaryVertex( ivtx )->GetPrimary( ipp );
      if ( primary_particle ) {
        primary_particle->SetUserInformation( new WASAPrimaryParticleInformation( 
          count_particles, primary_particle->GetPDGcode(), 
          primary_particle->GetMomentum() ) );
        count_particles++;              
      }
    } 
  }
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4ParticleGun* WASAPrimaryGeneratorAction::GetParticleGun() {
G4GeneralParticleSource* WASAPrimaryGeneratorAction::GetParticleGPS() {
  return fParticleGPS;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

