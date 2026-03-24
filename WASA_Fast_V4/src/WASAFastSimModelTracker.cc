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
/// \file WASAFastSimModelTracker.cc
/// Based on Geant4 Par02 parameterisations example
//  Andre Nepomuceno - February 2026


#include "WASAFastSimModelTracker.hh"
#include "WASAEventInformation.hh"
#include "WASAPrimaryParticleInformation.hh"
#include "WASASmearer.hh"
#include "WASAOutput.hh"

#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

#include "Randomize.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4TransportationManager.hh"
#include "G4PathFinder.hh"
#include "G4FieldTrack.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAFastSimModelTracker::WASAFastSimModelTracker( G4String aModelName, 
  G4Region* aEnvelope, WASADetectorParametrisation::Parametrisation aType ) :
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( aType ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAFastSimModelTracker::WASAFastSimModelTracker( G4String aModelName, 
                                                    G4Region* aEnvelope ) :
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( WASADetectorParametrisation::eWASA ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAFastSimModelTracker::WASAFastSimModelTracker( G4String aModelName ) :
  G4VFastSimulationModel( aModelName ), fCalculateParametrisation(),
  fParametrisation( WASADetectorParametrisation::eWASA ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAFastSimModelTracker::~WASAFastSimModelTracker() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WASAFastSimModelTracker::IsApplicable( const G4ParticleDefinition& 
                                                                   aParticleType ) {
  return aParticleType.GetPDGCharge() != 0;  // Applicable for all charged particles
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WASAFastSimModelTracker::ModelTrigger( const G4FastTrack& aFastTrack ) {
   return true;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WASAFastSimModelTracker::DoIt( const G4FastTrack& aFastTrack,
                                     G4FastStep& aFastStep ) {

//  G4cout << " ________Tracker model triggered _________" << G4endl;
// Kill particles with kinetic energy below 10 keV (no track in TPC in this case)
   G4double KE = aFastTrack.GetPrimaryTrack()->GetKineticEnergy()/MeV;
   if (KE < 0.01 ) {
     aFastStep.KillPrimaryTrack();
     return;
   }
  
   G4int pdgID = aFastTrack.GetPrimaryTrack()-> GetDefinition()->GetPDGEncoding();
   G4double time = aFastTrack.GetPrimaryTrack()->GetGlobalTime();

//Calculate the dE/dx for the given paticle with kinetic energy KE
   G4double dEdx_rho = WASASmearer::Instance()->
                       BetheBloch(aFastTrack.GetPrimaryTrack()-> GetDefinition(), 
                       KE, WASASmearer::eArCO2_8020);

// Calculate the final position (at the outer boundary of the tracking detector)
// of the particle with the momentum at the entrance of the tracking detector.
  G4Track track = * aFastTrack.GetPrimaryTrack();
  G4FieldTrack aFieldTrack( '0' );
  G4FieldTrackUpdator::Update( &aFieldTrack, &track );
  G4double retSafety = -1.0;
  ELimited retStepLimited;
  G4FieldTrack endTrack( 'a' );
  G4double currentMinimumStep = 10.0*m;  // not the real pathlengh
  G4PathFinder* fPathFinder = G4PathFinder::GetInstance();
  auto* nav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  
  fPathFinder->ComputeStep( aFieldTrack,
                            currentMinimumStep,
                            0,
                            aFastTrack.GetPrimaryTrack()->GetCurrentStepNumber(),
                            retSafety,
                            retStepLimited,
                            endTrack,
                            aFastTrack.GetPrimaryTrack()->GetVolume() );
 
  //Calculate the energy the particle *would* deposite given the path length, considering it crosses the whole TPC
  //If the particle energy is very small, it does not leave the TPC.
  G4ThreeVector posInitial = aFastTrack.GetPrimaryTrack()->GetPosition();
  G4ThreeVector dirInitial = aFastTrack.GetPrimaryTrack()->GetMomentumDirection();
  G4ThreeVector dir = dirInitial.unit();
  G4ThreeVector posFinal = endTrack.GetPosition();
  G4ThreeVector track_vec = (posFinal - posInitial);
  G4double pathLength = (posFinal - posInitial).mag()/cm;
  posFinal += dir * 1*um;
  nav->LocateGlobalPointAndSetup(posFinal, 0, true);                
  aFastStep.ProposePrimaryTrackFinalPosition( posFinal );         
  
   G4double res = fCalculateParametrisation->
        GetResolution( WASADetectorParametrisation::eTRACKER, 
                       fParametrisation, KE, pdgID );    

  G4double dEdx_sm = dEdx_rho*WASASmearer::Instance()->Gauss(1.0,res);
  G4double Edep = ( dEdx_sm*pathLength )/MeV;
  G4double DeltaKE = std::max( (KE - Edep), 0.0 ) ;
  aFastStep.ProposePrimaryTrackFinalKineticEnergy(DeltaKE);
  
   WASAOutput::Instance()->SaveTrack( WASAOutput::eSaveTracker,
                                             0,
                                             pdgID,
                                             dEdx_rho,
                                             posInitial/mm,
                                             track_vec/mm,
                                             res,
                                             1.0,
                                             dEdx_sm,                                         
                                             time/ps);
 
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

