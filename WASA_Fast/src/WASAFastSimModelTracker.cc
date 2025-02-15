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
/// \brief Implementation of the WASAFastSimModelTracker class

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

G4bool WASAFastSimModelTracker::ModelTrigger( const G4FastTrack& /*aFastTrack*/ ) {
  return true;  // No kinematical restrictions to apply the parametrisation
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WASAFastSimModelTracker::DoIt( const G4FastTrack& aFastTrack,
                                     G4FastStep& aFastStep ) {

  G4cout << " ________Tracker model triggered _________" << G4endl;

  // Calculate the final position (at the outer boundary of the tracking detector)
  // of the particle with the momentum at the entrance of the tracking detector.

  G4Track track = * aFastTrack.GetPrimaryTrack();
  G4FieldTrack aFieldTrack( '0' );
  G4FieldTrackUpdator::Update( &aFieldTrack, &track );

  G4double retSafety = -1.0;
  ELimited retStepLimited;
  G4FieldTrack endTrack( 'a' );
  G4double currentMinimumStep = 10.0*m;  // Temporary: change that to sth connected
                                         // to particle momentum.
  G4PathFinder* fPathFinder = G4PathFinder::GetInstance();
  /*G4double lengthAlongCurve = */ 
  fPathFinder->ComputeStep( aFieldTrack,
                            currentMinimumStep,
                            0,
                            aFastTrack.GetPrimaryTrack()->GetCurrentStepNumber(),
                            retSafety,
                            retStepLimited,
                            endTrack,
                            aFastTrack.GetPrimaryTrack()->GetVolume() );

  // Place the particle at the tracking detector exit 
  // (at the place it would reach without the change of its momentum).
  aFastStep.ProposePrimaryTrackFinalPosition( endTrack.GetPosition() );

  // Consider only primary tracks (do nothing else for secondary charged particles)
  G4ThreeVector Porg = aFastTrack.GetPrimaryTrack()->GetMomentum();
  if ( ! aFastTrack.GetPrimaryTrack()->GetParentID() ) {
    WASAEventInformation* info = (WASAEventInformation*) 
                            G4EventManager::GetEventManager()->GetUserInformation();
    if ( info->GetDoSmearing() ) {
      // Smearing according to the tracking detector resolution
      G4double res = fCalculateParametrisation->
        GetResolution( WASADetectorParametrisation::eTRACKER, 
                       fParametrisation, Porg.mag() );
      G4double eff = fCalculateParametrisation->
        GetEfficiency( WASADetectorParametrisation::eTRACKER,
                       fParametrisation, Porg.mag() );
      G4ThreeVector Psm;
      Psm = WASASmearer::Instance()->
                                 SmearMomentum( aFastTrack.GetPrimaryTrack(), res );
      WASAOutput::Instance()->FillHistogram( 0, ((Psm.mag()/MeV) / (Porg.mag()/MeV)) );
      // Setting the values of Psm, res and eff
      ( (WASAPrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetTrackerMomentum( Psm );
      ( (WASAPrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetTrackerResolution( res );
      ( (WASAPrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetTrackerEfficiency( eff );
    } else {
      // No smearing: simply setting the value of Porg
      ( (WASAPrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetTrackerMomentum( Porg );
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

