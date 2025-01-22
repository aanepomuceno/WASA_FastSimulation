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
// Based on Geant4 WASA parametrization example
// Andre Nepomuceno - Winter 2023

#include "WASAFastSimModelEMCal.hh"
#include "WASAEventInformation.hh"
#include "WASAPrimaryParticleInformation.hh"
#include "WASASmearer.hh"
#include "WASAOutput.hh"

#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

#include "G4PathFinder.hh"
#include "G4FieldTrack.hh"
#include "G4FieldTrackUpdator.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAFastSimModelEMCal::WASAFastSimModelEMCal( G4String aModelName, 
  G4Region* aEnvelope, WASADetectorParametrisation::Parametrisation aType ) :
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( aType ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAFastSimModelEMCal::WASAFastSimModelEMCal( G4String aModelName, 
                                                G4Region* aEnvelope ) : 
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( WASADetectorParametrisation::eWASA ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAFastSimModelEMCal::WASAFastSimModelEMCal( G4String aModelName ) :
  G4VFastSimulationModel( aModelName ), fCalculateParametrisation(), 
  fParametrisation( WASADetectorParametrisation::eWASA ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAFastSimModelEMCal::~WASAFastSimModelEMCal() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WASAFastSimModelEMCal::IsApplicable( 
  const G4ParticleDefinition& aParticleType ) {
  // Applicable for electrons, positrons, and gammas
  return &aParticleType == G4Electron::Definition()    ||
         &aParticleType == G4Positron::Definition()    ||
         &aParticleType == G4Gamma::Definition()       ||
         &aParticleType == G4PionMinus::Definition()   ||
         &aParticleType == G4PionPlus::Definition()    ||
         &aParticleType == G4Proton::Definition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WASAFastSimModelEMCal::ModelTrigger( const G4FastTrack& /*aFastTrack*/ ) {
  return true;  // No kinematical restrictions to apply the parametrisation
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WASAFastSimModelEMCal::DoIt( const G4FastTrack& aFastTrack,
                                   G4FastStep& aFastStep ) {
  //G4cout << " ________EMCal model triggered _________" << G4endl;

  G4int pdgID = 0;
  pdgID = aFastTrack.GetPrimaryTrack()-> GetDefinition()->GetPDGEncoding();
    
  if ( abs(pdgID) != 13) {
 // Kill the parameterised particle at the entrance of the electromagnetic calorimeter
  aFastStep.KillPrimaryTrack();
  aFastStep.ProposePrimaryTrackPathLength( 0.0 );
  }

//This part is specifi for cosmic muons crossing the detector
//--------------------------------------------------
  if ( abs(pdgID) == 13) {
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
 }
//End of cosmic muon specific part
//------------------------------------------------------------------------------------------

  G4double KE = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();
  G4ThreeVector Pos = aFastTrack.GetPrimaryTrack()->GetPosition();
  G4double time = aFastTrack.GetPrimaryTrack()->GetGlobalTime();
  
  WASAEventInformation* info = (WASAEventInformation*) 
                            G4EventManager::GetEventManager()->GetUserInformation();
                            
    if ( info->GetDoSmearing() ) {
      // Smearing according to the electromagnetic calorimeter resolution taken from DetectorParametrisation
      G4ThreeVector Porg = aFastTrack.GetPrimaryTrack()->GetMomentum();
      G4double res = fCalculateParametrisation->GetResolution( 
               WASADetectorParametrisation::eEMCAL, fParametrisation, KE ,pdgID ); //p->pdgID
    
      G4double med = fCalculateParametrisation->GetMedian( 
               WASADetectorParametrisation::eEMCAL, fParametrisation, KE,pdgID ); //p->pdgID
     // std::cout << med << std::endl;//

      G4double eff = fCalculateParametrisation->GetEfficiency( 
               WASADetectorParametrisation::eEMCAL, fParametrisation, Porg.mag() );

      G4double Esm;
      Esm = std::abs( WASASmearer::Instance()->
                        SmearEnergy( aFastTrack.GetPrimaryTrack(), res, med ) );


   //Save histogram and trees
      WASAOutput::Instance()->FillHistogram( 1, (Esm/MeV) / (KE/MeV) );
  
      //pdgID = aFastTrack.GetPrimaryTrack()-> GetDefinition()->GetPDGEncoding();
      WASAOutput::Instance()->SaveTrack( WASAOutput::eSaveEMCal,
                                         0,
                                         pdgID,
                                         Pos/mm,
                                         res,
                                         eff,
                                         Esm/MeV,
                                         KE/MeV,
                                         time/ns);

      // The (smeared) energy of the particle is deposited in the step
      // (which corresponds to the entrance of the electromagnetic calorimeter)
      aFastStep.ProposeTotalEnergyDeposited( Esm );
    } else {
      // No smearing: simply setting the value of KE
      // The (initial) energy of the particle is deposited in the step
      // (which corresponds to the entrance of the electromagnetic calorimeter)
      aFastStep.ProposeTotalEnergyDeposited( KE );
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

