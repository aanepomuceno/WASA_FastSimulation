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
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Gamma.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4Proton.hh"

#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"

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
         &aParticleType == G4MuonPlus::Definition()    ||
         &aParticleType == G4MuonMinus::Definition()   ||
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
  bool insideDetector = true;
  G4ThreeVector entrance_pos = G4ThreeVector(0,0,0);
  G4double dEdx_rho = 5.63079; //energy_loss/cm for muons
  G4double pathLength = 0;
  G4double Edep_muon = 0.0;
  
  pdgID = aFastTrack.GetPrimaryTrack()-> GetDefinition()->GetPDGEncoding();
  G4double KE = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();
  auto* nav = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  
  if ( abs(pdgID) != 13) {
 // Kill the parameterised particle at the initial point
        aFastStep.KillPrimaryTrack();
        aFastStep.ProposePrimaryTrackPathLength( 0.0 );
        G4ThreeVector posInitial = aFastTrack.GetPrimaryTrack()->GetPosition();
        G4ThreeVector dirInitial = aFastTrack.GetPrimaryTrack()->GetMomentumDirection();
        G4ThreeVector pos = posInitial;
        G4ThreeVector dir = dirInitial.unit();
        entrance_pos = pos;
        G4VPhysicalVolume* pv = nav->LocateGlobalPointAndSetup(pos);
       
        //The FastSim model is triggered when the particle enter the "Empty_Area", so we make it hit the detector and record
        // the entrance position 
        if ( pv->GetName() == "Empty_Area") {
           G4double retSafety = -1.0;
           G4double step = nav->ComputeStep(pos, dir, DBL_MAX, retSafety);
           aFastStep.ProposePrimaryTrackPathLength(step);
           pos += dir * step + dir * 1*um;
           nav->LocateGlobalPointAndSetup(pos, 0, true);
           aFastStep.ProposePrimaryTrackFinalPosition(pos);
           entrance_pos = pos;
        }
       
       //If the particles does not hit the detector, we do not record its position  
        pv = nav->LocateGlobalPointAndSetup(pos);
        if (pv->GetName() == "MOTHER_World_PV") insideDetector = false;
        
  }
 
//This part is specific for cosmic muons, that can cross the detector
//-------------------------------------------------------------------------------------------------------------------------
  if ( abs(pdgID) == 13) {
    
        G4ThreeVector posInitial = aFastTrack.GetPrimaryTrack()->GetPosition();
        G4ThreeVector dirInitial = aFastTrack.GetPrimaryTrack()->GetMomentumDirection();
            
        G4ThreeVector pos = posInitial;
        G4ThreeVector dir = dirInitial.unit();
        nav->LocateGlobalPointAndSetup(pos);
        //Due the WASA geometry, the muon will cross different layers of the detecor with "Empty_Area" between them, 
        //so we track the particle and restart the loop every time it hit an "Empty_Area", till it leave the detector
        while (true) {
    
           G4double retSafety = -1.0;
           G4double step = nav->ComputeStep(pos, dir, DBL_MAX, retSafety);
           aFastStep.ProposePrimaryTrackPathLength(step);
           pos += dir * step + dir * 1*um;
           nav->LocateGlobalPointAndSetup(pos, 0, true);
           G4VPhysicalVolume* pv = nav->LocateGlobalPointAndSetup(pos);
           //G4cout << "Position After Step: " << pos << G4endl;
           //G4cout << " Step Size : " << step << G4endl;
           //if (pv != NULL) G4cout << "volume= " << pv->GetName() << G4endl;
           pathLength = step;
           if (pv->GetName() == "Empty_Area") {
              entrance_pos = pos;
              break;
           }
          
           if (pv->GetName() == "MOTHER_World_PV") {
              aFastStep.ProposePrimaryTrackPathLength( 0.0 );
              aFastStep.ProposeTrackStatus(fStopAndKill);
              insideDetector = false;
              break;
           }
        }
   
        Edep_muon = ( dEdx_rho*(pathLength/cm) )/MeV;
  
        if (KE < Edep_muon) {
           aFastStep.KillPrimaryTrack();
           aFastStep.ProposePrimaryTrackPathLength( 0.0 );
           aFastStep.ProposePrimaryTrackFinalPosition(pos);
        }
    
        else {
           G4double DeltaKE = KE - Edep_muon;
           KE = Edep_muon;        
           aFastStep.ProposePrimaryTrackFinalPosition(pos);
           aFastStep.ProposePrimaryTrackFinalKineticEnergy(DeltaKE);
        }
 
  }
//End of cosmic muon specific part  
//-------------------------------------------------------------------------------------------------------------------------
 
  if (insideDetector ) {

        G4ThreeVector Pos = entrance_pos;
        G4double time = aFastTrack.GetPrimaryTrack()->GetGlobalTime();
          
        WASAEventInformation* info = (WASAEventInformation*) 
                            G4EventManager::GetEventManager()->GetUserInformation();
                            
        if ( info->GetDoSmearing() ) {
            // Smearing according to the electromagnetic calorimeter resolution taken from WASADetectorParametrisation
            G4ThreeVector Porg = aFastTrack.GetPrimaryTrack()->GetMomentum();
            G4double res = fCalculateParametrisation->GetResolution( 
               WASADetectorParametrisation::eEMCAL, fParametrisation, KE ,pdgID ); 
    
            G4double med = fCalculateParametrisation->GetMedian( 
               WASADetectorParametrisation::eEMCAL, fParametrisation, KE,pdgID ); 
     

            G4double eff = fCalculateParametrisation->GetEfficiency( 
               WASADetectorParametrisation::eEMCAL, fParametrisation, Porg.mag() );

            G4double Esm;
            Esm = std::abs( WASASmearer::Instance()->
                        SmearEnergy( aFastTrack.GetPrimaryTrack(), res, med, KE ) );


            //Save histogram and trees
            WASAOutput::Instance()->FillHistogram( 1, (Esm/MeV) / (KE/MeV) );
  
      
            WASAOutput::Instance()->SaveTrack( WASAOutput::eSaveEMCal,
                                             0,
                                             pdgID,
                                             KE/MeV,
                                             Pos/mm,
                                             res,
                                             eff,
                                             Esm/MeV,                                         
                                             time/ns);

      // The (smeared) energy of the particle is deposited in the step
            aFastStep.ProposeTotalEnergyDeposited( Esm );
        } else {
        // No smearing: simply setting the value of KE
        // The (initial) energy of the particle is deposited in the step
        // (which corresponds to the entrance of the electromagnetic calorimeter)
            aFastStep.ProposeTotalEnergyDeposited( KE );
          }
  }//if insideDetector
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

