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
/// \file WASAActionInitialization.cc
/// \brief Implementation of the WASAActionInitialization class

#include "WASAActionInitialization.hh"
#include "WASAPrimaryGeneratorAction.hh"
#include "WASARunAction.hh"
#include "WASAEventAction.hh"
#include "WASATrackingAction.hh"
#include "G4UIcommand.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAActionInitialization::WASAActionInitialization() : 
  G4VUserActionInitialization(), fFileName( "WASAFastOutput" ), fSmear( true ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAActionInitialization::WASAActionInitialization( const G4String aOutName, 
                                                      const G4String aSmear ) :
  G4VUserActionInitialization(), fFileName( aOutName ), 
  fSmear( G4UIcommand::ConvertToBool( aSmear ) ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAActionInitialization::WASAActionInitialization( const G4String aOutName ) :
  G4VUserActionInitialization(), fFileName( aOutName ), fSmear( true ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASAActionInitialization::~WASAActionInitialization() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WASAActionInitialization::BuildForMaster() const {
  SetUserAction( new WASARunAction( fFileName ) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WASAActionInitialization::Build() const {
  SetUserAction( new WASAPrimaryGeneratorAction );
  SetUserAction( new WASARunAction( fFileName ) );
  SetUserAction( new WASAEventAction( fSmear ) );
  SetUserAction( new WASATrackingAction );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

