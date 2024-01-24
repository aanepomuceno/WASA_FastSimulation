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
//----------------------------------------------------------------------------
/// \file WASADetectorParametrisation.cc
// Based on Geant4 WASA parametrization example
// Andre Nepomuceno - Winter 2023
////---------------------------------------------------------------------------

#include "WASADetectorParametrisation.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
G4double p1 = 1.0; //probability for sorting double gaussian

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASADetectorParametrisation::WASADetectorParametrisation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASADetectorParametrisation::~WASADetectorParametrisation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WASADetectorParametrisation::GetResolution( Detector aDetector, 
                                                      Parametrisation aParam, 
                                                      G4double aMomentum, G4double prob ) {

G4double res = 1.0;
//-------------------------------------------------------------------------- 
//Generic Detector (based on different detectors and test beam data)
  if ( aParam == eGENERIC ) {
    aMomentum /= GeV;  //aMomentum must be in GeV
    switch ( aDetector ) {
      case WASADetectorParametrisation::eEMCAL :
           res = 0.056/std::sqrt(aMomentum) + 0.011;
           break;
      case WASADetectorParametrisation::eHCAL :
           res = std::sqrt( std::pow( 0.51/std::sqrt( aMomentum ),2) + std::pow( 0.07, 2 ) );
           break;
      case WASADetectorParametrisation::eTRACKER :
           res = 0.013;
           break;
    }
   }
 //--------------------------------------------------------------------------
 //NNBAR detector (parametrization from Full Sim)
   else if ( aParam == eNNBAR ) { 
        aMomentum /= MeV;  //aMomentum must be in MeV
        if (aDetector == WASADetectorParametrisation::eEMCAL ) {

         if (aMomentum < 90.) p1 = 0.83;     
         else p1 = 0.73;
          
           if (prob < p1) {
             res = (0.01014472*aMomentum + 0.72389385)/aMomentum; 
             }
           else {
             res = (0.02066737*aMomentum + 1.11823963)/aMomentum;
             }
       }
   
      if (aDetector == WASADetectorParametrisation::eHCAL ) {
       res = 1.0/aMomentum;
       }
       
      if (aDetector == WASADetectorParametrisation::eTRACKER ) {
       res = 1.0;
       }

   }
 //--------------------------------------------------------------------------  
//WASA Detector
   else if ( aParam == eWASA ) {
    aMomentum /= GeV;  //aMomentum must be in GeV
    switch ( aDetector ) {
      case WASADetectorParametrisation::eEMCAL :
           res = 0.05/std::sqrt(aMomentum); 
           break;
      case WASADetectorParametrisation::eHCAL :
           res = 1.0/aMomentum;                   
           break;
      case WASADetectorParametrisation::eTRACKER :
           res = 1.0;
           break;
    }
  }

  return res;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WASADetectorParametrisation::GetMedian( Detector aDetector, 
                                                      Parametrisation aParam,
                                                      G4double aMomentum,G4double prob ) {
    G4double med = 1.0;
    if ( aParam == eNNBAR ) { 
       aMomentum /= MeV;  //aMomentum in MeV
      if (aDetector == WASADetectorParametrisation::eEMCAL ) {

         if (aMomentum < 90.) p1 = 0.83;
         else p1 = 0.73;

          if (prob < p1) {
             med =  (0.86014766*aMomentum + 0.40650026)/aMomentum;
             }
          else {
             med = (0.83779033*aMomentum - 0.67570678)/aMomentum;
             }
        }
   
     if (aDetector == WASADetectorParametrisation::eHCAL ) {
         med = 1.0;
       }

     if (aDetector == WASADetectorParametrisation::eTRACKER ) {
         med = 1.0;
       }

   }
  return med;
 }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WASADetectorParametrisation::GetEfficiency( Detector aDetector, 
                                                      Parametrisation /*aParam*/,
                                                      G4double /*aMomentum*/ ) {
  // For the time being, we set the efficiency to 1.0
  G4double eff = 1.0;
  switch ( aDetector ) {
    case WASADetectorParametrisation::eTRACKER :
      eff = 1.0;
      break;
    case WASADetectorParametrisation::eEMCAL :
      eff = 1.0;
      break;
    case WASADetectorParametrisation::eHCAL :
      eff = 1.0;
      break;
  }
  return eff;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

