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
/// \file WASASmearer.cc
// Based on Geant4 WASA parametrization example
// Andre Nepomuceno - Wintr  2026


#include "WASASmearer.hh"
#include "WASAPrimaryParticleInformation.hh"
#include "G4PrimaryParticle.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include <ctime>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASASmearer* WASASmearer::fWASASmearer = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASASmearer::WASASmearer() {
  time_t seed = time( NULL );
  fRandomEngine = new CLHEP::HepJamesRandom( static_cast< long >( seed ) );
  fRandomGauss = new CLHEP::RandGauss( fRandomEngine );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASASmearer::~WASASmearer() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASASmearer* WASASmearer::Instance() {
  if ( ! fWASASmearer ) {
    fWASASmearer = new WASASmearer();
  }
  return fWASASmearer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector WASASmearer::SmearMomentum( const G4Track* aTrackOriginal, 
                                           G4double aResolution ) {
  return SmearGaussian( aTrackOriginal, aResolution );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WASASmearer::SmearEnergy( const G4Track* aTrackOriginal, 
                                    G4double aResolution, G4double aMedian, G4double Kenergy ) {
  G4double newE = -1.0;
  while ( newE < 0.0 ) {  // To ensure that the resulting value is not negative
                          // (vital for energy smearing, does not change direction
                          // for momentum smearing)
    if ( aResolution != -1.0 ) {
        newE = Kenergy * Gauss( aMedian, aResolution );
    } else {
      newE = aTrackOriginal->GetKineticEnergy();
    }
  }
  return newE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector WASASmearer::SmearGaussian( const G4Track* aTrackOriginal, 
                                           G4double aResolution ) {
  G4ThreeVector originP = aTrackOriginal->GetMomentum();
  G4ThreeVector originPos = aTrackOriginal->GetPosition();
  G4double rdm = Gauss( 1.0, aResolution );
  G4ThreeVector smearedMom( originP.x()*rdm, originP.y()*rdm, originP.z()*rdm );
  return smearedMom;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WASASmearer::Gauss( G4double aMean, G4double aStandardDeviation ) {
  return fRandomGauss->fire( aMean, aStandardDeviation );
}

//Returns a dE/dx from Bethe-Bloch function - WVS - 2025-07-08
G4double WASASmearer::BetheBloch( G4ParticleDefinition* aParticle, G4double aKinectE, PotentialGas aPotentialGas) {
  const G4double c = 1;
  const G4double z = 1;
  const G4double Me = 0.511;
  const G4double PathLen = 10; //cm;
  const G4double cBetheBloch = 0.307075;
  
  //For Ar/CO2 90/10 by default
  G4double ExPotencialGas = 1.77e-4;
  G4double EffZA = 0.455;
  G4double GasDensity = 1.70e-3;

  switch (aPotentialGas) {
    case eArCO2_9010: {  //To initial composition 90/10
       ExPotencialGas = 1.77e-4;
       EffZA = 0.455;
       GasDensity = 1.70e-3;
       break;

    }
    case eArCO2_8020: {  //composition 80/20
       //ExPotencialGas = 1.674e-4;
       ExPotencialGas = 1.5e-4;
       EffZA = 0.46;
       GasDensity = 1.773e-3;
       break;

    }
    case eAir: {  //to confirm values for Air
       ExPotencialGas = 1.0;
       EffZA = 1.0;
       GasDensity = 1.0;
       break;

    }
  }

  //results nan for all
  G4double KinectE;
  G4double MassParticle;
  G4double TotalE;
  G4double FLorentz;
  G4double Beta;
  G4double Tmax;
  G4double dEdx;

  aKinectE /= MeV; //MeV 
  KinectE = aKinectE;   
  MassParticle = aParticle->GetPDGMass();
  TotalE = MassParticle + KinectE;
  FLorentz = TotalE/MassParticle;
  Beta = std::sqrt((1-(1/std::pow(FLorentz,2))));
  Tmax = (2 * Me * std::pow(c,2) * std::pow(Beta,2) * std::pow(FLorentz,2)) 
         / (1 + (2*FLorentz*(Me/MassParticle)) + std::pow(Me/MassParticle,2));
  dEdx = cBetheBloch * std::pow(z,2) * (EffZA/std::pow(Beta,2)) 
         * (0.5*std::log((2 * Me * std::pow(c,2) * std::pow(Beta,2) * std::pow(FLorentz,2) * Tmax) / std::pow(ExPotencialGas,2)) 
         - 2 * std::pow(Beta,2));
  dEdx = dEdx * GasDensity; //include gas density

  return dEdx;  //no pathlength applied
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

