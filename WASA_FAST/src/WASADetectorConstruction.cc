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
/// \brief Implementation of the WASADetectorConstruction class
/// Based on Geant4 Par02 parameterisations example
//  Adapted by A Nepomuceno - Winter 2023

#include "WASADetectorConstruction.hh"
#include "G4ProductionCuts.hh"
#include "G4SystemOfUnits.hh"
#include "G4RegionStore.hh"
#include "G4GDMLParser.hh"
#include "G4AutoDelete.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ProductionCuts.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4RegionStore.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASADetectorConstruction::WASADetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WASADetectorConstruction::~WASADetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WASADetectorConstruction::Construct() {
  
   G4double a;
   G4double z;
   G4double density;
   G4int fCheckOverlaps = 0;
   G4NistManager* nistManager = G4NistManager::Instance();
   

//-----------Vacuum-------------------------------------------
  
   new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
   auto defaultMaterial = G4Material::GetMaterial("Galactic");
   
//-----------SEC EMCAL----------------------------------
     nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
     auto csI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
   
/*
//-----------Scintilator HCAL---------------------------------
   //Elements for scintilator (Material will be named Scint)
   G4Element* elH = nistManager->FindOrBuildElement("H");
   G4Element* elC = nistManager->FindOrBuildElement("C");

   //Scintilator Material
   G4Material* Scint = new G4Material("Scint", 1.023*g/cm3, 2);
   Scint->AddElement(elH, 0.524573);
   Scint->AddElement(elC, 1 - 0.524573);
   auto scintMaterial = G4Material::GetMaterial("Scint");
  
  */  
   
    //SEC Dimensions
    G4double innerRadius = 32.5*cm;
    G4double outerRadius = 63.5*cm;
    G4double hz = 54.5*cm;  //HALF length in Z
    G4double startAngle = 0.*deg;
    G4double spanningAngle = 360.*deg;

    G4double worldSizeXY = 2.5 * outerRadius;
    G4double worldSizeZ  = 2.5 * hz;

//-----------Build World Volume--------------------------------------

   auto worldS = new G4Box("World", worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); 
   auto worldLV = new G4LogicalVolume(worldS,defaultMaterial,"World");
   auto worldPV  = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

//-----------Build SEC Calorimeter---------------------------------------------------

   auto secS = new G4Tubs("SECS",innerRadius,outerRadius,hz,startAngle,spanningAngle);
   auto secLV = new G4LogicalVolume(secS,csI,"SECLV");
   new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),secLV,"SECPV",worldLV,false,0,fCheckOverlaps);

//-------Makes SEC volume becoming a G4Region---------------------

   G4Region* caloRegion = new G4Region("EM_calo_region");
   caloRegion->AddRootLogicalVolume(secLV);
//   caloRegion->SetProductionCuts(new G4ProductionCuts());
//   caloRegion->GetProductionCuts()->SetProductionCut( 0.1*m, idxG4GammaCut );
//   caloRegion->GetProductionCuts()->SetProductionCut( 1.0*cm, idxG4GammaCut );

    return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WASADetectorConstruction::ConstructSDandField() {
	
  
   G4RegionStore* regionStore = G4RegionStore::GetInstance();
   G4Region* caloRegion = regionStore->GetRegion("EM_calo_region");
  
  WASAFastSimModelEMCal* fastSimModelEMCal
      = new WASAFastSimModelEMCal( "fastSimModelEMCal", caloRegion,
                                    WASADetectorParametrisation::eWASA );
  // Register the EM fast simulation model for deleting
    G4AutoDelete::Register(fastSimModelEMCal);
    
 // Add global magnetic field
  G4ThreeVector fieldValue = G4ThreeVector(0.,0.,0.);
  fMagFieldMessenger = new G4GlobalMagFieldMessenger( fieldValue );
  fMagFieldMessenger->SetVerboseLevel(1);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

