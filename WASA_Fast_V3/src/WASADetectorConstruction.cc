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
/// \brief Implementation of the WASADetectorConstruction class
/// Based on Geant4 Par02 parameterisations example
//  Adapted by B Meirose - Winter 2025

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
   
  // G4Material* Air = new G4Material("Air", density=1.290*mg/cm3, ncomponents=2);
  // G4Element* Nitrogen = FindOrBuildElement("N"); 
  // G4Element* Oxygen = FindOrBuildElement("O");   

   //Air->AddElement(Nitrogen, 70.0*perCent); // 70% Nitrogen by mass
   //Air->AddElement(Oxygen, 30.0*perCent);   // 30% Oxygen by mass
   G4NistManager* nist = G4NistManager::Instance();
   G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");

  
//Artifical box (for debbuging)
//-------------------------------------------------------------------------------------------------------------   
   G4Box* solidWorld = new G4Box("World", 1.0*m, 1.0*m, 1.0*m); // Huge world box
   G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Air, "World");
   G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, true);
//-------------------------------------------------------------------------------------------------------------

//Test air cilinder

    G4double innerRadius = 0.0*cm;
    G4double outerRadius = 33.4*cm;
    G4double hz = 52.8*cm;  //HALF length in Z
    G4double startAngle = 0.*deg;
    G4double spanningAngle = 360.*deg;

   
   //-----------Load the GDML file-------------------------------------
   G4GDMLParser parser;
   //parser.Read("my_detector.gdml");
    parser.Read("filtered_barrel.gdml"); // Replace with your actual GDML file path


//Geometry tests
//-------------------------------------------------------------------------------------------------------------
//   G4VPhysicalVolume* physDetector = parser.GetWorldVolume(); // This is the ROOT "world"
 //  G4LogicalVolume* logicGDMLDetector = physDetector->GetLogicalVolume();
  // new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicGDMLDetector, "Detector", logicWorld, false, 0, true);
//-----------------------------------------------------------------------------------------------------------------

//Get world from GDML
  G4VPhysicalVolume* worldPV = parser.GetWorldVolume();
  G4LogicalVolume* logicGDMLDetector = worldPV->GetLogicalVolume();
//Air cilinder
//------------------------------------------------------
  //auto tubeS = new G4Tubs("tubeS",innerRadius,outerRadius,hz,startAngle,spanningAngle);
  //auto tubeLV = new G4LogicalVolume(tubeS,Air,"tubeLV");
  //new G4PVPlacement(0,G4ThreeVector(0., 0., 3.*cm),tubeLV,"tubePV",logicGDMLDetector,false,0,true);
//------------------------------------------------------------

   // Retrieve the world volume from the GDML
   //auto worldPV = parser.GetWorldVolume();

   // Here, you could still set up additional parameters, regions, or cuts if necessary
   
 
   for (int i = 0; i <= 16; i++) {
        std::ostringstream sec_modulo;
        sec_modulo << "SECE" << i;  // construct SECEX, X=0,1,...,16.
        G4LogicalVolume* secLogVol = parser.GetVolume(sec_modulo.str());
        if (secLogVol) {
           fECalList.push_back( new G4Region( secLogVol->GetName() ) );
           fECalList.back()->AddRootLogicalVolume( secLogVol );
           G4cout << "Found volume: " << secLogVol->GetName() << G4endl;
           
        } else {
            G4cout << "Volume " << sec_modulo.str() << " not found!" << G4endl;
        }
   }
   
    G4LogicalVolume* secLogVol = parser.GetVolume("MOTHER1");
    fECalList.push_back( new G4Region( secLogVol->GetName() ) );
    fECalList.back()->AddRootLogicalVolume( secLogVol );
    G4cout << "Found volume: " << secLogVol->GetName() << G4endl;
   
   

//Try get whole worl volume
//  G4LogicalVolume* secLogVol = worldPV->GetLogicalVolume();
 // fECalList.push_back( new G4Region( secLogVol->GetName() ) );
 // fECalList.back()->AddRootLogicalVolume( secLogVol );


   //G4LogicalVolume* secLogVol = parser.GetVolume("SECE0");
   //if (secLogVol) {
    //G4cout << "Found SEC logical volume: " << secLogVol->GetName() << G4endl;
//} else {
  //  G4cout << "SEC logical volume not found!" << G4endl;
//}

   // return physWorld;
   return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WASADetectorConstruction::ConstructSDandField() {

  // G4RegionStore* regionStore = G4RegionStore::GetInstance();
   //G4Region* caloRegion = regionStore->GetRegion("EM_calo_region");

  for ( G4int iterECal = 0; iterECal < G4int( fECalList.size() ); iterECal++ ) {
    // Bound the fast simulation model for the electromagnetic calorimeter
    // to all the corresponding Geant4 regions
    WASAFastSimModelEMCal* fastSimModelEMCal
      = new WASAFastSimModelEMCal( "fastSimModelEMCal", fECalList[ iterECal ],
                                    WASADetectorParametrisation::eWASA );
                                     
    // Register the fast simulation model for deleting
    G4AutoDelete::Register(fastSimModelEMCal);
  }
  
  // if (caloRegion) {
  //     WASAFastSimModelEMCal* fastSimModelEMCal
      //     = new WASAFastSimModelEMCal("fastSimModelEMCal", caloRegion,
        //                                WASADetectorParametrisation::eWASA);
       //G4AutoDelete::Register(fastSimModelEMCal);
   //} else {
    //   G4cerr << "Warning: EM_calo_region not found! FastSimModel will not be applied." << G4endl;
   //}

   // Add global magnetic field
   G4ThreeVector fieldValue = G4ThreeVector(0., 0., 0.);
   fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
   fMagFieldMessenger->SetVerboseLevel(1);
}

