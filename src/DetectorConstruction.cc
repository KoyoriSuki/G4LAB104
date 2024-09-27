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
/// \file radioactivedecay/rdecay01/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SteppingAction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "FridgeConstruction.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
    messenger = new DetectorMessenger(this);
    std::cout << shieldingHeight << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    std::cout << "DetectorConstruction::Construct() Called" << std::endl;
  //
  // World volume
  //   
  
  // Material ---> Vacuum  
  G4Material* Vacuum =
  G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  
  // Full sphere shape
  G4double solidWorld_lengthmax = 1000 * cm;
  G4Box*
      solidWorld = new G4Box("World",                          // its name
          solidWorld_lengthmax, solidWorld_lengthmax, solidWorld_lengthmax);                // its size 

  G4LogicalVolume*                         
  logicWorld = new G4LogicalVolume(solidWorld,             // its solid
                                   Vacuum,                 // its material
                                   "World");               // its name
  G4VPhysicalVolume*                                   
  physiWorld = new G4PVPlacement(0,                        // no rotation
                                 G4ThreeVector(),          // at (0,0,0)
                                 logicWorld,               // its logical volume
                                 "World",                  // its name
                                 0,                        // its mother  volume
                                 false,                    // no boolean operation
                                 0);                       // copy number

  //
  // Crystal  
  //

  // 1. Define the materials
  G4String symbol;
  G4double a;                                              // atomic mass
  G4double z;                                              // atomic number 
  G4double density;
  G4int ncomponents, natoms;

  // Elements
  G4Element* Cd = new G4Element("Cadmium", symbol = "Cd", z = 48., a = 116.00 * g / mole);
  G4Element* Mo = new G4Element("Molybdenum", symbol = "Mo", z = 42., a = 100.0 * g / mole);
  G4Element* O = new G4Element("Oxygen", symbol = "O", z = 8., a = 16.0 * g / mole);
  G4Element* Li = new G4Element("Lithium", symbol = "Li", z = 3., a = 6.941 * g / mole);

  // Material ---> CdMoO_4
  G4Material* CdMoO4 = new G4Material("CdMoO4", density = 6.2 * g / cm3, ncomponents = 3);  //density should be checked!!!
  CdMoO4->AddElement(Cd, natoms = 1);
  CdMoO4->AddElement(Mo, natoms = 1);
  CdMoO4->AddElement(O, natoms = 4);

  // Material ---> Li_2MoO_4
  G4Material* Li2MoO4 = new G4Material("Li2MoO4", density = 2.66 * g / cm3, ncomponents = 3);  //density should be checked!!!
  Li2MoO4->AddElement(Li, natoms = 2);
  Li2MoO4->AddElement(Mo, natoms = 1);
  Li2MoO4->AddElement(O, natoms = 4);

  // 2. place the crystal
  if (crystalType == 0)
  {
      G4Tubs* SolidCrystal =
          new G4Tubs("Crystal", 0, 12.5 * mm, 22.5 * mm, 0, CLHEP::twopi);


      G4LogicalVolume*
          logicCrystal = new G4LogicalVolume(SolidCrystal,                                     // its solid
              CdMoO4,                                             // its material
              "Crystal");                                          // its name 

      //  G4VPhysicalVolume*  physiCrystal =
      new G4PVPlacement(0,                                                  // no rotation
          G4ThreeVector(0, 0, -391.5 * mm),                                     // at (0,0,0)
          logicCrystal,                                        // its logical volume
          "Crystal",                                           // its name
          logicWorld,                                           // its mother  volume
          false,                                               // no boolean operation
          0);                                                  // copy number
  }
  else if (crystalType == 1)
  {
      G4Box* SolidCrystal =
          new G4Box("Crystal", 22.5 * mm, 22.5 * mm, 22.5 * mm);


      G4LogicalVolume*
          logicCrystal = new G4LogicalVolume(SolidCrystal,                                     // its solid
              Li2MoO4,                                             // its material
              "Crystal");                                          // its name 

      //  G4VPhysicalVolume*  physiCrystal =
      new G4PVPlacement(0,                                                  // no rotation
          G4ThreeVector(0, 0, -450 * mm),
          logicCrystal,                                        // its logical volume
          "Crystal",                                           // its name
          logicWorld,                                           // its mother  volume
          false,                                               // no boolean operation
          0);                                                  // copy number
  }
  else if (crystalType == 2)
  {
      G4Box* SolidCrystal =
          new G4Box("Crystal", 5 * mm, 5 * mm, 5 * mm);


      G4LogicalVolume*
          logicCrystal = new G4LogicalVolume(SolidCrystal,                                     // its solid
              Li2MoO4,                                             // its material
              "Crystal");                                          // its name 

      //  G4VPhysicalVolume*  physiCrystal =
      new G4PVPlacement(0,                                                  // no rotation
          G4ThreeVector(0, 0, -460 * mm),
          logicCrystal,                                        // its logical volume
          "Crystal",                                           // its name
          logicWorld,                                           // its mother  volume
          false,                                               // no boolean operation
          0);                                                  // copy number
  }

  //
  // Refridgerator and shielding construction
  //

  Layer0InnerConstruction(logicWorld);
  Layer0OuterConstruction(logicWorld);
  Layer1Construction(logicWorld);
  Layer2Construction(logicWorld);
  Layer3Construction(logicWorld);
  Layer4Construction(logicWorld);
  Layer5Construction(logicWorld);
  Bucket1Construction(logicWorld);
  Bucket2Construction(logicWorld);
  Bucket3Construction(logicWorld);
  Shield4Construction(logicWorld);

  //
  // shielding construction
  // 
  
  //CopperHolderConstruction(logicWorld);
  if (shieldingFlag == 1)
  {
      //LeadShieldConstruction(logicWorld, shieldingHeight * cm);
      CuShieldConstruction(logicWorld, shieldingHeight * cm, shieldingHeight2 * cm);
      //CuInnerShieldConstruction(logicWorld);
      SpringPlateConstruction(logicWorld);
  }
  //CuTopShieldConstruction(logicWorld, shieldingHeight);
  //CuShieldConstruction(logicWorld);
  //WaterShieldConstruction(logicWorld);
  //SubstanceForCosmicRayConstruction(logicWorld);
  
  //
  //always return the physical World
  //  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
