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
#include "SteppingAction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PbBoxConstruction(G4LogicalVolume* logicWorld)
{
	//
// Pb shell
//

// Material ---> Pb
	G4NistManager* manager = G4NistManager::Instance();
	G4Material* matPb1 = manager->FindOrBuildMaterial("G4_Pb");

	G4Box* boxIn =
		new G4Box("Box", 20 * mm, 20 * mm, 75 * mm);

	G4Box* boxOut =
		new G4Box("Box", 70 * mm, 70 * mm, 125 * mm);
	//new G4Box("Box", 18 * mm, 18 * mm, 72 * mm);

	G4SubtractionSolid* solidPbShell = new G4SubtractionSolid("PbShell", boxOut, boxIn);

	G4LogicalVolume*
		logicPbShell = new G4LogicalVolume(solidPbShell, // its solid
			matPb1,              // its material
			"PbShell");           // its name 

	G4VPhysicalVolume* physiPbShell =
		new G4PVPlacement(0,                  // no rotation
			G4ThreeVector(0, 0, -81/*-36*/ * mm),    // at (0,0,0)
			logicPbShell,   // its logical volume
			"PbShell",              // its name
			logicWorld,         // its mother volume
			false,              // no boolean operation
			0);                 // copy number
}

void collimatorConstruction(G4LogicalVolume* logicWorld)
{
	// Material ---> Pb
	G4NistManager* manager = G4NistManager::Instance();
	G4Material* matPb1 = manager->FindOrBuildMaterial("G4_Pb");

	G4Box* PbBoxOut =
		new G4Box("Box", 38 * mm, 38 * mm, 45 * mm);

	G4Box* PbBoxIn =
		new G4Box("Box", 23 * mm, 23 * mm, 50 * mm);

	G4ThreeVector Pbtranslation = G4ThreeVector(0. * cm, 0. * cm, -20. * mm);

	G4SubtractionSolid* solidPbShellwithouthole = new G4SubtractionSolid("PbShell", PbBoxOut, PbBoxIn, 0, Pbtranslation);

	G4Tubs* PbHole =
		new G4Tubs("Crystal", 0, 5 * mm, 100 * mm, 0, CLHEP::twopi);

	G4SubtractionSolid* solidPbShell = new G4SubtractionSolid("PbShell", solidPbShellwithouthole, PbHole);

	G4LogicalVolume*
		logicPbShell = new G4LogicalVolume(solidPbShell, // its solid
			matPb1,              // its material
			"PbShell");           // its name 

	G4VPhysicalVolume* physiPbShell =
		new G4PVPlacement(0,                  // no rotation
			G4ThreeVector(0, 0, -34/*-36*/ * mm),    // at (0,0,0)
			logicPbShell,   // its logical volume
			"PbShell",              // its name
			logicWorld,         // its mother volume
			false,              // no boolean operation
			0);                 // copy number

	// Material ---> Al
	G4Material* matAl = manager->FindOrBuildMaterial("G4_Al");

	G4Box* AloutBoxOut =
		new G4Box("Box", 43 * mm, 43 * mm, 49 * mm);

	G4Box* AloutBoxIn =
		new G4Box("Box", 38 * mm, 38 * mm, 45 * mm);

	G4ThreeVector AloutTranslation = G4ThreeVector(0. * cm, 0. * cm, -1. * mm);

	G4SubtractionSolid* solidAloutShellwithouthole = new G4SubtractionSolid("AloutShell", AloutBoxOut, AloutBoxIn, 0, AloutTranslation);

	G4Tubs* AloutHole =
		new G4Tubs("Alout", 0, 5 * mm, 100 * mm, 0, CLHEP::twopi);
	
	G4Box* AloutBackHole =
		new G4Box("Box", 18.2 * mm, 18.2 * mm, 25 * mm);

	G4SubtractionSolid* solidAloutShell = new G4SubtractionSolid("Alout", solidAloutShellwithouthole, AloutHole);

	solidAloutShell = new G4SubtractionSolid("Alout", solidAloutShell, AloutBackHole, 0, G4ThreeVector(0. * cm, 0. * cm, -50. * mm));

	G4LogicalVolume*
		logicAloutShell = new G4LogicalVolume(solidAloutShell, // its solid
			matAl,              // its material
			"AloutShell");           // its name 

	G4VPhysicalVolume* physiAloutShell =
		new G4PVPlacement(0,                  // no rotation
			G4ThreeVector(0, 0, -34/*-36*/ * mm),    // at (0,0,0)
			logicAloutShell,   // its logical volume
			"AloutShell",              // its name
			logicWorld,         // its mother volume
			false,              // no boolean operation
			0);                 // copy number

	// Material ---> Al

	G4Box* AlinBoxOut =
		new G4Box("Box", 23 * mm, 23 * mm, 37.5 * mm);

	G4Box* AlinBoxIn =
		new G4Box("Box", 18 * mm, 18 * mm, 37.5 * mm);

	G4ThreeVector AlinTranslation = G4ThreeVector(0. * cm, 0. * cm, -5. * mm);

	G4SubtractionSolid* solidAlinShellwithouthole = new G4SubtractionSolid("AlinShell", AlinBoxOut, AlinBoxIn, 0, AlinTranslation);

	G4Tubs* AlinHole =
		new G4Tubs("Alin", 0, 5 * mm, 100 * mm, 0, CLHEP::twopi);

	G4SubtractionSolid* solidAlinShell = new G4SubtractionSolid("Alin", solidAlinShellwithouthole, AlinHole);

	G4LogicalVolume*
		logicAlinShell = new G4LogicalVolume(solidAlinShell, // its solid
			matAl,              // its material
			"AlinShell");           // its name 

	G4VPhysicalVolume* physiAlinShell =
		new G4PVPlacement(0,                  // no rotation
			G4ThreeVector(0, 0, -41.5/*-36*/ * mm),    // at (0,0,0)
			logicAlinShell,   // its logical volume
			"AlinShell",              // its name
			logicWorld,         // its mother volume
			false,              // no boolean operation
			0);                 // copy number

	// Material ---> Al

	G4Box* AlBackBoxOut =
		new G4Box("Box", 23.8 * mm, 23.8 * mm, 35.5 * mm);

	G4Box* AlBackBoxIn =
		new G4Box("Box", 17.8 * mm, 17.8 * mm, 100 * mm);

	G4ThreeVector AlBackTranslation = G4ThreeVector(0. * cm, 0. * cm, 0. * mm);

	G4SubtractionSolid* solidAlBackShell = new G4SubtractionSolid("AlBackShell", AlBackBoxOut, AlBackBoxIn, 0, AlBackTranslation);

	G4LogicalVolume*
		logicAlBackShell = new G4LogicalVolume(solidAlBackShell, // its solid
			matAl,              // its material
			"AlBackShell");           // its name 

	G4VPhysicalVolume* physiAlBackShell =
		new G4PVPlacement(0,                  // no rotation
			G4ThreeVector(0, 0, -118.5/*-36*/ * mm),    // at (0,0,0)
			logicAlBackShell,   // its logical volume
			"AlBackShell",              // its name
			logicWorld,         // its mother volume
			false,              // no boolean operation
			0);                 // copy number
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //
  // World volume
  //   
  
  // Material ---> Vacuum  
  G4Material* Vacuum =
  G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  
  // Full sphere shape
  G4double solidWorld_rmax = 100*cm;
  G4Orb*
  solidWorld = new G4Orb("World",                          // its name
		          solidWorld_rmax);                // its size 

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
  // Pb shell
  //

  // Material ---> Pb
  G4NistManager* manager = G4NistManager::Instance();
  //G4Element* Pb = manager->FindOrBuildElement("G4_Pb");
  //G4Material* matPb = manager->FindOrBuildMaterial("G4_Pb");
  //G4double density1 = 11.3437 * g / cm3;
  //G4double a1 = 207.2 * g / mole;
  //G4Material* Pb =
	 // new G4Material("Pb", 82 , a1, density1);

  //G4Box* box11 =
	 // new G4Box("Box", 38 * mm, 38 * mm, 10 * mm);
  //G4Tubs* cy1 =
	 // new G4Tubs("Cylinder", 0, 5 * mm, 10.1 * mm, 0, CLHEP::twopi);
  //G4SubtractionSolid* solidPbShell = new G4SubtractionSolid("PbShell", box11, cy1);
  //
  //G4LogicalVolume*
  //logicPbshell = new G4LogicalVolume(solidPbShell, // its solid
		//                         matPb,              // its material
		//			 "PbSehll");           // its name 
  
//  G4VPhysicalVolume* physiPbShell = 
        //             new G4PVPlacement(0,                  // no rotation
		      //                 G4ThreeVector(),    // at (0,0,0)
					   //logicPbshell,   // its logical volume
				    //   "PbShell",              // its name
				    //   logicWorld,         // its mother volume
				    //   false,              // no boolean operation
				    //   0);                 // copy number

//
// Fe box
//

// Material ---> Fe
	//G4NistManager* manager = G4NistManager::Instance();
  G4Material* matFe = manager->FindOrBuildMaterial("G4_Fe");

  G4Box* solidFeBox =
	  new G4Box("Box", 8 * mm, 5 * mm, 15 * mm);

  G4LogicalVolume*
	  logicFeBox = new G4LogicalVolume(solidFeBox, // its solid
		  matFe,              // its material
		  "FeBox");           // its name 

  //  G4VPhysicalVolume* physiPbShell = 
  //new G4PVPlacement(0,                  // no rotation
	 // G4ThreeVector(-5 * mm, -2.5 * mm, -91 * mm),    // at (0,0,0)
	 // logicFeBox,   // its logical volume
	 // "FeBox",              // its name
	 // logicWorld,         // its mother volume
	 // false,              // no boolean operation
	 // 0);                 // copy number

//
// Al shell
//

// Material ---> Al
	//G4NistManager* manager = G4NistManager::Instance();
	G4Material* matAl = manager->FindOrBuildMaterial("G4_Al");

	G4Box* box31 =
		new G4Box("Box", 15 * mm, 15 * mm, 69/*24*/ * mm);
		//new G4Box("Box", 15 * mm, 15 * mm, 69 * mm);
	G4Box* box32 =
		new G4Box("Box", 18 * mm, 18 * mm, 72/*27*/ * mm);
		//new G4Box("Box", 18 * mm, 18 * mm, 72 * mm);
	G4SubtractionSolid* solidAlShell = new G4SubtractionSolid("AlShell", box32, box31);

	G4LogicalVolume*
		logicAlShell = new G4LogicalVolume(solidAlShell, // its solid
			matAl,              // its material
			"AlShell");           // its name 

	//  G4VPhysicalVolume* physiPbShell = 
	new G4PVPlacement(0,                  // no rotation
		G4ThreeVector(0,0,-81/*-36*/ * mm),    // at (0,0,0)
		logicAlShell,   // its logical volume
		"AlShell",              // its name
		logicWorld,         // its mother volume
		false,              // no boolean operation
		0);                 // copy number

//
// Plastic shell
//

// Material ---> plastic
	G4Material* matPET = manager->FindOrBuildMaterial("G4_POLYCARBONATE");

	G4Box* box21 =
		new G4Box("Box", 5.1 * mm, 2.6 * mm, 10.1 * mm);
	G4Box* box22 =
		new G4Box("Box", 10 * mm, 5 * mm, 10 * mm);
	G4SubtractionSolid* solidPETShell = new G4SubtractionSolid("PETShell", box22, box21);
		
	G4LogicalVolume*
		logicPETShell = new G4LogicalVolume(solidPETShell, // its solid
			matPET,              // its material
			"PETShell");           // its name 

	//  G4VPhysicalVolume* physiPbShell = 
	new G4PVPlacement(0,                  // no rotation
		G4ThreeVector(0, 0, -34 * mm),    // at (0,0,0)
		logicPETShell,   // its logical volume
		"PETShell",              // its name
		logicWorld,         // its mother volume
		false,              // no boolean operation
		0);                 // copy number

	//second one
	G4Box* solidPETShell2 =
		new G4Box("Box", 14.1 * mm, 1 * mm, 51 * mm);

	G4LogicalVolume*
		logicPETShell2 = new G4LogicalVolume(solidPETShell2, // its solid
			matPET,              // its material
			"PETShell2");           // its name 

	//  G4VPhysicalVolume* physiPbShell = 
	new G4PVPlacement(0,                  // no rotation
		G4ThreeVector(0, -9 * mm, -71 * mm),    // at (0,0,0)
		logicPETShell2,   // its logical volume
		"PETShell2",              // its name
		logicWorld,         // its mother volume
		false,              // no boolean operation
		0);                 // copy number




  //
  // source
  //
  
  // Material ---> Na22Cl
  G4double density2 = 0.97 * g / cm3;
  G4double a2 = 21.994 * g / mole;
  G4Material* Na22 =
	  new G4Material("Pb", 11 , a2, density2);
  G4Tubs* SolidNa22Source =
	  new G4Tubs("Na22Source", 0, 1.5 * mm, 0.15 * mm, 0, CLHEP::twopi);

  G4Element* Na = new G4Element("Sodium", "Na", 11., 22 * g / mole);
  G4Element* Cl = new G4Element("Chlorine", "Cl", 17., 35.5 * g / mole);

  G4Material* NaCl = new G4Material("NaCl", 2.165 * g / cm3, 2);  //density should be checked!!!
  NaCl->AddElement(Na, 1);
  NaCl->AddElement(Cl, 1);
  
  G4LogicalVolume*
  logicNa22Source = new G4LogicalVolume(SolidNa22Source,               // its solid
		                  NaCl,                     // ite material
				 "Na22Source");                  // its name
  
//  G4VPhysicalVolume*  physiCell = 
	   //   new G4PVPlacement(0,                         // no rotation
		  //              G4ThreeVector(0*mm,0*mm,-5*mm),           // at(0,0,0)
				//logicNa22Source,                 // its logical volume
				//"Na22Source",                    // its name
			 // logicWorld,          // its mother volume
				//false,                     // no boolean operation
				//0);                        // copy number

  //
  // Crystal  
  //

  // Material ---> CdZnTe
  G4String symbol;
  G4double a;                                              // atomic mass
  G4double z;                                              // atomic number 
  G4double density;
  G4int ncomponents, natoms;	  

  G4Element* Cd = new G4Element("Cadmium",    symbol= "Cd", z= 48., a= 116.00*g/mole);
  G4Element* Zn = new G4Element("Zinc", symbol= "Zn", z= 30., a= 65.38*g/mole);
  G4Element* Te = new G4Element("Tellurium",     symbol= "Te" , z= 52. , a= 127.6*g/mole);

  G4Material* CdZnTe = new G4Material("CdZnTe",  density=4.78*g/cm3, ncomponents=3);  //density should be checked!!!
  CdZnTe->AddElement(Cd, natoms=9);
  CdZnTe->AddElement(Zn, natoms=2);
  CdZnTe->AddElement(Te, natoms=9);

  //G4Material* matCdZnTe = manager->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE"); 

  G4Box*
	  solidCrystal = new G4Box("Crystal",                                                  // its name
		  5 * mm, 2.5 * mm, 5 * mm);        // its size
		 

  G4LogicalVolume* 
  logicCrystal = new G4LogicalVolume(solidCrystal,                                     // its solid
		                   CdZnTe,                                             // its material
	                          "Crystal");                                          // its name 

//  G4VPhysicalVolume*  physiCrystal =
                new G4PVPlacement(0,                                                  // no rotation
                                  G4ThreeVector(0,0,-34*mm),                                     // at (0,0,0)
                                  logicCrystal,                                        // its logical volume
                                  "Crystal",                                           // its name
                                  logicWorld,                                           // its mother  volume
                                  false,                                               // no boolean operation
                                  0);                                                  // copy number
   

	//collimatorConstruction(logicWorld);

  //
  //always return the physical World
  //  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
