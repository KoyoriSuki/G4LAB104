#include "FridgeConstruction.hh"

void Layer0InnerConstruction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	//
	//Layer1 Plate
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4Tubs* SolidLayer0Plate =
		new G4Tubs("Layer0Plate", 0, 322.5 * mm, 10.0 * mm, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicLayer0Plate = new G4LogicalVolume(SolidLayer0Plate,               // its solid
			matCu,                     // ite material
			"Layer0Plate");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, 476 * mm),           // at(0,0,0)
		LogicLayer0Plate,                 // its logical volume
		"Layer0Plate",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number
}

void Layer0OuterConstruction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	//
	//Layer1 Plate
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4Tubs* SolidLayer0Plate =
		new G4Tubs("Layer0Plate", 322.5, 480 * mm, 20.0 * mm, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicLayer0Plate = new G4LogicalVolume(SolidLayer0Plate,               // its solid
			matCu,                     // ite material
			"Layer0Plate");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, 446 * mm),           // at(0,0,0)
		LogicLayer0Plate,                 // its logical volume
		"Layer0Plate",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number
}

void Layer1Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	//
	//Layer1 Plate
	// 
	
	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4Tubs* SolidLayer1Plate =
		new G4Tubs("Layer1Plate", 0, 238 * mm, 3.994 * mm, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicLayer1Plate = new G4LogicalVolume(SolidLayer1Plate,               // its solid
			matCu,                     // ite material
			"Layer1Plate");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, 351 * mm),           // at(0,0,0)
		LogicLayer1Plate,                 // its logical volume
		"Layer1Plate",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number
	
}

void Layer2Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	//
	//Layer2 Plate
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4Tubs* SolidLayer2Plate =
		new G4Tubs("Layer2Plate", 0, 215 * mm, 3.9135 * mm, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicLayer2Plate = new G4LogicalVolume(SolidLayer2Plate,               // its solid
			matCu,                     // ite material
			"Layer2Plate");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, 152 * mm),           // at(0,0,0)
		LogicLayer2Plate,                 // its logical volume
		"Layer2Plate",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number

}

void Layer3Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	//
	//Layer3 Plate
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4Tubs* SolidLayer3Plate =
		new G4Tubs("Layer3Plate", 0, 185 * mm, 2.4965 * mm, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicLayer3Plate = new G4LogicalVolume(SolidLayer3Plate,               // its solid
			matCu,                     // ite material
			"Layer3Plate");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, 0 * mm),           // at(0,0,0)
		LogicLayer3Plate,                 // its logical volume
		"Layer3Plate",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number

}

void Layer4Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	//
	//Layer4 Plate
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4Tubs* SolidLayer4Plate =
		new G4Tubs("Layer4Plate", 0, 165 * mm, 2.5085 * mm, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicLayer4Plate = new G4LogicalVolume(SolidLayer4Plate,               // its solid
			matCu,                     // ite material
			"Layer4Plate");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, -155 * mm),           // at(0,0,0)
		LogicLayer4Plate,                 // its logical volume
		"Layer4Plate",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number

}

void Layer5Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	//
	//Layer5 Plate
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Pb");

	G4Tubs* SolidLayer5Plate =
		new G4Tubs("Layer5Plate", 0, 145 * mm, 2.474 * mm, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicLayer5Plate = new G4LogicalVolume(SolidLayer5Plate,               // its solid
			matCu,                     // ite material
			"Layer5Plate");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, -354 * mm),           // at(0,0,0)
		LogicLayer5Plate,                 // its logical volume
		"Layer5Plate",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number

}

void Bucket1Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Al
	G4Material* matAl = manager->FindOrBuildMaterial("G4_Al");

	G4double halfinnerRadius = 270.077 * mm;
	G4double halfouterRadius = 275. * mm;
	G4double halfheight = 710. * mm;
	G4double halfBottomThickness =  4.  * mm;

	G4Tubs* solidCylinder = 
		new G4Tubs("Cylinder", halfinnerRadius, halfouterRadius, halfheight, 0, CLHEP::twopi);
	G4Tubs* solidDisk = 
		new G4Tubs("Bottom", 0, halfouterRadius, halfBottomThickness, 0, CLHEP::twopi);
	G4ThreeVector diskPosition(0., 0., -(halfheight + halfBottomThickness) * mm);
	G4UnionSolid* solidBucket1 = new G4UnionSolid("Bucket1", solidCylinder, solidDisk, nullptr, diskPosition);

	G4LogicalVolume*
		logicBucket1 = new G4LogicalVolume(solidBucket1,                                     // its solid
			matAl,                                             // its material
			"Bucket1");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -164 * mm - (halfheight - 630 * mm)),                                     // at (0,0,0)
		logicBucket1,                                        // its logical volume
		"Bucket1",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}

void Bucket2Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Al
	G4Material* matAl = manager->FindOrBuildMaterial("G4_Al");

	G4double halfinnerRadius = 235.976 * mm;
	G4double halfouterRadius = 238. * mm;
	G4double halfheight =  640. * mm;
	G4double halfBottomThickness = 2. * mm;

	G4Tubs* solidCylinder =
		new G4Tubs("Cylinder", halfinnerRadius, halfouterRadius, halfheight, 0, CLHEP::twopi);
	G4Tubs* solidDisk =
		new G4Tubs("Bottom", 0, halfouterRadius, halfBottomThickness, 0, CLHEP::twopi);
	G4ThreeVector diskPosition(0., 0., -(halfheight + halfBottomThickness) * mm);
	G4UnionSolid* solidBucket2 = new G4UnionSolid("Bucket2", solidCylinder, solidDisk, nullptr, diskPosition);

	G4LogicalVolume*
		logicBucket2 = new G4LogicalVolume(solidBucket2,                                     // its solid
			matAl,                                             // its material
			"Bucket2");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -217 * mm - (halfheight - 560 * mm)),                                     // at (0,0,0)
		logicBucket2,                                        // its logical volume
		"Bucket2",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}

void Bucket3Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Al
	G4Material* matAl = manager->FindOrBuildMaterial("G4_Al");

	G4double halfinnerRadius = 212.976 * mm;
	G4double halfouterRadius = 215. * mm;
	G4double halfheight = 530.5 * mm;
	G4double halfBottomThickness = 2. * mm;

	G4Tubs* solidCylinder =
		new G4Tubs("Cylinder", halfinnerRadius, halfouterRadius, halfheight, 0, CLHEP::twopi);
	G4Tubs* solidDisk =
		new G4Tubs("Bottom", 0, halfouterRadius, halfBottomThickness, 0, CLHEP::twopi);
	G4ThreeVector diskPosition(0., 0., -(halfheight + halfBottomThickness) * mm);
	G4UnionSolid* solidBucket3 = new G4UnionSolid("Bucket3", solidCylinder, solidDisk, nullptr, diskPosition);

	G4LogicalVolume*
		logicBucket3 = new G4LogicalVolume(solidBucket3,                                     // its solid
			matAl,                                             // its material
			"Bucket3");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -306.5 * mm - (halfheight - 450.5 * mm)),                                     // at (0,0,0)
		logicBucket3,                                        // its logical volume
		"Bucket3",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}

void Shield4Construction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Cu
	G4Material* matAl = manager->FindOrBuildMaterial("G4_Cu");

	G4double halfinnerRadius = 184.127 * mm;
	G4double halfouterRadius = 185 * mm;
	G4double halfheight = 367 * mm;
	G4double halfBottomThickness = 0.4365 * mm;

	G4Tubs* solidCylinder =
		new G4Tubs("Cylinder", halfinnerRadius, halfouterRadius, halfheight, 0, CLHEP::twopi);
	G4Tubs* solidDisk =
		new G4Tubs("Bottom", 0, halfouterRadius, halfBottomThickness, 0, CLHEP::twopi);
	G4ThreeVector diskPosition(0., 0., -(halfheight + halfBottomThickness) * mm);
	G4UnionSolid* solidShield4 = new G4UnionSolid("Shield4", solidCylinder, solidDisk, nullptr, diskPosition);

	G4LogicalVolume*
		logicShield4 = new G4LogicalVolume(solidShield4,                                     // its solid
			matAl,                                             // its material
			"Shield4");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -295 * mm - (halfheight - 287. * mm)),                                     // at (0,0,0)
		logicShield4,                                        // its logical volume
		"Shield4",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}

void CopperHolderConstruction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	//
	//Cu holder
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4double halfinnerRadius = 18.0 * mm;
	G4double halfouterRadius = 22.0 * mm;
	G4double halfheight = 28.0 * mm;
	G4double halfBottomThickness = 2 * mm;

	G4Tubs* solidCylinder =
		new G4Tubs("Cylinder", halfinnerRadius, halfouterRadius, halfheight = 28.0 * mm, 0, CLHEP::twopi);
	G4Tubs* solidDisk =
		new G4Tubs("Bottom", 10 * mm, halfouterRadius, halfBottomThickness, 0, CLHEP::twopi);
	G4ThreeVector diskPosition(0., 0., -(halfheight + halfBottomThickness) * mm);
	G4UnionSolid* SolidCopperHolder = new G4UnionSolid("CopperHolder", solidCylinder, solidDisk, nullptr, diskPosition);

	//G4Tubs* SolidCopperHolder =
	//	new G4Tubs("CopperHolder", 13.0 * mm, 17.0 * mm, 28.0 * mm, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicCopperHolder = new G4LogicalVolume(SolidCopperHolder,               // its solid
			matCu,                     // ite material
			"CopperHolder");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, -391.5 * mm),           // at(0,0,0)
		LogicCopperHolder,                 // its logical volume
		"CopperHolder",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number

}

void LeadShieldConstruction(G4LogicalVolume* logicWorld, G4double height) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Pb
	G4Material* matPb = manager->FindOrBuildMaterial("G4_Pb");

	G4double innerRadius = 400 * mm;
	G4double outerRadius = 500 * mm;
	G4double halfBottomThickness = 50 * mm;

	G4Tubs* solidCylinder =
		new G4Tubs("Cylinder", innerRadius, outerRadius, height / 2., 0, CLHEP::twopi);
	G4Tubs* solidDisk =
		new G4Tubs("Bottom", 0, outerRadius, halfBottomThickness, 0, CLHEP::twopi);
	G4ThreeVector diskPosition(0., 0., -(height / 2. + halfBottomThickness) * mm);
	G4UnionSolid* solidLeadShield = new G4UnionSolid("LeadShield", solidCylinder, solidDisk, nullptr, diskPosition);

	//G4Box* solidOut =
	//	new G4Box("solidOut", halfouterRadius, halfouterRadius, halfheight);
	//G4Box* solidIn =
	//	new G4Box("solidIn", halfinnerRadius, halfinnerRadius, halfheight);
	//G4ThreeVector diskPosition(0., 0., +100 * mm);
	//G4SubtractionSolid* solidLeadShield = new G4SubtractionSolid("LeadShield", solidOut, solidIn, nullptr, diskPosition);

	G4LogicalVolume*
		logicLeadShield = new G4LogicalVolume(solidLeadShield,                                     // its solid
			matPb,                                             // its material
			"LeadShield");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -1020. * mm + height / 2.),                                     // at (0,0,0)
		logicLeadShield,                                        // its logical volume
		"LeadShield",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}

void CuShieldConstruction(G4LogicalVolume* logicWorld, G4double height, G4double height2) {

	//height -= 20 * mm;
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4double innerRadius = 325 * mm;
	G4double outerRadius = 475 * mm;
	G4double halfBottomThickness = 75 * mm;

	G4Tubs* solidCylinder =
		new G4Tubs("Cylinder", innerRadius, outerRadius, height / 2., 0, CLHEP::twopi);
	G4Tubs* solidDisk =
		new G4Tubs("Bottom", 0, outerRadius, halfBottomThickness, 0, CLHEP::twopi);
	G4ThreeVector diskPosition(0., 0., -(height / 2. + halfBottomThickness) * mm);
	G4UnionSolid* solidCuShield = new G4UnionSolid("CuShield", solidCylinder, solidDisk, nullptr, diskPosition);

	G4LogicalVolume*
		logicCuShield = new G4LogicalVolume(solidCuShield,                                     // its solid
			matCu,                                             // its material
			"CuShield");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -1000. * mm + height / 2.),                                     // at (0,0,0)
		logicCuShield,                                        // its logical volume
		"CuShield",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number

	//thinner

	//G4double height2 = 20 * cm;

	G4Tubs* solidCuShield_up =
		new G4Tubs("CuShield_up", innerRadius, innerRadius + 10 * cm, height2 / 2., 0, CLHEP::twopi);

	G4LogicalVolume*
		logicCuShield_up = new G4LogicalVolume(solidCuShield_up,                                     // its solid
			matCu,                                             // its material
			"CuShield_up");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -1000. * mm + height + height2 / 2.),                                     // at (0,0,0)
		logicCuShield_up,                                        // its logical volume
		"CuShield_up",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}

/*
void CuShieldConstruction(G4LogicalVolume* logicWorld, G4double height) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4double halfinnerRadius = 380 * mm;
	G4double halfouterRadius = 400 * mm;
	G4double halfheight = 640 * mm;
	G4double BottomThickness = 30 * mm;

	//G4Tubs* solidCylinder =
	//	new G4Tubs("Cylinder", halfinnerRadius, halfouterRadius, halfheight, 0, CLHEP::twopi);
	//G4Tubs* solidDisk =
	//	new G4Tubs("Bottom", 0, halfouterRadius, halfBottomThickness, 0, CLHEP::twopi);
	//G4ThreeVector diskPosition(0., 0., -(halfheight + halfBottomThickness) * mm);
	//G4UnionSolid* solidLeadShield = new G4UnionSolid("LeadShield", solidCylinder, solidDisk, nullptr, diskPosition);

	G4Box* solidOut =
		new G4Box("solidOut", halfouterRadius, halfouterRadius, halfheight);
	G4Box* solidIn =
		new G4Box("solidIn", halfinnerRadius, halfinnerRadius, halfheight);
	G4ThreeVector diskPosition(0., 0., BottomThickness * mm);
	G4SubtractionSolid* solidCuShield = new G4SubtractionSolid("CuShield", solidOut, solidIn, nullptr, diskPosition);

	G4LogicalVolume*
		logicCuShield = new G4LogicalVolume(solidCuShield,                                     // its solid
			matCu,                                             // its material
			"CuShield");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -759 * mm + halfheight),                                     // at (0,0,0)
		logicCuShield,                                        // its logical volume
		"CuShield",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}
*/

void WaterShieldConstruction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Water
	G4Material* matWater = manager->FindOrBuildMaterial("G4_WATER");

	G4double halfinnerRadius = 400 * mm;
	G4double halfouterRadius = 900 * mm;
	G4double halfheight = 700 * mm;
	G4double halfBottomThickness = 250 * mm;

	G4Tubs* solidCylinder =
		new G4Tubs("Cylinder", halfinnerRadius, halfouterRadius, halfheight, 0, CLHEP::twopi);
	G4Tubs* solidDisk =
		new G4Tubs("Bottom", 0, halfouterRadius, halfBottomThickness, 0, CLHEP::twopi);
	G4ThreeVector diskPosition(0., 0., -(halfheight + halfBottomThickness) * mm);
	G4UnionSolid* solidWaterShield = new G4UnionSolid("WaterShield", solidCylinder, solidDisk, nullptr, diskPosition);

	G4LogicalVolume*
		logicWaterShield = new G4LogicalVolume(solidWaterShield,                                     // its solid
			matWater,                                             // its material
			"WaterShield");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, -689 * mm + halfheight),                                     // at (0,0,0)
		logicWaterShield,                                        // its logical volume
		"WaterShield",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}



void CuTopShieldConstruction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4double halfx = 1150 * mm;
	G4double halfy = 575 * mm;
	G4double halfz = 20 * mm;

	G4Box* solidCuTopShield = new G4Box("CuTopShield", halfx, halfy, halfz);

	G4LogicalVolume*
		logicCuTopShield = new G4LogicalVolume(solidCuTopShield,                                     // its solid
			matCu,                                             // its material
			"CuTopShield");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, 1246 * mm),                                     // at (0,0,0)
		logicCuTopShield,                                        // its logical volume
		"CuTopShield",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}

void SubstanceForCosmicRayConstruction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();
	// Material ---> Pb
	G4Material* matPb = manager->FindOrBuildMaterial("G4_Pb");

	G4double halfx = 2500 * mm;
	G4double halfy = 2500 * mm;
	G4double halfz = 200 * mm;

	G4Box* solidSubstanceForCosmicRayTopShield = new G4Box("SubstanceForCosmicRayTopShield", halfx, halfy, halfz);

	G4LogicalVolume*
		logicSubstanceForCosmicRayTopShield = new G4LogicalVolume(solidSubstanceForCosmicRayTopShield,                                     // its solid
			matPb,                                             // its material
			"SubstanceForCosmicRayTopShield");                                          // its name 

	//  G4VPhysicalVolume*  physiCrystal =
	new G4PVPlacement(0,                                                  // no rotation
		G4ThreeVector(0, 0, 1800 * mm),                                     // at (0,0,0)
		logicSubstanceForCosmicRayTopShield,                                        // its logical volume
		"SubstanceForCosmicRayTopShield",                                           // its name
		logicWorld,                                           // its mother  volume
		false,                                               // no boolean operation
		0);                                                  // copy number
}

void CuInnerShieldConstruction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	G4double height = 75 * mm;

	//
	//CuInnerShield
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4Tubs* SolidCuInnerShield =
		new G4Tubs("CuInnerShield", 0, 125.0 * mm, height, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicCuInnerShield = new G4LogicalVolume(SolidCuInnerShield,               // its solid
			matCu,                     // ite material
			"CuInnerShield");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, -450 * mm - (height - 75.0 * mm)),           // at(0,0,0)
		LogicCuInnerShield,                 // its logical volume
		"CuInnerShield",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number

}

void SpringPlateConstruction(G4LogicalVolume* logicWorld) {
	G4NistManager* manager = G4NistManager::Instance();

	G4double height = 5.0 * mm;

	//
	//SpringPlate
	// 

	// Material ---> Cu
	G4Material* matCu = manager->FindOrBuildMaterial("G4_Cu");

	G4Tubs* SolidSpringPlate =
		new G4Tubs("SpringPlate", 0, 50. * mm, height, 0, CLHEP::twopi);

	G4LogicalVolume*
		LogicSpringPlate = new G4LogicalVolume(SolidSpringPlate,               // its solid
			matCu,                     // ite material
			"SpringPlate");                  // its name

	//  G4VPhysicalVolume*  physiCell = 
	new G4PVPlacement(0,                         // no rotation
		G4ThreeVector(0 * mm, 0 * mm, -420 * mm),           // at(0,0,0)
		LogicSpringPlate,                 // its logical volume
		"SpringPlate",                    // its name
		logicWorld,          // its mother volume
		false,                     // no boolean operation
		0);                        // copy number

}