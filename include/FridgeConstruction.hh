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
#include "G4UnionSolid.hh"

//plate 1->5 up->down 0:top one
void Layer0OuterConstruction(G4LogicalVolume* logicWorld);
void Layer0InnerConstruction(G4LogicalVolume* logicWorld);
void Layer1Construction(G4LogicalVolume* logicWorld);
void Layer2Construction(G4LogicalVolume* logicWorld);
void Layer3Construction(G4LogicalVolume* logicWorld);
void Layer4Construction(G4LogicalVolume* logicWorld);
void Layer5Construction(G4LogicalVolume* logicWorld);

// bucket/shield 1->4 outter->inner
void Bucket1Construction(G4LogicalVolume* logicWorld);
void Bucket2Construction(G4LogicalVolume* logicWorld);
void Bucket3Construction(G4LogicalVolume* logicWorld);
void Shield4Construction(G4LogicalVolume* logicWorld);

//other component
void CopperHolderConstruction(G4LogicalVolume* logicWorld);
void LeadShieldConstruction(G4LogicalVolume* logicWorld, G4double height);
void WaterShieldConstruction(G4LogicalVolume* logicWorld);
void CuShieldConstruction(G4LogicalVolume* logicWorld, G4double height, G4double height2);
void CuTopShieldConstruction(G4LogicalVolume* logicWorld);
void SubstanceForCosmicRayConstruction(G4LogicalVolume* logicWorld);
void CuInnerShieldConstruction(G4LogicalVolume* logicWorld);
void SpringPlateConstruction(G4LogicalVolume* logicWorld);