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

#include "G4UImessenger.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//...oooOOO0OOOooo......................................................
DetectorMessenger::DetectorMessenger(DetectorConstruction* detConstruct)
	:myDetectorConstruction(detConstruct)
{
	usrDefinedGeoDirectory = new G4UIdirectory("/usrDefinedGeo/");
	usrDefinedGeoDirectory->SetGuidance("Geometry control commands.");

	crystalTypeCmd = new G4UIcmdWithAnInteger("/usrDefinedGeo/setCrystalType", this);
	crystalTypeCmd->SetGuidance("Set Crystal Type.");
	crystalTypeCmd->SetParameterName("crystalType", true);
	crystalTypeCmd->SetDefaultValue(0);

	shieldingHeightCmd = new G4UIcmdWithADouble("/usrDefinedGeo/setShieldingHeight", this);
	shieldingHeightCmd->SetGuidance("Set Shielding Height.");
	shieldingHeightCmd->SetParameterName("shieldingHeight", true);
	shieldingHeightCmd->SetDefaultValue(50.0);

	shieldingHeight2Cmd = new G4UIcmdWithADouble("/usrDefinedGeo/setShieldingHeight2", this);
	shieldingHeight2Cmd->SetGuidance("Set Shielding Height2.");
	shieldingHeight2Cmd->SetParameterName("shieldingHeight2", true);
	shieldingHeight2Cmd->SetDefaultValue(50.0);

	shieldingFlagCmd = new G4UIcmdWithAnInteger("/usrDefinedGeo/setShieldingFlag", this);
	shieldingFlagCmd->SetGuidance("Set shielding flag.");
	shieldingFlagCmd->SetParameterName("shieldingFlag", true);
	shieldingFlagCmd->SetDefaultValue(0);
}

//...oooOOO0OOOooo.......................................................
DetectorMessenger::~DetectorMessenger()
{
	delete crystalTypeCmd;
	delete shieldingHeightCmd;
}

//...oooOOO0OOOooo.......................................................
void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if (command == crystalTypeCmd) {
		myDetectorConstruction->SetCrystalType(crystalTypeCmd->GetNewIntValue(newValue));
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
		std::cout << "DetectorMessenger Called" << std::endl;
	}
	else if (command == shieldingHeightCmd) {
		myDetectorConstruction->SetShieldingHeight(shieldingHeightCmd->GetNewDoubleValue(newValue));
		std::cout << "Setting new shielding height: " << newValue << std::endl;
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
	}
	else if (command == shieldingHeight2Cmd) {
		myDetectorConstruction->SetShieldingHeight2(shieldingHeight2Cmd->GetNewDoubleValue(newValue));
		std::cout << "Setting new shielding height2: " << newValue << std::endl;
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
	}
	else if (command == shieldingFlagCmd) {
		myDetectorConstruction->SetShieldingFlag(shieldingFlagCmd->GetNewIntValue(newValue));
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
		std::cout << "DetectorMessenger Called" << std::endl;
	}
}