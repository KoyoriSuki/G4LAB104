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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "TrackingAction.hh"

#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
#include <math.h>
#include <iostream>
#include "VoltageFieldData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
SteppingAction* SteppingAction::fgInstance =0;
SteppingAction* SteppingAction::Instance()
{
  return fgInstance;
}
//
SteppingAction::SteppingAction(DetectorConstruction* det,
			       TrackingAction*       tracking
    ,EventAction* CurrentEventAction
) 
  : G4UserSteppingAction(),
    fDetector(det), 
    fTrackingAction(tracking),
    fEventAction(CurrentEventAction)
{
  fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ 
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // get volume of the current step
  G4String preVolumeName 
    = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if(preVolumeName != "Crystal" && preVolumeName != "PbShell") 
    return;
  
  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4ThreeVector positionPost = aStep->GetPostStepPoint()->GetPosition();
  if(edep > 1e-10 && preVolumeName == "Crystal"){
      fTrackingAction->AddEdep_Crystal(edep);
      fEventAction->EffectiveEventFlag = 1;
  }
  if(edep > 1e-10 && preVolumeName == "PbShell"){
    //fTrackingAction->AddEdep_PbShell(EdepTailCorrection(edep, positionPost));
  } 

  // if any particle passes through the detector, set the flag to 1
  //if (preVolumeName == "Crystal")
  //{
  //    
  //}
}
 
void SteppingAction::Reset()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
