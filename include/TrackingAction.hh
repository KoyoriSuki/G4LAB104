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
/// \file DBDecay/include/TrackingAction.hh
/// \brief Definition of the TrackingAction class
//
//
// $Id: TrackingAction.hh 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include "HistoManager.hh"
#include "SteppingAction.hh"
#include <iomanip>
//class SteppingVerbose;
class RunAction;
class EventAction;
class TrackingMessenger;
class TrackInfo;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackingAction : public G4UserTrackingAction {

  public:  
  TrackingAction(RunAction*,EventAction*, HistoManager*);
   ~TrackingAction();
   
    virtual void  PreUserTrackingAction(const G4Track*);
    virtual void  PostUserTrackingAction(const G4Track*);
    
    void SetFullChain(G4bool flag) { fFullChain = flag;};
  
    void AddEdep_Crystal(G4double edep){
      fTrackEdepInCrystal += edep;
      return;
    }
    void AddEdep_PbShell(G4double edep){
      fTrackEdepInPbShell += edep;
      return;
    }
  
  private:
  RunAction* fRun;
  EventAction*        fEvent;
  TrackingMessenger*  fTrackMessenger;
  HistoManager* fHistoManager;
    
  G4double fCharge, fBaryonNo, fMass;
  G4double  fParticleEnCode;        
  G4double fTrackEdepInCrystal;
  G4double fTrackEdepInPbShell;
  G4bool   fFullChain;
  //SteppingAction fSteppingAction_Tracking;
 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
