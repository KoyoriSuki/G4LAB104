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
/// \file DBDecay/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 68030 2015-03-13 13:51:27Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
//#include "EventMessenger.hh"

#include "G4Event.hh"
#include <iomanip>
#include <HistoManager.hh>
//#include <TTree.h>
#include "RunAction.hh"
#include "G4UnitsTable.hh"
#include "PrimaryGeneratorAction.hh"
#include <iostream>
#include <stdio.h>
#include <random>
#include <chrono>
#include <math.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(PrimaryGeneratorAction* kin,HistoManager* histo)
:G4UserEventAction(),
 fPrintModulo(1000),fDecayChain(), fPrimary(kin), fHistoManager_Event(histo)
{
  //fHistoManager_Event = new HistoManager();
  //fEventMessenger = new EventMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  //delete fHistoManager_Event;
  //delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
 fDecayChain = " ";
 fHistoManager_Event->fTrackInfo.reset();
 EffectiveEventFlag = 0;
// G4cout<<"begin of event"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID(); 
 //printing survey
 //
 if (evtNb%fPrintModulo == 0) 
   G4cout << "\n end of event " << std::setw(6) << evtNb 
          << " :" + fDecayChain << G4endl;

 //
 // calculate energy deposition of the current event
 //
 if (EffectiveEventFlag == 1)
 {
     int i = 0;
     double EventEdepInCrystal = 0;
     double EventEdepPbShell = 0;
     for (auto first = fHistoManager_Event->fTrackInfo.fEnergyDepositionInCrystal.begin(); first < fHistoManager_Event->fTrackInfo.fEnergyDepositionInCrystal.end(); ++first)
     {
         EventEdepInCrystal += double(fHistoManager_Event->fTrackInfo.fEnergyDepositionInCrystal[i]);
         EventEdepPbShell += double(fHistoManager_Event->fTrackInfo.fEnergyDepositionInShielding[i]);
         ++i;
     }
     fHistoManager_Event->fEventInfo.fEnergyDepositionInCrystal.push_back(EventEdepInCrystal);
     fHistoManager_Event->fEventInfo.fEnergyDepositionInShielding.push_back(EventEdepPbShell);
     //std::cout << EventEdepInCrystal << std::endl;

     //
     // calculate initial energy of the current event
     //
     G4double eprimary = fPrimary->GetParticleGun()->GetCurrentSource()->GetParticleEnergy();
     fHistoManager_Event->fEventInfo.fInitialEnergy.push_back(eprimary);

     //
     // calculate the zenith angle of the current event
     //
     double zenithAngle = 0;
     G4PrimaryVertex* primaryVertex = evt->GetPrimaryVertex(0);
     if (primaryVertex) {
         G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
         if (primaryParticle) {
             G4ThreeVector momentumDirection = primaryParticle->GetMomentumDirection();
             zenithAngle = std::acos(momentumDirection.z());
             //std::cout << zenithAngle << std::endl;
         }
     }
     fHistoManager_Event->fEventInfo.fStartPositionTheta.push_back(-cos(zenithAngle));
 }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


