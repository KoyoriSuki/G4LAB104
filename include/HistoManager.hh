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
// $Id: HistoManager.hh 68017 2013-03-13 13:29:53Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef HistoManager_h
#define HistoManager_h 1

#include "HistoMessenger.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <G4ThreeVector.hh>
#include <TTree.h>
#include <TFile.h>
#include <sys/stat.h>

//const G4int kMAXTrack=5000;//should be checked!!!
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class EventInfo
{
public:
  std::vector<double> fInitialEnergy;
  std::vector<double> fEnergyDepositionInCrystal;
  std::vector<double> fEnergyDepositionInShielding;
  std::vector<double> fStartPositionTheta;
  
  void reset()
  {
      fInitialEnergy.clear();
      fEnergyDepositionInCrystal.clear();
      fEnergyDepositionInShielding.clear();
      fStartPositionTheta.clear();
  };

  EventInfo()
  {
      fInitialEnergy.clear();
      fEnergyDepositionInCrystal.clear();
      fEnergyDepositionInShielding.clear();
      fStartPositionTheta.clear();
   }
};

class TrackInfo
{
public:
    std::vector<double> fInitialEnergy;
    std::vector<double> fEnergyDepositionInCrystal;
    std::vector<double> fEnergyDepositionInShielding;
    std::vector<double> fTime;

    void reset()
    {
        fInitialEnergy.clear();
        fEnergyDepositionInCrystal.clear();
        fEnergyDepositionInShielding.clear();
        fTime.clear();
    };

    TrackInfo()
    {
        fInitialEnergy.clear();
        fEnergyDepositionInCrystal.clear();
        fEnergyDepositionInShielding.clear();
        fTime.clear();
    }
};

class HistoManager
{
public:
   HistoManager(G4String);
  ~HistoManager();
  void save();
  void book();
  EventInfo fEventInfo;
  TrackInfo fTrackInfo;
  //std::fstream ResponseData;
  //FILE* ResponseData_c;
  //FILE* countsRecord;
  G4String processNumber;

private:
    HistoMessenger* messenger;
  void Book();
  
public:  
  TFile* fRootFile;
  TTree* fNtuple;
  inline void SetFileName(G4String f)
  {
      fFileName = f;
  }

public:
  G4String fFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

