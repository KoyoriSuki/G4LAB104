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
// $Id: HistoManager.cc 72249 2013-07-12 08:57:49Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include <stdio.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager(G4String processNo)
  : fFileName("shieldingHeightResult"),
    processNumber(processNo)
{
    messenger = new HistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  //delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{
    //
    // Initialize the histograms
    //
    std::string folderPath = "./" + fFileName + "_results"; // 定义文件夹路径
    mkdir(folderPath.c_str(), 0777);

  G4cout<<"------>create rootfile"<<G4endl;
  std::string fRootFileName = "./" + fFileName + "_results" + "/" + fFileName + processNumber + ".root";
  fRootFile = new TFile(fRootFileName.c_str(), "RECREATE");
  fNtuple = new TTree("TreeOfData","data");
  fNtuple->Branch("Initial_Energy", &fEventInfo.fInitialEnergy);
  fNtuple->Branch("Energy_Deposition_In_Crystal", &fEventInfo.fEnergyDepositionInCrystal);
  fNtuple->Branch("Energy_Deposition_In_Shielding", &fEventInfo.fEnergyDepositionInShielding);
  fNtuple->Branch("Direction_Zenith_Angle", &fEventInfo.fStartPositionTheta);
}



void HistoManager::save()
{
    //for (double elem : fEventInfo.fEnergyDepositionInShielding) {
    //    std::cout << elem << " ";
    //}
    fNtuple->Fill();
    fRootFile->Write();
    fRootFile->Close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
