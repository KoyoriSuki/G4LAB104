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
/// \file DBDecay/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "PrimaryGeneratorAction.hh"
#include <iostream>
#include <iomanip>
#include "PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
//#include "G4HEPEvtInterface.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4Gamma.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4UnionSolid.hh"
#include "G4NeutrinoE.hh"
#include "G4UImanager.hh"

//#include "CRYSetup.h"
//#include "CRYGenerator.h"
//#include "CRYParticle.h"
//#include "CRYUtils.h"
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector RndPointInTub(double radiusInner, double radiusOuter, double height, G4ThreeVector center) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disR2(radiusInner * radiusInner, radiusOuter * radiusOuter);
    std::uniform_real_distribution<> disTheta(0, 2 * M_PI);
    std::uniform_real_distribution<> disZ(center.z() - height / 2, center.z() + height / 2);
    double r = sqrt(disR2(gen));
    double theta = disTheta(gen);
    double z = disZ(gen);
    G4ThreeVector p;
    p.setX(center.x() + r * cos(theta));
    p.setY(center.y() + r * sin(theta));
    p.setZ(z);
    return p;
}
double rndSeedCRY()
{
    static std::random_device rd;
    static std::mt19937 generator(rd());
    static std::uniform_real_distribution<double> seedRnd(0.0, 1.0);
    return seedRnd(generator);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
 : G4VUserPrimaryGeneratorAction(),
   fGParticleSource(),
   fDetector(det)
{
  fGParticleSource  = new G4GeneralParticleSource();

  // Where the UI commands works
  messenger = new PrimaryGeneratorMessenger(this);

  // If using code for GPS settings
  if (sourceType == 1)
  {
      //
      //Energy Distribution Initialization
      //
      std::ifstream histFile("unfoldedSpectra_centerValue.txt");
      std::vector<double> histData;
      double freq;
      double totalFreq = 0.0;
      while (histFile >> freq) {
          histData.push_back(freq);
          totalFreq += freq;
      }
      std::transform(histData.begin(), histData.end(), histData.begin(),
          [totalFreq](double freq1) { return freq1 / totalFreq; });
      EnergyDistribution = std::discrete_distribution<int>(histData.begin(), histData.end());
      generator.seed(std::random_device()());
  }

  // If using CRY
  if (sourceType == 2)
  {
      //// Read the cry input file
      //std::ifstream inputFile;
      //inputFile.open(inputFileName, std::ios::in);
      //char buffer[1000];
      //if (inputFile.fail()) {
      //    if (inputFileName != "")  //....only complain if a filename was given
      //        G4cout << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputFileName << G4endl;
      //    exit(-1);
      //}
      //else
      //{
      //    std::string setupString("");
      //    while (!inputFile.getline(buffer, 1000).eof())
      //    {
      //        setupString.append(buffer);
      //        setupString.append(" ");
      //    }
      //    CRYSetup* setup = new CRYSetup(setupString, "../data");
      //    gen = new CRYGenerator(setup);
      //    // set random number generator
      //    setup->setRandomFunction(rndSeedCRY);
      //}
      //// create a vector to store the CRY particle properties
      //vect = new std::vector<CRYParticle*>;
      //// Create the table containing all particle names
      //particleTable = G4ParticleTable::GetParticleTable();
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  //delete HEPEvt;
  delete fGParticleSource;
  delete messenger;

  //fclose(countsRecord);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    if (sourceType == 1) // Random position in shielding Copper
    {
        G4ThreeVector rndPointInCuShielding(0, 0, 0);
        double radiusInner = 0.375 * m, radiusOuter = 0.525 * m, wallHeight = 1.400 * m, diskHeight = 0.150 * m;
        G4ThreeVector centerPointWall(0, 0, -1.000 * m + wallHeight / 2.);
        G4ThreeVector centerPointDisk(0, 0, -1.000 * m - diskHeight / 2.);
        double wallVolume = 3.14159 * (radiusOuter * radiusOuter - radiusInner * radiusInner) * wallHeight;
        double diskVolume = 3.14159 * (radiusOuter * radiusOuter) * diskHeight;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> volumeChoice(0, wallVolume + diskVolume);
        double a = volumeChoice(gen);
        if (a < wallVolume)
        {
            rndPointInCuShielding = RndPointInTub(radiusInner,
                radiusOuter,
                wallHeight,
                centerPointWall);
        }
        else
        {
            rndPointInCuShielding = RndPointInTub(0,
                radiusOuter,
                diskHeight,
                centerPointDisk);
        }
        std::cout << rndPointInCuShielding.x() << " " << rndPointInCuShielding.y() << " " << rndPointInCuShielding.z() << std::endl;
        fGParticleSource->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
        fGParticleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(rndPointInCuShielding);
        //generate
        fGParticleSource->GeneratePrimaryVertex(anEvent);
    }
    else if (sourceType == 2)
    {
        //G4String particleName;
        //vect->clear();
        //gen->genEvent(vect);

        //for (unsigned j = 0; j < vect->size(); j++) {
        //    particleName = CRYUtils::partName((*vect)[j]->id());

        //    //....debug output  
        //    //cout << "  " << particleName << " "
        //    //    << "charge=" << (*vect)[j]->charge() << " "
        //    //    << setprecision(4)
        //    //    << "energy (MeV)=" << (*vect)[j]->ke() * MeV << " "
        //    //    << "pos (m)"
        //    //    << G4ThreeVector((*vect)[j]->x(), (*vect)[j]->y(), (*vect)[j]->z())
        //    //    << " " << "direction cosines "
        //    //    << G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w())
        //    //    << " " << endl;

        //    fGParticleSource->SetParticleDefinition(particleTable->FindParticle((*vect)[j]->PDGid()));
        //    fGParticleSource->GetCurrentSource()->GetEneDist()->SetMonoEnergy((*vect)[j]->ke() * MeV);
        //    fGParticleSource->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
        //    fGParticleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector((*vect)[j]->x() * m, (*vect)[j]->y() * m, (*vect)[j]->z() * m));
        //    fGParticleSource->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector((*vect)[j]->u(), (*vect)[j]->v(), (*vect)[j]->w()));
        //    fGParticleSource->GetCurrentSource()->SetParticleTime((*vect)[j]->t());
        //    //generate
        //    fGParticleSource->GeneratePrimaryVertex(anEvent);
        //    delete (*vect)[j];
        //}
    }
    else
    {
        //generate
        fGParticleSource->GeneratePrimaryVertex(anEvent);
    }


   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
