#include "HistoMessenger.hh"
#include "HistoManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

HistoMessenger::HistoMessenger(HistoManager* histomanager)
    : G4UImessenger(),
    fHistoManager(histomanager),
    fDirectory(0),
    fSetFileNameCmd(0)
{
    fDirectory = new G4UIdirectory("/histo/");
    fDirectory->SetGuidance("Histogram manager control");

    fSetFileNameCmd = new G4UIcmdWithAString("/histo/setFileName", this);
    fSetFileNameCmd->SetGuidance("Set the name of the output file.");
    fSetFileNameCmd->SetParameterName("fileName", false);
    fSetFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

HistoMessenger::~HistoMessenger()
{
    delete fSetFileNameCmd;
    delete fDirectory;
}

void HistoMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    if (command == fSetFileNameCmd) {
        fHistoManager->SetFileName(newValue);
    }
}
