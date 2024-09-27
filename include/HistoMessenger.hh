// HistoMessenger.hh
#ifndef HistoMessenger_h
#define HistoMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class HistoManager;

class HistoMessenger : public G4UImessenger
{
public:
    HistoMessenger(HistoManager* histomanager);
    virtual ~HistoMessenger();

    virtual void SetNewValue(G4UIcommand* command, G4String newValue);

private:
    HistoManager* fHistoManager;
    G4UIdirectory* fDirectory;
    G4UIcmdWithAString* fSetFileNameCmd;
};

#endif
