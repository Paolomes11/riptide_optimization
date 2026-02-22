#ifndef RIPTIDE_RUN_ACTION_HPP
#define RIPTIDE_RUN_ACTION_HPP

#include <G4UserRunAction.hh>

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    ~RunAction();

    virtual void BeginOfRunAction(const G4Run *);
    virtual void EndOfRunAction(const G4Run *);
};

#endif // RIPTIDE_RUN_ACTION_HPP