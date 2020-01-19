#include "QBBC.hh"
#include "globals.hh"

class B1PhysicsList : public QBBC
{
public:

  explicit B1PhysicsList(G4int ver = 1, const G4String& type = "B1PhysicsList");

  virtual ~B1PhysicsList();

private:

  // copy constructor and hide assignment operator
  B1PhysicsList(B1PhysicsList &);
  B1PhysicsList & operator=(const B1PhysicsList &right);

};

