#include "superResidue.h"

superResidue::superResidue() 
{
  id = 0;
  rot = 0;
}

superResidue::superResidue(int _identity, int _rotamer)
{
  id = _identity;
  rot = _rotamer;
}

superResidue::~superResidue()
{
}

int superResidue::operator!=(const superResidue & res) const
{
  if (id != res.getIdentity() || 
      rot != res.getRotamer()) {
    return 1;
  }
  return 0;
}

int superResidue::operator==(const superResidue & rhs) const
{
  if (getIdentity() == rhs.getIdentity() && 
      getRotamer() == rhs.getRotamer()) {
    return 1;
  }
  return 0;
}
  
superResidue superResidue::operator=(const superResidue & res)
{
  id = res.getIdentity();
  rot = res.getRotamer();
  return res;
}

ostream& operator<<(ostream& theStream, const superResidue& res)
{
  theStream << res.getIdentity()<< " " << res.getRotamer();
  return theStream;
}
      
void superResidue::printRes() const
{
  cout << "id: " << id << "  rot: " << rot << "\n";
}
