#include "assert.h"
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <molecule.h>
#include <ga/ga.h>

#ifndef SRESIDUE_H
#define SRESIDUE_H

//contains the identity and rotamer of an amino acid.  
//also contains some other friend functions so that the GA can use it.
class superResidue
{
 public:
  superResidue();
  superResidue(int _id, int _rotamer);
  ~superResidue();
  friend ostream& operator<<(ostream &, const superResidue &);
  int operator!=(const superResidue &) const;
  int operator==(const superResidue &) const;
  superResidue  operator=(const superResidue &);

  int getIdentity() const { return id; }
  int getRotamer() const { return rot; }
  void printRes() const;

 protected:
  int id;
  int rot;
};

#endif
