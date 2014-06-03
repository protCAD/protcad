//filename MyElement.h contains a header file for myelement
#include "assert.h"
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include "rules.h"


#ifndef MYELEMENT_H
#define MYELEMENT_H

//is the element contained by the GA.  just contains an integer.
//ummmmm...at one time it actually did something.  now i guess it doesn't.
class myElement
{
 public:
  myElement();
  myElement(int _val);
  ~myElement();
  friend ostream& operator<<(ostream &, const myElement &);
  int operator!=(const myElement &) const;
  int operator==(const myElement &) const;
  myElement operator=(const myElement &);
  int getValue() const { return val; }
  void print() const { cout << val << "\n"; }
  void setValue(int _val);

 protected:
  int val;
};

#endif
