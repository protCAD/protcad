#include "assert.h"
#include <iostream.h>
#include <ga/ga.h>
#include "superResidue.h"
#include <chainPosition.h>

#ifndef RULES_H
#define RULES_H

//this class constructs an array which holds all of the possible
//residues which can exist at a chain position.  if parses through the 
//information given by the position, and sets up the array.
//this allows the GA to just hold ints, and can pick random ints within
//a range given by rules.  this random number can then be interpreted by rules
//as needed so that it can be converted into a residue or a float value.
//
//this helps lessen the load on the ga, since it just creates a bunch of random
//integers.
class rules
{
 public:
  //constructors
  rules();
  rules(chainPosition *_pos, int _chain, int _residue);
  ~rules();
  
  //accessors
  int getSize() const { return size; }
  int getNewValue();
  int getFlag() const { return flag; }
  chainPosition* getPosition() const { return pos; }
  int getChain() const { return chain; }
  double getMax();
  double getStep();
  int getResidue() const { return residue; }
  vector<superResidue> getData() const { return data; }

  //other functions
  rules operator=(const rules &);
  void print();
  void print(int n);
  superResidue translateRes(int _orig);
  
 protected:
  int size;
  int flag;
  double max;
  double step;
  chainPosition * pos;
  vector<superResidue> data;
  int chain;
  int residue;
};
#endif
