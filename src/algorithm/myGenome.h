#include "assert.h"
#include <stdio.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <ga/ga.h>
#include "rules.h"
#include "superResidue.h"
#include "myElement.h"
#include "protein.h"

#ifndef MY_GENOME_H
#define MY_GENOME_H

class myGenome
{
 public:

  //function declarations:
  
  myGenome(molecule* _pmol);                 //adds vars to the global environment
  void run();                                //instantiates the GA with the global variables
  void setSeed(const UInt _seed);            // set the random seed from the default
  static void ArrayInitializer(GAGenome &);  //initiates the GA as specified
  UInt getSeed() const {return seed;}
  UInt getSize() const { return size; }      //find out the size of the GA
  UInt getPopSize() { return popsize; }
  UInt getNumGen() { return ngen; }
  float getProbCross() { return pcross; }
  float getProbMut() { return pmut; }


  //the objective function is the most imprortant function.
  //It determines the fitness of individual genes.
  static float Objective(GAGenome &);

  static int MyMut(GAGenome &, float pmut);  //mutation function


  //set the population, number of generations, probability of mutation,
  //and probability of cross-over.
  void setPopSize(const UInt _pop);
  void setNumGen(const UInt _ngen);
  void setProbMut(const float _pmut);
  void setProbCross(const float _pcross);
  
  static molecule* pmol;

 protected:
  void makeIndex();

 protected:
  UInt popsize;
  UInt ngen;
  float pmut;
  float pcross;
  rules* pindex;
  UInt seed;
  UInt size;
};

#endif
