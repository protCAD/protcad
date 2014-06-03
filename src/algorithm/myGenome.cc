#include "myGenome.h"

molecule* myGenome::pmol = 0;
myGenome::myGenome(molecule* _pmol) {
  popsize = 50;
  ngen = 400;
  pmut = 0.01;
  pcross = 0.6;
  pmol = _pmol;
  seed = 12;
  
  makeIndex();
}

void myGenome::run() {
  GA1DArrayGenome<myElement> genome(size, myGenome::Objective, (void *)pindex);
  genome.initializer(ArrayInitializer);
  genome.mutator(MyMut);

  double energy;
  int val, id, rot, chn, resPos;
  superResidue res;
  rules * temp;

  GASimpleGA ga(genome);
  GASigmaTruncationScaling scaling; 
  ga.populationSize(popsize);
  ga.nGenerations(ngen);
  ga.pMutation(pmut);
  ga.pCrossover(pcross);
  ga.scaling(scaling);
  ga.scoreFilename("bog.dat");
  ga.scoreFrequency(10);
  ga.flushFrequency(50);
  ga.evolve(seed);

  genome = ga.statistics().bestIndividual();
  

  cout << "the results are: \n";
  for (int i = 0; i < genome.length(); i++)
    {
      temp = &pindex[i];
      val = genome.gene(i).getValue();
      temp->print(val);
      res = temp->translateRes(val);
      id = res.getIdentity();
      rot = res.getRotamer();
      chn = temp->getChain();
      resPos = temp->getResidue();
      static_cast<protein*>(pmol)->mutateWBC((UInt)chn, (UInt)resPos, (UInt)id);
      static_cast<protein*>(pmol)->setRotamer((UInt)chn, (UInt)resPos, 0, rot);
    }
  energy = ((protein*)(pmol))->intraEnergy();
  cout << "the energy is: " << energy << "\n"; 
}


//this function takes the integer values, interprets them into residues,
//applies those changes to the protein pointer, and then computes
//the energy of the new protein.  it then assigns a fitness to this new protein.
float myGenome::Objective(GAGenome& g)
{
  GA1DArrayGenome<myElement> & gen = (GA1DArrayGenome<myElement> &)g;
  float value = gen.length();
  superResidue res;
  int val, id, rot, chn, resPos;
  double energy;
  int sz = gen.length();
  rules * index = (rules *)g.userData();
  rules * temp = new rules;

  for (int i = 0; i < sz; i++)
    {
      temp = &index[i];
      val = gen.gene(i).getValue();
      res = temp->translateRes(val);
      id = res.getIdentity();
      rot = res.getRotamer();
      chn = temp->getChain();
      resPos = temp->getResidue();
      static_cast<protein*>(pmol)->mutateWBC((UInt)chn, (UInt)resPos, (UInt)id);
      static_cast<protein*>(pmol)->setRotamer((UInt)chn, (UInt)resPos, 0, rot);
    }
  energy = ((protein*)(pmol))->intraEnergy();
#ifdef MYGENOME_DEBUG
  cout << "the energy is: " << energy << "\n"; 
#endif
  value = -energy;
  return value;
}

//initializes each array to contain a random set of integers based on the 
//rules that are specified.
void myGenome::ArrayInitializer(GAGenome & c) {
  GA1DArrayGenome<myElement> & genome = (GA1DArrayGenome<myElement> &)c;
  rules *index = (rules *)c.userData();
  myElement mine;
  rules * temp = new rules;
  int n = genome.length();
  int val;

  for (int i = 0; i < n; i++) {
    temp = &index[i];
    val = temp->getNewValue();
    mine = myElement(val);
    genome.gene(i, mine);
  }
}

//picks a new random value for the integer, within the range spefified by rules.
int myGenome::MyMut(GAGenome & c, float pmut) {
  GA1DArrayGenome<myElement> & genome = (GA1DArrayGenome<myElement> &)c;
  rules * index = (rules *)c.userData();
  rules * rul = new rules;
  if (pmut <= 0.0) return (0);
  int nMut = 0;
  int val;
  myElement temp;
  for (int i = genome.length() - 1; i >= 0; i--) 
    {
      if (GAFlipCoin(pmut)) 
	{
	  rul = &index[i];
	  val = rul->getNewValue();
	  temp = myElement(val);
	  genome.gene(i, temp);
	}
    }
  return nMut;
}

void myGenome::setSeed(const UInt _seed) 
{
  seed = _seed;
}

void myGenome::setPopSize(const UInt _pop) 
{
  popsize = _pop;
}

void myGenome::setNumGen(const UInt _ngen)
{
  ngen = _ngen;
}

void myGenome::setProbMut(const float _pmut)
{
  pmut = _pmut;
}

void myGenome::setProbCross(const float _pcross)
{
  pcross = _pcross;
}

//should initialize the index matrix with the indeces to all of the
//residues.
//the index matrix contains pointers to each chain position that is being
//modified, so that the objective function can apply any changes.
void myGenome::makeIndex()
{
  int count = 0;
  vector<chainPosition*> actives;
  UInt n;

  UInt number = static_cast<protein*>(pmol)->getNumChains();
  for (UInt i = 0; i < number; i++) 
    {
      actives = static_cast<protein*>(pmol)->getChainPositionVector(i);
      n = actives.size();
      for (UInt j = 0; j < n; j++)
	{
	  if (actives[j]) 
	    {
	      count++;
	      cout << "the number of residues: " << n << " the number of actives: " << count << "\n";
	    }
	}
    }
  size = count;
  pindex = new rules[count];
  count = 0;
  for (UInt i = 0; i < number; i++)
    {
      actives = static_cast<protein*>(pmol)->getChainPositionVector(i);
      n = actives.size();
      for (UInt j = 0; j < n; j++)
	if (actives[j])
	  {
	    pindex[count] = rules(actives[j], (int)i, (int)j);
#ifdef MYGENOME_DEBUG
	    cout << "numbers being fed in:  " << i << "  " << j << "\n";
#endif
	    count++;
	  }
    }

}
