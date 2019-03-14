#include <pdbReader.h>
#include <pdbWriter.h>
#include <atomIterator.h>
#include <ran.h>
#include <time.h>
#include <myGenome.h>

int menu();

int main()
{
  time_t startTime,endTime,runTime;
  time(&startTime);
  string filename;
  int _x = 0;
  UInt n;
  float nn;

  cout << "enter the name of a pdb file: ";
  cin >> filename;
  
  molecule* pTheProtein1 = pdbReader(filename);
  if (pTheProtein1 == 0)
    {	
      return 1;
    }
  cout << "number of residues generated: " << residue::getHowMany() << endl;
  cout << "number of chains generated:" << chain::getHowMany() << endl;
  while (_x >= 0)
    {
      cout << "enter a chain position to activate, or -1 to end: ";
      cin >> _x;
      if (_x >= 0) 
	{
	  static_cast<protein*>(pTheProtein1)->activateChainPosition(0,_x);
	  static_cast<protein*>(pTheProtein1)->setOnlyROCHydrophobic(0,_x);
	}
    }
  static_cast<protein*>(pTheProtein1)->setAllAlphaAminoAcids();

  myGenome theGA = myGenome(pTheProtein1);
  _x = 0;
  while (_x != 7) 
    {
      _x = menu();
      switch (_x)
	{
	case (1):
	  cout << "enter the new poputation: ";
	  cin >> n;
	  theGA.setPopSize(n);
	  break;
	case (2):
	  cout << "enter the new number of generations: ";
	  cin >> n;
	  theGA.setNumGen(n);
	  break;
	case (3):
	  cout << "enter new probability of mutation: ";
	  cin >> nn;
	  theGA.setProbMut(nn);
	  break;
	case (4):
	  cout << "enter new probability of crossover: ";
	  cin >> nn;
	  theGA.setProbCross(nn);
	  break;
	case (5):
	  cout << "enter new random seed: ";
	  cin >> n;
	  theGA.setSeed(n);
	  break;
	case (6):
	  cout << "the size is: " << theGA.getPopSize() << "\n";
	  cout << "the number of generations is: " << theGA.getNumGen() << "\n";
	  cout << "the probability of crossover and mutation is: ";
	  cout << theGA.getProbCross() << " and " << theGA.getProbMut() << "\n";
	  cout << "the seed is: " << theGA.getSeed() << "\n";
	  break;
	case (7):
	  break;
	default:
	  cout << "please enter a valid option.\n";
	  break;
	}
    }

  
  theGA.run();
  
  string outfilename = "JoshsBaby.pdb";
  pdbWriter(static_cast<protein*>(pTheProtein1),outfilename);
  
  delete pTheProtein1;
  
  time(&endTime);
  outfilename = "gaTime.dat";
  ofstream oFile;
  oFile.open(outfilename.c_str());
  oFile << "startTime = " << startTime << "\n";
  oFile << "endTime = " << endTime << "\n";
  runTime = endTime - startTime;
  oFile << "runTime = " << runTime << "\n";
  oFile.close();
  return 0; 
  
} 

int menu()
{
  int n;
  cout << "enter a menu option:\n";
  cout << "1: set population size\n";
  cout << "2: set number of generations\n";
  cout << "3: set the probability of mutation\n";
  cout << "4: set the probability of a crossover\n";
  cout << "5: specify a random seed\n";
  cout << "6: view current settings\n";
  cout << "7: run the genetic algorithm\n";
  cout << ": ";
  cin >> n;
  return n;
}
