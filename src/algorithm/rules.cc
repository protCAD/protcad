#include "rules.h"

rules::rules()
{
  size = 0;
}

rules::rules(chainPosition * _pos, int _chain, int _residue)
{
  flag = 2;
  pos = _pos;
  vector< vector< vector< UInt > > > DB;
  DB = _pos->getAllowedDB();
  UInt maxRes = DB.size();
  UInt maxRots;
  UInt rotamer;
  superResidue temp;
  chain = _chain;
  residue = _residue;

  for (UInt i = 0; i < maxRes; i++)
    {
      if (DB[i].size() > 0)
	{
	  maxRots = DB[i][0].size();
	  for (UInt j = 0; j < maxRots; j++)
	    {
	      rotamer = (int) DB[i][0][j];
	      temp = superResidue(i, rotamer);
	      //cout << "made a residue\n";
	      //temp.printRes();
	      data.push_back(temp);
	    }
	}
    }
  size = data.size();
}

rules::~rules()
{
}

rules rules::operator=(const rules &rhs)
{
  if (rhs.getFlag() == 2)
    {
      flag = rhs.getFlag();
      size = rhs.getSize();
      pos = rhs.getPosition();
      data = rhs.getData();
      chain = rhs.getChain();
      residue = rhs.getResidue();
    }
  return rhs;
}

int rules::getNewValue()
{
  int n = GARandomInt(0, size);
  return n;
}

void rules::print()
{
  if (flag == 2) 
    {
      for (int i = 0; i < size; i++)
	{
	  data[i].printRes();
	}
    }
}
void rules::print(int n)
{
	if (flag == 2)
		{
			int rot, id;
			superResidue sres = data[n];
			id = sres.getIdentity();
			rot = sres.getRotamer();
			cout << "in chain " << chain << " at residue ";
			cout << residue << " the id is: " << id;
			cout << "  and the rotamer is: " << rot << "\n";
		}
}
superResidue rules::translateRes(int _orig)
{
  if (flag == 2) 
    {
      superResidue ret;
      ret = data[_orig];
      return ret;
    }
  return superResidue();
}
