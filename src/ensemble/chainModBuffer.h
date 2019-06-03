#include "assert.h"
#include <string.h>
#include <vector>
#include "typedef.h"
#include "enums.h"
#include "chainPosition.h"
#include "residue.h"


#ifndef CHAIN_MOD_BUFFER
#define CHAIN_MOD_BUFFER
class chainModBuffer
{
	public:
		chainModBuffer();
		~chainModBuffer();
		// deep copy constructor
		chainModBuffer(const chainModBuffer& _rhs);

//  Undo and buffering accessors

		void printAll();

		void setIndexInChainBuffer(int _index);
		void resetIndexInChainBuffer();
		int  getIndexInChainBuffer() const;

		void setRotamerIndexBuffer(vector<UInt> _rotamer);
		void resetRotamerIndexBuffer();
		vector<UInt> getRotamerIndexBuffer() const;

		void setResidueIdentityBuffer(int _type);
		void resetResidueIdentityBuffer();
		int getResidueIdentityBuffer() const;

		void setSidechainDihedralAngleBuffer(vector< vector< double > > theAngles);
		void resetSidechainDihedralAngleBuffer();
		vector< vector< double> > getSidechainDihedralAngleBuffer() const;
		
		void setBackboneDihedralAngleBuffer(vector< double > theAngles);
		void resetBackboneDihedralAngleBuffer();
		vector< double> getBackboneDihedralAngleBuffer() const;
		
		void resetAllBuffers();
		bool containsData() { return setFlag;}

	private:
		int indexInChainBuffer;
		vector< UInt> rotamerIndexBuffer;
		vector< vector<double> > sidechainDihedralAngleBuffer;
		vector<double> backboneDihedralAngleBuffer;
		int residueIdentityBuffer;
		bool setFlag;
};
#endif
