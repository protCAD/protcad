#include "chainModBuffer.h"

chainModBuffer::chainModBuffer()
{
	resetAllBuffers();
}

chainModBuffer::~chainModBuffer()
{
}

chainModBuffer::chainModBuffer(const chainModBuffer& _rhs)
{	indexInChainBuffer = _rhs.indexInChainBuffer;
	residueIdentityBuffer = _rhs.residueIdentityBuffer;
	rotamerIndexBuffer = _rhs.rotamerIndexBuffer;
	sidechainDihedralAngleBuffer = _rhs.sidechainDihedralAngleBuffer;
	setFlag = _rhs.setFlag;
}

void chainModBuffer::resetAllBuffers()
{
	resetIndexInChainBuffer();
	resetRotamerIndexBuffer();
	resetResidueIdentityBuffer();
	resetSidechainDihedralAngleBuffer();
	setFlag = false;
}

void chainModBuffer::printAll()
{
	cout << "Identity: " << residueIdentityBuffer << endl;
	cout << "IndexInChain: " << indexInChainBuffer << endl;
	if (rotamerIndexBuffer.size() != 0)
	{
		for (UInt i=0; i< rotamerIndexBuffer.size(); i++)
		{
			cout << "Rotamer Index " << i << ": " << rotamerIndexBuffer[i] << endl;
		}
	}
	else
	{	cout << "Rotamer Index : null" << endl;
	}

	if (sidechainDihedralAngleBuffer.size() != 0)
	{
		for (UInt i=0; i< sidechainDihedralAngleBuffer.size(); i++)
		{
			for (UInt j=0; j< sidechainDihedralAngleBuffer[i].size(); j++)
			{ cout << "Dihedral " << j << " bpt " << i << ": " << sidechainDihedralAngleBuffer[i][j] << endl;
			}
		}
	}
	else
	{
		cout << "Dihedral : null" << endl;
	}
	cout << endl;
}
	
void chainModBuffer::setIndexInChainBuffer(int _index)
{
	indexInChainBuffer = _index;
	setFlag = true;
}

void chainModBuffer::resetIndexInChainBuffer()
{
	indexInChainBuffer = -1;
}

int chainModBuffer::getIndexInChainBuffer() const
{
	return indexInChainBuffer;
}

void chainModBuffer::setRotamerIndexBuffer( vector< UInt> _rotamer)
{
	rotamerIndexBuffer = _rotamer;
	setFlag = true;
}

void chainModBuffer::resetRotamerIndexBuffer()
{
	rotamerIndexBuffer.resize(0);
}

vector <UInt> chainModBuffer::getRotamerIndexBuffer() const
{	
	return rotamerIndexBuffer;
}

void chainModBuffer::setResidueIdentityBuffer(int _id)
{
	residueIdentityBuffer = _id;
	setFlag = true;
}

void chainModBuffer::resetResidueIdentityBuffer()
{
	residueIdentityBuffer = -1;
}

int chainModBuffer::getResidueIdentityBuffer() const
{
	return residueIdentityBuffer;
}
void chainModBuffer::setSidechainDihedralAngleBuffer(vector< vector< double > > theAngles)
{
	sidechainDihedralAngleBuffer = theAngles;
	setFlag = true;
}

void chainModBuffer::resetSidechainDihedralAngleBuffer()
{
	sidechainDihedralAngleBuffer.resize(0);
}

vector< vector< double> > chainModBuffer::getSidechainDihedralAngleBuffer() const
{
	return sidechainDihedralAngleBuffer;
}
