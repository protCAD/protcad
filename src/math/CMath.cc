// filename: CMath.cpp
// contents: CMath.h defined functions' implementation

#include "CMath.h"
#include <stdio.h>

double CMath::dotProduct(const dblVec& _vec1, const dblVec& _vec2)
{
	int N = _vec1.dim();
	ASSERT (N == _vec2.dim() );
	double sum = 0.0;
	for (int i=0; i<N; i++)
	{
		sum += _vec1[i] * _vec2[i];
	}
	return sum;
}

dblVec CMath::dotProduct(const dblMat& _mat1, const dblVec& _vec1)
{
	int N = _vec1.dim();
	ASSERT (N == _mat1.num_cols() );
	dblVec result(N);
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			result[i] += _vec1[j] * _mat1[i][j];
		}
	}
	return result;
}

dblVec CMath::translate(const dblVec* _pdv1,const dblVec* _pdv2)
{	dblVec dv3 = *(_pdv1) + *(_pdv2);
	return dv3;
}

dblVec CMath::translate(const dblVec* _pdv1,const dblVec& _dv2)
{	dblVec dv3 = *(_pdv1) + _dv2;
	return dv3;
}

dblVec CMath::translate(const dblVec& _dv1, const dblVec* _pdv2)
{	dblVec dv3 = _dv1 + *(_pdv2);
	return dv3;
}

dblVec CMath::translate(const dblVec& _dv1,const dblVec& _dv2)
{	dblVec dv3 = _dv1 + _dv2;
	return dv3;
}

dblVec CMath::transform(const dblVec* _pdv,const dblMat* _pdm)
{
	dblVec dv2 = dotProduct( (*_pdm) , (*_pdv) );
	return dv2;
}

dblVec CMath::transform(const dblVec* _pdv,const dblMat& _dm)
{
	dblVec dv2 = dotProduct( _dm, (*_pdv) );
	return dv2;
}

dblVec CMath::transform(const dblVec& _dv,const dblMat* _pdm)
{
	dblVec dv2 = dotProduct( (*_pdm) , _dv );
	return dv2;
}

dblVec CMath::transform(const dblVec& _dv,const dblMat& _dm)
{
	dblVec dv2 = dotProduct(_dm, _dv);
	return dv2;
}

double CMath::distance( const dblVec* _pdv1, const dblVec* _pdv2)
{
/*
	return	sqrt(   pow( (*_pdv1)[0] - (*_pdv2)[0], 2) +
			pow( (*_pdv1)[1] - (*_pdv2)[1], 2) +
			pow( (*_pdv1)[2] - (*_pdv2)[2], 2));
*/
	return	sqrt(   ( ((*_pdv1)[0] - (*_pdv2)[0]) * ((*_pdv1)[0] - (*_pdv2)[0]) ) +
			( ((*_pdv1)[1] - (*_pdv2)[1]) * ((*_pdv1)[1] - (*_pdv2)[1]) ) +
			( ((*_pdv1)[2] - (*_pdv2)[2]) * ((*_pdv1)[2] - (*_pdv2)[2]) )  );
}

double CMath::distance(const dblVec& _dv1,const dblVec* _pdv2)
{
/*
	return	sqrt(   pow( _dv1[0] - (*_pdv2)[0], 2) +
			pow( _dv1[1] - (*_pdv2)[1], 2) +
			pow( _dv1[2] - (*_pdv2)[2], 2));
*/
	return	sqrt(   ( (_dv1[0] - (*_pdv2)[0]) * (_dv1[0] - (*_pdv2)[0]) ) +
			( (_dv1[1] - (*_pdv2)[1]) * (_dv1[1] - (*_pdv2)[1]) ) +
			( (_dv1[2] - (*_pdv2)[2]) * (_dv1[2] - (*_pdv2)[2]) )  );
}

double CMath::distance(const dblVec* _pdv1, const dblVec& _dv2)
{
/*
	return	sqrt(   pow( (*_pdv1)[0] - _dv2[0], 2) +
			pow( (*_pdv1)[1] - _dv2[1], 2) +
			pow( (*_pdv1)[2] - _dv2[2], 2));
*/
	return	sqrt(   ( ((*_pdv1)[0] - _dv2[0]) * ((*_pdv1)[0] - _dv2[0]) ) +
			( ((*_pdv1)[1] - _dv2[1]) * ((*_pdv1)[1] - _dv2[1]) ) +
			( ((*_pdv1)[2] - _dv2[2]) * ((*_pdv1)[2] - _dv2[2]) )  );
}


double CMath::distance(const dblVec& _dv1, const dblVec& _dv2)
{
/*
	return	sqrt(   pow( _dv1[0] - _dv2[0], 2) +
			pow( _dv1[1] - _dv2[1], 2) +
			pow( _dv1[2] - _dv2[2], 2));
*/
	return	sqrt(   ( (_dv1[0] - _dv2[0]) * (_dv1[0] - _dv2[0]) ) +
			( (_dv1[1] - _dv2[1]) * (_dv1[1] - _dv2[1]) ) +
			( (_dv1[2] - _dv2[2]) * (_dv1[2] - _dv2[2]) )  );
}

double CMath::distanceSquared(const dblVec& _dv1, const dblVec& _dv2)
{
/*
	return	pow( _dv1[0] - _dv2[0], 2) +
		pow( _dv1[1] - _dv2[1], 2) +
		pow( _dv1[2] - _dv2[2], 2);
*/
	return	   ( (_dv1[0] - _dv2[0]) * (_dv1[0] - _dv2[0]) ) +
			( (_dv1[1] - _dv2[1]) * (_dv1[1] - _dv2[1]) ) +
			( (_dv1[2] - _dv2[2]) * (_dv1[2] - _dv2[2]) ) ;
}

double CMath::dihedral(const dblVec& _V_a, const dblVec& _V_b,
			const dblVec& _V_c, const dblVec& _V_d)

// Adapted from F77 program Propak by J.W. Ponder
{	dblVec V_ab(3);
	dblVec V_bc(3);
	dblVec V_cd(3);

	V_ab = _V_b - _V_a;
	V_bc = _V_c - _V_b;
	V_cd = _V_d - _V_c;

	// Calculate Cross Products

	dblVec V_abXbc(3);
	dblVec V_bcXcd(3);

	V_abXbc = CMath::cross( V_ab, V_bc );
	V_bcXcd = CMath::cross( V_bc, V_cd );

	double normV_abXbc = dotProduct( V_abXbc, V_abXbc );
	double normV_bcXcd = dotProduct( V_bcXcd, V_bcXcd );

	double angle;
	double radToDeg = 57.29577951308;
	double angleCos;
	if ( normV_abXbc > 0.0 && fabs(normV_bcXcd) > 0.0)
	{
		angleCos = dotProduct( V_abXbc, V_bcXcd )/sqrt( normV_abXbc * normV_bcXcd );
	}
	else
	{
		angle = 1000.00;
		cout << "no dihedral" << endl;
		return angle;
	}

	if (fabs(angleCos) > 1.0) angleCos = 1.0 * (angleCos/fabs(angleCos));

	double angleAbs;
	angleAbs  = radToDeg * acos(angleCos);

	double signTest;
	signTest = dotProduct( V_bcXcd, CMath::cross( V_bc, V_abXbc ) );

	if (signTest > 0.0)
	{	angle = angleAbs;
	}
	else
	{	angle = -angleAbs;
	}

	return angle;
}

double CMath::angle(const dblVec& _one, const dblVec& _two, const dblVec& _three)
{
	dblVec first = _one - _two;
	dblVec second = _three - _two;

	double dotProd;
	double norm1=0.0;
	double norm2=0.0;

	dotProd = dotProduct(first,second);
	norm1 = sqrt(dotProduct(first,first));
	norm2 = sqrt(dotProduct(second,second));

	double angle;
	double radToDeg = 57.29577951308;
	double angleCos;


	if ( fabs(norm1) >= 0.01 && fabs(norm2) >= 0.01 )
	{
		angleCos = dotProd / (norm1 * norm2);
	}
	else
	{
		angle = 1000.00;
		cout << "Angle undefined" << endl;
		return angle;
	}

	if (fabs(angleCos) >= 1.0) angleCos = 1.0 * (angleCos/fabs(angleCos));

	double angleAbs;
	angleAbs  = radToDeg * acos(angleCos);

	double signTest;
	signTest = 1.0;

	angle = 0.0;

	if (signTest == 0.0)
	{
		angle = angleAbs;
	}
	else
	{
		angle = angleAbs * (signTest/fabs(signTest));
	}

	return angle;
}

dblMat CMath::rotationMatrix(const dblVec& _dv,const double& _theta)
{	dblMat R(3,3,0.0);

// Theta should be passed in as an angle in degrees!!!
// Here we convert _theta to radians for upcoming calculation....

	double degToRad = 0.017453293;
	double theta = degToRad * _theta;
	double sinTheta = sin(theta);
	double cosTheta = cos(theta);
	double norm = sqrt(dotProduct(_dv,_dv));
	double n1 = _dv[0]/norm;
	double n2 = _dv[1]/norm;
	double n3 = _dv[2]/norm;
	double n11 = n1*n1;
	double n12 = n1*n2;
	double n13 = n1*n3;
	double n22 = n2*n2;
	double n23 = n2*n3;
	double n33 = n3*n3;

	R[0][0] = n11 + (1-n11)*cosTheta;
	R[0][1] = n12*(1-cosTheta) - n3*sinTheta;
	R[0][2] = n13*(1-cosTheta) + n2*sinTheta;

	R[1][0] = n12*(1-cosTheta) + n3*sinTheta;
	R[1][1] = n22 + (1-n22)*cosTheta;
	R[1][2] = n23*(1-cosTheta) - n1*sinTheta;

	R[2][0] = n13*(1-cosTheta) - n2*sinTheta;
	R[2][1] = n23*(1-cosTheta) + n1*sinTheta;
	R[2][2] = n33 + (1-n33)*cosTheta;

	return R;
}

dblVec CMath::centroid(const vector<dblVec>* thePoints, const vector<double>* theWeights)
{
	unsigned int vsize = ((*thePoints)[0]).dim();
	dblVec rv(vsize,0.0);

	if (vsize != 3)
	{
			cout << "You've passed function centroid some ";
			cout << "nonstandard vector... extent " << vsize << endl;
			return rv;
	}

	unsigned int np = thePoints->size();
	if (np != theWeights->size())
	{
			cout << "Error in correspondence between Points (";
			cout << np << ") and Weights (" << theWeights->size();
			cout << ") in function centroid." << endl;
			return rv;
	}

	double totalWeight = 0.0;
	for (unsigned int i=0; i < np; i++)
	{
			totalWeight += (*theWeights)[i];
			for (unsigned int j=0; j<vsize; j++)
			{
					rv[j] += ((*thePoints)[i])[j] * (*theWeights)[i] ;
			}
	}

	for (unsigned int j=0; j<vsize; j++)
	{
			rv[j] /= totalWeight;
	}
	return rv;

}

dblVec CMath::cross(const dblVec& _one, const dblVec& _two)
{
	dblVec result(3);
	result[0] =   (_one[1] * _two[2]) - (_one[2] * _two[1]);
	result[1] =   (_one[2] * _two[0]) - (_one[0] * _two[2]);
	result[2] =   (_one[0] * _two[1]) - (_one[1] * _two[0]);
	return result;
}

double CMath::determinant(const dblMat& _mat)
{
	double determinant;
	determinant = _mat[0][0]* (_mat[1][1]*_mat[2][2] - _mat[2][1]*_mat[1][2]) +
                      _mat[1][0]* (_mat[2][1]*_mat[0][2] - _mat[0][1]*_mat[2][2]) +
                      _mat[2][0]* (_mat[0][1]*_mat[1][2] - _mat[1][1]*_mat[0][2]);
        return determinant;
}

double CMath::determinant(const dblMat* _pMat)
{
	double determinant;
	determinant = (*_pMat)[0][0]* ((*_pMat)[1][1]*(*_pMat)[2][2] - (*_pMat)[2][1]*(*_pMat)[1][2]) +
                       (*_pMat)[1][0]* ((*_pMat)[2][1]*(*_pMat)[0][2] - (*_pMat)[0][1]*(*_pMat)[2][2]) +
                       (*_pMat)[2][0]* ((*_pMat)[0][1]*(*_pMat)[1][2] - (*_pMat)[1][1]*(*_pMat)[0][2]);
        return determinant;
}

double CMath::linearInterpolate(const double _m, const double _b, const
			double _x)
{

	return _m*_x + _b;
}

