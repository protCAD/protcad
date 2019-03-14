	#include "typedef.h"
	#include "CMath.h"

int main()
{
	dblVec vec(3);
	vec[0] = 1.0;
	vec[1] = 2.0;
	vec[2] = 3.0;
	dblMat mat(3,3);

	double counter = 3.0;
	for (int i=0; i<3; i++)
	{	for (int j=0; j<3; j++)
		{	mat[i][j] = double(i*j) * counter;
			counter = counter + 2.0;
		}
	}

	cout << vec << endl;

	cout << mat << endl;

	dblVec vec2(3);
	//using namespace CMath;
	vec2 = CMath::dotProduct(mat,vec);
	cout << vec2 << endl;

	return 1;
}	
