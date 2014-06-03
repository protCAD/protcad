void residue::rotateLocal(atom* _pAtom1, atom* _pAtom2, double deltaTheta, double _theta, int distance, int direction)
{
#ifdef __RES_DEBUG
	cout << _pAtom1->getName() << " " << _pAtom1->getCoords() << endl;
	cout << _pAtom2->getName() << " " << _pAtom2->getCoords() << endl;
#endif
	if (direction == 0)
	{	
		#ifdef __RES_DEBUG
		_pAtom2->queryChildrensCoords();
		#endif
		dblVec toOrigin = _pAtom1->getCoords() * (-1.0);
		dblVec backHome = _pAtom1->getCoords();

		_pAtom1->translate(toOrigin);
		_pAtom2->translate(toOrigin);
		_pAtom2->translateChildren(toOrigin);

		if (pItsNextRes)
		{	pItsNextRes->recursiveTranslateLocal(toOrigin, direction);
		}

		#ifdef __RES_DEBUG
		_pAtom2->queryChildrensCoords();
		#endif
		dblVec atomCoords = _pAtom2->getCoords();
		dblMat R(3,3,0.0);
		R = CMath::rotationMatrix(atomCoords, deltaTheta);
		_pAtom2->transformChildren(R);

		if (distance == 0)
		{
			if (pItsNextRes)
			{	
				pItsNextRes->recursiveTransformLocal(atomCoords, deltaTheta, direction);
			}

		}
		if (distance == 1)
		{
			if (pItsNextRes)
			{	
				pItsNextRes->recursiveTransform(R);
			}
		}
		#ifdef __RES_DEBUG
		_pAtom2->queryChildrensCoords();
		#endif
		_pAtom1->translate(backHome);
		_pAtom2->translate(backHome);
		_pAtom2->translateChildren(backHome);

		if (pItsNextRes)
		{	pItsNextRes->recursiveTranslateLocal(backHome, direction);
		}

		#ifdef __RES_DEBUG
		_pAtom2->queryChildrensCoords();
		#endif
	}

	if (direction == 1)
	{
		#ifdef __RES_DEBUG
		_pAtom2->queryChildrensCoords();
		#endif
		dblVec toOrigin = _pAtom1->getCoords() * (-1.0);
		dblVec backHome = _pAtom1->getCoords();

		_pAtom1->translate(toOrigin);
		_pAtom2->translate(toOrigin);
		_pAtom2->translateChildren(toOrigin);

		if (pItsPrevRes)
		{	pItsPrevRes->recursiveTranslateLocal(toOrigin, direction);
		}

		#ifdef __RES_DEBUG
		_pAtom2->queryChildrensCoords();
		#endif
		dblVec atomCoords = _pAtom2->getCoords();
		dblMat R(3,3,0.0);
		R = CMath::rotationMatrix(atomCoords, deltaTheta);
		//_pAtom2->transformChildren(R);

		if (distance == 0)
		{
			if (pItsPrevRes)
			{	
				pItsPrevRes->recursiveTransformLocal(atomCoords, deltaTheta, direction);
			}
		}

		if (distance == 1)
		{
			if (pItsPrevRes)
			{	
				pItsPrevRes->recursiveTransformR(R);
			}
		}
		#ifdef __RES_DEBUG
		_pAtom2->queryChildrensCoords();
		#endif
		_pAtom1->translate(backHome);
		_pAtom2->translate(backHome);
		_pAtom2->translateChildren(backHome);

		if (pItsPrevRes)
		{	pItsPrevRes->recursiveTranslateLocal(backHome, direction);
		}

		#ifdef __RES_DEBUG
		_pAtom2->queryChildrensCoords();
		#endif
	}
}
