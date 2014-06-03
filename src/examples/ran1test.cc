#include "ran1.cpp"
#include <iostream>

int main()
{	
	int histogram1[100];
	for (int i = 0;i<100;i++)
		histogram1[i] = 0;
	ran1* jason = new ran1(33234);
	for(unsigned int round =0;round < 5000; round ++)
	{	
		double temp;
		temp = jason->getNext();
		double bin;
		bin = 0.0;
		for (int i = 0;i<100;i++)
		{
			if (temp >= bin &&
			    temp < bin + 0.01)
				histogram1[i]++;
			bin += 0.01;
		}
//		cout << temp << endl;
	}

	cout << endl;
	for (int i=0;i<100;i++)
		cout << histogram1[i] << endl;

//	int histogram2[11];
//	for (int i = 0;i<11;i++)
//		histogram2[i] = 0;
/*
	//jason->initialize(432);
	for(unsigned int round =0;round < 10000; round ++)
	{	
		int temp;
		temp = jason->getNext(5,7);
		double bin;
		bin = 0.0;
		for (int i = 0;i<11;i++)
		{
			if (temp == i)
				histogram2[i]++;
		}
		cout << temp << endl;
	}
	cout << endl;
	for (int i=0;i<11;i++)
		cout << histogram2[i] << endl;
*/
	for (int i=0;i<100;i++)
	cout << jason->getNext(-20.0,20.0) << endl;
        delete jason;
	return 0;
}
