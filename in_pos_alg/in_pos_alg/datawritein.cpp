#include <iostream>
#include "datawritein.h"
using namespace std;
datareadinfo::datareadinfo(int nrnum,int rownum,int column)
{
	nr = nrnum;
#ifdef debug
	cout << nr << endl;
#endif
	measuredata = new float[nrnum];
	//cout << "measuredata:\t" << measuredata << endl;
	basestation_coor = new float[nrnum * 3];
	//cout << "basestation_coor:\t" << basestation_coor << endl;
	planes = new float *[column];
	for (int i = 0; i < column; i++)
		planes[i] = new float[rownum];
	//cout << "palanes:\t" << planes << endl;
}

datareadinfo::~datareadinfo()
{
	delete[] measuredata;
	delete[] basestation_coor;

	int col = (sizeof(planes) / sizeof(float)) / (sizeof(planes[0]) / sizeof(float));
	for (int j = 0; j < col; j++)
		delete[] planes[j];
	delete[] planes;
}
void datareadinfo::Set_Measure_Data(float *measuredatatemp)
{
	measuredata = measuredatatemp;
#ifdef debug
	cout << *measuredata << endl;
#endif
}
void datareadinfo::Set_Basestation_Coor(float *basestation_coor_temp)
{
	basestation_coor = basestation_coor_temp;
#ifdef debug
	cout << *basestation_coor << endl;
#endif // debug
}
void datareadinfo::Read_1dTo2d(float *temp1d, int temph, int templ, float &taget,int maxcolumn)
{
	taget = temp1d[temph*maxcolumn + templ];
}