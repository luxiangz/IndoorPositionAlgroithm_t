#pragma once
// analysis the data.txt
//struct dataread
//{
//	struct info
//	{
//		int nr;//the number of basestation
//		//int txpoint
//	};
//
//
//};
//#define debug
class datareadinfo
{
public:
	datareadinfo(int nrnum, int rownum, int column);
	~datareadinfo();
	void Set_Measure_Data(float *measuredatatemp);
	void Set_Basestation_Coor(float *basestation_coor_temp);
	void Read_1dTo2d(float *temp1d, int temph, int templ, float &taget, int maxcolumn);
	void After_R_W2Planes();
//private:
	int nr;
	float *measuredata;
	float *basestation_coor;
	float **planes;
};

