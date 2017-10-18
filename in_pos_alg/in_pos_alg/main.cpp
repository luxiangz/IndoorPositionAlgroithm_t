
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>

#include "stdlib.h"
#include "file_operation.h"
#include "datawritein.h"
#include "total_reflection.h"
using namespace std;
struct potcollection
{
	int firstp;
	int secondp;
	int thirdp;
	int fourthp;
};
struct coordinatedata
{
public:
	float x = {};
	float y = {};
	float z = {};
	/*potcollection *index = {};
	coordinatedata(int pointnum0, int formsurface0)
	{
		float *x = new float[pointnum0];
		float *y = new float[pointnum0];
		float *z = new float[pointnum0];
		potcollection *index = new potcollection[pointnum0];
		cout << "struct构造函数" << endl;
	}
	~coordinatedata()
	{
		delete[] x;
		delete[] y;
		delete[] z;
		delete[] index;
		cout << "struct析构函数" << endl;
	}*/

};
int formsurface[] = {
4,5,6,1,
5,8,7,6,
8,14,11,7,
16,14,10,11,
15,16,11,12,
13,15,12,9,
3,13,9,2,
3,4,1,2,
21,22,19,20,
24,21,20,17,
24,23,18,17,
22,23,18,19,
20,19,18,14
};

#define pointnum 24
#define dimension 3
float xtemp[pointnum] = { 000.0f, 000.0f, 000.0f, 000.0f, 491.5f, 491.5f, 491.5f, 491.5f, 089.2f, 108.3f, 108.3f, 089.2f, 089.2f, 108.3f, 089.2f, 108.3f, 200.0f, 350.0f, 350.0f, 200.0f, 200.0f, 350.0f, 350.0f, 200.0f };
float ytemp[pointnum] = { 000.0f, 483.7f, 483.7f, 000.0f, 000.0f, 000.0f, 483.7f, 483.7f, 483.7f, 483.7f, 458.1f, 458.1f, 483.7f, 483.7f, 458.1f, 458.1f, 260.0f, 260.0f, 200.0f, 200.0f, 200.0f, 200.0f, 260.0f, 260.0f };
float ztemp[pointnum] = { 200.0f, 200.0f, 000.0f, 000.0f, 000.0f, 200.0f, 200.0f, 000.0f, 200.0f, 200.0f, 200.0f, 200.0f, 000.0f, 000.0f, 000.0f, 000.0f, 050.0f, 050.0f, 050.0f, 050.0f, 000.0f, 000.0f, 000.0f, 000.0f };

float Cal_Pol(int clie, int chang);
void Perfect_planes(int maxmatrixl, int maxmatrixh, float **planes);
int Write_Into_Txt(char* filename, bool writeflag);
float maxmax(float temp1, float temp2, float temp3, float temp4);
float minmin(float temp1, float temp2, float temp3, float temp4);
void findborder(int maxmatrixl, float **planes, int *line, float & min3, float & max3);
void Set_vect(float * vect, int nb, float min4, float max4);
int main()
{
	//char *a=NULL;
	//cout << "a的长度   "<<sizeof(a) << endl;
	char *filename = "./data.txt";//name of file(文件名称地址)
	bool writeflag = false;//是否写入数据标志位
	if (Write_Into_Txt(filename, writeflag) == 0)
	{
		return -1;
	}

	fstream datafile;
	int pli_norm_num = 1;
	int distance_origin_num = 1;
	int plpoi_maxmin_num = 2;
	int planepointnum = 4;//the number of points consituted a planes
	int column = (planepointnum + plpoi_maxmin_num + pli_norm_num)*dimension + distance_origin_num;
	int row = CountDataRow(filename);//13
	const int basestationnum = 4;
	datareadinfo datatxt(basestationnum,row,column);
#ifdef debug
	cout << "row " << row << endl;
#endif // debug
	float distanceme[basestationnum] = { 7.5f, 895.4f, 56.9f, 78.5f};
	float basecoor_temp[basestationnum][dimension] = {
		{ 0.0f, 12.0f, 16.0f },
		{ 1.0f, 13.0f, 17.0f },
		{ 2.0f, 14.0f, 18.0f },
		{ 3.0f, 15.0f, 19.0f } };
	for (int h = 0; h < basestationnum; h++)
		for (int l = 0; l < dimension;l++)
			datatxt.basestation_coor[h][l] = basecoor_temp[h][l];
	//cout << datatxt.basestation_coor[2][0] << " ceshi30" << endl;
	datatxt.Set_Measure_Data(distanceme);
	//datatxt.Set_Basestation_Coor(basecoor);
	float tagetcoor = 0.0f;
#ifdef debug
	datatxt.Read_1dTo2d(basecoor, 3, 2, tagetcoor,dimension);//prevent exceeding the array size
	cout << "tagetcoor: "<<tagetcoor << endl;
#endif // debug

	datafile.open(filename, ios::in);//ios::in represents read,ios::out represents write,and once this step is done,the contents of the file wiil be cleared.
	if (!datafile)
	{
		cerr << "Error opening file,read!" << endl;
		return -2;
	}
	//int num_sort = 0;
	if (!datafile.fail())
	{
		int columntemp = 0;
		for (int readhighline = 0 ; readhighline < row; readhighline++)
		{
			string linetemp = ReadDataRow(filename, readhighline+1);
			string floatemp;
			//cout << linetemp << endl;
			for (auto c : linetemp)
			{
				if (isspace(c))
				{//equivalent to the transfer of data in the file assigned to the array转置
					datatxt.planes[columntemp][readhighline] = (float)atof(floatemp.c_str());
					//float atest = (float)atof(floatemp.c_str());
					//cout << atest << endl;
					//cout << "datatex.p" << datatxt.planes << endl; //[columntemp][readhighline]
					//cout << "datatex.b" << datatxt.basestation_coor << endl;
					//cout << "datatex.m" << datatxt.measuredata << endl;
					floatemp = "";
					columntemp++;
					//cout << floatemp << endl;
				}
				else
				{
					floatemp += c;
				}
			}
			columntemp = 0;
		}
		//After_R_W2Planes();	
	}
	datafile.close();//read the data in the file
#ifdef debug
	for (int i = 0; i < column; ++i)
	{
		for (int j = 0; j < row; ++j)
		{
			cout << datatxt.planes[i][j] << "\t";
		}
		cout << endl;
	}
	cout << "maxmax" << maxmax(17.0f, 14.0f, 18.0f, 9.0f) << endl;
	cout << "minmin" << minmin(12.0f, 14.0f, 13.0f, 8.0f) << endl;
#endif // debug
	Perfect_planes(row, column, datatxt.planes);
#ifdef debug
	for (int i = 0; i < column; ++i)
	{
		for (int j = 0; j < row; ++j)
		{
			cout << datatxt.planes[i][j] << "\t";
		}
		cout << endl;
	}
#endif // debug
	//load_flp^|
	//建立测试点程序
	//float *temp;
	float xmin = 0.0f, xmax = 0.0f, ymin = 0.0f, ymax = 0.0f, zmin = 0.0f, zmax = 0.0f;
	int lineofx[3] = { 1, 4, 7 };
	findborder(row,datatxt.planes,lineofx, xmin, xmax);
	//cout << "xmin " << xmin << endl;
	//cout << "xmax " << xmax << endl;
	int lineofy[3] = { 2, 5, 8 };
	findborder(row, datatxt.planes, lineofy, ymin, ymax);
	int lineofz[3] = { 3, 6, 9 };
	findborder(row, datatxt.planes, lineofz, zmin, zmax);
	
	int nb_pts = 100;
	float *vectx = new float[nb_pts];
	Set_vect(vectx, nb_pts, xmin, xmax);
	//for (int i = 0; i < nb_pts; i++)
	//	cout << vectx[i] << endl;
	float *vecty = new float[nb_pts];
	Set_vect(vecty, nb_pts, ymin, ymax);
	float t_FlatH = (zmax - zmin) / 2;

	float **meshgrid = new float *[dimension+1];
	for (int len = 0; len < (dimension+1); len++)
		meshgrid[len] = new float[nb_pts*nb_pts];
	for (int xnum = 0; xnum < nb_pts; xnum++)
	{
		for (int ynum = 0; ynum < nb_pts; ynum++)
		{
			meshgrid[0][100 * xnum + ynum] = vectx[xnum];
			meshgrid[1][100 * xnum + ynum] = vecty[ynum];
			meshgrid[2][100 * xnum + ynum] = t_FlatH;
		}
	}

#ifdef debug
	for (int i = 0; i < 201; i++)
	{
		cout << meshgrid[0][i] << " " << meshgrid[1][i] << " " << meshgrid[2][i] << endl;
	}
#endif
	//setup three reflections
	int nr = 3;//the times of reflection
	int Nr = nb_pts*nb_pts;
	float *rs_amp = new float[basestationnum*Nr];//total power Nt*Nr

	Mfct_Ref(dimension, datatxt.planes, datatxt.basestation_coor, meshgrid, basestationnum, Nr, row, nr, rs_amp);//material//fc//flp_scale
		//(getdimensions,datatxt.planes，**basecoor_temp0**meshgrid, basestationnum,nb_pts, row,     nr)




	for (int len = 0; len < dimension; len++)
		delete[] meshgrid[len];
	delete[] meshgrid;
	delete[] vectx;
	delete[] vecty;
	delete[] rs_amp;
 	system("pause");
	return 0;
	
}
void Set_vect(float * vect, int nb, float min4, float max4)
{
	float step = (max4 - min4) / (nb - 1);
	vect[0] = min4;
	vect[nb - 1] = max4;
	for (int i = 1; i < nb - 1; i++)
		vect[i] = min4 + step*i;
	
}
void findborder(int maxmatrixl, float **planes, int *line, float & min3, float & max3)
{
	int elecut = maxmatrixl*dimension;
	float *temp=new float[elecut];
	int current = 0;
	for (int j = 0; j < dimension; j++)
	{
		for (int i = 0; i < maxmatrixl; i++)
		{
			//cout << "line" << line[j];
			temp[current] = planes[line[j]][i];
			//cout << *temp << " " << endl;
			current++;
		}
	}
	min3 = *min_element(temp, temp + (elecut - 1));
	max3 = *max_element(temp, temp + (elecut - 1));
	delete[] temp;
}
float maxmax(float temp1, float temp2, float temp3, float temp4)
{
	if (temp1 < temp2)
		temp1 = temp2;
	if (temp1 < temp3)
		temp1 = temp3;
	if (temp1 < temp4)
		temp1 = temp4;
	return temp1;
}
float minmin(float temp1, float temp2, float temp3, float temp4)
{
	if (temp1 > temp2)
		temp1 = temp2;
	if (temp1 > temp3)
		temp1 = temp3;
	if (temp1 > temp4)
		temp1 = temp4;
	return temp1;
}

void Perfect_planes(int maxmatrixl, int maxmatrixh,float **planes)
{
	for (int i = 0; i < maxmatrixl; ++i)
	{
		if (maxmatrixh==22)
		{
			float pli_x1 = planes[0][i];
			float pli_y1 = planes[1][i];
			float pli_z1 = planes[2][i];

			float pli_x2 = planes[3][i];
			float pli_y2 = planes[4][i];
			float pli_z2 = planes[5][i];

			float pli_x3 = planes[6][i];
			float pli_y3 = planes[7][i];
			float pli_z3 = planes[8][i];

			float pli_x4 = planes[9][i];
			float pli_y4 = planes[10][i];
			float pli_z4 = planes[11][i];

			float pli_norm_x =
				(pli_y1 - pli_y2)*(pli_z1 + pli_z2) +
				(pli_y2 - pli_y3)*(pli_z2 + pli_z3) +
				(pli_y3 - pli_y1)*(pli_z3 + pli_z1);
			float pli_norm_y =
				(pli_z1 - pli_z2)*(pli_x1 + pli_x2) +
				(pli_z2 - pli_z3)*(pli_x2 + pli_x3) +
				(pli_z3 - pli_z1)*(pli_x3 + pli_x1);
			float pli_norm_z =
				(pli_x1 - pli_x2)*(pli_y1 + pli_y2) +
				(pli_x2 - pli_x3)*(pli_y2 + pli_y3) +
				(pli_x3 - pli_x1)*(pli_y3 + pli_y1);
			
			float mag_vec = (float)sqrt(pli_norm_x*pli_norm_x + pli_norm_y*pli_norm_y + pli_norm_z*pli_norm_z);
			pli_norm_x = pli_norm_x / mag_vec;
			pli_norm_y = pli_norm_y / mag_vec;
			pli_norm_z = pli_norm_z / mag_vec;
			float distance_origin = -(pli_x1*pli_norm_x + pli_y1*pli_norm_y + pli_z1*pli_norm_z);
			
			planes[12][i] = pli_norm_x;
			planes[13][i] = pli_norm_y;
			planes[14][i] = pli_norm_z;
			planes[15][i] = distance_origin;

			float xmaxt = maxmax(pli_x1, pli_x2, pli_x3, pli_x4);
			float xmint = minmin(pli_x1, pli_x2, pli_x3, pli_x4);
			float ymaxt = maxmax(pli_y1, pli_y2, pli_y3, pli_y4);
			float ymint = minmin(pli_y1, pli_y2, pli_y3, pli_y4);
			float zmaxt = maxmax(pli_z1, pli_z2, pli_z3, pli_z4);
			float zmint = minmin(pli_z1, pli_z2, pli_z3, pli_z4);
			planes[16][i] = xmint;
			planes[17][i] = xmaxt;
			planes[18][i] = ymint;
			planes[19][i] = ymaxt;
			planes[20][i] = zmint;
			planes[21][i] = zmaxt;

		} 
		else if (maxmatrixh==2)
		{
		}
		
	}

}
int Write_Into_Txt(char* filename, bool writeflag)
{
	fstream datafile0;
	int surfacenum = sizeof(formsurface) / sizeof(formsurface[0]) / 4;
#ifdef debug
	cout << "surfacenum" << surfacenum << endl;
#endif
	coordinatedata datac[pointnum];
	
	if (writeflag)
	{
		datafile0.open(filename, ios::in | ios::out);
		if (!datafile0)
		{
			cerr << "Can't open the file,write!" << endl;
			return 0;
		}
		if (!datafile0.fail())
		{//select the four point clockwise to form a face
			for (int i = 0; i < pointnum; i++)
			{
				datac[i].x = xtemp[i];
				datac[i].y = ytemp[i];
				datac[i].z = ztemp[i];
#ifdef debug
				cout << "data[" << i << "]\t" << datac[i].x << endl;
				//datafile << setw(2) << datac[i].x << "\t" << datac[i].y << "\t" << datac[i].z << "\n";
#endif
			}
			for (int j = 0; j < surfacenum; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					int tempindex = formsurface[j * 4 + k];
					datafile0 << datac[tempindex - 1].x << "\t" << datac[tempindex - 1].y << "\t" << datac[tempindex - 1].z << "\t";
				}
				datafile0 << "\n";
			}
			//cout << datac[9].x << "\t" << datac[9].y << "\t" << datac[9].z << "\n";
			//datafile0 << "what happend?" << "\n";
			datafile0.close();
		}
	}
	return 1;

}
