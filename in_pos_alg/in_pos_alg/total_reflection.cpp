#include <cmath>
#include "total_reflection.h"
#include <cstdlib>
#define PI 3.14159265358979323846f
float c = 3e8, Ant = 1.0, loss_factor = 1.0, t_tol = 1.0;
/*--------------------------------------------------------------------------------------*/
void Mfct_Ref(int getdimensions,  float **planes, float **TXpoint, float **RXpoint, int Nt, int Nr, int nplanes, int nr,float *rs_amp)//material//fc//flp_scale
//			 (getdimensions,           datatxt.planes，     **basecoor_temp0,    **meshgrid, basestationnum,     nb_pts,        row,     nr)
{
	float *TX = new float[getdimensions];
	float *RX = new float[getdimensions];
	float *TX_i = new float[getdimensions];
	float *TX_j = new float[getdimensions];
	float *TX_k = new float[getdimensions];
	float *Phiti = new float[getdimensions];
	float *Phitj = new float[getdimensions];
	float *Phitij = new float[getdimensions];
	float *Phitk = new float[getdimensions];
	float *Phitjk = new float[getdimensions];
	float *Phitijk = new float[getdimensions];
	float *Phit_nt = new float[getdimensions];

	//int Nr = nb_pts0*nb_pts0;
	int npath1, npath2, npath3, npath, currentpath, id, jd, kd, nb_pts, indice,ind,ld;

	int tflag0, tflag11, tflag12, tflag21, tflag22, tflag23, tflag31, tflag32, tflag33, tflag34;
	bool rflagi = false, rflagj = false, rflagij = false, rflagk = false, rflagjk = false, rflagijk = false;
	float *rpath1 = NULL, *rpath2 = NULL, *rpath3 = NULL;
	float *path0 = NULL, *path1 = NULL, *path2 = NULL, *path3 = NULL;

	float *Thits0 = NULL;
	float *Thits11 = NULL, *Thits12 = NULL;
	float *Thits21 = NULL, *Thits22 = NULL, *Thits23 = NULL;
	float *Thits31 = NULL, *Thits32 = NULL, *Thits33 = NULL, *Thits34 = NULL;
	float *distance_mfct, *lossfac_mfct;
	//xfloat *rs_amp;
	float temp;
	//output parament
	float **FlatDis = new float *[Nt];
	for (int i = 0; i < Nt; i++)
		FlatDis[i] = new float[Nr];
	//the outermost loop 
	for (int r = 0; r < Nr; r++)
	{
		ld = r*Nt;
		for (int i = 0; i < getdimensions; i++)
			RX[i] = RXpoint[i][r];
		//inner loop of the firing point
		for (int n = 0; n < Nt; n++)
		{
			for (int j = 0; j < getdimensions; j++)
				TX[j] = TXpoint[j][n];
			npath1 = npath2 = npath3 = 0;
			currentpath = 0;
			if (nr > 0)
			{
				for (int cur_nplanes = 0; cur_nplanes < nplanes; cur_nplanes++)
				{
					fct_calimage(TX, planes, cur_nplanes, TX_i);
					rflagi = fct_calrpt(TX_i, RX, planes, cur_nplanes, Phiti);
					if (rflagi == true)
					{
						id = npath1 * 5;
						npath1++;
						rpath1 = (float *)realloc((float *)rpath1, (npath1)* 5 * sizeof(float));
						//if(!rpath1){Error handing}//reallocate memory without changing the original allocation
						rpath1[0 + id] = Phiti[0];
						rpath1[1 + id] = Phiti[1];
						rpath1[2 + id] = Phiti[2];//交点坐标
						rpath1[3 + id] = (float)cur_nplanes;//反射平面序号
						rpath1[4 + id] = 1.0f;
					}
					if (nr > 1)
					{
						for (int nplanes_tempj = 0; nplanes_tempj < nplanes; nplanes_tempj++)
						{
							if (cur_nplanes != nplanes_tempj)//除去一次反射面，剩下所有面未发生反射的发生二次反射
							{
								fct_calimage(TX_i, planes, nplanes_tempj, TX_j);
								rflagj = fct_calrpt(TX_j, RX, planes, nplanes_tempj, Phitj);
								if (rflagj == true)
								{
									rflagij = fct_calrpt(TX_i, Phitj, planes, cur_nplanes, Phitij);
									if (rflagij == true)
									{
										jd = npath2 * 10;
										npath2++;
										rpath2 = (float *)realloc((float *)rpath2, (npath2)* 10 * sizeof(float));
										rpath2[0 + jd] = Phitij[0];
										rpath2[1 + jd] = Phitij[1];
										rpath2[2 + jd] = Phitij[2];
										rpath2[3 + jd] = (float)cur_nplanes;
										rpath2[4 + jd] = 1.0f;
										rpath2[5 + jd] = Phitj[0];
										rpath2[6 + jd] = Phitj[1];
										rpath2[7 + jd] = Phitj[2];
										rpath2[8 + jd] = (float)nplanes_tempj;
										rpath2[9 + jd] = 1.0f;
									}
								}
								if (nr > 2)
								{
									for (int nplanes_tempk = 0; nplanes_tempk < nplanes_tempk; nplanes_tempk++)
									{
										if (nplanes_tempk != nplanes_tempj)
										{
											fct_calimage(TX_j, planes, nplanes_tempk, TX_k);
											rflagk = fct_calrpt(TX_k, RX, planes, nplanes_tempk, Phitk);
											if (rflagk = true)
											{
												rflagjk = fct_calrpt(TX_j, Phitk, planes, nplanes_tempj, Phitjk);
												if (rflagjk == true)
												{
													rflagijk = fct_calrpt(TX_i, Phitjk, planes, cur_nplanes, Phitijk);
													if (rflagijk = true)
													{
														kd = npath3 * 15;
														npath3++;
														rpath3 = (float *)realloc((float *)rpath3, (npath3)* 15 * sizeof(float));
														rpath3[0 + kd] = Phitijk[0];
														rpath3[1 + kd] = Phitijk[1];
														rpath3[2 + kd] = Phitijk[2];
														rpath3[3 + kd] = (float)cur_nplanes;
														rpath3[4 + kd] = 1.0;
														rpath3[5 + kd] = Phitjk[0];
														rpath3[6 + kd] = Phitjk[1];
														rpath3[7 + kd] = Phitjk[2];
														rpath3[8 + kd] = (float)nplanes_tempj;
														rpath3[9 + kd] = 1.0;
														rpath3[10 + kd] = Phitk[0];
														rpath3[11 + kd] = Phitk[1];
														rpath3[12 + kd] = Phitk[2];
														rpath3[13 + kd] = (float)nplanes_tempk;
														rpath3[14 + kd] = 1.0;
													}

												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			npath = npath1 + npath2 + npath3 + 1;
			distance_mfct = new float[npath];
			lossfac_mfct = new float[npath];

			//direct path (no reflection)
			Thits0 = fct_caltpts(TX, RX, planes, nplanes, &tflag0, Phit_nt);
			nb_pts = 2 + tflag0;//产生交点的数目
			path0 = new float[nb_pts * 5];
			path0[0] = TX[0];
			path0[1] = TX[1];
			path0[2] = TX[2];
			path0[3] = 0.0f;
			path0[4] = 0.0f;
			indice = (tflag0 + 1) * 5;
			for (int indice_i = 5; indice_i < indice; indice_i++)
			{
				path0[indice_i] = Thits0[indice_i - 5];
			}
			path0[0 + indice] = RX[0];
			path0[1 + indice] = RX[1];
			path0[2 + indice] = RX[2];
			path0[3 + indice] = 0.0f;
			path0[4 + indice] = 0.0f;
			compute_distloss(path0, nb_pts, planes, currentpath, distance_mfct, lossfac_mfct);

			delete[] Thits0;
			Thits0 = NULL;
			delete[] path0;
			path0 = NULL;//与free不同之处

			//path with 1 reflection
			for (int npathi = 0; npathi < npath1; npathi++)
			{
				currentpath++;
				id = npathi * 5;
				Phiti[0] = rpath1[0 + id];
				Phiti[1] = rpath1[1 + id];
				Phiti[2] = rpath1[2 + id];
				Thits11 = fct_caltpts(TX, Phiti, planes, nplanes, &tflag11, Phit_nt);
				Thits12 = fct_caltpts(Phiti, RX, planes, nplanes, &tflag12, Phit_nt);
				nb_pts = tflag11 + tflag12 + 3;
				path1 = new float[nb_pts * 5];
				path1[0] = TX[0];
				path1[1] = TX[1];
				path1[2] = TX[2];
				path1[3] = 0.0f;
				path1[4] = 0.0f;
				indice = (tflag11 + 1) * 5;
				for (int indice_j = 5; indice_j < indice; indice_j++)
				{
					path1[indice_j] = Thits11[indice_j - 5];
				}
				path1[0 + indice] = Phiti[0];
				path1[1 + indice] = Phiti[1];
				path1[2 + indice] = Phiti[2];
				path1[3 + indice] = rpath1[3+id];
				path1[4 + indice] = rpath1[4+id];
				ind = (tflag11 + 2) * 5;
				indice = (tflag11 + tflag12 + 2) * 5;
				for (int ind_j = ind; ind_j < indice; ind_j++)
				{
					path1[ind_j] = Thits12[ind_j - ind];
				}
				path1[0 + indice] = RX[0];
				path1[1 + indice] = RX[1];
				path1[2 + indice] = RX[2];
				path1[3 + indice] = 0.0f;
				path1[4 + indice] = 0.0f;
				//compute_distloss(path0, nb_pts, planes, currentpath, distance_mfct, lossfac_mfct);
				compute_distloss(path1, nb_pts, planes, currentpath, distance_mfct, lossfac_mfct);
				delete[] Thits11;
				Thits11 = NULL;
				delete[] Thits12;
				Thits12 = NULL;
				delete[] path1;
				path1 = NULL;

			}
			//path with 2 reflection
			for (int npath2i = 0; npath2i < npath2; npath2i++)
			{
				currentpath++;
				id = npath2i * 10;
				Phitij[0] = rpath2[0 + id];
				Phitij[1] = rpath2[1 + id];
				Phitij[2] = rpath2[2 + id];
				Phitj[0] = rpath2[5 + id];
				Phitj[1] = rpath2[6 + id];
				Phitj[2] = rpath2[7 + id];
				Thits21 = fct_caltpts(TX, Phitij, planes, nplanes, &tflag21, Phit_nt);
				Thits22 = fct_caltpts(Phitij, Phitj, planes, nplanes, &tflag22, Phit_nt);
				Thits23 = fct_caltpts(Phitj, RX, planes, nplanes, &tflag23, Phit_nt);
				nb_pts = tflag21 + tflag22 + tflag23 + 4;
				path2 = new float[nb_pts * 5];
				path2[0] = TX[0];
				path2[1] = TX[1];
				path2[2] = TX[2];
				path2[3] = 0.0f;
				path2[4] = 0.0f;
				indice = (tflag21 + 1) * 5;
				for (int indice_j = 5; indice_j < indice; indice_j++)
				{
					path2[indice_j] = Thits21[indice_j - 5];
				}
				path2[0 + indice] = Phitij[0];
				path2[1 + indice] = Phitij[1];
				path2[2 + indice] = Phitij[2];
				path2[3 + indice] = rpath2[3 + id];
				path2[4 + indice] = rpath2[4 + id];
				ind = (tflag21 + 2) * 5;
				indice = (tflag21 + tflag22 + 2) * 5;
				for (int indice_j = ind; indice_j < indice; indice_j++)
				{
					path2[indice_j] = Thits22[indice_j - ind];
				}
				path2[0 + indice] = Phitj[0];
				path2[1 + indice] = Phitj[1];
				path2[2 + indice] = Phitj[2];
				path2[3 + indice] = rpath2[8 + id];
				path2[4 + indice] = rpath2[9 + id];
				ind = (tflag21 + tflag22 + 3) * 5;
				indice = (tflag21 + tflag22 + tflag23 + 3) * 5;
				for (int indice_j = ind; indice_j < indice; indice_j++)
				{
					path2[indice_j] = Thits23[indice_j - ind];
				}
				path2[0 + indice] = RX[0];
				path2[1 + indice] = RX[1];
				path2[2 + indice] = RX[2];
				path2[3 + indice] = 0.0f;
				path2[4 + indice] = 0.0f;
				compute_distloss(path2, nb_pts, planes, currentpath, distance_mfct, lossfac_mfct);
				delete[] Thits21;
				Thits21 = NULL;
				delete[] Thits22;
				Thits22 = NULL;
				delete[] Thits23;
				Thits23 = NULL;
				delete[] path2;
				path2 = NULL;
			}
			//path with 3 reflection
			for (int npath3i = 0; npath3i < npath3; npath3i++)
			{
				currentpath++;
				id = npath3i * 15;
				Phitijk[0] = rpath3[0 + id];
				Phitijk[1] = rpath3[1 + id];
				Phitijk[2] = rpath3[2 + id];
				Phitjk[0] = rpath3[5 + id];
				Phitjk[1] = rpath3[6 + id];
				Phitjk[2] = rpath3[7 + id];
				Phitk[0] = rpath3[10 + id];
				Phitk[1] = rpath3[11 + id];
				Phitk[2] = rpath3[12 + id];

				Thits31 = fct_caltpts(TX, Phitijk, planes, nplanes, &tflag31, Phit_nt);
				Thits32 = fct_caltpts(Phitijk, Phitjk, planes, nplanes, &tflag32, Phit_nt);
				Thits33 = fct_caltpts(Phitjk, Phitk, planes, nplanes, &tflag33, Phit_nt);
				Thits34 = fct_caltpts(Phitk, RX, planes, nplanes, &tflag34, Phit_nt);
				nb_pts = (tflag31 + tflag32 + tflag33 + tflag34 + 5);

				path3 = new float[nb_pts * 5];
				path3[0] = TX[0];
				path3[1] = TX[1];
				path3[2] = TX[2];
				path3[3] = 0.0f;
				path3[4] = 0.0f;
				indice = (tflag31 + 1) * 5;
				for (int indice_j = 5; indice_j < indice; indice_j++)
				{
					path3[indice_j] = Thits32[indice_j - 5];
				}
				path3[0 + indice] = Phitijk[0];
				path3[1 + indice] = Phitijk[1];
				path3[2 + indice] = Phitijk[2];
				path3[3 + indice] = rpath3[3 + id];
				path3[4 + indice] = rpath3[4 + id];
				ind = (tflag31 + tflag32 + 3) * 5;
				indice = (tflag31 + tflag32 + tflag33 + 3) * 5;
				for (int indice_j = ind; indice_j < indice; indice_j++)
				{
					path3[indice_j] = Thits33[indice_j - ind];
				}
				path3[0 + indice] = Phitjk[0];
				path3[1 + indice] = Phitjk[1];
				path3[2 + indice] = Phitjk[2];
				path3[3 + indice] = rpath3[8 + id];
				path3[4 + indice] = rpath3[9 + id];
				ind = (tflag31 + tflag32 + 3) * 5;
				indice = (tflag31 + tflag32 + tflag33 + 3) * 5;
				for (int indice_j = ind; indice_j < indice; indice_j++)
				{
					path3[indice_j] = Thits33[indice_j - ind];
				}
				path3[0 + indice] = Phitk[0];
				path3[1 + indice] = Phitk[1];
				path3[2 + indice] = Phitk[2];
				path3[3 + indice] = rpath3[13 + id];
				path3[4 + indice] = rpath3[14 + id];

				ind = (tflag31 + tflag32 + tflag33 + 4) * 5;
				indice = (tflag31 + tflag32 + tflag33 + tflag34 + 4) * 5;
				for (int indice_j = ind; indice_j < indice; indice_j++)
				{
					path3[indice_j] = Thits34[indice_j - ind];
				}

				path3[0 + indice] = RX[0];
				path3[1 + indice] = RX[1];
				path3[2 + indice] = RX[2];
				path3[3 + indice] = 0.0;
				path3[4 + indice] = 0.0;
				compute_distloss(path3, nb_pts, planes, currentpath, distance_mfct, lossfac_mfct);
				delete[] Thits31;
				Thits31 = NULL;
				delete[] Thits32;
				Thits32 = NULL;
				delete[] Thits33;
				Thits33 = NULL;
				delete[] Thits34;
				Thits34 = NULL;
				delete[] path3;
				path3 = NULL;
			}
			//compute final value
			temp = 0.0f;
			for (int temp_i = 0; temp_i < npath; temp_i++)
			{
				temp += (lossfac_mfct[temp_i] * lossfac_mfct[temp_i]);
			}
			rs_amp[n + ld] = sqrtf(temp);
			delete[] rpath1;
			rpath1 = NULL;
			delete[] rpath2;
			rpath2 = NULL;
			delete[] rpath3;
			rpath3 = NULL;
			if (nr == 0)
			{
				delete[] path1;
				delete[] path2;
				delete[] path3;
				path1 = path2 = path3 = NULL;
			}
			if (nr == 1)
			{
				delete[] path2;
				delete[] path3;
				path2 = path3 = NULL;
			}
			if (nr == 2)
			{
				delete[] path3;
				path3 = NULL;
			}
			delete[] distance_mfct;
			distance_mfct = NULL;
			delete[] lossfac_mfct;
			lossfac_mfct = NULL;
		}
	}
	delete[] TX_i;
	delete[] TX_j;
	delete[] TX_k;
	delete[] Phiti;
	delete[] Phitj;
	delete[] Phitij;
	delete[] Phitk;
	delete[] Phitjk;
	delete[] Phitijk;
	delete[] Phit_nt;
	delete[] TX;
	delete[] RX;


	
}
//get a reflection symmetry point
void fct_calimage(float *TX0, float **planes, int i, float *TX_i)
{
	register float t_imi, ux = planes[12][i], uy = planes[13][i], uz = planes[14][i];
	register float TX0_x = TX0[0], TX0_y = TX0[1], TX0_z = TX0[2];
	t_imi = -2.0f*(TX0_x*ux + TX0_y*uy + TX0_z*uz + planes[15][i]);
	TX_i[0] = TX0_x + ux*t_imi;
	TX_i[1] = TX0_y + uy*t_imi;
	TX_i[2] = TX0_z + uz*t_imi;//slove the symmetry point
}
bool fct_calrpt(float *P1, float *P2, float **planes, int i, float *Phit)
{
	float ux = planes[12][i], uy = planes[13][i], uz = planes[14][i];
	register float P1_x = P1[0], P1_y = P1[1], P1_z = P1[2];
	float P1P2_x, P1P2_y, P1P2_z;
	//the coordinates of the receiving point and the virtual base station coordinate different
	P1P2_x = P2[0] - P1_x;
	P1P2_y = P2[1] - P1_y;
	P1P2_z = P2[2] - P1_z;
	//calculate the distance
	float ti_max = sqrtf(P1P2_x*P1P2_x + P1P2_y*P1P2_y + P1P2_z*P1P2_z);//虚像点和接收点之间的距离
	float ti_maxtol = ti_max - t_tol;
	
	float cte,Rdx,Rdy,Rdz;
	float Vdp;
	float ti;
	if (ti_max > t_tol)
	{
		cte = 1.0f / ti_max;
		Rdx = P1P2_x*cte;
		Rdy = P1P2_y*cte;
		Rdz = P1P2_z*cte;//the cos value of the vector
		//dot product to get the angle
		//both vectors are normailizied, so the parament is the cosine of the angle//half correct half wrong
		Vdp = ux*Rdx + uy*Rdy + uz*Rdz;
		if (Vdp != 0.0)//not a vertical reflection
		{
			ti = -(ux*P1_x + uy*P1_y + uz*P1_z + planes[15][i]) / Vdp;
			if ((t_tol < ti) && (ti < ti_maxtol))
			{//coordinate value of the reflection intersection point
				Phit[0] = P1_x + Rdx*ti;
				Phit[1] = P1_y + Rdy*ti;
				Phit[2] = P1_z + Rdz*ti;
				if ((planes[16][i] <= Phit[0]) && (Phit[0] <= planes[17][i]) && (planes[18][i] <= Phit[1]) && (Phit[1] <= planes[19][i]) && (planes[20][i] <= Phit[2]) && (Phit[2] <= planes[21][i]))
					return true;
			}
		}
	}
	return false;
}
float * fct_caltpts(float *P1, float *P2, float **planes, int nplanes, int *nthits, float *Phit_nt)
{
	float P1_x = P1[0], P1_y = P1[1], P1_z = P1[2];
	float P1P2_x, P1P2_y, P1P2_z,tmp,Rdx_nt,Rdy_nt,Rdz_nt,t_nt;
	float ux_nt, uy_nt, uz_nt, Vd_nt;
	int nt, number = 0, jd, id, i, ind;//ntd;
	float *result_caltpts=NULL;
	float *tab_pts = NULL, *sorted_tab_pts = NULL;
	int *index = NULL;
	bool alread_there = false;
	//the coordinates of the receiving point and the virtual base station coordinate different
	P1P2_x = P2[0] - P1_x;
	P1P2_y = P2[1] - P1_y;
	P1P2_z = P2[2] - P1_z;
	float tmax_nt = sqrtf(P1P2_x*P1P2_x + P1P2_y*P1P2_y + P1P2_z*P1P2_z);
	float tmax_ntol = tmax_nt - t_tol;
	if (tmax_nt > t_tol)
	{
		tmp = 1.0f / tmax_nt;
		Rdx_nt = P1P2_x*tmp;
		Rdy_nt = P1P2_y*tmp;
		Rdz_nt = P1P2_z*tmp;
		for (nt = 0; nt < nplanes; nt++)
		{
			//ntd =
			ux_nt = planes[12][nt];
			uy_nt = planes[13][nt];
			uz_nt = planes[14][nt];
			Vd_nt = ux_nt*Rdx_nt + uy_nt*Rdy_nt + uz_nt*Rdz_nt;
			if (Vd_nt != 0.0f)
			{
				////发射点到墙面第一个点的向量和墙面法向向量的内积除以上述余弦值 
				t_nt = -(ux_nt*P1_x + uy_nt*P1_y + uz_nt*P1_z + planes[15][nt]) / Vd_nt;
				if ((t_tol < t_nt) && (t_nt < tmax_ntol))
				{
					Phit_nt[0] = P1_x + Rdx_nt*t_nt;
					Phit_nt[1] = P1_y + Rdy_nt*t_nt;
					Phit_nt[2] = P1_z + Rdz_nt*t_nt;
					//if ((planes[16][nt]<=Phit_nt[0])&&(Phit_nt[0]<=planes))
					if ((planes[16][nt] <= Phit_nt[0]) && (Phit_nt[0] <= planes[17][nt]) && (planes[18][nt] <= Phit_nt[1]) && (Phit_nt[1] <= planes[19][nt]) && (planes[20][nt] <= Phit_nt[2]) && (Phit_nt[2] <= planes[21][nt]))
					{
						if (number==0)
						{
							jd = number * 6;
							number++;
							tab_pts = (float *)realloc((float *)tab_pts, (jd+ 6) * sizeof(float));
							tab_pts[0 + jd] = Phit_nt[0];
							tab_pts[1 + jd] = Phit_nt[1];
							tab_pts[2 + jd] = Phit_nt[2];
							tab_pts[3 + jd] = (float)nt;
							tab_pts[4 + jd] = 2.0;
							tab_pts[5 + jd] = t_nt;
						}
						else
						{
							alread_there = false;
							for (i = 0; i < number; i++)
							{
								id = 6 * i;
								if ((fabsf(tab_pts[0 + id] - Phit_nt[0]) < t_tol) && (fabsf(tab_pts[1 + id] - Phit_nt[1]) < t_tol) && (fabsf(tab_pts[2 + id] - Phit_nt[2]) < t_tol))
								{
									alread_there = true;
								}
							}
							if (alread_there == false)
							{
								jd = number * 6;
								number++;
								tab_pts = (float *)realloc((float *)tab_pts, (jd+6) * sizeof(float));
								tab_pts[0 + jd] = Phit_nt[0];
								tab_pts[1 + jd] = Phit_nt[1];
								tab_pts[2 + jd] = Phit_nt[2];
								tab_pts[3 + jd] = (float)nt;
								tab_pts[4 + jd] = 2.0;
								tab_pts[5 + jd] = t_nt;
							}
						}
					}
				}
			}
		}
		if (number > 0)
		{
			sorted_tab_pts = new float[number];
			index = new int[number];
			for (i = 0; i < number;i++)
			{
				index[i] = i;
				sorted_tab_pts[i] = tab_pts[5 + i * 6];
			}
			qsindex(sorted_tab_pts, index, 0, number - 1);
			result_caltpts = new float[5 * number];
			for (i = 0; i < number; i++)
			{
				ind = 6 * index[i];
				id = 5 * i;
				result_caltpts[0 + id] = tab_pts[0 + ind];
				result_caltpts[1 + id] = tab_pts[1 + ind];
				result_caltpts[2 + id] = tab_pts[2 + ind];
				result_caltpts[3 + id] = tab_pts[3 + ind];
				result_caltpts[4 + id] = tab_pts[4 + ind];
			}
		}
	}
	delete[]tab_pts;
	delete[]sorted_tab_pts;
	delete[]index;
	nthits[0] = number;

	return result_caltpts;
}
void compute_distloss(float *path, int num_pts, float **planes, int currentpath, float *distance_distloss, float *lossfac)
{
	//change parament
	float flp_scale = 100.0f;
	float fc = 1000000000;//Hz
	float material[] = { 0.0f, 0.0f, 0.0f, 0.0f, 5.0f, 0.04f };

	float lossreal = 1.0f, lossimag = 0.0f, tmp1, tmp2, tmp;
	float Vd, vs, vs2;
	float x1, x2, y1, y2, x3, y3, x4, y4, x5, y5, r1, r2;
	float P1P2_x, P1P2_y, P1P2_z, uy, dj, loss_jreal, loss_jimag, epsreal, epsimag;
	int j, jd, flagP2, ind, ind2;
	float lambda, Ao, cte = 1.0f / flp_scale;

	lambda = c / fc;
	Ao = lambda / (4.0f*PI);
	distance_distloss[currentpath] = 0.0f;
	lossfac[currentpath] = 1.0f;
	for (j = 0; j < num_pts; j++)
	{
		jd = j * 5;
		P1P2_x = path[5 + jd] - path[0 + jd];
		P1P2_y = path[6 + jd] - path[1 + jd];
		P1P2_z = path[7 + jd] - path[2 + jd];
		
		jd += 5;
		flagP2 = (int)path[4 + jd];
		ind = (int)path[3 + jd];
		//ind1=ind
		uy = planes[13][ind];
		dj = sqrtf(P1P2_x*P1P2_x + P1P2_y*P1P2_y + P1P2_z*P1P2_z);
		loss_jreal = 1.0f;
		loss_jimag = 0.0f;
		if (flagP2 == 1)
		{
			epsreal = material[4]; 
			epsimag = -60.0f*lambda*material[5];
			tmp = 1.0f / dj;
			Vd = planes[12][ind] * P1P2_x*tmp + uy*P1P2_y*tmp + planes[14][ind] * P1P2_z*tmp;
			vs = fabsf(Vd);
			vs2 = vs*vs;
			x1 = epsreal - (1.0f - vs2);
			y1 = epsimag;
			r1 = sqrt(x1*x1 + y1*y1);
			x2 = sqrtf((r1 + x1)*0.5f);
			y2 = y1 / (2.0f*x2);
			if (fabsf(uy) > 0.5)
			{
				x3 = epsreal*vs - x2;
				y3 = epsimag*vs - y2;
				x4 = epsreal*vs + x2;
				y4 = epsimag*vs + y2;
				r2 = 1.0f / (x4*x4 + y4*y4);
				loss_jreal = (x3*x4 + y3*y4)*r2;
				loss_jimag = (y3*x4 - x3*y4)*r2;
			}
			else
			{
				x3 = vs - x2;
				y3 = -y2;
				x4 = vs + x2;
				y4 = y2;
				r2 = 1.0f / (x4*x4 + y4*y4);
				loss_jreal = (x3*x4 + y3*y4)*r2;
				loss_jimag = (y3*x4 - x3*y4)*r2;
			}
		}
		if ((flagP2 == 2))
		{
			ind2 = ind * 6;
			epsreal = material[4];
			epsimag = -60.0f*lambda*material[5];
			tmp = 1.0f / dj;
			Vd = planes[12][ind] * P1P2_x*tmp + uy*P1P2_y*tmp + planes[14][ind] * P1P2_z*tmp;
			vs = fabsf(Vd);
			vs2 = vs*vs;
			x1 = epsreal - (1.0f - vs2);
			y1 = epsimag;
			r1 = sqrtf(x1*x1 + y1*y1);
			x2 = sqrtf((r1 + x1)*0.5f);
			y2 = y1 / (2.0f*x2);
			if (fabsf(uy) > 0.5f)
			{
				x3 = epsreal*vs - x2;
				y3 = epsimag*vs - y2;
				x4 = epsreal*vs + x2;
				y4 = epsimag*vs + y2;
				r2 = 1.0f / (x4*x4 + y4*y4);
				x5 = (x3*x4 + y3*y4)*r2;
				y5 = (y3*x4 - x3*y4)*r2;
				loss_jreal = sqrtf(loss_factor*(1.0f - (x5*x5 + y5*y5)));
			}
			else
			{
				x3 = vs - x2;
				y3 = -y2;
				x4 = vs + x2;
				y4 = y2;
				r2 = 1.0f / (x4*x4 + y4*y4);
				x5 = (x3*x4 + y3*y4)*r2;
				y5 = (y3*x4 - x3*y4)*r2;
				loss_jreal = sqrtf(loss_factor*(1.0f - (x5*x5 + y5*y5)));
			}
		}
		distance_distloss[currentpath] += dj;
		tmp1 = (lossreal*loss_jreal - lossimag*loss_jimag);
		tmp2 = (lossreal*loss_jimag + lossimag*loss_jreal);
		lossreal = tmp1;
		lossimag = tmp2;
	}
	lossfac[currentpath] = (Ao*sqrtf(lossreal*lossreal + lossimag*lossimag)) / (distance_distloss[currentpath] * cte);
}
void qsindex(float *a, int *index, int lo, int hi)//sort from small to large
{
	//  lo is the lower index, hi is the upper index of the region of array a that is to be sorted
	int i = lo, j = hi, ind;
	float x = a[(lo + hi) / 2], h;
	do 
	{
		while (a[i] < x) i++;
		while (a[j] > x) j--;
		if (i <= j)
		{
			h = a[i];
			a[i] = a[j];
			a[j] = h;
			ind = index[i];
			index[i] = index[j];
			index[j] = ind;
			i++;
			j--;
		}
	} while (i<=j);
	/*  recursion */
	if (lo < j) qsindex(a, index, lo, j);
	if (i < hi) qsindex(a, index, i, hi);
}