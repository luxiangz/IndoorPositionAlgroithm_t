#pragma once
void Mfct_Ref(int getdimensions, float **planes, float **TXpoint, float **RXpoint, int Nt, int Nr, int nplanes, int nr, float *rs_amp);/*material//fc//flp_scale */ // (getdimensions, datatxt.planes£¬ **basecoor_temp0, **meshgrid, basestationnum, nb_pts, row, nr);
void fct_calimage(float *TX0, float **planes, int i, float *TX_i);
bool fct_calrpt(float *P1, float *P2, float **planes, int i, float *Phit);
float * fct_caltpts(float *P1, float *P2, float **planes, int nplanes, int *nthits, float *Phit_nt);
void compute_distloss(float *path, int num_pts, float **planes, int currentpath, float *distance_distloss, float *lossfac);
void qsindex(float *, int *, int, int);
