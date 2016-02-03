#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>


#include "splinesurf.h"

/*
  Program: csection

  Author: Vinicius Vaz da Cruz - viniciusvcruz@gmail.com 
	
  History: Based on Raman written by Vinicius

  Purpose: Compute the ABSORPTION and RIXS cross section in strong laser fields

  Stockholm, 19th of January of 2016
*/

int center_fft(fftw_complex *out,int N);

int csection(double ti, double stept, int n,double *re_wr12, double *im_wr12,double *re_wr23, double *im_wr23, double *G12, double *G23, int n_fourier,double *detun){
  int i,j;

  //--- fft variables
  int nE;
  double *T,*E,Ei,aE,stepE,steptspl;
  fftw_complex *rho12_t,*rho23_t,*rho12_v,*rho23_v,*G12_t, *G23_t,*G12_v, *G23_v;

  //--- spline variables
  int ntg,kx;
  double stept_spl;
  double *bcoefre12,*bcoefim12,*bcoefre23,*bcoefim23,*tknot;
  double *bcoefG12,*bcoefG23;
  
  //--default values--------
  kx = 6.0e+0;
  //------------------------

  ntg = pow(2,n_fourier); // number of points to be used in the fourier transform


  //---- allocate vectors

  T = malloc(n*sizeof(double));
  tknot = malloc((n+kx)*sizeof(double));

  bcoefre12 = malloc(n*sizeof(double));
  bcoefim12 = malloc(n*sizeof(double));

  bcoefre23 = malloc(n*sizeof(double));
  bcoefim23 = malloc(n*sizeof(double));

  bcoefG12  = malloc(n*sizeof(double));

  bcoefG23  = malloc(n*sizeof(double));

  rho12_t = fftw_malloc(ntg * sizeof(fftw_complex));
  rho12_v = fftw_malloc(ntg * sizeof(fftw_complex));

  rho23_t = fftw_malloc(ntg * sizeof(fftw_complex));
  rho23_v = fftw_malloc(ntg * sizeof(fftw_complex));

  G12_t = fftw_malloc(ntg * sizeof(fftw_complex));
  G12_v = fftw_malloc(ntg * sizeof(fftw_complex));

  G23_t = fftw_malloc(ntg * sizeof(fftw_complex));
  G23_v = fftw_malloc(ntg * sizeof(fftw_complex));
  

  //------ do all splines ---------------------------
  

  for(i=0;i<n;i++) T[i] = ti + i*stept;

  dbsnak_ (&n, T, &kx, tknot);

  // rho_12
  dbsint_ (&n,T,re_wr12,&kx,tknot,bcoefre12);
  dbsint_ (&n,T,im_wr12,&kx,tknot,bcoefim12);

  // rho_23
  dbsint_ (&n,T,re_wr12,&kx,tknot,bcoefre23);
  dbsint_ (&n,T,im_wr12,&kx,tknot,bcoefim23);

  //G12
  dbsint_ (&n,T,G12,&kx,tknot,bcoefG12);

  //G23
  dbsint_ (&n,T,G23,&kx,tknot,bcoefG23);

  //--------------------------------------------------

  //--- generate data in ntg = 2 ^ M points using splines

  gen_data(n,ntg,kx,ti,tknot,bcoefre12,bcoefim12,rho12_t); // rho_12
  gen_data(n,ntg,kx,ti,tknot,bcoefre23,bcoefim23,rho23_t); // rho_23
  gen_data_real(n,ntg,kx,ti,tknot,bcoefG12,G12_t);         // G12
  gen_data_real(n,ntg,kx,ti,tknot,bcoefG23,G23_t);         // G23

  // free spline coefficient matrices
  free(bcoefre12);free(bcoefim12);
  free(bcoefre23);free(bcoefim23);
  free(bcoefG12);free(bcoefG23);free(tknot);

  //--- do fourier transforms ---------------------------
  steptspl=(2.0*T[n-1])/ntg;
  stepE = 2*M_PI/(ntg*steptspl);
  E = malloc(ntg*sizeof(double));

  do_fft(ntg,rho12_t,rho12_v);
  do_fft(ntg,rho23_t,rho23_v);
  do_fft(ntg,G12_t,G12_v);
  do_fft(ntg,G23_t,G23_v);
  

  //--- compute cross-sections -----------------------

  //--------------------------------------------------


  // free fourier vectors
  fftw_free(rho12_t);fftw_free(rho23_t);
  fftw_free(rho12_v);fftw_free(rho23_v);
  fftw_free(G12_t);fftw_free(G23_t);
  fftw_free(G12_v);fftw_free(G23_v);

  return 0;
}


//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//


//routine that generates the fftw using splines
// in this case we can't mirror the data (due to equations) but fftw periodicity can still be enforced since rho_ij(0) = 0
int gen_data(int n,int NTG,int kx,double ti,double *tknot, double *bcoefre, double *bcoefim, fftw_complex *workin){
  int i;
  double t,steptspl;

  //generate data  
  for(i=0;i<NTG;i++){
    t = ti + i*steptspl;
    workin[i][0] =  dbsval_ (&t,&kx,tknot,&n,bcoefre);
    workin[i][1] =  dbsval_ (&t,&kx,tknot,&n,bcoefim);
  }
  
  return 0;
}

//routine that generates the fftw using splines for real data
int gen_data_real(int n,int NTG,int kx,double ti,double *tknot,double *bcoefre, fftw_complex *workin){
  int i;
  double t,steptspl;

  //generate data  
  for(i=0;i<NTG;i++){
    t = ti + i*steptspl;
    workin[i][0] =  dbsval_ (&t,&kx,tknot,&n,bcoefre);
    workin[i][1] =  0.0e+0;
  }
  
  return 0;
}

//----------------------------

int do_fft(int NTG,fftw_complex *workin, fftw_complex *workout){
  // subroutine to perform the fourier transform using fftw
  fftw_plan p;

  //centers t=0 at the zero frequency position of the fft input vector(first) due to periodicity requirement of fft
  //center_fft(workin,NTG); not necessary because we're not inputting data from [-t,t] but rather [0,t]

  // FFTW_BACKWARD -> sign in the exponent = 1.
  p = fftw_plan_dft_1d(NTG,workin,workout,FFTW_BACKWARD,FFTW_ESTIMATE);

  fftw_execute(p);  

  //shifts the fft result back to the center of the array
  center_fft(workout,NTG);

  fftw_destroy_plan(p);

  //--------------------------------------------------
  return 0;
}


/*
 shifts zero frequency to the center of the array
*/

int center_fft(fftw_complex *out,int N){
  int i;
  double work;

  //centering the fourier transform
  for(i=0;i<N/2;i++){
    work= out[N/2 +i][0];
    out[N/2 +i][0] = out[i][0];
    out[i][0] = work;

    work= out[N/2 +i][1];
    out[N/2 +i][1] = out[i][1];
    out[i][1] = work;  
  }
}
