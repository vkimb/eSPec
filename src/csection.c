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



int gen_data(int n,double stept_spl,int NTG,int kx,double ti,double *tknot, double *bcoefre, double *bcoefim, fftw_complex *workin);
int gen_data_real(int n,double stept_spl,int NTG,int kx,double ti,double *tknot,double *bcoefre, fftw_complex *workin);
int mult_detun(int ntg,double ti,double stept_spl,double detun,fftw_complex *rho_t, int s);
int do_fft(int NTG,fftw_complex *workin, fftw_complex *workout);
int center_fft(fftw_complex *out,int N);



int csection(double ti, double stept, int n,double *re_wr12, double *im_wr12,double *re_wr23, double *im_wr23, double *G12, double *G23, double *detun, int n_fourier){

  int i,j,verbose;

  double *sigma_xas, *sigma_rixs;
  
  //--- fft variables
  int nE,ierr;
  double *T,*E,Ei,aE,stepE,work,tol,wk[2];
  fftw_complex *rho12_t,*rho23_t,*rho12_v,*rho23_v,*G12_t, *G23_t,*G12_v, *G23_v;

  //--- spline variables
  int ntg,kx;
  double stept_spl;
  double *bcoefre12,*bcoefim12,*bcoefre23,*bcoefim23,*tknot;
  double *bcoefG12,*bcoefG23;

  //--- debug file
  FILE *deb=fopen("debug_csection.dat","w");
  FILE *debf=fopen("debug_fourier.dat","w");
  
  //--default values--------
  kx = 6.0e+0;
  tol = 1.0e-15;
  verbose = 2; //2 -> debug
  //------------------------

  if(verbose=2) printf("\n Starting cross-section routine \n");

  ntg = pow(2,n_fourier); // number of points to be used in the fourier transform
  
  if(verbose=2){
    printf(" n = %d \n",n);
    printf(" n_fourier = %d \n",n_fourier);
    printf(" number of points in the fourier transform: %d \n",ntg);
    printf(" initial time : %lf \n",ti);
    printf(" time step : %lf \n",stept);
    printf("detunings %E %E \n",detun[0],detun[1]);

    /* for(i=0;i<n;i++){ */
    /*   T[0] = ti + i*stept; */
    /*   fprintf(deb,"%lf %lf %lf %lf %lf \n",T[0],re_wr12[i],im_wr12[i],re_wr23[i],im_wr23[i]); */
    /* } */
  }

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

  if(verbose=2) printf("\n arrays allocated \n");

  //------ do all splines ---------------------------
  

  for(i=0;i<n;i++) T[i] = ti + i*stept;

  for(i=0;i<n;i++)fprintf(deb,"%lf %E %E \n",T[i],re_wr12[i],im_wr12[i]);
  fprintf(deb,"\n\n\n");

  dbsnak_ (&n, T, &kx, tknot);

  // rho_12
  dbsint_ (&n,T,re_wr12,&kx,tknot,bcoefre12);
  dbsint_ (&n,T,im_wr12,&kx,tknot,bcoefim12);

  // rho_23
  dbsint_ (&n,T,re_wr23,&kx,tknot,bcoefre23);
  dbsint_ (&n,T,im_wr23,&kx,tknot,bcoefim23);

  //G12
  dbsint_ (&n,T,G12,&kx,tknot,bcoefG12);

  //G23
  dbsint_ (&n,T,G23,&kx,tknot,bcoefG23);

  
   //check if spline is correct
  ierr=chk_spl(n,kx,T,tknot,re_wr12,bcoefre12,tol);if(ierr!=0)return 1;
  ierr=chk_spl(n,kx,T,tknot,im_wr12,bcoefim12,tol);if(ierr!=0)return 1;
  ierr=chk_spl(n,kx,T,tknot,re_wr23,bcoefre23,tol);if(ierr!=0)return 1;
  ierr=chk_spl(n,kx,T,tknot,im_wr23,bcoefim23,tol);if(ierr!=0)return 1;
  ierr=chk_spl(n,kx,T,tknot,G12,bcoefG12,tol);if(ierr!=0)return 1;
  ierr=chk_spl(n,kx,T,tknot,G23,bcoefG23,tol);if(ierr!=0)return 1;
  
  //--------------------------------------------------

  if(verbose=2) printf("\n spline coefficients generated! \n");

  //--- generate data in ntg = 2 ^ M points using splines
  stept_spl = (T[n-1] - T[0])/ntg;

  gen_data(n,stept_spl,ntg,kx,ti,tknot,bcoefre12,bcoefim12,rho12_t); // rho_12

  for(i=0;i<ntg;i++){
    T[0] = ti + i*stept_spl;
    fprintf(deb,"%lf %E %E \n",T[0],rho12_t[i][0],rho12_t[i][1] );
  }
  fprintf(deb,"\n\n\n");

  gen_data(n,stept_spl,ntg,kx,ti,tknot,bcoefre23,bcoefim23,rho23_t); // rho_23
  gen_data_real(n,stept_spl,ntg,kx,ti,tknot,bcoefG12,G12_t);         // G12
  gen_data_real(n,stept_spl,ntg,kx,ti,tknot,bcoefG23,G23_t);         // G23

  // multiply by exp(+s *i * omega * t) ==> R_ij(t) = rho_ij(t) e^{+s *i Omega t}
  mult_detun(ntg,ti,stept_spl,detun[0],rho12_t,-1);
  mult_detun(ntg,ti,stept_spl,detun[1],rho23_t,+1);
  

  // free spline coefficient matrices
  free(bcoefre12);free(bcoefim12);
  free(bcoefre23);free(bcoefim23);
  free(bcoefG12);free(bcoefG23);free(tknot);

  //--- do fourier transforms ---------------------------
  stepE = 2*M_PI/(ntg*stept_spl*41.3411);
  Ei = -ntg*stepE/(2.0E+0);
  E = malloc(ntg*sizeof(double));
  for(i=0;i<ntg;i++)E[i] = Ei + i*stepE;

  do_fft(ntg,rho12_t,rho12_v);
  do_fft(ntg,rho23_t,rho23_v);
  do_fft(ntg,G12_t,G12_v);
  do_fft(ntg,G23_t,G23_v);

  

  sigma_xas = malloc(ntg * sizeof(double));
  sigma_rixs = malloc(ntg * sizeof(double));

  //--- compute cross-sections -----------------------

  
  for(i=0;i<ntg;i++){
    T[0] = ti + i*stept_spl;

    //wk[0] = rho12_v[i][0] * G12_v[i][0] + rho12_v[i][1] * G12_v[i][1];
    //wk[1] = rho12_v[i][1] * G12_v[i][0] - rho12_v[i][0] * G12_v[i][1];
    sigma_xas[i] = rho12_v[i][1] * G12_v[i][0] - rho12_v[i][0] * G12_v[i][1];
    sigma_rixs[i] = rho23_v[i][1] * G23_v[i][0] - rho23_v[i][0] * G23_v[i][1];

    fprintf(deb,"%lf %E %E %E %E \n",T[0],rho12_t[i][0],rho12_t[i][1],G12_t[i][0],G12_t[i][1] );
    fprintf(debf,"%lf %E %E %E %E %E %E \n",E[i],rho12_v[i][0],rho12_v[i][1],G12_v[i][0],G12_v[i][1],sigma_xas[i],sigma_rixs[i]);
  }

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
int gen_data(int n,double stept_spl,int NTG,int kx,double ti,double *tknot, double *bcoefre, double *bcoefim, fftw_complex *workin){
  int i;
  double t;

  //generate data  
  for(i=0;i<NTG;i++){
    t = ti + i*stept_spl;
    workin[i][0] =  dbsval_ (&t,&kx,tknot,&n,bcoefre);
    workin[i][1] =  dbsval_ (&t,&kx,tknot,&n,bcoefim);
  }
  
  return 0;
}

//routine that generates the fftw using splines for real data
int gen_data_real(int n,double stept_spl,int NTG,int kx,double ti,double *tknot,double *bcoefre, fftw_complex *workin){
  int i;
  double t;

  //generate data  
  for(i=0;i<NTG;i++){
    t = ti + i*stept_spl;
    workin[i][0] =  dbsval_ (&t,&kx,tknot,&n,bcoefre);
    workin[i][1] =  0.0e+0;
  }
  
  return 0;
}

//----------------------------

int mult_detun(int ntg,double ti,double stept_spl,double detun,fftw_complex *rho_t, int s){
  int i;
  double t;

  for(i=0;i<ntg;i++){
    t = ti + i*stept_spl;
    t = t * 41.3411;
    rho_t[i][0] = rho_t[i][0] * cos(detun * t) - s * rho_t[i][1] * sin(detun * t);
    rho_t[i][1] = rho_t[i][1] * cos(detun * t) + s * rho_t[i][0] * sin(detun * t);
  }

  return 0;
}

//---------------------------

int do_fft(int NTG,fftw_complex *workin, fftw_complex *workout){
  // subroutine to perform the fourier transform using fftw
  fftw_plan p;

  //centers t=0 at the zero frequency position of the fft input vector(first) due to periodicity requirement of fft
  //center_fft(workin,NTG); not necessary because we're not inputting data from [-t,t] but rather [0,t]

  // FFTW_BACKWARD -> sign in the exponent = 1.
  p = fftw_plan_dft_1d(NTG,workin,workout,FFTW_FORWARD,FFTW_ESTIMATE);

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


/* subroutine to check if the spline routine produced correct values */
int chk_spl(int n, int kx, double *x, double *xknot, double *y, double *bcoef, double tol){
  int i,m;
  double chk_val,diff;

  m=n/4;

  for(i=0;i<n;i=i+m){
    
    chk_val = dbsval_ (&x[i],&kx,xknot,&n,bcoef);
    diff = fabs(chk_val - y[i]);
    
    if(diff > tol){
      printf("error!! spline routine error is above the %E tolerance! \n\n",tol);
      return 1;
    }
    
  }

  return 0;
}
