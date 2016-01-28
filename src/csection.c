#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>


#include "splinesurf.h"
#include "fourier.h"

/*
  Program: csection

  Author: Vinicius Vaz da Cruz - viniciusvcruz@gmail.com 
	
  History: Based on Raman written by Vinicius

  Purpose: Compute the ABSORPTION and RIXS cross section in strong laser fields

  Stockholm, 19th of January of 2016
*/



int csection(double ti, double stept, int n,double *re_wr12, double *im_wr12,double *re_wr23, double *im_wr23, double *G12, double *G23, int n_fourier,double *detun){
  int i,j;

  //--- fft variables
  int nE;
  double *T,*E,Ei,aE,stepE;
  fftw_complex *rho12_t,*rho23_t,*rho12_v,*rho23_v,*G12_t, *G23_t,*G12_v, *G23_v;

  //--- spline variables
  int ntg,kx;
  double stept_spl;
  double *bcoefre12,*bcoefim12,*bcoefre23,*bcoefim23,*tknot;
  double *bcoefG12,*bcoefG23;
  
  //--default values--------
  kx = 6.0e+0;

  

  //------ do all splines ---------------------------

  T = malloc(n*sizeof(double));
  tknot = malloc((n+kx)*sizeof(double));

  bcoefre12 = malloc(n*sizeof(double));
  bcoefim12 = malloc(n*sizeof(double));
  bcoefre23 = malloc(n*sizeof(double));
  bcoefim23 = malloc(n*sizeof(double));
  bcoefG12  = malloc(n*sizeof(double));
  bcoefG23  = malloc(n*sizeof(double));
  

  for(i=0;i<n;i++) T[i] = ti + i*stept;

  dbsnak_ (&n, T, &kx, tknot);

  // rho_12
  dbsint_ (&n,T,re_wr12,&kx,tknot,bcoefre12);
  dbsint_ (&n,T,im_wr12,&kx,tknot,bcoefim12);

  // rho_23
  dbsint_ (&n,T,re_wr12,&kx,tknot,bcoefre12);
  dbsint_ (&n,T,im_wr12,&kx,tknot,bcoefim12);

  //G12
  dbsint_ (&n,T,G12,&kx,tknot,bcoefG12);

  //G23
  dbsint_ (&n,T,G23,&kx,tknot,bcoefG23);

  //--------------------------------------------------


  //--- fourier transforms ---------------------------
  ntg = pow(2,twopow);
  steptspl=(2.0*T[nf-1])/ntg;
  stepE = 2*M_PI/(ntg*steptspl);
  E = malloc(ntg*sizeof(double));


  
  //--- compute cross-sections -----------------------

  //--------------------------------------------------

  return 0;
}


//create routine that generates the fftw mirrored vectos using splines
int gen_data(int NTG,double tf,double *bcoefre, double *bcoefim, fftw_complex *Y){
  int i,s;
  double t,steptspl;

  steptspl = (2.0e+0 * tf)/NTG;

  //generate mirrored data
  for(i=0;i<NTG;i++){
    t = -T[nf-1] + i*steptspl;
    fprintf(deb,"%lf ",t);
    if(t < 0) s = -1.0e+0;
    else s = 1.0e+0;

    t = fabs(t);

    if(strncasecmp(windtype,".SGAUSS",7)==0){
      //TAUX = TMAX**2/(LOG(ONE/WP))**2 espec
      //taux = pow(T[nf-1],2)/(log(1.000/width)/log(M_E));
      taux = pow(T[nf-1],2)/pow(log(1.000/width),2);
      window = exp(-pow(t,2)/taux);
    }else if(strncasecmp(windtype,".EXPDEC",7)==0){
      window = exp(-width*t);
    }else if(strncasecmp(windtype,".NONE",7)==0){
      window = 1.0e0;
    }

    //debug
    //window=1.000;
    workin[i][0] =   window*dbsval_ (&t,&kx,tknot,&nf,bcoefre);
    workin[i][1] = s*window*dbsval_ (&t,&kx,tknot,&nf,bcoefim);
    fprintf(deb,"%lf %lf \n",workin[i][0], workin[i][1]);
  }


}

  //dbsval_ (&omega,&kx,xas_knot,&nxas,xas_bcoef);

}

//----------------------------

int do_fft(int NTG,fftw_complex *workin, fftw_complex *workout){
  // subroutine to perform the fourier transform using fftw
  fftw_plan p;

  //centers t=0 at the zero frequency position of the fft input vector(first) due to periodicity requirement of fft
  center_fft(workin,NTG);

  // FFTW_BACKWARD -> sign in the exponent = 1.
  p = fftw_plan_dft_1d(NTG,workin,workout,FFTW_BACKWARD,FFTW_ESTIMATE);

  fftw_execute(p);  

  //shifts the fft result back to the center of the array
  center_fft(workout,NTG);

  fftw_destroy_plan(p);

  //--------------------------------------------------
  return 0;
}
