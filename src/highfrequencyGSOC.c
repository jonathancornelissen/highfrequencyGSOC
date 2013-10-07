#include <stdio.h> 
#include <stdlib.h> 
#include <R.h> 
#include "highfrequencyGSOC.h" 
#include <math.h> 

void heavy_likelihoodR(double *parameters, double *data, int *T1, int *K1, double *means, int *p, int *q, int *pMax1, int *qMax1, 
      double *backcast, double *LB, double *UB, int *compconst1, double *h, double *lls, double *llRM, double *ll) 
{
   
   int i,j,t,l,k;
   double sum=0.0;
   double sum1,sum2;
   double lll = 0.0;
   double htemp;
  
   int K,T,pMax,qMax,compconst;
   K = K1[0];
   T = T1[0];
   
   pMax = pMax1[0]; qMax = qMax1[0];
   compconst = compconst1[0];
   
   double likConst = K*log(2*M_PI);
  
   double O[K];
   double A[K*K*pMax];
   double B[K*K*qMax];
   double temp[K]; 

   //---- initialize -----------------
   for(i=0;i<K;i++)
   {
     O[i] = 0.0;
     for(k=0;k<K;k++)
     {
       for(j=0;j<pMax;j++)
       {A[K*K*j + K*i+k] = 0.0;}
       for(j=0;j<qMax;j++)
       {B[K*K*j + K*i+k] = 0.0;}              
     } 
     llRM[i] = 0.0;
   }   
   
   if(compconst == 1)
   {heavy_parameter_transformR(parameters, K, p, q, &O[0], &A[0], &B[0]);}
   else
   {heavy_parameter_transform_RetrackR(parameters, K, p, q, means, &O[0], &A[0], &B[0]);}
   
   
   for(i=0;i<K;i++) if(O[i] < 0) {O[i] = .000001;} //?? O(O<0) = realmin();
     
   for(t = 0; t < T; t++) //for t=1:T
   {
     
     //for(i=0;i<K;i++) {SETM(model->h,t,i,O[i]);} //h(:,t) = O
     
     for(i=0;i<K;i++) {h[t*K + i] = O[i];}
     
     //------------------------------------------------

     for(j=0;j<pMax;j++) // j=1:pMax
     {
        if(t-j > 0)
	{        
	  //h(:,t) = h(:,t) + A(:,:,j)*data(:,t-j);
	  for(l=0;l<K;l++) 
	  {
	    sum = 0.0;
	    for(k=0;k<K;k++)
	    {sum = sum + A[K*K*j + K*l + k]*data[(t-j-1)*K+k];} //{sum = sum + A[K*K*j + K*l + k]*GETM(model->data,t-j-1,k);}

	    htemp = h[t*K+l]; // htemp = GETM(model->h,t,l); 
	    h[t*K+l] = htemp+sum; //SETM(model->h,t,l,htemp + sum);
	  }
	}    
	else
	{
          //h(:,t) = h(:,t) + A(:,:,j)*backCast';
	  for(l=0;l<K;l++) 
	  {
	    sum = 0.0;
	    for(k=0;k<K;k++)
	    {sum = sum + A[K*K*j + K*l + k]*backcast[k];} //{sum = sum + A[K*K*j + K*l + k]*GET(model->backcast,k);}
	    //htemp = h[t*K+l]; //htemp = GETM(model->h,t,l); 
	    h[t*K + l] = h[t*K+l]+sum; //SETM(model->h,t,l,htemp + sum);
	  } 
	}

    }
    
    
    for(j=0;j<qMax;j++) // j=1:pMax
    {
        if(t-j > 0)
	{        
	  //h(:,t) = h(:,t) + B(:,:,j)*h(:,t-j);
	  for(l=0;l<K;l++) 
	  {
	    sum = 0.0;
	    for(k=0;k<K;k++)
	    {sum = sum + B[K*K*j + K*l + k]*h[(t-j-1)*K+k];} //{sum = sum + B[K*K*j + K*l + k]*GETM(model->h,t-j-1,k);}
	    //htemp = GETM(h,l,t); 
	    temp[l] = sum;
	    //printf("%lf\n", temp[l]);
	  }
	  for(l=0;l<K;l++) 
	  {h[t*K + l] = h[t*K+l]+sum;}	  
	}    
	else
	{
          //h(:,t) = h(:,t) + B(:,:,j)*backCast';
	  for(l=0;l<K;l++) 
	  {
	    sum = 0.0;
	    for(k=0;k<K;k++)
	    {sum = sum + B[K*K*j + K*l + k]*backcast[k];}

	    h[t*K + l] = h[t*K+l]+sum;
	  } 
	}
    }    
    
   // printf("%lf\n", GETM(model->h,t,1));
    
    for(j=0;j<K;j++) //j=1:K
    {
       
       htemp = h[t*K + j];
       
       if(htemp < LB[j])       //h(j,t)<lb(j)
       { 
	 //printf("%lf\n", GETM(model->h,t,j)); 
	 //h(j,t) = lb(j) * 1./(1-(h(j,t)-lb(j)));
         //htemp = GETM(h,j,t);
	 htemp = LB[j]/(1.0 - (htemp - LB[j]));
	 h[t*K+j] = htemp;
       }
       else if(htemp > UB[j])
       {
        h[t*K+j] = UB[j] + log(htemp - UB[j]);
       }
    }
    
    //lls(t)  = 0.5*(likConst + sum(log(h(:,t))) + sum(data(:,t)./h(:,t)));
    sum1 = 0.0; sum2 = 0.0;
    for(l=0;l<K;l++)
    {
     sum1 = sum1 + log(h[t*K+l]);
     sum2 = sum2 + data[t*K+l]/h[t*K+l]; 
     
     llRM[l] = llRM[l] - .5*(log(h[t*K+l]) + data[t*K+l]/h[t*K+l]);     
    }
    
    
    
    //SET(model->lls,t, 0.5*(likConst + sum1 + sum2));
    lls[t] = 0.5*(likConst + sum1 + sum2);
    //printf("%lf\n", GET(model->lls,t));
    
    lll = lll + lls[t];
  }
  
  //plotData(model->lls->data, model->n_obs);
  
  //ll = gsl_blas_dasum(model->lls);

  *ll = lll;
  
  //printf("likelihood = %lf\n", ll);
  
  if(isnan(*ll) || isinf(*ll)) 
  {*ll = 1000000000000.0;}

}


void heavy_parameter_transformR(double *parameters, int K, int *p, int *q, double *O, double *A, double *B)
{

  int i,j,k;
  
  for(k=0;k<K;k++) {O[k] = parameters[k];} //O(:) = parameters(1:K);

  int pc,qc;
  int offset = K;

  for(i=0;i<K;i++)
  {
    for(j=0;j<K;j++)
    {
      pc = p[i*K+j];
      for(k=0;k<pc;k++)
      {A[K*K*k + K*i + j] = parameters[offset + k];}
      offset = offset + pc;
    }
  }
  
  for(i=0;i<K;i++)
  {
    for(j=0;j<K;j++)
    {
      qc = q[i*K+j];
      for(k=0;k<qc;k++)
      {B[K*K*k + K*i + j] = parameters[offset + k];}
      offset = offset + qc;
    }
  } 
  
}


void heavy_parameter_transform_RetrackR(double *parameters, int K,  int *p,  int *q,  double* means, double *O, double *A, double *B)
{
  
  int i,j,k;
  double kappa;

  int pc,qc;
  int offset = 0;

  for(i=0;i<K;i++)
  {
    for(j=0;j<K;j++)
    {
      pc = p[i*K+j];
      for(k=0;k<pc;k++)
      {A[K*K*k + K*i + j] = parameters[offset + k];}
      offset = offset + pc;
    }
  }
  
  for(i=0;i<K;i++)
  {
    for(j=0;j<K;j++)
    {
      qc = q[i*K+j];
      for(k=0;k<qc;k++)
      {B[K*K*k + K*i + j] = parameters[offset + k];}
      offset = offset + qc;
    }
  }
  
  if(K==0)
  {O[0] = means[0]*(1.0 - A[0] - B[0]);}  
  else if(K==1)
  {
    kappa = means[1]/means[0];  
    O[0] = means[0]*(1.0 - A[0] - kappa*A[1] - B[0]);
    O[1] = means[1]*(1.0 - A[3] - B[3]);
  }
  else if(K==3)
  {
     O[1] = means[1]*(1.0 - A[4] - A[5] - B[4]);
     O[2] = means[2]*(1.0 - A[8] - B[8]);
  }     

}