
/**************************************************************************
% Fuzzy Random Impulse Noise Reduction Method  (Second Detection Method)
%
%  The FRINR method was proposed in: 
%
%  Stefan Schulte, Valérie De Witte, Mike Nachtegael, 
%  Dietrich Van der Weken and  Etienne E. Kerre 
%  Fuzzy Sets and Systems 158(3), 2007 pp. 270-283  
%
% Stefan Schulte (stefan.schulte@Ugent.be):
% Last modified: 15/01/06
%**************************************************************************/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double LARGE (double x, double p1, double p2) {
   double res = 0.0;
   if ((x > p1) && (x < p2))   res = (x-p1)/(p2-p1);
   else if (x > p2)            res = 1.0;
   else                        res = 0.0;
   return res;
}

double absol(double a) {
   double b;
   if(a<0)   b=-a;
   else   b=a;
   return b;
}

/*************************************************************************
*  The main Function 
**************************************************************************/
void callMem(double **A1,double *mem1,int M, int N,int W) { 
   int i,j,k,l;   
   double **gem;
   int **xy, **adj;
   double som, tel, min;
   double a = 0, b = 0, noise, no_noise;
   double diff1, diff2, diff3, l1, l2, l3;
   
   int rand1a = 0;
   int rand1b = 0;
   int rand2a = 0;
   int rand2b = 0;
               
   gem = malloc(M*sizeof(double));
   for(i=0;i<M;i++)
      gem[i] = malloc(N*sizeof(double));

   xy = malloc(9*sizeof(int));
   adj = malloc(9*sizeof(int));
   for(i=0;i<9;i++){
      xy[i] = malloc(2*sizeof(int));
      adj[i] = malloc(2*sizeof(int));
   }
   
   xy[0][0] = -1; xy[0][1] = -1; xy[1][0] = -1; xy[1][1] = 0; 
   xy[2][0] = -1; xy[2][1] =  1; xy[3][0] =  0; xy[3][1] = -1; 
   xy[4][0] =  0; xy[4][1] =  0; xy[5][0] =  0; xy[5][1] = 1; 
   xy[6][0] =  1; xy[6][1] = -1; xy[7][0] =  1; xy[7][1] = 0; 
   xy[8][0] =  1; xy[8][1] =  1; 

   adj[0][0] = 1; adj[0][1] = -1; adj[1][0] =  0; adj[1][1] = -1; 
   adj[2][0] = 1; adj[2][1] =  1; adj[3][0] =  1; adj[3][1] =  0; 
   adj[4][0] = 0; adj[4][1] =  0; adj[5][0] =  1; adj[5][1] =  0; 
   adj[6][0] = 1; adj[6][1] =  1; adj[7][0] =  0; adj[7][1] = -1; 
   adj[8][0] = 1; adj[8][1] = -1; 

   for(i=0; i<M; i++){
      for(j=0; j<N; j++){
         /* step 1. Determine the local window*/      
         if(i < W) {
            rand1a = i;
            rand1b = W;
         }
         else {
              if (i>M-W-1){
               rand1a = W;
               rand1b = M-i-1;
            }
            else{
               rand1a = W;
               rand1b = W;
            }
         }

         
         if(j < W) {
            rand2a = j;
            rand2b = W;
         }
         else {
            if (j > N-W-1){
               rand2a = W;
               rand2b = N-j-1;
            }
            else{
               rand2a = W;
               rand2b = W;
            }
         }
         /* end step 1. */      
         
        
         
         /* step 2. Calculation of the membership degrees */      
         som = 0.0; tel = 0.0;
         for (k=i-rand1a; k<=i+rand1b; k++){
	        for (l=j-rand2a; l<=j+rand2b; l++){
	           if ((k==i)&&(l==j)) continue;
	           som +=  absol(A1[k][l] - A1[i][j]);
               tel += 1;
            }
         }
         gem[i][j] = som/tel;
      }
   }
         
       
   for(i=2; i<M-2; i++){
      for(j=2; j<N-2; j++){
         /* step 1. Determine the local window*/      
         if(i < W) {
            rand1a = i;
            rand1b = W;
         }
         else {
              if (i>M-W-1){
               rand1a = W;
               rand1b = M-i-1;
            }
            else{
               rand1a = W;
               rand1b = W;
            }
         }

         
         if(j < W) {
            rand2a = j;
            rand2b = W;
         }
         else {
            if (j > N-W-1){
               rand2a = W;
               rand2b = N-j-1;
            }
            else{
               rand2a = W;
               rand2b = W;
            }
         }
         /* end step 1. */      
         
         som = 0.0; tel = 0.0;
         min = 25555;
         for (k = i-rand1a; k<=i+rand1b; k++){
	        for (l = j-rand2a; l<=j+rand2b; l++){
	           if ((i==k)&&(j==l)) continue;
	           som +=  gem[k][l];
               tel += 1;
            }
         }
         a = som/tel;
         b = 1.1*a;
         
         no_noise = 0;
         noise = 0;
         for (k = 0; k<9; k++){
            if (k==4) continue;
            diff1 = A1[i+xy[k][1]][j+xy[k][0]] - A1[i][j];
            diff2 = A1[i+xy[k][1]+adj[k][1]][j+xy[k][0]+adj[k][0]] - A1[i+adj[k][1]][j+adj[k][0]];
            diff3 = A1[i+xy[k][1]-adj[k][1]][j+xy[k][0]-adj[k][0]] - A1[i-adj[k][1]][j-adj[k][0]];
            
            l1 = LARGE(absol(diff1),a,b);
            l2 = LARGE(absol(diff2),a,b);
            l3 = LARGE(absol(diff3),a,b);

            no_noise += l1*l2*l3 + (1-l1)*(1-l2)*(1-l3) - l1*l2*l3*(1-l1)*(1-l2)*(1-l3);
            noise    += (1-l1)*l2*l3 + l1*(1-l2)*(1-l3) + l1*l2*(1-l3) + l1*(1-l2)*l3 - 
                         (1-l1)*l2*l3 * l1*(1-l2)*(1-l3)*l1*l2*(1-l3)*l1*(1-l2)*l3;
         }
         if (noise == 0) mem1[i+j*M] = 8;
         else  mem1[i+j*M] = no_noise/noise;   // if mem >= 1 no_noise else noise      
      }   
   }         
   
   for(i=0;i<9;i++)  free(xy[i]);
   free(xy); 
   for(i=0;i<9;i++)  free(adj[i]);
   free(adj); 
   
   for(i=0;i<M;i++)  free(gem[i]);
   free(gem); 
   
   
}  /* End of callFuzzyShrink */


#define Im1      prhs[0]
#define WSIZE    prhs[1]

#define Mem1 plhs[0]



/**
*  The interaction with Matlab (mex):
*        nlhs = amount of output arguments (= 1)
*        nrhs = amount of input arguments (= 2)
*     *plhs[] = link to the output 
*     *prhs[] = link to the input 
*
**/
void mexFunction( int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[] ) {
    int row, col, i, M, N, K, W;
    double **mem1, **A1;
    
    if (nlhs!=1)
        mexErrMsgTxt("It requires three output arguments only [M1].");
    if (nrhs!=2)
       mexErrMsgTxt("this method requires 4 input argument [Im1 Wsize]");

    /* Get input values */  
    M = mxGetM(Im1);
    N = mxGetN(Im1);
    W = mxGetScalar(WSIZE);

    /**
    * Allocate memory for return matrices 
    **/
    Mem1 = mxCreateDoubleMatrix(M, N, mxREAL);  
    mem1 = mxGetPr(Mem1);

    /**
    * Dynamic allocation of memory for the input array
    **/
    A1 = malloc(M*sizeof(int));
    for(i=0;i<M;i++)
      A1[i] = malloc(N*sizeof(double));

     /**
     * Convert ARRAY_IN and INPUT_MASK to 2x2 C arrays (MATLAB stores a two-dimensional matrix 
     * in memory as a one-dimensional array) 
     ***/
     for (col=0; col < N; col++)
         for (row=0; row < M; row++) {
             A1[row][col] = (mxGetPr(Im1))[row+col*M];
	      }
	
    callMem(A1,mem1,M,N,W);

    for(i=0;i<M;i++)  free(A1[i]);
    free(A1); 
}
/* end mexFcn*/