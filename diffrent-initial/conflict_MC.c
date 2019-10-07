#include <math.h>
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{
 
    int i, j, k, r1, r2;
    int n, b, t;
    double *pInp, *pW, *pT;
    double tri_bf, tri_af;
    
    n = mxGetN(prhs[0]);
    n = (int)sqrt(n);
    
    pInp = mxGetPr(prhs[0]);
    pT = mxGetPr(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(n*n, 1, mxREAL);
    
    pW = mxGetPr(plhs[0]);

    for(i=0; i<n*n; i++)
    {
        *(pW+i) = *(pInp+i);
    }

    for(j=0; j<50*1225; j++)
    {       

        r1 = (rand()%n);        
        r2 = (rand()%n);
        tri_bf = 0;
        tri_af = 0;
        if(r1!=r2)
        {
            for(i=0; i<n; i++)
            {            
                if((r1!=i)&&(r2!=i))
                {
                    tri_bf += *(pW+(r1*n)+i) * *(pW+(r2*n)+i) * *(pW+(r1*n)+r2);
                    tri_af -= *(pW+(r1*n)+i) * *(pW+(r2*n)+i) * *(pW+(r1*n)+r2);                
                }
            }
            if(tri_af>tri_bf||(1.0*rand()/RAND_MAX)<exp((tri_af-tri_bf)/ *pT))
            {
                *(pW+(r1*n)+r2) *= -1;
                *(pW+(r2*n)+r1) *= -1;
            }
        }
            
    }

    tri_af = 0;
    for(i=0; i<n; i++)
    {
        for(j=i+1; j<n; j++)
        {
            for(k=j+1; k<n; k++)
            {
                tri_af += *(pW+(i*n)+j) * *(pW+(j*n)+k) * *(pW+(k*n)+i)/n/(n-1)/(n-2)*6;
            }
        }
    }
    
    *(pW+0) = tri_af;
    
    return;
                   
}