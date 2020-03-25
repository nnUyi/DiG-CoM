#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "mex.h"


/*  A method for Anisotropic TV splitBregman derain */

void derain(double **D, double **B, double **R, double **E, double **dBx, double **dBy,
            double **dRtheta, double **dEx, double **dEy, double **bBx, double **bBy, double **bRtheta, double **bEx,
            double **bEy, double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double w1, double w2, double alpha1, double alpha2,
            double alpha3, int nGS, int nBreg, int width, int height, int flag);

  /*************** Minimization Methods **************/

void gs_an(double **D, double **B, double **R, double **E, double **dBx, double **dBy,
           double **dRtheta, double **dEx, double **dEy, double **bBx, double **bBy, double **bRtheta, double **bEx, double **bEy,
           double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double w1, double w2, double alpha1, double alpha2,
           double alpha3, int width, int height, int iter, int flag);

  /*************** Relaxation Methods ***************/
void gsB(double **D, double **B, double **R, double **E, double **dBx, double **dBy,
         double **bBx, double **bBy, double alpha1, double lambda, int width, int height); // update B
// void gsB(double **D, double **B, double **R, double **dBy, double **bBy, double alpha1, int width, int height);

void gsBX(double **B, double **dBx, double **bBx, double lambda1, double alpha1,
          int width, int height);  // update dBx

void gsBY(double **B, double **dBy, double **bBy, double lambda1, double alpha1,
          int width, int height);  // update dBy

void gsRT(double **D, double **B, double **R, double **E, double **dRtheta, double **bRtheta,
              double lambda2, double alpha2, double w1, double w2, int width, int height);

void gsRX(double **R, double **dRtheta, double **bRtheta, double lambda3, double alpha2,
          int width, int height);  //update dRy

void gsRtheta(double **R, double **dRtheta, double **bRtheta, double lambda3, double alpha2,
              double w1, double w2, int width, int height);

void updateE(double **E, double **R, double **D, double w1, double w2, int width, int height);

void gsEX(double **E, double **dEx, double **bEx, double lambda4, double alpha3,
          int width, int height);  // update dEx

void gsEY(double **E, double **dEy, double **bEy, double lambda4, double alpha3,
          int width, int height);  // update dEy


  /**************** Split Bregman Methods *************/
void bregmanBX(double **B, double **dBx, double **bBx, int width, int height);  //update bBx
void bregmanBY(double **B, double **dBy, double **bBy, int width, int height);  //update bBy
void bregmanRX(double **R, double **dRx, double **bRx, int width, int height);  //update bRy
void bregmanRtheta(double **R, double **dRtheta, double **bRtheta, double w1, double w2, int width, int height);
  /**************** Memory ***************************/
double **newMatrix(int rows, int cols);
void deleteMatrix(double **a);
double **get2dArray(const mxArray *mx, int isCopy);
double copy(double **source, double **dest, int rows, int cols);


/****************** The MEX Interface to derain ***************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
       /* get the size of the image */
    int rows = mxGetN(prhs[0]);
    int cols = mxGetM(prhs[0]);
    mexPrintf("rows: %d, cols: %d\n", rows, cols);
       /* get the fidelity and convergence parameters */
    double lambda1 = (double) (mxGetScalar(prhs[3]));
    double lambda2 = (double) (mxGetScalar(prhs[4]));
    double lambda3 = (double) (mxGetScalar(prhs[5]));
    double lambda4 = (double) (mxGetScalar(prhs[6]));
    double lambda5 = (double) (mxGetScalar(prhs[7]));
    double w1 = (double) (mxGetScalar(prhs[8]));
    double w2 = (double) (mxGetScalar(prhs[9]));
    double tol = (double) (mxGetScalar(prhs[10]));

    //double w1 = 1.0;
    //double w2 = 0.0;
    double alpha1 = 0.1;
    double alpha2 = 0.2;
    double alpha3 = 0.1;
    mexPrintf("lambda1: %.4f\n", lambda1);
    mexPrintf("lambda2: %.4f\n", lambda2);
    mexPrintf("lambda3: %.4f\n", lambda3);
    mexPrintf("w1: %.4f\n", w1);
    mexPrintf("w2: %.4f\n", w2);
    mexPrintf("tol: %.4f\n", tol);

       /* get the image, and declare memory */
    double **D = get2dArray(prhs[0], 0);
    double **B = get2dArray(prhs[1], 0);
    double **R = get2dArray(prhs[2], 0);
    double **E = get2dArray(prhs[11], 0);
    double **dBx = newMatrix(rows, cols);
    double **dBy = newMatrix(rows, cols);
    double **bBx = newMatrix(rows, cols);
    double **bBy = newMatrix(rows, cols);
    //double **dRx = newMatrix(rows, cols);
    //double **bRx = newMatrix(rows, cols);
    double **dRtheta = newMatrix(rows, cols);
    double **bRtheta = newMatrix(rows, cols);
    double **dEx = newMatrix(rows, cols);
    double **dEy = newMatrix(rows, cols);
    double **bEx = newMatrix(rows, cols);
    double **bEy = newMatrix(rows, cols);

    double **uOld;
    double *outArrayB;
    double *outArrayR;
    double *outArrayE;
    double diff;
    int count;
    int i, j;

       /* Check conditions ***************/
    if(nrhs != 12) {mexErrMsgTxt("Twelve input arguments required.");}
    if(nlhs != 3) {mexErrMsgTxt("Too many output arguments.");}
    if(!(mxIsDouble(prhs[0]))) {mexErrMsgTxt("Input array must be type of double.");}
    /* denoised the image */
    uOld = newMatrix(rows, cols);
    count = 0;
    int flag = 0;
    for(int i=1; i<rows-1; i++)
    {
        for(int j=0; j<cols; j++)
        {
            dEx[i][j] = (E[i][j]-E[i-1][j]);
            dBx[i][j] = (B[i][j]-B[i-1][j]);
        }
    }
    for(int j = 0; j<cols; j++)
    {
        dEx[0][j] = E[0][j];
        dBx[0][j] = B[0][j];
    }
    
    for(int i=0; i<rows; i++)
    {
        for(int j=1; j<cols; j++)
        {
            dEy[i][j] = (E[i][j] - E[i][j-1]);
            dBy[i][j] = (B[i][j] - B[i][j-1]);
        }
        
    }
    for(int i=0; i<rows; i++)
    {
        dEy[i][0] = E[i][0];
        dBy[i][0] = B[i][0];
    }
    
    for(int j=0; j<cols-1; j++)
    {
        dRtheta[0][j] = R[0][j] - R[0][j+1];
    }
    for(int i=1; i<rows; i++)
    {
        for(int j=0; j<cols-1; j++)
            dRtheta[i][j] = R[i][j] - w1*R[i][j+1] - w2*R[i-1][j+1];
    }
    for(int i=1; i<rows; i++)
        dRtheta[i][cols-1] = w2*(R[i][cols-1]-R[i-1][cols-1]);
    
    do
    {
        if(count%2 == 1) {flag=1;}
        derain(D, B, R, E, dBx, dBy, dRtheta, dEx, dEy, bBx, bBy, bRtheta, bEx, bEy, lambda1, lambda2, lambda3,
               lambda4, lambda5, w1, w2, alpha1, alpha2, alpha3, 1, 10, rows, cols, flag);
        diff = copy(B, uOld, rows, cols);
        count ++;
        printf("------------ tol: %.4f-------------\n", diff);
    }while( ((diff>tol) && count < 100) || count <5);

    /* copy denoised image to output vector */
    plhs[0] = mxCreateDoubleMatrix(cols, rows, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(cols, rows, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(cols, rows, mxREAL);
    outArrayB = mxGetPr(plhs[0]);
    outArrayR = mxGetPr(plhs[1]);
    outArrayE = mxGetPr(plhs[2]);
    for(i=0; i<rows; i++)
    {
        for(j=0; j<cols; j++)
            outArrayB[(i*cols+j)] = B[i][j];
    }
    for(i=0; i<rows; i++)
    {
        for(j=0; j<cols; j++)
            outArrayR[(i*cols+j)] = R[i][j];
    }
    for(i=0; i<rows; i++)
    {
        for(j=0; j<cols; j++)
            outArrayE[(i*cols+j)] = E[i][j];
    }

    /* Free Memory */
    deleteMatrix(dBx);
    deleteMatrix(dBy);
    deleteMatrix(bBx);
    deleteMatrix(bBy);
    //deleteMatrix(dRx);
    //deleteMatrix(bRx);
    deleteMatrix(dRtheta);
    deleteMatrix(bRtheta);
    deleteMatrix(dEx);
    deleteMatrix(dEy);
    deleteMatrix(bEx);
    deleteMatrix(bEy);
    deleteMatrix(uOld);

}

double **get2dArray(const mxArray *mx, int isCopy)
{
    double *oned = mxGetPr(mx);
    int rows = mxGetN(mx);
    int cols = mxGetM(mx);
    double **rval = (double **) malloc(rows*sizeof(double *));
    int r;
    if(isCopy)
    {
        double *copy = (double *) malloc(rows*cols*sizeof(double));
        int i, sent = rows*cols;
        for(i=0; i<sent; i++)
            copy[i] = oned[i];
        oned = copy;
    }
    for(r=0; r<rows; r++)
        rval[r] = &oned[cols*r];
    return rval;
}



/*                 IMPLEMENTATION BELOW THIS LINE        */
void derain(double **D, double **B, double **R, double **E, double **dBx, double **dBy,
            double **dRtheta, double **dEx, double **dEy, double **bBx, double **bBy, double **bRtheta,
            double **bEx, double **bEy, double lambda1, double lambda2, double lambda3, double lambda4,double lambda5,
            double w1, double w2, double alpha1, double alpha2, double alpha3, int nGS, int nBreg, int width, int height, int flag)
{
    int breg;
    for(breg=0; breg<nBreg; breg++)
    {
        gs_an(D, B, R, E, dBx, dBy, dRtheta, dEx, dEy, bBx, bBy, bRtheta, bEx, bEy, lambda1, lambda2, lambda3, lambda4, lambda5, w1, w2, alpha1,
              alpha2, alpha3, width, height, nGS, flag);
        bregmanBX(B, dBx, bBx, width, height);
        bregmanBY(B, dBy, bBy, width, height);
        // bregmanRX(R, dRx, bRx, width, height);  // matlab column first, so u(i, j)-u(i-1, j) becomes u(i, j)-u(i, j-1)
        bregmanRtheta(R, dRtheta, bRtheta, w1, w2, width, height);
        bregmanBX(E, dEx, bEx, width, height); // update bEx the same as bBx
        bregmanBY(E, dEy, bEy, width, height); // update bEy the same as bBy
    }
}

void gs_an(double **D, double **B, double **R, double **E, double **dBx, double **dBy,
           double **dRtheta, double **dEx, double **dEy, double **bBx, double **bBy, double **bRtheta, double **bEx, double **bEy,
           double lambda1, double lambda2, double lambda3, double lambda4, double lambda5, double w1, double w2, double alpha1, double alpha2,
           double alpha3, int width, int height, int iter, int flag)
{
    int j;
    for(j=0; j<iter; j++)
    {
        // gsB(D, B, R, dBy, bBy, alpha1, width, height);
        gsB(D, B, R, E, dBx, dBy, bBx, bBy, alpha1, 0, width, height);
    }
   
    for(j=0; j<iter; j++)
    {
       gsB(D, E, R, B, dEx, dEy, bEx, bEy, alpha3, lambda5, width, height); // update E the same as B
    }
    
    for(j=0; j<iter; j++)
    {
      gsRT(D, B, R, E, dRtheta, bRtheta, lambda2, alpha2, w1, w2, width, height);
      for(int w=0; w<width; w++)
        {
            for(int h=0; h<height; h++)
            {
               if(R[w][h]<0) {R[w][h]=0;}
               if(R[w][h]>D[w][h]) {R[w][h]=D[w][h];}

            }
        }
    }
    /*
     for(j=0; j<iter; j++)
    {
       gsB(D, E, R, B, dEx, dEy, bEx, bEy, alpha3, lambda5, width, height); // update E the same as B
    }
    */
        gsBX(B, dBx, bBx, lambda1, alpha1, width, height);
        gsBY(B, dBy, bBy, lambda1, alpha1, width, height);
        // gsRX(R, dRx, bRx, lambda3, alpha2, width, height);
        gsRtheta(R, dRtheta, bRtheta, lambda3, alpha2, w1, w2, width, height);
        gsEX(E, dEx, bEx, lambda4, alpha3, width, height); // update dEx
        gsEY(E, dEy, bEy, lambda4, alpha3, width, height); // update dEy
        
        // updateE(E, R, D, w1, w2, width, height);
        
/*

        for(int w=0; w<width; w++)
        {
            for(int h=0; h<height; h++)
            {
               if(R[w][h]<0) {R[w][h]=0;}
               if(R[w][h]>D[w][h]) {R[w][h]=D[w][h];}

            }
        }
*/
       

        //gsBX(E, dEx, bEx, -lambda4, alpha3, width, height);
        //gsBY(E, dEy, bEy, -lambda4, alpha3, width, height);

 }


/*** Relaxation operators****/

void gsB(double **D, double **B, double **R, double **E, double **dBx, double **dBy,
         double **bBx, double **bBy, double alpha1, double lambda, int width, int height)
{
    int w, h;
    double sum;
    double normConst = 1.0/(2+4*alpha1+2*lambda);
    int wSent = width - 1, hSent = height - 1;
    for(w=1; w<wSent; w++)                     /* do the central pixels */
    {
        for(h=1; h<hSent; h++)
        {
            sum = B[w][h+1] + B[w][h-1] + B[w+1][h] + B[w-1][h];
            sum += (dBx[w][h]-dBx[w+1][h]+dBy[w][h]-dBy[w][h+1]);
            sum += (-bBx[w][h]+bBx[w+1][h]-bBy[w][h]+bBy[w][h+1]);
            sum *= alpha1;
            sum += (2*(D[w][h]-R[w][h]-E[w][h]));
            sum *= normConst;
            B[w][h] = sum;
        }
    }
    w = 0;
    for(h=1; h<hSent; h++)                   /* do the top pixels */
    {
        sum = B[w+1][h] + B[w][h+1] + B[w][h-1];
        sum += (-dBx[w+1][h]+dBy[w][h]-dBy[w][h+1]);
        sum += (bBx[w+1][h]-bBy[w][h]+bBy[w][h+1]);
        sum *= alpha1;
        sum += (2*(D[w][h]-R[w][h]-E[w][h]));
        sum /= (2+3*alpha1+2*lambda);
        B[w][h] = sum;
    }
    w = wSent;
    for(h=1; h<hSent; h++)                 /* do the bottom pixels */
    {
        sum = B[w][h+1] + B[w-1][h] + B[w][h-1];
        sum += (dBx[w][h]+dBy[w][h]-dBy[w][h+1]);
        sum += (-bBx[w][h]-bBy[w][h]+bBy[w][h+1]);
        sum *= alpha1;
        sum += (2*(D[w][h]-R[w][h]-E[w][h]));
        sum /= (2+3*alpha1+2*lambda);
        B[w][h] = sum;
    }
    h = 0;
    for(w=1; w<wSent; w++)               /* do the left pixels */
    {
        sum = B[w+1][h] + B[w-1][h] + B[w][h+1];
        sum += (dBx[w][h] - dBx[w+1][h] - dBy[w][h+1]);
        sum += (-bBx[w][h] + bBx[w+1][h] + bBy[w][h+1]);
        sum *= alpha1;
        sum += (2*(D[w][h]-R[w][h]-E[w][h]));
        sum /= (2+3*alpha1+2*lambda);
        B[w][h] = sum;
    }
    h = hSent;
    for(w=1; w<wSent; w++)              /* do the right pixels */
    {
        sum = B[w][h-1]+B[w-1][h]+B[w+1][h];
        sum += (dBx[w][h]-dBx[w+1][h]+dBy[w][h]);
        sum += (-bBx[w][h] + bBx[w+1][h] - bBy[w][h]);
        sum *= alpha1;
        sum += (2*(D[w][h]-R[w][h]-E[w][h]));
        sum /= (2+3*alpha1+2*lambda);
        B[w][h] = sum;
    }
    w = h = 0;                         /* do the up-left pixel */
    sum = B[w+1][h] + B[w][h+1];
    sum += (-dBx[w+1][h] - dBy[w][h+1]);
    sum += (bBx[w+1][h] + bBy[w][h+1]);
    sum *= alpha1;
    sum += (2*(D[w][h]-R[w][h]-E[w][h]));
    sum /= (2+2*alpha1+2*lambda);
    B[w][h] = sum;
    w = 0, h = hSent;                 /* do the up-right pixel */
    sum = B[w][h-1]+B[w+1][h];
    sum += (-dBx[w+1][h] + dBy[w][h+1]);
    sum += (bBx[w+1][h] - bBy[w][h+1]);
    sum *= alpha1;
    sum += (2*(D[w][h]-R[w][h]-E[w][h]));
    sum /= (2+2*alpha1+2*lambda);
    B[w][h] = sum;
    w = wSent, h = 0;                /* do the bottom-left pixel */
    sum = B[w][h+1]+B[w-1][h];
    sum += (dBx[w][h] - dBy[w][h+1]);
    sum += (-bBx[w][h] + bBy[w][h+1]);
    sum *= alpha1;
    sum += (2*(D[w][h]-R[w][h]-E[w][h]));
    sum /= (2+2*alpha1+2*lambda);
    B[w][h] = sum;
    w = wSent, h = hSent;           /* do the bottom-right pixel */
    sum = B[w][h-1] + B[w-1][h];
    sum += (dBx[w][h]+dBy[w][h]);
    sum += (-bBx[w][h]-bBy[w][h]);
    sum *= alpha1;
    sum += (2*(D[w][h]-R[w][h]-E[w][h]));
    sum /= (2+2*alpha1+2*lambda);
    B[w][h] = sum;
}

void gsBX(double **B, double **dBx, double **bBx, double lambda1, double alpha1,
          int width, int height)
{
    int w, h;
    double base;
    const double flux = lambda1/alpha1;
    const double mflux = -lambda1/alpha1;
    for(w=1; w<width; w++)
    {
        for(h=0; h<height; h++)
        {
            base = B[w][h]-B[w-1][h]+bBx[w][h];
            if(base > flux) {dBx[w][h] = base - flux; continue;}
            if(base < mflux) {dBx[w][h] = base + flux; continue;}
            dBx[w][h] = 0;
        }
    }
    w = 0;
    for(h=0; h<height; h++)
    {
        base = bBx[w][h];
        if(base > flux) {dBx[w][h]=base - flux;continue;}
        if(base < mflux) {dBx[w][h] = base + flux;continue;}
        dBx[w][h]=0;
    }
}

void gsBY(double **B, double **dBy, double **bBy, double lambda1, double alpha1,
          int width, int height)
{
    int w, h;
    double base;
    const double flux = lambda1/alpha1;
    const double mflux = -lambda1/alpha1;
    for(w=0; w<width; w++)
    {
        for(h=1; h<height; h++)
        {
            base = B[w][h]-B[w][h-1]+bBy[w][h];
            if(base > flux) {dBy[w][h]=base - flux;continue;}
            if(base < mflux) {dBy[w][h] = base + flux; continue;}
            dBy[w][h]=0;
        }
    }
    h = 0;
    for(w=0; w<width; w++)
    {
        base = bBy[w][h];
        if(base > flux) {dBy[w][h] = base - flux; continue;}
        if(base < mflux) {dBy[w][h] = base + flux; continue;}
        dBy[w][h] = 0;
    }
}

void bregmanBX(double **B, double **dBx, double **bBx, int width, int height)
{
    int w, h;
    for(w=1; w<width; w++)
    {
        for(h=0; h<height; h++)
        {
            bBx[w][h] += (B[w][h]-B[w-1][h]-dBx[w][h]);
        }
    }
    w = 0;
    for(h=0; h<height; h++)
        bBx[w][h] += (-dBx[w][h]);
}

void bregmanBY(double **B, double **dBy, double **bBy, int width, int height)
{
    int w, h;
    for(w=0; w<width; w++)
    {
        for(h=1; h<height; h++)
            bBy[w][h] += (B[w][h]-B[w][h-1]-dBy[w][h]);
    }
    h = 0;
    for(w=0; w<width; w++)
        bBy[w][h] += (-dBy[w][h]);
}

/*********************** memory *********************/
double **newMatrix(int rows, int cols)
{
    double *a = (double *) malloc(rows*cols*sizeof(double));
    double **rval = (double **) malloc(rows*sizeof(double *));
    if((a==NULL) || (rval == NULL)) {printf("memory failed.\n");}
    int j, g;
    rval[0] = a;
    for(j=1; j<rows; j++)
        rval[j] = &a[j*cols];
    for(j=0; j<rows; j++)
    {
        for(g=0; g<cols; g++)
            rval[j][g] = 0;
    }
    return rval;
}

void deleteMatrix(double **a)
{
    free(a[0]);
    free(a);
}

double copy(double **source, double **dest, int rows, int cols)
{
    int r, c;
    double temp, sumDiff=0.0, sum=0.0;
    for(r=0; r<rows; r++)
    {
        for(c=0; c<cols; c++)
        {
            temp = dest[r][c];
            sum += temp*temp;
            temp -= source[r][c];
            sumDiff += temp*temp;
            dest[r][c] = source[r][c];
        }
    }
    // printf("---- sumDiff -------- %.5lf", sumDiff);
    // printf("---- sum ------------ %.5lf", sum);
    return sqrt(sumDiff/sum);
}

void gsRT(double **D, double **B, double **R, double **E, double **dRtheta, double **bRtheta,
              double lambda2, double alpha2, double w1, double w2, int width, int height)
{
    int w, h;
    double sum;
    double normConst;
    int wSent=width-1, hSent=height-1;
    for(w=1; w<wSent; w++)                  /* do the central pixels */
    {
        for(h=1; h<hSent; h++)
        {
            sum = (-w1*w2*R[w-1][h]+w2*R[w-1][h+1]+w1*R[w][h+1]-w1*w2*R[w+1][h]+w2*R[w+1][h-1]+w1*R[w][h-1]);
            sum += (dRtheta[w][h]-w1*dRtheta[w][h-1]-w2*dRtheta[w+1][h-1]);
            sum += (-bRtheta[w][h] + w1*bRtheta[w][h-1]+w2*bRtheta[w+1][h-1]);
            sum *= alpha2;
            sum += (2*(D[w][h]-B[w][h]-E[w][h]));
            normConst = 2+2*lambda2+(1+w1*w1+w2*w2)*alpha2;
            sum /= normConst;
            R[w][h] = sum;
        }
    }

    w = 0;
    for(h=1; h<hSent; h++)               /* do the top pixels */
    {
        sum = R[w][h-1]+R[w][h+1]+w2*R[w+1][h-1]-w1*w2*R[w+1][h];
        sum +=(dRtheta[w][h]-dRtheta[w][h-1]-w2*dRtheta[w+1][h-1]);
        sum +=(-bRtheta[w][h]+bRtheta[w][h-1]+w2*bRtheta[w+1][h-1]);
        sum *= alpha2;
        sum += (2*(D[w][h]-B[w][h]-E[w][h]));
        normConst = 2+2*lambda2+(2+w2*w2)*alpha2;
        sum /= normConst;
        R[w][h]=sum;
    }

    w = wSent;
    for(h=1; h<hSent; h++)              /* do the bottom pixels */
    {
        sum = w1*R[w][h-1]-w1*w2*R[w-1][h]+w2*R[w-1][h+1]+w1*R[w][h+1];
        sum += (dRtheta[w][h]-w1*dRtheta[w][h-1]);
        sum += (-bRtheta[w][h] + w1*bRtheta[w][h-1]);
        sum *= alpha2;
        sum += (2*(D[w][h]-B[w][h]-E[w][h]));
        normConst = 2+2*lambda2+(1+w1*w1)*alpha2;
        sum /= normConst;
        R[w][h] = sum;
    }

    h = 0;
    for(w=1; w<wSent; w++)             /* do the left pixels */
    {
        sum = w1*R[w][h+1] + w2*R[w-1][h+1];
        sum += (dRtheta[w][h]);
        sum += (-bRtheta[w][h]);
        sum *= alpha2;
        sum += (2*(D[w][h]-B[w][h]-E[w][h]));
        normConst = 2+2*lambda2+1*alpha2;
        sum /= normConst;
        R[w][h] = sum;
    }

    h = hSent;
    for(w=1; w<wSent; w++)              /* do the right pixels */
    {
        sum = w1*R[w][h-1]+(w2*w2-w1*w2)*R[w-1][h]+(w2*w2-w1*w2)*R[w+1][h]+w2*R[w+1][h-1];
        sum += (w2*dRtheta[w][h]-w1*dRtheta[w][h-1]-w2*dRtheta[w+1][h-1]-w2*dRtheta[w+1][h]);
        sum += (-w2*bRtheta[w][h] + w1*bRtheta[w][h-1] + w2*bRtheta[w+1][h-1] + w2*bRtheta[w+1][h]);
        sum *= alpha2;
        sum += (2*(D[w][h]-B[w][h]-E[w][h]));
        normConst = 2+2*lambda2+(w1*w1+3*w2*w2)*alpha2;
        sum /= normConst;
        R[w][h] = sum;
    }

    w = h = 0;                           /* do the up-left pixel */
    sum = R[w][h+1];
    sum += (dRtheta[w][h]);
    sum += (-bRtheta[w][h]);
    sum *= alpha2;
    sum += (2*(D[w][h]-B[w][h]-E[w][h]));
    normConst = 2+2*lambda2+1*alpha2;
    sum /= normConst;
    R[w][h] = sum;

    w = 0, h = hSent;                   /* do the up-right pixel*/
    sum = R[w][h-1]+w2*R[w+1][h-1]+(w2*w2-w1*w2)*R[w+1][h];
    sum += (-dRtheta[w][h-1]-w2*dRtheta[w+1][h-1]-w2*dRtheta[w+1][h]);
    sum += (bRtheta[w][h-1]+w2*bRtheta[w+1][h-1]+w2*bRtheta[w+1][h]);
    sum *= alpha2;
    sum += (2*(D[w][h]-B[w][h]-E[w][h]));
    normConst = 2+2*lambda2+(1+w2*w2*2)*alpha2;
    sum /= normConst;
    R[w][h] = sum;

    w = wSent, h = 0;                    /* do the bottom-left pixel */
    sum = w1*R[w][h+1]+w2*R[w-1][h+1];
    sum += (dRtheta[w][h]);
    sum += (-bRtheta[w][h]);
    sum *= alpha2;
    sum += (2*(D[w][h]-B[w][h]-E[w][h]));
    normConst = 2+2*lambda2+1*alpha2;
    sum /= normConst;
    R[w][h] = sum;

    w = wSent, h = hSent;               /* do the bottom-right pixel */
    sum = w1*R[w][h-1] + (w2*w2-w1*w2)*R[w-1][h];
    sum += (w2*dRtheta[w][h]-w1*dRtheta[w][h-1]);
    sum += (-w2*bRtheta[w][h] + w1*bRtheta[w][h-1]);
    sum *= alpha2;
    sum += (2*(D[w][h]-B[w][h]-E[w][h]));
    normConst = 2+2*lambda2+(w1*w1+w2*w2)*alpha2;
    sum /= normConst;
    R[w][h] = sum;
}

void gsRtheta(double **R, double **dRtheta, double **bRtheta, double lambda3, double alpha2,
              double w1, double w2, int width, int height)
{
    int w, h;
    double base;
    double flux = lambda3/alpha2;
    double mflux = -lambda3/alpha2;
    int hSent = height-1;
    for(w=1; w<width; w++)                             /* do the central pixels */
    {
        for(h=0; h<hSent; h++)
        {
            base = R[w][h]-w1*R[w][h+1]-w2*R[w-1][h+1]+bRtheta[w][h];
            if(base > flux) {dRtheta[w][h] = base - flux; continue;}
            if(base < mflux) {dRtheta[w][h] = base +flux; continue;}
            dRtheta[w][h]=0;
        }
    }
    w = 0;
    for(h=0; h<hSent; h++)              /* do the top pixels */
    {
        base = R[w][h]-R[w][h+1]+bRtheta[w][h];
        if(base > flux) {dRtheta[w][h] = base - flux; continue;}
        if(base < mflux) {dRtheta[w][h] = base + flux; continue;}
        dRtheta[w][h]=0;
    }
    h = hSent;
    for(w=1; w<width; w++)
    {
        base = w2*R[w][h]-w2*R[w-1][h]+bRtheta[w][h];
        if(base > flux) {dRtheta[w][h] = base - flux; continue;}
        if(base < mflux) {dRtheta[w][h] = base + flux; continue;}
        dRtheta[w][h]=0;
    }
    w = 0, h = hSent;                 /* do the top-right pixel */
    base = bRtheta[w][h];
    if(base > flux)
    {
        dRtheta[w][h] = base - flux;
    }
    else if(base < mflux)
    {
        dRtheta[w][h] = base + flux;
    }
    else
    {
        dRtheta[w][h] = 0;
    }
}

void bregmanRtheta(double **R, double **dRtheta, double **bRtheta, double w1, double w2, int width, int height)
{
    int w, h, hSent = height-1;
    for(w=1; w<width; w++)           /* do central pixels */
    {
        for(h=0; h<hSent; h++)
            bRtheta[w][h] += (R[w][h]-w1*R[w][h+1]-w2*R[w-1][h+1]-dRtheta[w][h]);
    }
    w = 0;
    for(h=0; h<hSent; h++)           /* do top pixels */
    {
        bRtheta[w][h] += (R[w][h]-R[w][h+1]-dRtheta[w][h]);
    }
    h = hSent;
    for(w=1; w<width; w++)
    {
        bRtheta[w][h] += (w2*R[w][h]-w2*R[w-1][h]-dRtheta[w][h]);
    }
    w = 0, h = hSent;               /* do the top-right pixel */
    bRtheta[w][h] += (-dRtheta[w][h]);
}

void updateE(double **E, double **R, double **D, double w1, double w2, int width, int height)
{
    int w, h;
	double tmp;
	double threshold = 0.05*255;
	double rho=0;
    for(w=1; w<width-1; w++)  // do the central pixels
    {
        for(h=1; h<height-1; h++)
        {
			// note that in C language, 1/9 = 0
            // E[w][h] = 1.0/9*(-D[w-1][h-1]-D[w-1][h]-D[w-1][h+1]-D[w][h-1]-D[w][h+1]-D[w+1][h-1]-D[w+1][h]-D[w+1][h+1]+8*D[w][h]);
			// E[w][h] -= 1.0/9*(-R[w-1][h-1]-R[w-1][h]-R[w-1][h+1]-R[w][h-1]-R[w][h+1]-R[w+1][h-1]-R[w+1][h]-R[w+1][h+1]+8*R[w][h]);
			// E[w][h] = E[w][h];// + (R[w][h] - w1*R[w][h+1]-w2*R[w-1][h+1]);
			// double tmp1 = 1.0/9*(-B[w-1][h-1]-B[w-1][h]-B[w-1][h+1]-B[w][h-1]-B[w][h+1]-B[w+1][h-1]-B[w+1][h]-B[w+1][h+1]+8*B[w][h]);
            E[w][h] = 0.5*(D[w][h] - w1*D[w][h+1]-w2*D[w-1][h+1]);
            E[w][h] -= 0.5*(R[w][h] - w1*R[w][h+1]-w2*R[w-1][h+1]);
           // E[w][h] = rho*E[w][h] + (1-rho)*(R[w][h] - w1*R[w][h+1]-w2*R[w-1][h+1]);
			if (abs(E[w][h]) < threshold) E[w][h] = 0;
        }
    }
    w = 0; // do the top pixels
    for(h=1; h<height-1; h++)
    {
        // E[w][h] = 1.0/9*(-2*D[w][h-1]-2*D[w][h+1]-D[w+1][h-1]-D[w+1][h]-D[w+1][h+1]+7*D[w][h]);
        // E[w][h] -= 1.0/9*(-2*R[w][h-1]-2*R[w][h+1]-R[w+1][h-1]-R[w+1][h]-R[w+1][h+1]+7*R[w][h]);
		//E[w][h] = E[w][h];// + (R[w][h]-R[w][h+1]);
       E[w][h] = 0.5*(D[w][h]-D[w][h+1]);
       E[w][h] -= 0.5*(R[w][h]-R[w][h+1]);
        //E[w][h] = rho*E[w][h] +(1-rho)*(R[w][h]-R[w][h+1]);
		if (abs(E[w][h]) < threshold) E[w][h] = 0;
    }
    w = width-1; // do the bottom pixels
    for(h=1; h<height-1; h++)
    {
        // E[w][h] = 1.0/9*(-2*D[w][h-1]-2*D[w][h+1]-D[w-1][h-1]-D[w-1][h]-D[w-1][h+1]+7*D[w][h]);
        // E[w][h] -= 1.0/9*(-2*R[w][h-1]-2*R[w][h+1]-R[w-1][h-1]-R[w-1][h]-R[w-1][h+1]+7*R[w][h]);
		// E[w][h] = E[w][h];// + (R[w][h]-w1*R[w][h+1]-w2*R[w-1][h+1]);
      E[w][h] = 0.5*(D[w][h]-w1*D[w][h+1]-w2*D[w-1][h+1]);
      E[w][h] -= 0.5*(R[w][h]-w1*R[w][h+1]-w2*R[w-1][h+1]);
       // E[w][h] = rho*E[w][h] + (1-rho)*(R[w][h]-w1*R[w][h+1]-w2*R[w-1][h+1]);
		if (abs(E[w][h]) < threshold) E[w][h] = 0;
    }
    h = 0; // do the left pixels
    for(w=1; w<width-1; w++)
    {
        // E[w][h] = 1.0/9*(-2*D[w-1][h]-2*D[w+1][h]-D[w-1][h+1]-D[w][h+1]-D[w+1][h+1]+7*D[w][h]);
        // E[w][h] -= 1.0/9*(-2*R[w-1][h]-2*R[w+1][h]-R[w-1][h+1]-R[w][h+1]-R[w+1][h+1]+7*R[w][h]);
		// E[w][h] = E[w][h];// + (R[w][h] - w1 * R[w][h + 1] - w2 * R[w - 1][h + 1]);
      E[w][h] = 0.5*(D[w][h] - w1 * D[w][h + 1] - w2 * D[w - 1][h + 1]);
      E[w][h] -= 0.5*(R[w][h] - w1 * R[w][h + 1] - w2 * R[w - 1][h + 1]);
        //E[w][h] = rho*E[w][h] + (1-rho)*(R[w][h] - w1 * R[w][h + 1] - w2 * R[w - 1][h + 1]);
		if (abs(E[w][h]) < threshold) E[w][h] = 0;
    }
    h = height-1; // do the right pixels
    for(w=1; w<width-1; w++)
    {
        // E[w][h] = 1.0/9*(-2*D[w-1][h]-2*D[w+1][h]-D[w-1][h-1]-D[w][h-1]-D[w+1][h-1]+7*D[w][h]);
        // E[w][h] -= 1.0/9*(-2*R[w-1][h]-2*R[w+1][h]-R[w-1][h-1]-R[w][h-1]-R[w+1][h-1]+7*R[w][h]);
		// E[w][h] =  E[w][h];// + (w2*R[w][h] - w2 * R[w-1][h]);
       E[w][h] = 0.5*(w2*D[w][h] - w2 * D[w-1][h]);
       E[w][h] -= 0.5*(w2*R[w][h] - w2 * R[w-1][h]);
        //E[w][h] = E[w][h]*rho + (1-rho)*(w2*R[w][h] - w2 * R[w-1][h]);
		if (abs(E[w][h]) < threshold) E[w][h] = 0;
    }
    w = 0; h = 0; // top-left
    //E[w][h] = 1.0/9*(-2*D[w][h+1]-2*D[w+1][h]-D[w+1][h+1] + 5*D[w][h]);
    //E[w][h] -= 1.0/9*(-2*R[w][h+1]-2*R[w+1][h]-R[w+1][h+1] + 5*R[w][h]);
	//E[w][h] =  E[w][h];// + (R[w][h] - R[w][h + 1]);
   E[w][h] = 0.5*(D[w][h] - D[w][h + 1]);
   E[w][h] -= 0.5*(R[w][h] - R[w][h + 1]);
    //E[w][h] = rho*E[w][h] + (1-rho)*(R[w][h] - R[w][h + 1]);
	if (abs(E[w][h]) < threshold) E[w][h] = 0;
    w = 0; h = height-1; // top-right
    //E[w][h] = 1.0/9*(-2*D[w][h-1]-2*D[w+1][h]-D[w+1][h-1] + 5*D[w][h]);
    //E[w][h] -= 1.0/9*(-2*R[w][h-1]-2*R[w+1][h]-R[w+1][h-1] + 5*R[w][h]);
   	E[w][h] =  0.5*D[w][h];
   	E[w][h] -=  0.5*R[w][h];
    //E[w][h] = E[w][h]*rho + (1-rho)*R[w][h];
	if (abs(E[w][h]) < threshold) E[w][h] = 0;
    w = width-1; h = 0; // bottom-left
    // E[w][h] = 1.0/9*(-2*D[w-1][h]-2*D[w][h+1]-D[w-1][h+1] + 5*D[w][h]);
    // E[w][h] -= 1.0/9*(-2*R[w-1][h]-2*R[w][h+1]-R[w-1][h+1] + 5*R[w][h]);
	// E[w][h] = E[w][h]; //+ (R[w][h] - w1 * R[w][h + 1] - w2 * R[w - 1][h + 1]);
   E[w][h] = 0.5*(D[w][h] - w1 * D[w][h + 1] - w2 * D[w - 1][h + 1]);
   E[w][h] -= 0.5*(R[w][h] - w1 * R[w][h + 1] - w2 * R[w - 1][h + 1]);
    //E[w][h] = rho*E[w][h] + (1-rho)*(R[w][h] - w1 * R[w][h + 1] - w2 * R[w - 1][h + 1]);
	if (abs(E[w][h]) < threshold) E[w][h] = 0;
    w = width-1; h = height-1; // bottom-right
    //E[w][h] = 1.0/9*(-2*D[w][h-1]-2*D[w-1][h]-D[w-1][h-1] + 5*D[w][h]);
    //E[w][h] -= 1.0/9*(-2*R[w][h-1]-2*R[w-1][h]-R[w-1][h-1] + 5*R[w][h]);
	//E[w][h] =  E[w][h]; //+ (w2*R[w][h] - w2 * R[w-1][h]);
  E[w][h] = 0.5*(w2*D[w][h] - w2 * D[w-1][h]);
  E[w][h] -= 0.5*(w2*R[w][h] - w2 * R[w-1][h]);
    //E[w][h] = E[w][h]*rho + (1-rho)*(w2*R[w][h] - w2 * R[w-1][h]);
	if (abs(E[w][h]) < threshold) E[w][h] = 0;
}

void gsEX(double **E, double **dEx, double **bEx, double lambda4, double alpha3, int width, int height)
{
    int w, h;
    double base;
    const double flux = lambda4/alpha3;
    for(w=1; w<width; w++)
    {
        for(h=0; h<height; h++)
        {
            base = E[w][h]-E[w-1][h]+bEx[w][h];
            if(base > flux) {dEx[w][h] = base + flux; continue;}
            if(base < -flux) {dEx[w][h] = base - flux; continue;}
            dEx[w][h] = 0;
        }
    }
    w = 0;
    for(h=0; h<height; h++)
    {
        base = bEx[w][h];
        if(base > flux) {dEx[w][h]=base + flux;continue;}
        if(base < -flux) {dEx[w][h] = base - flux;continue;}
        dEx[w][h]=0;
    }
}

void gsEY(double **E, double **dEy, double **bEy, double lambda4, double alpha3,
          int width, int height)
{
    int w, h;
    double base;
    const double flux = lambda4/alpha3;
    for(w=0; w<width; w++)
    {
        for(h=1; h<height; h++)
        {
            base = E[w][h]-E[w][h-1]+bEy[w][h];
            if(base > flux) {dEy[w][h]=base + flux;continue;}
            if(base < -flux) {dEy[w][h] = base - flux; continue;}
            dEy[w][h]=0;
        }
    }
    h = 0;
    for(w=0; w<width; w++)
    {
        base = bEy[w][h];
        if(base > flux) {dEy[w][h] = base + flux; continue;}
        if(base < -flux) {dEy[w][h] = base - flux; continue;}
        dEy[w][h] = 0;
    }
}
