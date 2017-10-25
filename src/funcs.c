#include "mppade.h"


static long int N, Jay, d;
static mpfr_t *u, *y, n0, n1, n2, n3, *U;
static mpfr_t *M; // *wght;
static mpfc_t *P, *Q, *W, *X, *Y, *z, *res, *zsc, *ressc;
static mpfc_t **g, *eye;
static mpfr_t stk1, stk2, stk3, stk4;
static mpfc_t stk1c, stk2c, stk3c, stk4c;
static mpfc_t *A, *B, *C, *D;
static mpfc_t *overQ, **c;
static mpfr_t pie, R, L;
static mpfr_t *absQ2;
//static FILE* fhlog;
static mpfr_t tolGS, tolEA;
static double S, T, G, Rd; //, Ld;
static int prec, tolNEWTON;
static char resname[480];
static int *abFLAG, phaseFLAG;
static int zetaFLAG, eFlag;


void init_program(int *argc, char **argv){ 
 char str[480], line[4096], value[480], param[480], fname[480];
 char v1[480], v2[480], v3[480], v4[480], v5[480];
 char v6[480], v7[480], v8[480], v9[480], v10[480];
 int iFlag = 0;
 mpfr_t aux;

 sprintf(str, "Usage:\n\t%s input", argv[0]);
 if (*argc != 2) err_msg(str);
 sprintf(str, "%s.cfg", argv[1]);
 FILE *fh = fopen(str,"r"); 
 if (fh == NULL) err_msg("Cannot open input file");
 while (fgets(line, 4096, fh)!=NULL) {
    sscanf(line, "%s\t%s", param, value);
    if (strcmp(param,"#resname=") == 0) sprintf(resname,"%s", value);
    if (strcmp(param,"#npoints=") == 0) N = atoi(value);
    if (strcmp(param,"#nrefine=") == 0) Jay = atoi(value);
    if (strcmp(param,"#precisn=") == 0) {
	prec = atoi(value);
 	mpfr_set_default_prec(prec);
 	prec = mpfr_get_default_prec();
	printf("Precision set to %d bits\n", prec);
    }
    if (strcmp(param,"#posit_x=") == 0) T = atof(value);
    if (strcmp(param,"#posit_y=") == 0) S = atof(value);
    if (strcmp(param,"#magnitu=") == 0) G = atof(value);
    if (strcmp(param,"#n_poles=") == 0) d = atoi(value);
    if (strcmp(param,"#zradius=") == 0) Rd = atof(value);
    if (strcmp(param,"#sortphs=") == 0) phaseFLAG = atoi(value);
    if (strcmp(param,"#newttol=") == 0) tolNEWTON = atoi(value);
    if (strcmp(param,"#mulzeta=") == 0) zetaFLAG = atoi(value);
    if (strcmp(param,"#refin_l=") == 0) {
	mpfr_init_set_ui(L, 1, RMODE);
	mpfr_set_str(L, value, 10, RMODE);
	//printf("L = %s\n", value);
	//mpfr_printf("L = %.200Re\n", L);
    }
    if (strcmp(param,"#invertl=") == 0) iFlag = atoi(value);
    if (strcmp(param,"#evlpade=") == 0) {
	eFlag = atoi(value);
	if (eFlag == 1) printf("Just Evaluate Pade on grid and Exit.\n");
    }
    if (strcmp(param,"#padetxt=") == 0) {
	sprintf(fname, "%s", value);
        if (eFlag == 1) printf("Reading two column data from %s\n", fname);	
    }
 }
 fclose(fh);
 if (iFlag == 1) mpfr_ui_div(L, 1, L, RMODE);
 mpfr_init_set_d(tolGS, 0., RMODE); // tolerance of GS coefficients
 mpfr_init_set_d(tolEA, 0.1, RMODE); // tolerance of root finding
 //mpfr_init_set_d(L, Ld, RMODE); // set L in multiprecision
 mpfr_init(aux);
 mpfr_pow_ui(tolEA, tolEA, tolNEWTON, RMODE); 
 mpfr_init(pie); 
 mpfr_init_set_d(R, Rd, RMODE);
 mpfr_mul(R, R, L, RMODE);
 mpfr_const_pi(pie, RMODE);
 u = malloc(N*sizeof(mpfr_t));
 U = malloc(N*sizeof(mpfr_t));
 y = malloc(N*sizeof(mpfr_t));
 eye = malloc(N*sizeof(mpfc_t)); 
 M = malloc(N*sizeof(mpfr_t));
 P = malloc(N*sizeof(mpfc_t));
 Q = malloc(N*sizeof(mpfc_t));
 X = malloc(N*sizeof(mpfc_t));
 Y = malloc(N*sizeof(mpfc_t));
 W = malloc(N*sizeof(mpfc_t));
 overQ = malloc(N*sizeof(mpfc_t));
 absQ2 = malloc(N*sizeof(mpfr_t));
 g = malloc((2*d+1)*sizeof(mpfc_t *));
 c = malloc((2*d+1)*sizeof(mpfc_t *));
 A = malloc(5*sizeof(mpfc_t));
 B = malloc(5*sizeof(mpfc_t));
 C = malloc(5*sizeof(mpfc_t));
 D = malloc(5*sizeof(mpfc_t));
 //stk1c = malloc(sizeof(mpfc_t));
 //stk2c = malloc(sizeof(mpfc_t));
 //stk3c = malloc(sizeof(mpfc_t));
 //stk4c = malloc(sizeof(mpfc_t));
 z = malloc(d*sizeof(mpfc_t));
 zsc = malloc(d*sizeof(mpfc_t));
 res = malloc(d*sizeof(mpfc_t));
 ressc = malloc(d*sizeof(mpfc_t));
 abFLAG = malloc(d*sizeof(int));
 mpfr_t y0; 
 mpfr_init(stk1);
 mpfr_init(stk2);
 mpfr_init(stk3);
 mpfr_init(stk4);
 mpfr_init(stk1c.re);
 mpfr_init(stk1c.im);
 mpfr_init(stk2c.re);
 mpfr_init(stk2c.im);
 mpfr_init(stk3c.re);
 mpfr_init(stk3c.im);
 mpfr_init(stk4c.re);
 mpfr_init(stk4c.im);
 mpfr_init(y0);
 mpfr_init(n0);
 mpfr_init(n1);
 mpfr_init(n2);
 mpfr_init(n3);
 for (int jj = 0; jj < 5; jj++) {
   mpfr_init(A[jj].re);
   mpfr_init(A[jj].im);
   mpfr_init(B[jj].re);
   mpfr_init(B[jj].im);
   mpfr_init(C[jj].re);
   mpfr_init(C[jj].im);
   mpfr_init(D[jj].re);
   mpfr_init(D[jj].im);
 }
 for (int jj = 0; jj < 2*d+1; jj++) {
   g[jj] = malloc(N*sizeof(mpfc_t)); 
   c[jj] = malloc(4*sizeof(mpfc_t));
   for (int ll = 0; ll < 4; ll++) {
     mpfr_init_set_ui(c[jj][ll].re, 0, RMODE);
     mpfr_init_set_ui(c[jj][ll].im, 0, RMODE);
   }
   for (int kk = 0; kk < N; kk++) {
   mpfr_init_set_ui(g[jj][kk].re, 0, RMODE);
   mpfr_init_set_ui(g[jj][kk].im, 0, RMODE);
   }
 }
 for (int jj = 0; jj < d; jj++) {
   mpfr_init(z[jj].re);		mpfr_init(z[jj].im);
   mpfr_init(zsc[jj].re);	mpfr_init(zsc[jj].im);
   mpfr_init(res[jj].re);	mpfr_init(res[jj].im);
   mpfr_init(ressc[jj].re);	mpfr_init(ressc[jj].im);
   abFLAG[jj] = 0;
 }
 for (int jj = 0; jj < N; jj++) {
   mpfr_init(u[jj]);
   mpfr_init_set_si(y[jj], jj-N/2, RMODE);
   mpfr_div_ui(y[jj], y[jj], N, RMODE);
   mpfr_mul(y[jj], y[jj], pie, RMODE);

   mpfr_init_set(U[jj], y[jj], RMODE);
   mpfr_mul_ui(U[jj], U[jj], 2, RMODE);
   
   mpfr_tan(y[jj], y[jj], RMODE);  
   mpfr_div(u[jj], y[jj], L, RMODE); // modified with L


   mpfr_atan(u[jj], u[jj], RMODE);   
   mpfr_mul_ui(u[jj], u[jj], 2, RMODE);
   if (jj == 0) mpfr_mul_si(u[jj], u[jj], -1, RMODE);

   mpfr_init(absQ2[jj]);
   mpfr_init(M[jj]);
   mpfr_init(Q[jj].re);
   mpfr_init(Q[jj].im);
   mpfr_init(P[jj].re);
   mpfr_init(P[jj].im);
   mpfr_init(overQ[jj].re);
   mpfr_init(overQ[jj].im);
   mpfr_init(X[jj].re);
   mpfr_init(X[jj].im);
   mpfr_init(Y[jj].re);
   mpfr_init(Y[jj].im);
   mpfr_init_set_ui(eye[jj].re, 1, RMODE);
   mpfr_init_set_ui(eye[jj].im, 0, RMODE);
   mpfr_init(W[jj].re);
   mpfr_init(W[jj].im);
   mpfr_set_ui(M[jj], 1, RMODE);
 }
 fh = fopen(resname,"r"); 
 //FILE *fh2 = fopen("copy.txt","w");
 int jj = 0;
 int counter = 0;
 if (eFlag == 1) {
   fh = fopen(fname,"r"); 
   if (fh == NULL) err_msg("Cannot open Pade text data.\n"); 
   if (fh != NULL) {
     while (fgets(line, 2048, fh)!=NULL) {
       if (jj == 0) {
         fgets(line, 2048, fh); //dummy = 
         fgets(line, 2048, fh); //dummy =        
         fgets(line, 2048, fh); //dummy = 
       }
       sscanf(line, "%s\t%s\n", v1, v2);
       mpfr_set_str(zsc[jj].re, v1, 10, RMODE);   // replace v2
       mpfr_set_str(ressc[jj].re, v2, 10, RMODE);   // replace v3
       counter++;
       jj++;
     }
   }
   fclose(fh);
   printf("Read %d poles out of %ld expected.\n", counter, d);
   fh = fopen("recontruction.txt","w"); 
   fprintf(fh, "# 1. Re 2. Im\n# From %ld poles in %s\n\n", d, fname);
   for (int jj = 0; jj < N; jj++) {
     mpfr_set_ui(W[jj].re, 0, RMODE);
     mpfr_set_ui(W[jj].im, 0, RMODE);
     for (int dd = 0; dd < d; dd++) {
       mpfr_mul(aux, zsc[dd].re, zsc[dd].re, RMODE);
       mpfr_fma(aux, y[jj], y[jj], aux, RMODE);
       mpfr_div(aux, ressc[dd].re, aux, RMODE);
       mpfr_fma(W[jj].re, y[jj], aux, W[jj].re, RMODE);
       mpfr_fma(W[jj].im, zsc[dd].re, aux, W[jj].im, RMODE);
     }
     mpfr_fprintf(fh, "%.200Re\t%.200Re\t%.200Re\n", u[jj], W[jj].re, W[jj].im);
   }
   fclose(fh);
   //mpfr_ui_div(L, 1, L, RMODE);
   mpfr_printf("%.100Re\n", L);
   mpfr_printf("%.40Re\n", U[1]);
   mpfr_printf("%.40Re\n", y[1]);
   mpfr_printf("%.40Re\n", u[1]);
   exit(1);
 }
 jj = 0;

 mpfr_t Tf, Lf, uf, qf;
 mpfr_inits(Tf, Lf, uf, qf, (mpfr_ptr) NULL);
//# 1. q 2. u 3.-4. Q 5.-6. V 7.-8. Z 9.-10. Phi
//# Time = 1.22706191184548E-01	L = 4.00000000000000E-02	u* = 0.00000000000000E+00	q* = 0.00000000000000E+00
 if (fh == NULL) err_msg("Cannot open restart file");
 if (fh != NULL) {
   while (fgets(line, 4096, fh)!=NULL) {
     if (jj == 0) {
       fgets(line, 4096, fh);//dummy = 
       sscanf(line, "# Time = %s\tL = %s\tu* = %s\tq* = %s", v1, v2, v3, v4);
       fgets(line, 4096, fh);//dummy =        
       //sscanf(line, "# Time = %s\tL = %s\tu* = %s\tq* = %s", v1, v2, v3, v4);
       mpfr_set_str(Tf, v1, 10, RMODE);
       mpfr_set_str(Lf, v2, 10, RMODE);
       mpfr_set_str(uf, v3, 10, RMODE);
       mpfr_set_str(qf, v4, 10, RMODE);
       mpfr_printf("Time = %.12Re:\nConformal Map parameters:\nL = %.64Re\nu* = %.64Re\nq* = %.64Re\n\n", Tf, Lf, uf, qf);
       fgets(line, 4096, fh);//dummy = 
       printf("Last header Line:\t%s", line);
     }
     //sscanf(line, "%s\t%s\t%s\t%s\n", v1, v2, v3, v4);
     //sscanf(line, "%s\t%s\t%s\n", v1, v2, v3);
     sscanf(line, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10);	
     //printf("I = %d\tv3 = %s\tv4 = %s\n", jj, v3, v4);
     // setting 1/z_u
     mpfr_set_str(W[jj].re, v3, 10, RMODE);   // replace v2
     mpfr_set_str(W[jj].im, v4, 10, RMODE);   // replace v3
     mpfr_mul(aux, W[jj].im, W[jj].im, RMODE);
     mpfr_fms(W[jj].re, W[jj].re, W[jj].re, aux, RMODE);

     mpfr_set_str(aux, v3, 10, RMODE);   
     mpfr_mul(W[jj].im, aux, W[jj].im, RMODE);
     mpfr_mul_ui(W[jj].im, W[jj].im, 2, RMODE);
     // now cook-up z_u
     mpfr_mul(aux, W[jj].re, W[jj].re, RMODE);
     mpfr_fma(aux, W[jj].im, W[jj].im, aux, RMODE);

     mpfr_div(W[jj].re, W[jj].re, aux, RMODE);
     mpfr_div(W[jj].im, W[jj].im, aux, RMODE);
     mpfr_neg(W[jj].im, W[jj].im, RMODE);
     
     //mpfr_div_ui(y[jj], u[jj], 2, RMODE);
     //mpfr_tan(y[jj], y[jj], RMODE);
     //mpfr_sub_ui(stk1, y[jj], 1, RMODE);
     //mpfr_sqr(stk2, stk1, RMODE);
     //mpfr_add_ui(stk3, stk2, 1, RMODE);
     //mpfr_div(W[jj].re, stk1, stk3, RMODE);
     //mpfr_ui_div(W[jj].im, 1, stk3, RMODE);
     //mpfr_cosh(W[jj].re, stk1, RMODE);	         // nls soliton test
     //mpfr_ui_div(W[jj].re, 1, W[jj].re, RMODE);  //   
     //mpfr_set_ui(W[jj].im, 0, RMODE);
     //printf("Here\n%s\n", line);
     jj++;   
   }
   for (int j = N-1; j > -1; j--) {
     mpfr_sub(W[j].re, W[j].re, W[0].re, RMODE);
     mpfr_sub(W[j].im, W[j].im, W[0].im, RMODE);
   }

   write_cmplx(W, "W.txt");
   printf("Lines read = %d of %ld expected\n", jj, N);

   if ( jj != N) err_msg("I/O error\n");
   fclose(fh);
 }

 if (zetaFLAG == 1) {
   for (int jj = 0; jj < N; jj++) {
     mpfr_fma(aux, y[jj], y[jj], eye[jj].re, RMODE);
     mpfr_div_ui(aux, aux, 2, RMODE);
     mpfr_div(W[jj].re, W[jj].re, aux, RMODE);
     mpfr_div(W[jj].im, W[jj].im, aux, RMODE);
   }
 }
}

void generate_Q0() {
  //printf("T,S = %f\t%f\n", T, S);
  for (int jj = 0; jj < N; jj++) {
    mpfr_set_ui(Q[jj].re, 1, RMODE);
    mpfr_set_zero(Q[jj].im, 1);
    for (int j = 0; j < d; j++) {
      mpfr_set_d(B[0].re, T, RMODE);
      mpfr_set_d(B[0].im, S, RMODE);
      mpfr_mul_si(B[0].im, B[0].im, 2*j-d+1, RMODE);
      mpfr_div_ui(B[0].im, B[0].im, 2*d, RMODE);
      mpfr_add(B[0].im, B[0].im, y[jj], RMODE);
      mpfr_mul(stk1, Q[jj].im, B[0].im, RMODE);
      mpfr_mul(stk2, Q[jj].re, B[0].im, RMODE);
      mpfr_fms(Q[jj].re, Q[jj].re, B[0].re, stk1, RMODE);
      mpfr_fma(Q[jj].im, Q[jj].im, B[0].re, stk2, RMODE);
    }
  }
}

void set_weight() {
  invert(overQ, Q);
  absolute2(absQ2, overQ);
  for (int jj = 0; jj < N; jj++) {
    
    mpfr_div(M[jj], y[jj], L, RMODE);
    mpfr_fma(M[jj], M[jj], M[jj], eye[jj].re, RMODE);
    mpfr_div_ui(M[jj], M[jj], 2, RMODE);
    mpfr_mul(M[jj], M[jj], L, RMODE);

    //mpfr_fma(M[jj], y[jj], y[jj], eye[jj].re, RMODE);  // modified to include L
    //mpfr_div_ui(M[jj], M[jj], 2, RMODE);
    mpfr_mul(M[jj], M[jj], absQ2[jj], RMODE);
  }
  //write_cmplx(overQ, "overQ.dat");
}

void gram_schmidt() {
  assign_cmplx(g[0], W);
  dotpr(&n0, g[0], g[0]);
  dotpc(&stk2c, eye, g[0]);
  //dotpr(&stk2, eye, g[0]);
  mpfr_div(c[0][0].re, stk2c.re, n0, RMODE);
  mpfr_div(c[0][0].im, stk2c.im, n0, RMODE);
  for (int jj = 0; jj < N; jj++) {
    //mpfr_mul(Y[jj].re, c[0][0], g[0][jj].re, RMODE);
    mpfr_mul(Y[jj].re, c[0][0].im, g[0][jj].im, RMODE);
    mpfr_fms(Y[jj].re, c[0][0].re, g[0][jj].re, Y[jj].re, RMODE);

    mpfr_mul(Y[jj].im, c[0][0].re, g[0][jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[0][0].im, g[0][jj].re, Y[jj].im, RMODE);

    mpfr_sub(g[1][jj].re, eye[jj].re, Y[jj].re, RMODE);
    mpfr_sub(g[1][jj].im, eye[jj].im, Y[jj].im, RMODE);
    mpfr_mul(X[jj].re, g[0][jj].im, y[jj], RMODE);  
    mpfr_mul(X[jj].im, g[0][jj].re, y[jj], RMODE);
    mpfr_neg(X[jj].re, X[jj].re, RMODE);
  }
  dotpr(&n1, g[1], g[1]);
  //dotpr(&stk2, X, g[0]);
  dotpc(&stk2c, X, g[0]);
  //mpfr_div(c[1][1], stk2, n0, RMODE);
  mpfr_div(c[1][1].re, stk2c.re, n0, RMODE);
  mpfr_div(c[1][1].im, stk2c.im, n0, RMODE);
  if (mpfr_cmp(n0, tolGS) < 0) {
    mpfr_set_zero(c[1][1].re, 1);
    mpfr_set_zero(c[1][1].im, 1);
  }
  //dotpr(&stk2, X, g[1]);
  dotpc(&stk2c, X, g[1]);
  //mpfr_div(c[1][0], stk2, n1, RMODE);
  mpfr_div(c[1][0].re, stk2c.re, n1, RMODE);
  mpfr_div(c[1][0].im, stk2c.im, n1, RMODE);
  if (mpfr_cmp(n0, tolGS) < 0) {
    mpfr_set_zero(c[1][0].re, 1);
    mpfr_set_zero(c[1][0].im, 1);
  }
  for (int jj = 0; jj < N; jj++) {
    /*
    mpfr_mul(Y[jj].re, c[1][0], g[1][jj].re, RMODE);
    mpfr_mul(Y[jj].re, c[1][0], g[1][jj].re, RMODE);

    mpfr_mul(Y[jj].im, c[1][0], g[1][jj].im, RMODE);
    mpfr_mul(Y[jj].im, c[1][0], g[1][jj].im, RMODE);
    mpfr_fma(Y[jj].re, c[1][1], g[0][jj].re, Y[jj].re, RMODE);
    mpfr_fma(Y[jj].re, c[1][1], g[0][jj].re, Y[jj].re, RMODE);
    mpfr_fma(Y[jj].im, c[1][1], g[0][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[1][1], g[0][jj].im, Y[jj].im, RMODE);
    mpfr_sub(g[2][jj].re, X[jj].re, Y[jj].re, RMODE);
    mpfr_sub(g[2][jj].re, X[jj].re, Y[jj].re, RMODE);
    mpfr_sub(g[2][jj].im, X[jj].im, Y[jj].im, RMODE);
    mpfr_sub(g[2][jj].im, X[jj].im, Y[jj].im, RMODE);
    mpfr_mul(X[jj].re, g[1][jj].im, y[jj], RMODE);  
    mpfr_mul(X[jj].re, g[1][jj].im, y[jj], RMODE);  
    mpfr_mul(X[jj].im, g[1][jj].re, y[jj], RMODE);
    mpfr_mul(X[jj].im, g[1][jj].re, y[jj], RMODE);
    mpfr_neg(X[jj].re, X[jj].re, RMODE);
    mpfr_neg(X[jj].re, X[jj].re, RMODE);
    */




    mpfr_mul(Y[jj].re, c[1][0].im, g[1][jj].im, RMODE);
    mpfr_fms(Y[jj].re, c[1][0].re, g[1][jj].re, Y[jj].re, RMODE);
    mpfr_mul(Y[jj].im, c[1][0].re, g[1][jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[1][0].im, g[1][jj].re, Y[jj].im, RMODE);


    // mpfr_fma(Y[jj].re, c[1][1], g[0][jj].re, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[1][1].im, g[0][jj].im, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[1][1].re, g[0][jj].re, Y[jj].re, RMODE);    

    // mpfr_fma(Y[jj].im, c[1][1], g[0][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[1][1].re, g[0][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[1][1].im, g[0][jj].re, Y[jj].im, RMODE);
    
    mpfr_sub(g[2][jj].re, X[jj].re, Y[jj].re, RMODE);
    mpfr_sub(g[2][jj].im, X[jj].im, Y[jj].im, RMODE);
    mpfr_mul(X[jj].re, g[1][jj].im, y[jj], RMODE);  
    mpfr_mul(X[jj].im, g[1][jj].re, y[jj], RMODE);
    mpfr_neg(X[jj].re, X[jj].re, RMODE);
    
  }
  dotpr(&n2, g[2], g[2]);
  //dotpr(&stk2, X, g[0]);
  dotpc(&stk2c, X, g[0]);
  //mpfr_div(c[2][2], stk2, n0, RMODE);
  mpfr_div(c[2][2].re, stk2c.re, n0, RMODE);
  mpfr_div(c[2][2].im, stk2c.im, n0, RMODE);        
  if (mpfr_cmp(n0, tolGS) < 0) {
    mpfr_set_zero(c[2][2].re, 1);
    mpfr_set_zero(c[2][2].im, 1);
  }
  //dotpr(&stk2, X, g[1]);
  dotpc(&stk2c, X, g[1]);
  //mpfr_div(c[2][1], stk2, n1, RMODE);
  mpfr_div(c[2][1].re, stk2c.re, n1, RMODE);
  mpfr_div(c[2][1].im, stk2c.im, n1, RMODE);
  if (mpfr_cmp(n1, tolGS) < 0) {
    mpfr_set_zero(c[2][1].re, 1);
    mpfr_set_zero(c[2][1].im, 1);
  }
  //dotpr(&stk2, X, g[2]);
  dotpc(&stk2c, X, g[2]);
  //mpfr_div(c[2][0], stk2, n2, RMODE);
  mpfr_div(c[2][0].re, stk2c.re, n2, RMODE);
  mpfr_div(c[2][0].im, stk2c.im, n2, RMODE);
  if (mpfr_cmp(n2, tolGS) < 0) {
    mpfr_set_zero(c[2][0].re, 1);
    mpfr_set_zero(c[2][0].im, 1);
  }
  for (int jj = 0; jj < N; jj++) {
    mpfr_mul(Y[jj].re, c[2][0].im, g[2][jj].im, RMODE);
    mpfr_fms(Y[jj].re, c[2][0].re, g[2][jj].re, Y[jj].re, RMODE);

    mpfr_mul(Y[jj].im, c[2][0].re, g[2][jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[2][0].im, g[2][jj].re, Y[jj].im, RMODE);

    //mpfr_fma(Y[jj].re, c[2][1], g[1][jj].re, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[2][1].im, g[1][jj].im, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[2][1].re, g[1][jj].re, Y[jj].re, RMODE);
   
    //mpfr_fma(Y[jj].im, c[2][1], g[1][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[2][1].re, g[1][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[2][1].im, g[1][jj].re, Y[jj].im, RMODE);

    //mpfr_fma(Y[jj].re, c[2][2], g[0][jj].re, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[2][2].im, g[0][jj].im, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[2][2].re, g[0][jj].re, Y[jj].re, RMODE);

    //mpfr_fma(Y[jj].im, c[2][2], g[0][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[2][2].re, g[0][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[2][2].im, g[0][jj].re, Y[jj].im, RMODE);

    if ( d > 1) {
      mpfr_sub(g[3][jj].re, X[jj].re, Y[jj].re, RMODE);
      mpfr_sub(g[3][jj].im, X[jj].im, Y[jj].im, RMODE);
    }
    mpfr_mul(X[jj].re, g[2][jj].im, y[jj], RMODE);  
    mpfr_mul(X[jj].im, g[2][jj].re, y[jj], RMODE);
    mpfr_neg(X[jj].re, X[jj].re, RMODE);
  }
  if (d == 1) {
   //mpfr_printf("c00 = %.10Re + I%.10Re\tc10 = %.10Re+ I%.10Re\tc11 = %.10Re+ I%.10Re\n", c[0][0].re,c[0][0].im , c[1][0].re,c[1][0].im, c[1][1].re, c[1][1].im);
   return;
  }
  dotpr(&n3, g[3], g[3]);
  //dotpr(&stk2, X, g[0]);
  dotpc(&stk2c, X, g[0]);
  //mpfr_div(c[3][3], stk2, n0, RMODE);
  mpfr_div(c[3][3].re, stk2c.re, n0, RMODE);
  mpfr_div(c[3][3].im, stk2c.im, n0, RMODE);
  if (mpfr_cmp(n0, tolGS) < 0) {
    mpfr_set_zero(c[3][3].re, 1);
    mpfr_set_zero(c[3][3].im, 1);
  }
  //dotpr(&stk2, X, g[1]);
  dotpc(&stk2c, X, g[1]);
  //mpfr_div(c[3][2], stk2, n1, RMODE);
  mpfr_div(c[3][2].re, stk2c.re, n1, RMODE);
  mpfr_div(c[3][2].im, stk2c.im, n1, RMODE);
  if (mpfr_cmp(n1, tolGS) < 0) {
    mpfr_set_zero(c[3][2].re, 1);
    mpfr_set_zero(c[3][2].im, 1);
  }
  //dotpr(&stk2, X, g[2]);
  dotpc(&stk2c, X, g[2]);
  //mpfr_div(c[3][1], stk2, n2, RMODE);
  mpfr_div(c[3][1].re, stk2c.re, n2, RMODE);
  mpfr_div(c[3][1].im, stk2c.im, n2, RMODE);
  if (mpfr_cmp(n2, tolGS) < 0) {
    mpfr_set_zero(c[3][1].re, 1);
    mpfr_set_zero(c[3][1].im, 1);
  }
  //dotpr(&stk2, X, g[3]);
  dotpc(&stk2c, X, g[3]);
  //mpfr_div(c[3][0], stk2, n3, RMODE);
  mpfr_div(c[3][0].re, stk2c.re, n3, RMODE);
  mpfr_div(c[3][0].im, stk2c.im, n3, RMODE);
  if (mpfr_cmp(n3, tolGS) < 0) {
    mpfr_set_zero(c[3][0].re, 1);
    mpfr_set_zero(c[3][0].im, 1);
  }
  for (int jj = 0; jj < N; jj++) {
    //mpfr_mul(Y[jj].re, c[3][0], g[3][jj].re, RMODE);
    mpfr_mul(Y[jj].re, c[3][0].im, g[3][jj].im, RMODE);
    mpfr_fms(Y[jj].re, c[3][0].re, g[3][jj].re, Y[jj].re, RMODE);
    //mpfr_mul(Y[jj].im, c[3][0], g[3][jj].im, RMODE);
    mpfr_mul(Y[jj].im, c[3][0].re, g[3][jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[3][0].im, g[3][jj].re, Y[jj].im, RMODE);

    //mpfr_fma(Y[jj].re, c[3][1], g[2][jj].re, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[3][1].im, g[2][jj].im, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[3][1].re, g[2][jj].re, Y[jj].re, RMODE);

    //mpfr_fma(Y[jj].im, c[3][1], g[2][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[3][1].re, g[2][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[3][1].im, g[2][jj].re, Y[jj].im, RMODE);

    //mpfr_fma(Y[jj].re, c[3][2], g[1][jj].re, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[3][2].im, g[1][jj].im, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[3][2].re, g[1][jj].re, Y[jj].re, RMODE);

    //mpfr_fma(Y[jj].im, c[3][2], g[1][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[3][2].re, g[1][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[3][2].im, g[1][jj].re, Y[jj].im, RMODE);

    //mpfr_fma(Y[jj].re, c[3][3], g[0][jj].re, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[3][3].im, g[0][jj].im, Y[jj].re, RMODE);
    mpfr_fms(Y[jj].re, c[3][3].re, g[0][jj].re, Y[jj].re, RMODE);

    //mpfr_fma(Y[jj].im, c[3][3], g[0][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[3][3].re, g[0][jj].im, Y[jj].im, RMODE);
    mpfr_fma(Y[jj].im, c[3][3].im, g[0][jj].re, Y[jj].im, RMODE);

    mpfr_sub(g[4][jj].re, X[jj].re, Y[jj].re, RMODE);
    mpfr_sub(g[4][jj].im, X[jj].im, Y[jj].im, RMODE);
    mpfr_mul(X[jj].re, g[3][jj].im, y[jj], RMODE);  
    mpfr_mul(X[jj].im, g[3][jj].re, y[jj], RMODE);
    mpfr_neg(X[jj].re, X[jj].re, RMODE);
  }
  for (int k = 4; k < 2*d; k++) {
    mpfr_set(n0, n1, RMODE);
    mpfr_set(n1, n2, RMODE); 
    mpfr_set(n2, n3, RMODE); 
    dotpr(&n3, g[k], g[k]);
    //dotpr(&stk2, X, g[k-3]);
    dotpc(&stk2c, X, g[k-3]);
    //mpfr_div(c[k][3], stk2, n0, RMODE);
    mpfr_div(c[k][3].re, stk2c.re, n0, RMODE);
    mpfr_div(c[k][3].im, stk2c.im, n0, RMODE);
    if (mpfr_cmp(n0, tolGS) < 0) {
      mpfr_set_zero(c[k][3].re, 1);
      mpfr_set_zero(c[k][3].im, 1);
    } 
    //dotpr(&stk2, X, g[k-2]);
    dotpc(&stk2c, X, g[k-2]);
    //mpfr_div(c[k][2], stk2, n1, RMODE);
    mpfr_div(c[k][2].re, stk2c.re, n1, RMODE);
    mpfr_div(c[k][2].im, stk2c.im, n1, RMODE);
    if (mpfr_cmp(n1, tolGS) < 0) {
      mpfr_set_zero(c[k][2].re, 1);
      mpfr_set_zero(c[k][2].im, 1);
    }
    dotpc(&stk2c, X, g[k-1]);
    mpfr_div(c[k][1].re, stk2c.re, n2, RMODE);
    mpfr_div(c[k][1].im, stk2c.im, n2, RMODE);
    if (mpfr_cmp(n2, tolGS) < 0) {
      mpfr_set_zero(c[k][1].re, 1);
      mpfr_set_zero(c[k][1].im, 1);
    }
    dotpc(&stk2c, X, g[k]);
    mpfr_div(c[k][0].re, stk2c.re, n3, RMODE);
    mpfr_div(c[k][0].im, stk2c.im, n3, RMODE);
    if (mpfr_cmp(n3, tolGS) < 0) {
      mpfr_set_zero(c[k][0].re, 1);
      mpfr_set_zero(c[k][0].im, 1);
    }
    for (int jj = 0; jj < N; jj++) {
      //mpfr_mul(Y[jj].re, c[k][0], g[k-0][jj].re, RMODE);
      mpfr_mul(Y[jj].re, c[k][0].im, g[k-0][jj].im, RMODE);
      mpfr_fms(Y[jj].re, c[k][0].re, g[k-0][jj].re, Y[jj].re, RMODE);

      //mpfr_fma(Y[jj].re, c[k][1], g[k-1][jj].re, Y[jj].re, RMODE);
      mpfr_fms(Y[jj].re, c[k][1].im, g[k-1][jj].im, Y[jj].re, RMODE);
      mpfr_fms(Y[jj].re, c[k][1].re, g[k-1][jj].re, Y[jj].re, RMODE);

      //mpfr_fma(Y[jj].re, c[k][2], g[k-2][jj].re, Y[jj].re, RMODE);
      mpfr_fms(Y[jj].re, c[k][2].im, g[k-2][jj].im, Y[jj].re, RMODE);
      mpfr_fms(Y[jj].re, c[k][2].re, g[k-2][jj].re, Y[jj].re, RMODE);

      //mpfr_fma(Y[jj].re, c[k][3], g[k-3][jj].re, Y[jj].re, RMODE);
      mpfr_fms(Y[jj].re, c[k][3].im, g[k-3][jj].im, Y[jj].re, RMODE);
      mpfr_fms(Y[jj].re, c[k][3].re, g[k-3][jj].re, Y[jj].re, RMODE);

      //mpfr_mul(Y[jj].im, c[k][0], g[k-0][jj].im, RMODE);
      mpfr_mul(Y[jj].im, c[k][0].re, g[k-0][jj].im, RMODE);
      mpfr_fma(Y[jj].im, c[k][0].im, g[k-0][jj].re, Y[jj].im, RMODE);

      //mpfr_fma(Y[jj].im, c[k][1], g[k-1][jj].im, Y[jj].im, RMODE);
      mpfr_fma(Y[jj].im, c[k][1].re, g[k-1][jj].im, Y[jj].im, RMODE);
      mpfr_fma(Y[jj].im, c[k][1].im, g[k-1][jj].re, Y[jj].im, RMODE);

      //mpfr_fma(Y[jj].im, c[k][2], g[k-2][jj].im, Y[jj].im, RMODE);
      mpfr_fma(Y[jj].im, c[k][2].re, g[k-2][jj].im, Y[jj].im, RMODE);
      mpfr_fma(Y[jj].im, c[k][2].im, g[k-2][jj].re, Y[jj].im, RMODE);

      //mpfr_fma(Y[jj].im, c[k][3], g[k-3][jj].im, Y[jj].im, RMODE);
      mpfr_fma(Y[jj].im, c[k][3].re, g[k-3][jj].im, Y[jj].im, RMODE);
      mpfr_fma(Y[jj].im, c[k][3].im, g[k-3][jj].re, Y[jj].im, RMODE);

      mpfr_sub(g[k+1][jj].re, X[jj].re, Y[jj].re, RMODE);
      mpfr_sub(g[k+1][jj].im, X[jj].im, Y[jj].im, RMODE);
      mpfr_mul(X[jj].re, g[k][jj].im, y[jj], RMODE);  
      mpfr_mul(X[jj].im, g[k][jj].re, y[jj], RMODE);
      mpfr_neg(X[jj].re, X[jj].re, RMODE);
    }
  }
}

void update_Q(){
  mpfr_set_ui(B[0].re, 1, RMODE);
  mpfr_set_ui(B[0].im, 0, RMODE);
  //mpfr_neg(B[1].re, c[0][0], RMODE);
  //mpfr_set_ui(B[1].im, 0, RMODE);
  mpfr_neg(B[1].re, c[0][0].re, RMODE);
  mpfr_neg(B[1].im, c[0][0].im, RMODE);
  if (d == 1) {
    for (int jj = 0; jj < N; jj++) {
      //mpfr_mul(Q[jj].re, c[1][0], B[1].re, RMODE);
      mpfr_mul(Q[jj].re, c[1][0].im, B[1].im, RMODE);
      mpfr_fms(Q[jj].re, c[1][0].re, B[1].re, Q[jj].re, RMODE);

      //mpfr_fma(Q[jj].re, c[1][1], B[0].re, Q[jj].re, RMODE);
      mpfr_fms(Q[jj].re, c[1][1].im, B[0].im, Q[jj].re, RMODE);
      mpfr_fms(Q[jj].re, c[1][1].re, B[0].re, Q[jj].re, RMODE);

      mpfr_fma(Q[jj].re, y[jj], B[0].im, Q[jj].re, RMODE);
      mpfr_neg(Q[jj].re, Q[jj].re, RMODE);                    // update Q

      //mpfr_mul(Q[jj].im, c[1][0], B[1].im, RMODE);
      mpfr_mul(Q[jj].im, c[1][0].re, B[1].im, RMODE);
      mpfr_fma(Q[jj].im, c[1][0].im, B[1].re, Q[jj].im, RMODE);

      //mpfr_fma(Q[jj].im, c[1][1], B[0].im, Q[jj].im, RMODE);
      mpfr_fma(Q[jj].im, c[1][1].re, B[0].im, Q[jj].im, RMODE);
      mpfr_fma(Q[jj].im, c[1][1].im, B[0].re, Q[jj].im, RMODE);

      mpfr_fms(Q[jj].im, y[jj], B[0].re, Q[jj].im, RMODE);
    }
    return; 
  } else if (d == 2) { 
    for (int jj = 0; jj < N; jj++) {
      //mpfr_mul(B[2].re, c[1][0], B[1].re, RMODE);
      mpfr_mul(B[2].re, c[1][0].im, B[1].im, RMODE);
      mpfr_fms(B[2].re, c[1][0].re, B[1].re, B[2].re, RMODE);
      
      //mpfr_fma(B[2].re, c[1][1], B[0].re, B[2].re, RMODE);
      mpfr_fms(B[2].re, c[1][1].im, B[0].im, B[2].re, RMODE);
      mpfr_fms(B[2].re, c[1][1].re, B[0].re, B[2].re, RMODE);
      mpfr_fma(B[2].re, y[jj], B[0].im, B[2].re, RMODE);
      mpfr_neg(B[2].re, B[2].re, RMODE);

      //mpfr_mul(B[2].im, c[1][0], B[1].im, RMODE);
      mpfr_mul(B[2].im, c[1][0].re, B[1].im, RMODE);
      mpfr_fma(B[2].im, c[1][0].im, B[1].re, B[2].im, RMODE);

      //mpfr_fma(B[2].im, c[1][1], B[0].im, B[2].im, RMODE);  // b2
      mpfr_fma(B[2].im, c[1][1].re, B[0].im, B[2].im, RMODE);
      mpfr_fma(B[2].im, c[1][1].im, B[0].re, B[2].im, RMODE);

      mpfr_fms(B[2].im, y[jj], B[0].re, B[2].im, RMODE);
      //------------------------------------------------
      //mpfr_mul(B[3].re, c[2][0], B[2].re, RMODE);
      mpfr_mul(B[3].re, c[2][0].im, B[2].im, RMODE);
      mpfr_fms(B[3].re, c[2][0].re, B[2].re, B[3].re, RMODE);

      //mpfr_fma(B[3].re, c[2][1], B[1].re, B[3].re, RMODE);
      mpfr_fms(B[3].re, c[2][1].im, B[1].im, B[3].re, RMODE);
      mpfr_fms(B[3].re, c[2][1].re, B[1].re, B[3].re, RMODE);

      //mpfr_fma(B[3].re, c[2][2], B[0].re, B[3].re, RMODE);
      mpfr_fms(B[3].re, c[2][2].im, B[0].im, B[3].re, RMODE);
      mpfr_fms(B[3].re, c[2][2].re, B[0].re, B[3].re, RMODE);

      mpfr_fma(B[3].re, y[jj], B[1].im, B[3].re, RMODE);
      mpfr_neg(B[3].re, B[3].re, RMODE);                    // b3

      //mpfr_mul(B[3].im, c[2][0], B[2].im, RMODE);
      mpfr_mul(B[3].im, c[2][0].re, B[2].im, RMODE);
      mpfr_fma(B[3].im, c[2][0].im, B[2].re, B[3].im, RMODE);

      //mpfr_fma(B[3].im, c[2][1], B[1].im, B[3].im, RMODE);
      mpfr_fma(B[3].im, c[2][1].re, B[1].im, B[3].im, RMODE);
      mpfr_fma(B[3].im, c[2][1].im, B[1].re, B[3].im, RMODE);

      //mpfr_fma(B[3].im, c[2][2], B[0].im, B[3].im, RMODE);
      mpfr_fma(B[3].im, c[2][2].re, B[0].im, B[3].im, RMODE);
      mpfr_fma(B[3].im, c[2][2].im, B[0].re, B[3].im, RMODE);

      mpfr_fms(B[3].im, y[jj], B[1].re, B[3].im, RMODE);
      //-----------------------------------------------
      //mpfr_mul(Q[jj].re, c[3][0], B[3].re, RMODE);
      mpfr_mul(Q[jj].re, c[3][0].im, B[3].im, RMODE);
      mpfr_fms(Q[jj].re, c[3][0].re, B[3].re, Q[jj].re, RMODE);

      //mpfr_fma(Q[jj].re, c[3][1], B[2].re, Q[jj].re, RMODE);
      mpfr_fms(Q[jj].re, c[3][1].im, B[2].im, Q[jj].re, RMODE);
      mpfr_fms(Q[jj].re, c[3][1].re, B[2].re, Q[jj].re, RMODE);

      //mpfr_fma(Q[jj].re, c[3][2], B[1].re, Q[jj].re, RMODE);
      mpfr_fms(Q[jj].re, c[3][2].im, B[1].im, Q[jj].re, RMODE);
      mpfr_fms(Q[jj].re, c[3][2].re, B[1].re, Q[jj].re, RMODE);

      //mpfr_fma(Q[jj].re, c[3][3], B[0].re, Q[jj].re, RMODE);
      mpfr_fms(Q[jj].re, c[3][3].im, B[0].im, Q[jj].re, RMODE);
      mpfr_fms(Q[jj].re, c[3][3].re, B[0].re, Q[jj].re, RMODE);

      mpfr_fma(Q[jj].re, y[jj], B[2].im, Q[jj].re, RMODE);
      mpfr_neg(Q[jj].re, Q[jj].re, RMODE);                    // update Q
      //mpfr_mul(Q[jj].im, c[3][0], B[3].im, RMODE);
      mpfr_mul(Q[jj].im, c[3][0].re, B[3].im, RMODE);
      mpfr_fma(Q[jj].im, c[3][0].im, B[3].re, Q[jj].im, RMODE);

      //mpfr_fma(Q[jj].im, c[3][1], B[2].im, Q[jj].im, RMODE);
      mpfr_fma(Q[jj].im, c[3][1].re, B[2].im, Q[jj].im, RMODE);
      mpfr_fma(Q[jj].im, c[3][1].im, B[2].re, Q[jj].im, RMODE);

      //mpfr_fma(Q[jj].im, c[3][2], B[1].im, Q[jj].im, RMODE);
      mpfr_fma(Q[jj].im, c[3][2].re, B[1].im, Q[jj].im, RMODE);
      mpfr_fma(Q[jj].im, c[3][2].im, B[1].re, Q[jj].im, RMODE);

      //mpfr_fma(Q[jj].im, c[3][3], B[0].im, Q[jj].im, RMODE);
      mpfr_fma(Q[jj].im, c[3][3].re, B[0].im, Q[jj].im, RMODE);
      mpfr_fma(Q[jj].im, c[3][3].im, B[0].re, Q[jj].im, RMODE);

      mpfr_fms(Q[jj].im, y[jj], B[2].re, Q[jj].im, RMODE);
    }
  } else {
    for (int jj = 0; jj < N; jj++) {
      mpfr_set_ui(B[0].re, 1, RMODE);
      mpfr_set_ui(B[0].im, 0, RMODE);
      //mpfr_neg(B[1].re, c[0][0], RMODE);
      //mpfr_set_ui(B[1].im, 0, RMODE);
      mpfr_neg(B[1].re, c[0][0].re, RMODE);
      mpfr_neg(B[1].im, c[0][0].im, RMODE);

      //mpfr_mul(B[2].re, c[1][0], B[1].re, RMODE);
      mpfr_mul(B[2].re, c[1][0].im, B[1].im, RMODE);
      mpfr_fms(B[2].re, c[1][0].re, B[1].re, B[2].re, RMODE);

      //mpfr_fma(B[2].re, c[1][1], B[0].re, B[2].re, RMODE);
      mpfr_fms(B[2].re, c[1][1].im, B[0].im, B[2].re, RMODE);
      mpfr_fms(B[2].re, c[1][1].re, B[0].re, B[2].re, RMODE);

      mpfr_fma(B[2].re, y[jj], B[0].im, B[2].re, RMODE);
      mpfr_neg(B[2].re, B[2].re, RMODE);

      //mpfr_mul(B[2].im, c[1][0], B[1].im, RMODE);
      mpfr_mul(B[2].im, c[1][0].re, B[1].im, RMODE);
      mpfr_fma(B[2].im, c[1][0].im, B[1].re, B[2].im, RMODE);

      //mpfr_fma(B[2].im, c[1][1], B[0].im, B[2].im, RMODE);  // b2
      mpfr_fma(B[2].im, c[1][1].re, B[0].im, B[2].im, RMODE);  // b2
      mpfr_fma(B[2].im, c[1][1].im, B[0].re, B[2].im, RMODE);  // b2

      mpfr_fms(B[2].im, y[jj], B[0].re, B[2].im, RMODE);
      //------------------------------------------------
      //mpfr_mul(B[3].re, c[2][0], B[2].re, RMODE);
      mpfr_mul(B[3].re, c[2][0].im, B[2].im, RMODE);
      mpfr_fms(B[3].re, c[2][0].re, B[2].re, B[3].re, RMODE);

      //mpfr_fma(B[3].re, c[2][1], B[1].re, B[3].re, RMODE);
      mpfr_fms(B[3].re, c[2][1].im, B[1].im, B[3].re, RMODE);
      mpfr_fms(B[3].re, c[2][1].re, B[1].re, B[3].re, RMODE);

      //mpfr_fma(B[3].re, c[2][2], B[0].re, B[3].re, RMODE);
      mpfr_fms(B[3].re, c[2][2].im, B[0].im, B[3].re, RMODE);
      mpfr_fms(B[3].re, c[2][2].re, B[0].re, B[3].re, RMODE);

      mpfr_fma(B[3].re, y[jj], B[1].im, B[3].re, RMODE);
      mpfr_neg(B[3].re, B[3].re, RMODE);                    // b3
//      mpfr_mul(B[3].im, c[2][0], B[2].im, RMODE);
      mpfr_mul(B[3].im, c[2][0].re, B[2].im, RMODE);
      mpfr_fma(B[3].im, c[2][0].im, B[2].re, B[3].im, RMODE);

//      mpfr_fma(B[3].im, c[2][1], B[1].im, B[3].im, RMODE);
      mpfr_fma(B[3].im, c[2][1].re, B[1].im, B[3].im, RMODE);
      mpfr_fma(B[3].im, c[2][1].im, B[1].re, B[3].im, RMODE);

//      mpfr_fma(B[3].im, c[2][2], B[0].im, B[3].im, RMODE);
      mpfr_fma(B[3].im, c[2][2].re, B[0].im, B[3].im, RMODE);
      mpfr_fma(B[3].im, c[2][2].im, B[0].re, B[3].im, RMODE);

      mpfr_fms(B[3].im, y[jj], B[1].re, B[3].im, RMODE);
      //-----------------------------------------------
      //mpfr_mul(B[4].re, c[3][0], B[3].re, RMODE);
      mpfr_mul(B[4].re, c[3][0].im, B[3].im, RMODE);
      mpfr_fms(B[4].re, c[3][0].re, B[3].re, B[4].re, RMODE);

//      mpfr_fma(B[4].re, c[3][1], B[2].re, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[3][1].im, B[2].im, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[3][1].re, B[2].re, B[4].re, RMODE);

//      mpfr_fma(B[4].re, c[3][2], B[1].re, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[3][2].im, B[1].im, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[3][2].re, B[1].re, B[4].re, RMODE);

//      mpfr_fma(B[4].re, c[3][3], B[0].re, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[3][3].im, B[0].im, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[3][3].re, B[0].re, B[4].re, RMODE);

      mpfr_fma(B[4].re, y[jj], B[2].im, B[4].re, RMODE);
      mpfr_neg(B[4].re, B[4].re, RMODE);                    // b4


//      mpfr_mul(B[4].im, c[3][0], B[3].im, RMODE);
      mpfr_mul(B[4].im, c[3][0].re, B[3].im, RMODE);
      mpfr_fma(B[4].im, c[3][0].im, B[3].re, B[4].im, RMODE);

//      mpfr_fma(B[4].im, c[3][1], B[2].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[3][1].re, B[2].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[3][1].im, B[2].re, B[4].im, RMODE);

//      mpfr_fma(B[4].im, c[3][2], B[1].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[3][2].re, B[1].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[3][2].im, B[1].re, B[4].im, RMODE);

//      mpfr_fma(B[4].im, c[3][3], B[0].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[3][3].re, B[0].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[3][3].im, B[0].re, B[4].im, RMODE);

      mpfr_fms(B[4].im, y[jj], B[2].re, B[4].im, RMODE);
      for (int kk = 4; kk < 2*d; kk++) {
        mpfr_set(B[0].re, B[1].re, RMODE); mpfr_set(B[0].im, B[1].im, RMODE);
        mpfr_set(B[1].re, B[2].re, RMODE); mpfr_set(B[1].im, B[2].im, RMODE);
        mpfr_set(B[2].re, B[3].re, RMODE); mpfr_set(B[2].im, B[3].im, RMODE);
        mpfr_set(B[3].re, B[4].re, RMODE); mpfr_set(B[3].im, B[4].im, RMODE);
//        mpfr_mul(B[4].re, c[kk][0], B[3].re, RMODE);
        mpfr_mul(B[4].re, c[kk][0].im, B[3].im, RMODE);
        mpfr_fms(B[4].re, c[kk][0].re, B[3].re, B[4].re, RMODE);

//        mpfr_fma(B[4].re, c[kk][1], B[2].re, B[4].re, RMODE);
        mpfr_fms(B[4].re, c[kk][1].im, B[2].im, B[4].re, RMODE);
        mpfr_fms(B[4].re, c[kk][1].re, B[2].re, B[4].re, RMODE);

//        mpfr_fma(B[4].re, c[kk][2], B[1].re, B[4].re, RMODE);
        mpfr_fms(B[4].re, c[kk][2].im, B[1].im, B[4].re, RMODE);
        mpfr_fms(B[4].re, c[kk][2].re, B[1].re, B[4].re, RMODE);

//        mpfr_fma(B[4].re, c[kk][3], B[0].re, B[4].re, RMODE);
        mpfr_fms(B[4].re, c[kk][3].im, B[0].im, B[4].re, RMODE);
        mpfr_fms(B[4].re, c[kk][3].re, B[0].re, B[4].re, RMODE);

        mpfr_fma(B[4].re, y[jj], B[2].im, B[4].re, RMODE);
        mpfr_neg(B[4].re, B[4].re, RMODE);                    // b(k+1)
//        mpfr_mul(B[4].im, c[kk][0], B[3].im, RMODE);
        mpfr_mul(B[4].im, c[kk][0].re, B[3].im, RMODE);
        mpfr_fma(B[4].im, c[kk][0].im, B[3].re, B[4].im, RMODE);

//        mpfr_fma(B[4].im, c[kk][1], B[2].im, B[4].im, RMODE);
        mpfr_fma(B[4].im, c[kk][1].re, B[2].im, B[4].im, RMODE);
        mpfr_fma(B[4].im, c[kk][1].im, B[2].re, B[4].im, RMODE);

//        mpfr_fma(B[4].im, c[kk][2], B[1].im, B[4].im, RMODE);
        mpfr_fma(B[4].im, c[kk][2].re, B[1].im, B[4].im, RMODE);
        mpfr_fma(B[4].im, c[kk][2].im, B[1].re, B[4].im, RMODE);

//        mpfr_fma(B[4].im, c[kk][3], B[0].im, B[4].im, RMODE);
        mpfr_fma(B[4].im, c[kk][3].re, B[0].im, B[4].im, RMODE);
        mpfr_fma(B[4].im, c[kk][3].im, B[0].re, B[4].im, RMODE);

        mpfr_fms(B[4].im, y[jj], B[2].re, B[4].im, RMODE);	
      }
      mpfr_set(Q[jj].re, B[4].re, RMODE);
      mpfr_set(Q[jj].im, B[4].im, RMODE);
    }
  }
}

void find_P(){
  mpfr_set_ui(A[0].re, 0, RMODE);
  mpfr_set_ui(A[0].im, 0, RMODE);
  mpfr_set_si(A[1].re, -1, RMODE);
  mpfr_set_ui(A[1].im, 0, RMODE);
  if (d == 1) {
    for (int jj = 0; jj < N; jj++) {
//      mpfr_mul(P[jj].re, c[1][0], A[1].re, RMODE);
      mpfr_mul(P[jj].re, c[1][0].im, A[1].im, RMODE);
      mpfr_fms(P[jj].re, c[1][0].re, A[1].re, P[jj].re, RMODE);

//      mpfr_fma(P[jj].re, c[1][1], A[0].re, P[jj].re, RMODE);
      mpfr_fms(P[jj].re, c[1][1].im, A[0].im, P[jj].re, RMODE);
      mpfr_fms(P[jj].re, c[1][1].re, A[0].re, P[jj].re, RMODE);

      mpfr_fma(P[jj].re, y[jj], A[0].im, P[jj].re, RMODE);
      mpfr_neg(P[jj].re, P[jj].re, RMODE);                    // update P

//      mpfr_mul(P[jj].im, c[1][0], A[1].im, RMODE);
      mpfr_mul(P[jj].im, c[1][0].re, A[1].im, RMODE);
      mpfr_fma(P[jj].im, c[1][0].im, A[1].re, P[jj].im, RMODE);

//      mpfr_fma(P[jj].im, c[1][1], A[0].im, P[jj].im, RMODE);
      mpfr_fma(P[jj].im, c[1][1].re, A[0].im, P[jj].im, RMODE);
      mpfr_fma(P[jj].im, c[1][1].im, A[0].re, P[jj].im, RMODE);

      mpfr_fms(P[jj].im, y[jj], A[0].re, P[jj].im, RMODE);
    }
  } else if (d == 2) { 
    for (int jj = 0; jj < N; jj++) {
      mpfr_mul(A[2].re, c[1][0].im, A[1].im, RMODE);
      mpfr_fms(A[2].re, c[1][0].re, A[1].re, A[2].re, RMODE);
      mpfr_fms(A[2].re, c[1][1].im, A[0].im, A[2].re, RMODE);
      mpfr_fms(A[2].re, c[1][1].re, A[0].re, A[2].re, RMODE);
      mpfr_fma(A[2].re, y[jj], A[0].im, A[2].re, RMODE);
      mpfr_neg(A[2].re, A[2].re, RMODE);
      mpfr_mul(A[2].im, c[1][0].re, A[1].im, RMODE);
      mpfr_fma(A[2].im, c[1][0].im, A[1].re, A[2].im, RMODE);
      mpfr_fma(A[2].im, c[1][1].re, A[0].im, A[2].im, RMODE);  // a2
      mpfr_fma(A[2].im, c[1][1].im, A[0].re, A[2].im, RMODE);  // a2
      mpfr_fms(A[2].im, y[jj], A[0].re, A[2].im, RMODE);
      //------------------------------------------------

      mpfr_mul(A[3].re, c[2][0].im, A[2].im, RMODE);
      mpfr_fms(A[3].re, c[2][0].re, A[2].re, A[3].re, RMODE);
      mpfr_fms(A[3].re, c[2][1].im, A[1].im, A[3].re, RMODE);
      mpfr_fms(A[3].re, c[2][1].re, A[1].re, A[3].re, RMODE);
      mpfr_fms(A[3].re, c[2][2].im, A[0].im, A[3].re, RMODE);
      mpfr_fms(A[3].re, c[2][2].re, A[0].re, A[3].re, RMODE);
      mpfr_fma(A[3].re, y[jj], A[1].im, A[3].re, RMODE);
      mpfr_neg(A[3].re, A[3].re, RMODE);                    // a3
      mpfr_mul(A[3].im, c[2][0].re, A[2].im, RMODE);
      mpfr_fma(A[3].im, c[2][0].im, A[2].re, A[3].im, RMODE);
      mpfr_fma(A[3].im, c[2][1].re, A[1].im, A[3].im, RMODE);
      mpfr_fma(A[3].im, c[2][1].im, A[1].re, A[3].im, RMODE);
      mpfr_fma(A[3].im, c[2][2].re, A[0].im, A[3].im, RMODE);
      mpfr_fma(A[3].im, c[2][2].im, A[0].re, A[3].im, RMODE);
      mpfr_fms(A[3].im, y[jj], A[1].re, A[3].im, RMODE);
      //-----------------------------------------------
      mpfr_mul(P[jj].re, c[3][0].im, A[3].im, RMODE);
      mpfr_fms(P[jj].re, c[3][0].re, A[3].re, P[jj].re, RMODE);
      mpfr_fms(P[jj].re, c[3][1].im, A[2].im, P[jj].re, RMODE);
      mpfr_fms(P[jj].re, c[3][1].re, A[2].re, P[jj].re, RMODE);
      mpfr_fms(P[jj].re, c[3][2].im, A[1].im, P[jj].re, RMODE);
      mpfr_fms(P[jj].re, c[3][2].re, A[1].re, P[jj].re, RMODE);
      mpfr_fms(P[jj].re, c[3][3].im, A[0].im, P[jj].re, RMODE);
      mpfr_fms(P[jj].re, c[3][3].re, A[0].re, P[jj].re, RMODE);
      mpfr_fma(P[jj].re, y[jj], A[2].im, P[jj].re, RMODE);
      mpfr_neg(P[jj].re, P[jj].re, RMODE);                    // update P
      mpfr_mul(P[jj].im, c[3][0].re, A[3].im, RMODE);
      mpfr_fma(P[jj].im, c[3][0].im, A[3].re, P[jj].im, RMODE);
      mpfr_fma(P[jj].im, c[3][1].re, A[2].im, P[jj].im, RMODE);
      mpfr_fma(P[jj].im, c[3][1].im, A[2].re, P[jj].im, RMODE);
      mpfr_fma(P[jj].im, c[3][2].re, A[1].im, P[jj].im, RMODE);
      mpfr_fma(P[jj].im, c[3][2].im, A[1].re, P[jj].im, RMODE);
      mpfr_fma(P[jj].im, c[3][3].re, A[0].im, P[jj].im, RMODE);
      mpfr_fma(P[jj].im, c[3][3].im, A[0].re, P[jj].im, RMODE);
      mpfr_fms(P[jj].im, y[jj], A[2].re, P[jj].im, RMODE);
      /*if (jj == N/2) {
        mpfr_printf("A(0) = %.10Re + %.10ReI\tA(1) = %.10Re + %.10ReI\tA(2) = %.10Re + %.10ReI\tA(3) = %.10Re + %.10ReI\t\n\n", A[0].re, A[0].im, A[1].re, A[1].im, A[2].re, A[2].im, A[3].re, A[3].im);
        mpfr_printf("c(10) = %.10Re + %.10ReI\tc(11) = %.10Re + %.10ReI\n", c[1][0].re, c[1][0].im, c[1][1].re, c[1][1].im);
        mpfr_printf("c(20) = %.10Re + %.10ReI\tc(21) = %.10Re + %.10ReI\tc(22) = %.10Re + %.10ReI\n", c[2][0].re, c[2][0].im, c[2][1].re, c[2][1].im, c[2][2].re, c[2][2].im); 
        mpfr_printf("c(30) = %.10Re + %.10ReI\tc(31) = %.10Re + %.10ReI\tc(32) = %.10Re + %.10ReI\tc(33) = %.10Re + %.10ReI\n\n\n", c[3][0].re, c[3][0].im, c[3][1].re, c[3][1].im, c[3][2].re, c[3][2].im, c[3][3].re, c[3][3].im); 
        mpfr_printf("P(0) = %.10Re + %.10ReI\n\n\n", P[jj].re, P[jj].im);
        mpfr_printf("c33 = %.10Re + I%.10Re\ta0 =  %.10Re + I%.10Re\n", c[3][3].re, c[3][3].im, A[0].re, A[0].im);
        mpfr_printf("c32 = %.10Re + I%.10Re\ta1 =  %.10Re + I%.10Re\n", c[3][2].re, c[3][2].im, A[1].re, A[1].im);
        mpfr_printf("c31 = %.10Re + I%.10Re\ta2 =  %.10Re + I%.10Re\n", c[3][1].re, c[3][1].im, A[2].re, A[2].im);
        mpfr_printf("c30 = %.10Re + I%.10Re\ta3 =  %.10Re + I%.10Re\n", c[3][0].re, c[3][0].im, A[3].re, A[3].im);
        exit(1);
      }*/
    }
    //mpfr_printf("c%.10Re\n", *err); 
  } else {
    for (int jj = 0; jj < N; jj++) {
      mpfr_set_ui(A[0].re, 0, RMODE);
      mpfr_set_ui(A[0].im, 0, RMODE);
      mpfr_set_si(A[1].re, -1, RMODE);
      mpfr_set_ui(A[1].im, 0, RMODE);
      //-----------------------------------------------
//      mpfr_mul(A[2].re, c[1][0], A[1].re, RMODE);
      mpfr_mul(A[2].re, c[1][0].im, A[1].im, RMODE);
      mpfr_fms(A[2].re, c[1][0].re, A[1].re, A[2].re, RMODE);

//      mpfr_fma(A[2].re, c[1][1], A[0].re, A[2].re, RMODE);
      mpfr_fms(A[2].re, c[1][1].im, A[0].im, A[2].re, RMODE);
      mpfr_fms(A[2].re, c[1][1].re, A[0].re, A[2].re, RMODE);

      mpfr_fma(A[2].re, y[jj], A[0].im, A[2].re, RMODE);
      mpfr_neg(A[2].re, A[2].re, RMODE);

//      mpfr_mul(A[2].im, c[1][0], A[1].im, RMODE);
      mpfr_mul(A[2].im, c[1][0].re, A[1].im, RMODE);
      mpfr_fma(A[2].im, c[1][0].im, A[1].re, A[2].im, RMODE);

      //mpfr_fma(A[2].im, c[1][1], A[0].im, A[2].im, RMODE);  // a2
      mpfr_fma(A[2].im, c[1][1].re, A[0].im, A[2].im, RMODE);  // a2
      mpfr_fma(A[2].im, c[1][1].im, A[0].re, A[2].im, RMODE);  // a2

      mpfr_fms(A[2].im, y[jj], A[0].re, A[2].im, RMODE);
      //------------------------------------------------
//      mpfr_mul(A[3].re, c[2][0], A[2].re, RMODE);
      mpfr_mul(A[3].re, c[2][0].im, A[2].im, RMODE);
      mpfr_fms(A[3].re, c[2][0].re, A[2].re, A[3].re, RMODE);

//      mpfr_fma(A[3].re, c[2][1], A[1].re, A[3].re, RMODE);
      mpfr_fms(A[3].re, c[2][1].im, A[1].im, A[3].re, RMODE);
      mpfr_fms(A[3].re, c[2][1].re, A[1].re, A[3].re, RMODE);

//      mpfr_fma(A[3].re, c[2][2], A[0].re, A[3].re, RMODE);
      mpfr_fms(A[3].re, c[2][2].im, A[0].im, A[3].re, RMODE);
      mpfr_fms(A[3].re, c[2][2].re, A[0].re, A[3].re, RMODE);

      mpfr_fma(A[3].re, y[jj], A[1].im, A[3].re, RMODE);
      mpfr_neg(A[3].re, A[3].re, RMODE);                    // a3
//      mpfr_mul(A[3].im, c[2][0], A[2].im, RMODE);
      mpfr_mul(A[3].im, c[2][0].re, A[2].im, RMODE);
      mpfr_fma(A[3].im, c[2][0].im, A[2].re, A[3].im, RMODE);

//      mpfr_fma(A[3].im, c[2][1], A[1].im, A[3].im, RMODE);
      mpfr_fma(A[3].im, c[2][1].re, A[1].im, A[3].im, RMODE);
      mpfr_fma(A[3].im, c[2][1].im, A[1].re, A[3].im, RMODE);

//      mpfr_fma(A[3].im, c[2][2], A[0].im, A[3].im, RMODE);
      mpfr_fma(A[3].im, c[2][2].re, A[0].im, A[3].im, RMODE);
      mpfr_fma(A[3].im, c[2][2].im, A[0].re, A[3].im, RMODE);

      mpfr_fms(A[3].im, y[jj], A[1].re, A[3].im, RMODE);
      //-----------------------------------------------
//      mpfr_mul(A[4].re, c[3][0], A[3].re, RMODE);
      mpfr_mul(A[4].re, c[3][0].im, A[3].im, RMODE);
      mpfr_fms(A[4].re, c[3][0].re, A[3].re, A[4].re, RMODE);

//      mpfr_fma(A[4].re, c[3][1], A[2].re, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[3][1].im, A[2].im, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[3][1].re, A[2].re, A[4].re, RMODE);

//      mpfr_fma(A[4].re, c[3][2], A[1].re, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[3][2].im, A[1].im, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[3][2].re, A[1].re, A[4].re, RMODE);

//      mpfr_fma(A[4].re, c[3][3], A[0].re, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[3][3].im, A[0].im, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[3][3].re, A[0].re, A[4].re, RMODE);

      mpfr_fma(A[4].re, y[jj], A[2].im, A[4].re, RMODE);
      mpfr_neg(A[4].re, A[4].re, RMODE);                    // a4

//      mpfr_mul(A[4].im, c[3][0], A[3].im, RMODE);
      mpfr_mul(A[4].im, c[3][0].re, A[3].im, RMODE);
      mpfr_fma(A[4].im, c[3][0].im, A[3].re, A[4].im, RMODE);

//      mpfr_fma(A[4].im, c[3][1], A[2].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[3][1].re, A[2].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[3][1].im, A[2].re, A[4].im, RMODE);

//      mpfr_fma(A[4].im, c[3][2], A[1].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[3][2].re, A[1].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[3][2].im, A[1].re, A[4].im, RMODE);

//      mpfr_fma(A[4].im, c[3][3], A[0].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[3][3].re, A[0].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[3][3].im, A[0].re, A[4].im, RMODE);

      mpfr_fms(A[4].im, y[jj], A[2].re, A[4].im, RMODE);
      for (int kk = 4; kk < 2*d; kk++) {
        mpfr_set(A[0].re, A[1].re, RMODE); mpfr_set(A[0].im, A[1].im, RMODE);
        mpfr_set(A[1].re, A[2].re, RMODE); mpfr_set(A[1].im, A[2].im, RMODE);
        mpfr_set(A[2].re, A[3].re, RMODE); mpfr_set(A[2].im, A[3].im, RMODE);
        mpfr_set(A[3].re, A[4].re, RMODE); mpfr_set(A[3].im, A[4].im, RMODE);

//        mpfr_mul(A[4].re, c[kk][0], A[3].re, RMODE);
        mpfr_mul(A[4].re, c[kk][0].im, A[3].im, RMODE);
        mpfr_fms(A[4].re, c[kk][0].re, A[3].re, A[4].re, RMODE);

//        mpfr_fma(A[4].re, c[kk][1], A[2].re, A[4].re, RMODE);
        mpfr_fms(A[4].re, c[kk][1].im, A[2].im, A[4].re, RMODE);
        mpfr_fms(A[4].re, c[kk][1].re, A[2].re, A[4].re, RMODE);

//        mpfr_fma(A[4].re, c[kk][2], A[1].re, A[4].re, RMODE);
        mpfr_fms(A[4].re, c[kk][2].im, A[1].im, A[4].re, RMODE);
        mpfr_fms(A[4].re, c[kk][2].re, A[1].re, A[4].re, RMODE);

//        mpfr_fma(A[4].re, c[kk][3], A[0].re, A[4].re, RMODE);
        mpfr_fms(A[4].re, c[kk][3].im, A[0].im, A[4].re, RMODE);
        mpfr_fms(A[4].re, c[kk][3].re, A[0].re, A[4].re, RMODE);

        mpfr_fma(A[4].re, y[jj], A[2].im, A[4].re, RMODE);
        mpfr_neg(A[4].re, A[4].re, RMODE);                    // ak
//        mpfr_mul(A[4].im, c[kk][0], A[3].im, RMODE);
        mpfr_mul(A[4].im, c[kk][0].re, A[3].im, RMODE);
        mpfr_fma(A[4].im, c[kk][0].im, A[3].re, A[4].im, RMODE);

//        mpfr_fma(A[4].im, c[kk][1], A[2].im, A[4].im, RMODE);
        mpfr_fma(A[4].im, c[kk][1].re, A[2].im, A[4].im, RMODE);
        mpfr_fma(A[4].im, c[kk][1].im, A[2].re, A[4].im, RMODE);

//        mpfr_fma(A[4].im, c[kk][2], A[1].im, A[4].im, RMODE);
        mpfr_fma(A[4].im, c[kk][2].re, A[1].im, A[4].im, RMODE);
        mpfr_fma(A[4].im, c[kk][2].im, A[1].re, A[4].im, RMODE);

//        mpfr_fma(A[4].im, c[kk][3], A[0].im, A[4].im, RMODE);
        mpfr_fma(A[4].im, c[kk][3].re, A[0].im, A[4].im, RMODE);
        mpfr_fma(A[4].im, c[kk][3].im, A[0].re, A[4].im, RMODE);

        mpfr_fms(A[4].im, y[jj], A[2].re, A[4].im, RMODE);	
      }
      mpfr_set(P[jj].re, A[4].re, RMODE);
      mpfr_set(P[jj].im, A[4].im, RMODE);
    }
  }
}

void compute_approximation() {
  find_P();
  invert(overQ, Q);
  for (int jj = 0; jj < N; jj++) {
    mpfr_mul(X[jj].re, P[jj].im, overQ[jj].im, RMODE);
    mpfr_fms(X[jj].re, P[jj].re, overQ[jj].re, X[jj].re, RMODE);
    mpfr_mul(X[jj].im, P[jj].re, overQ[jj].im, RMODE);
    mpfr_fma(X[jj].im, P[jj].im, overQ[jj].re, X[jj].im, RMODE);
  } 
}

void polyevalN(mpfc_t *Pout, mpfc_t *Ppout, mpfc_t *Qout, mpfc_t *Qpout, mpfc_t *in, int NN){
  mpfr_set_ui(A[0].re, 0, RMODE);	  mpfr_set_zero(C[0].re, 1);
  mpfr_set_ui(A[0].im, 0, RMODE);	  mpfr_set_zero(C[0].im, 1);
  mpfr_set_si(A[1].re, -1, RMODE);	  mpfr_set_zero(C[1].re, 1);
  mpfr_set_ui(A[1].im, 0, RMODE);	  mpfr_set_zero(C[1].im, 1);
  mpfr_set_ui(B[0].re, 1, RMODE);	  mpfr_set_zero(D[0].re, 1);
  mpfr_set_ui(B[0].im, 0, RMODE);	  mpfr_set_zero(D[0].im, 1);
  /*mpfr_neg(B[1].re, c[0][0], RMODE);*/  mpfr_set_zero(D[1].re, 1);
  /*mpfr_set_ui(B[1].im, 0, RMODE);*/	  mpfr_set_zero(D[1].im, 1);
  mpfr_neg(B[1].re, c[0][0].re, RMODE);
  mpfr_neg(B[1].im, c[0][0].im, RMODE);
  if (d == 1) {
    for (int jj = 0; jj < NN; jj++) {
//     mpfr_mul(Pout[jj].re, c[1][0], A[1].re, RMODE);
     mpfr_mul(Pout[jj].re, c[1][0].im, A[1].im, RMODE);
     mpfr_fms(Pout[jj].re, c[1][0].re, A[1].re, Pout[jj].re, RMODE);

//     mpfr_fma(Pout[jj].re, c[1][1], A[0].re, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, c[1][1].im, A[0].im, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, c[1][1].re, A[0].re, Pout[jj].re, RMODE);

     mpfr_fma(Pout[jj].re, in[jj].im, A[0].im, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, in[jj].re, A[0].re, Pout[jj].re, RMODE);                    // find P(x + iy)
//     mpfr_mul(Pout[jj].im, c[1][0], A[1].im, RMODE);
     mpfr_mul(Pout[jj].im, c[1][0].re, A[1].im, RMODE);
     mpfr_fma(Pout[jj].im, c[1][0].im, A[1].re, Pout[jj].im, RMODE);

//     mpfr_fma(Pout[jj].im, c[1][1], A[0].im, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, c[1][1].re, A[0].im, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, c[1][1].im, A[0].re, Pout[jj].im, RMODE);

     mpfr_fms(Pout[jj].im, in[jj].im, A[0].re, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, in[jj].re, A[0].im, Pout[jj].im, RMODE);
     //-----------------------------------------------------
//     mpfr_mul(Qout[jj].re, c[1][0], B[1].re, RMODE);
     mpfr_mul(Qout[jj].re, c[1][0].im, B[1].im, RMODE);
     mpfr_fms(Qout[jj].re, c[1][0].re, B[1].re, Qout[jj].re, RMODE);

//     mpfr_fma(Qout[jj].re, c[1][1], B[0].re, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, c[1][1].im, B[0].im, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, c[1][1].re, B[0].re, Qout[jj].re, RMODE);

     mpfr_fma(Qout[jj].re, in[jj].im, B[0].im, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, in[jj].re, B[0].re, Qout[jj].re, RMODE);                    // find Q(x + iy)

//     mpfr_mul(Qout[jj].im, c[1][0], B[1].im, RMODE);
     mpfr_mul(Qout[jj].im, c[1][0].re, B[1].im, RMODE);
     mpfr_fma(Qout[jj].im, c[1][0].im, B[1].re, Qout[jj].im, RMODE);

//     mpfr_fma(Qout[jj].im, c[1][1], B[0].im, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, c[1][1].re, B[0].im, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, c[1][1].im, B[0].re, Qout[jj].im, RMODE);

     mpfr_fms(Qout[jj].im, in[jj].im, B[0].re, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, in[jj].re, B[0].im, Qout[jj].im, RMODE);
     //-----------------------------------------------------
//     mpfr_mul(Ppout[jj].re, c[1][0], C[1].re, RMODE);
     mpfr_mul(Ppout[jj].re, c[1][0].im, C[1].im, RMODE);
     mpfr_fms(Ppout[jj].re, c[1][0].re, C[1].re, Ppout[jj].re, RMODE);

//     mpfr_fma(Ppout[jj].re, c[1][1], C[0].re, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, c[1][1].im, C[0].im, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, c[1][1].re, C[0].re, Ppout[jj].re, RMODE);

     mpfr_fma(Ppout[jj].re, in[jj].im, C[0].im, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, in[jj].re, C[0].re, Ppout[jj].re, RMODE);                    // find P'(x + iy)
     mpfr_add(Ppout[jj].re, A[0].re, Ppout[jj].re, RMODE);
//     mpfr_mul(Ppout[jj].im, c[1][0], C[1].im, RMODE);
     mpfr_mul(Ppout[jj].im, c[1][0].re, C[1].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[1][0].im, C[1].re, Ppout[jj].im, RMODE);

//     mpfr_fma(Ppout[jj].im, c[1][1], C[0].im, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[1][1].re, C[0].im, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[1][1].im, C[0].re, Ppout[jj].im, RMODE);

     mpfr_fms(Ppout[jj].im, in[jj].im, C[0].re, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, in[jj].re, C[0].im, Ppout[jj].im, RMODE);
     mpfr_add(Ppout[jj].im, A[0].im, Ppout[jj].im, RMODE);
     //-----------------------------------------------------
//     mpfr_mul(Qpout[jj].re, c[1][0], D[1].re, RMODE);
     mpfr_mul(Qpout[jj].re, c[1][0].im, D[1].im, RMODE);
     mpfr_fms(Qpout[jj].re, c[1][0].re, D[1].re, Qpout[jj].re, RMODE);

//     mpfr_fma(Qpout[jj].re, c[1][1], D[0].re, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, c[1][1].im, D[0].im, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, c[1][1].re, D[0].re, Qpout[jj].re, RMODE);

     mpfr_fma(Qpout[jj].re, in[jj].im, D[0].im, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, in[jj].re, D[0].re, Qpout[jj].re, RMODE);                    // find Q'(x + iy)
     mpfr_add(Qpout[jj].re, B[0].re, Qpout[jj].re, RMODE);
//     mpfr_mul(Qpout[jj].im, c[1][0], D[1].im, RMODE);
     mpfr_mul(Qpout[jj].im, c[1][0].re, D[1].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[1][0].im, D[1].re, Qpout[jj].im, RMODE);

//     mpfr_fma(Qpout[jj].im, c[1][1], D[0].im, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[1][1].re, D[0].im, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[1][1].im, D[0].re, Qpout[jj].im, RMODE);

     mpfr_fms(Qpout[jj].im, in[jj].im, D[0].re, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, in[jj].re, D[0].im, Qpout[jj].im, RMODE);
     mpfr_add(Qpout[jj].im, B[0].im, Qpout[jj].im, RMODE);
    }
  } else if (d == 2) { 
    for (int jj = 0; jj < NN; jj++) {
//     mpfr_mul(A[2].re, c[1][0], A[1].re, RMODE);
     mpfr_mul(A[2].re, c[1][0].im, A[1].im, RMODE);
     mpfr_fms(A[2].re, c[1][0].re, A[1].re, A[2].re, RMODE);

//     mpfr_fma(A[2].re, c[1][1], A[0].re, A[2].re, RMODE);
     mpfr_fms(A[2].re, c[1][1].im, A[0].im, A[2].re, RMODE);
     mpfr_fms(A[2].re, c[1][1].re, A[0].re, A[2].re, RMODE);

     mpfr_fma(A[2].re, in[jj].im, A[0].im, A[2].re, RMODE);
     mpfr_fms(A[2].re, in[jj].re, A[0].re, A[2].re, RMODE);        
//     mpfr_mul(A[2].im, c[1][0], A[1].im, RMODE);
     mpfr_mul(A[2].im, c[1][0].re, A[1].im, RMODE);
     mpfr_fma(A[2].im, c[1][0].im, A[1].re, A[2].im, RMODE);

//     mpfr_fma(A[2].im, c[1][1], A[0].im, A[2].im, RMODE);  // a2
     mpfr_fma(A[2].im, c[1][1].re, A[0].im, A[2].im, RMODE); 
     mpfr_fma(A[2].im, c[1][1].im, A[0].re, A[2].im, RMODE);  

     mpfr_fms(A[2].im, in[jj].im, A[0].re, A[2].im, RMODE);
     mpfr_fma(A[2].im, in[jj].re, A[0].im, A[2].im, RMODE);
     //------------------------------------------------
//     mpfr_mul(B[2].re, c[1][0], B[1].re, RMODE);
     mpfr_mul(B[2].re, c[1][0].im, B[1].im, RMODE);
     mpfr_fms(B[2].re, c[1][0].re, B[1].re, B[2].re, RMODE);

//     mpfr_fma(B[2].re, c[1][1], B[0].re, B[2].re, RMODE);
     mpfr_fms(B[2].re, c[1][1].im, B[0].im, B[2].re, RMODE);
     mpfr_fms(B[2].re, c[1][1].re, B[0].re, B[2].re, RMODE);

     mpfr_fma(B[2].re, in[jj].im, B[0].im, B[2].re, RMODE);
     mpfr_fms(B[2].re, in[jj].re, B[0].re, B[2].re, RMODE);        
//     mpfr_mul(B[2].im, c[1][0], B[1].im, RMODE);
     mpfr_mul(B[2].im, c[1][0].re, B[1].im, RMODE);
     mpfr_fma(B[2].im, c[1][0].im, B[1].re, B[2].im, RMODE);

//     mpfr_fma(B[2].im, c[1][1], B[0].im, B[2].im, RMODE);  // b2
     mpfr_fma(B[2].im, c[1][1].re, B[0].im, B[2].im, RMODE);  
     mpfr_fma(B[2].im, c[1][1].im, B[0].re, B[2].im, RMODE);  

     mpfr_fms(B[2].im, in[jj].im, B[0].re, B[2].im, RMODE);
     mpfr_fma(B[2].im, in[jj].re, B[0].im, B[2].im, RMODE);
     //------------------------------------------------
//     mpfr_mul(C[2].re, c[1][0], C[1].re, RMODE);
     mpfr_mul(C[2].re, c[1][0].im, C[1].im, RMODE);
     mpfr_fms(C[2].re, c[1][0].re, C[1].re, C[2].re, RMODE);

//     mpfr_fma(C[2].re, c[1][1], C[0].re, C[2].re, RMODE);
     mpfr_fms(C[2].re, c[1][1].im, C[0].im, C[2].re, RMODE);
     mpfr_fms(C[2].re, c[1][1].re, C[0].re, C[2].re, RMODE);

     mpfr_fma(C[2].re, in[jj].im, C[0].im, C[2].re, RMODE);
     mpfr_fms(C[2].re, in[jj].re, C[0].re, C[2].re, RMODE);  // C2
     mpfr_add(C[2].re, A[0].re, C[2].re, RMODE);
//     mpfr_mul(C[2].im, c[1][0], C[1].im, RMODE);
     mpfr_mul(C[2].im, c[1][0].re, C[1].im, RMODE);
     mpfr_fma(C[2].im, c[1][0].im, C[1].re, C[2].im, RMODE);

//     mpfr_fma(C[2].im, c[1][1], C[0].im, C[2].im, RMODE);
     mpfr_fma(C[2].im, c[1][1].re, C[0].im, C[2].im, RMODE);
     mpfr_fma(C[2].im, c[1][1].im, C[0].re, C[2].im, RMODE);

     mpfr_fms(C[2].im, in[jj].im, C[0].re, C[2].im, RMODE);
     mpfr_fma(C[2].im, in[jj].re, C[0].im, C[2].im, RMODE);
     mpfr_add(C[2].im, A[0].im, C[2].im, RMODE);
     //-----------------------------------------------------
//     mpfr_mul(D[2].re, c[1][0], D[1].re, RMODE);
     mpfr_mul(D[2].re, c[1][0].im, D[1].im, RMODE);
     mpfr_fms(D[2].re, c[1][0].re, D[1].re, D[2].re, RMODE);

//     mpfr_fma(D[2].re, c[1][1], D[0].re, D[2].re, RMODE);
     mpfr_fms(D[2].re, c[1][1].im, D[0].im, D[2].re, RMODE);
     mpfr_fms(D[2].re, c[1][1].re, D[0].re, D[2].re, RMODE);

     mpfr_fma(D[2].re, in[jj].im, D[0].im, D[2].re, RMODE);
     mpfr_fms(D[2].re, in[jj].re, D[0].re, D[2].re, RMODE);   // D2
     mpfr_add(D[2].re, B[0].re, D[2].re, RMODE);
//     mpfr_mul(D[2].im, c[1][0], D[1].im, RMODE);
     mpfr_mul(D[2].im, c[1][0].re, D[1].im, RMODE);
     mpfr_fma(D[2].im, c[1][0].im, D[1].re, D[2].im, RMODE);

//     mpfr_fma(D[2].im, c[1][1], D[0].im, D[2].im, RMODE);
     mpfr_fma(D[2].im, c[1][1].re, D[0].im, D[2].im, RMODE);
     mpfr_fma(D[2].im, c[1][1].im, D[0].re, D[2].im, RMODE);

     mpfr_fms(D[2].im, in[jj].im, D[0].re, D[2].im, RMODE);
     mpfr_fma(D[2].im, in[jj].re, D[0].im, D[2].im, RMODE);
     mpfr_add(D[2].im, B[0].im, D[2].im, RMODE);
     //------------------------------------------------
     //------------------------------------------------
//     mpfr_mul(A[3].re, c[2][0], A[2].re, RMODE);
     mpfr_mul(A[3].re, c[2][0].im, A[2].im, RMODE);
     mpfr_fms(A[3].re, c[2][0].re, A[2].re, A[3].re, RMODE);

//     mpfr_fma(A[3].re, c[2][1], A[1].re, A[3].re, RMODE);
     mpfr_fms(A[3].re, c[2][1].im, A[1].im, A[3].re, RMODE);
     mpfr_fms(A[3].re, c[2][1].re, A[1].re, A[3].re, RMODE);

//     mpfr_fma(A[3].re, c[2][2], A[0].re, A[3].re, RMODE);
     mpfr_fms(A[3].re, c[2][2].im, A[0].im, A[3].re, RMODE);
     mpfr_fms(A[3].re, c[2][2].re, A[0].re, A[3].re, RMODE);

     mpfr_fma(A[3].re, in[jj].im, A[1].im, A[3].re, RMODE);
     mpfr_fms(A[3].re, in[jj].re, A[1].re, A[3].re, RMODE);                    // a3

//     mpfr_mul(A[3].im, c[2][0], A[2].im, RMODE);
     mpfr_mul(A[3].im, c[2][0].re, A[2].im, RMODE);
     mpfr_fma(A[3].im, c[2][0].im, A[2].re, A[3].im, RMODE);

//     mpfr_fma(A[3].im, c[2][1], A[1].im, A[3].im, RMODE);
     mpfr_fma(A[3].im, c[2][1].re, A[1].im, A[3].im, RMODE);
     mpfr_fma(A[3].im, c[2][1].im, A[1].re, A[3].im, RMODE);

//     mpfr_fma(A[3].im, c[2][2], A[0].im, A[3].im, RMODE);
     mpfr_fma(A[3].im, c[2][2].re, A[0].im, A[3].im, RMODE);
     mpfr_fma(A[3].im, c[2][2].im, A[0].re, A[3].im, RMODE);

     mpfr_fms(A[3].im, in[jj].im, A[1].re, A[3].im, RMODE);
     mpfr_fma(A[3].im, in[jj].re, A[1].im, A[3].im, RMODE);
     //------------------------------------------------
//     mpfr_mul(B[3].re, c[2][0], B[2].re, RMODE);
     mpfr_mul(B[3].re, c[2][0].im, B[2].im, RMODE);
     mpfr_fms(B[3].re, c[2][0].re, B[2].re, B[3].re, RMODE);

//     mpfr_fma(B[3].re, c[2][1], B[1].re, B[3].re, RMODE);
     mpfr_fms(B[3].re, c[2][1].im, B[1].im, B[3].re, RMODE);
     mpfr_fms(B[3].re, c[2][1].re, B[1].re, B[3].re, RMODE);

//     mpfr_fma(B[3].re, c[2][2], B[0].re, B[3].re, RMODE);
     mpfr_fms(B[3].re, c[2][2].im, B[0].im, B[3].re, RMODE);
     mpfr_fms(B[3].re, c[2][2].re, B[0].re, B[3].re, RMODE);

     mpfr_fma(B[3].re, in[jj].im, B[1].im, B[3].re, RMODE);
     mpfr_fms(B[3].re, in[jj].re, B[1].re, B[3].re, RMODE);                    // b3
//     mpfr_mul(B[3].im, c[2][0], B[2].im, RMODE);
     mpfr_mul(B[3].im, c[2][0].re, B[2].im, RMODE);
     mpfr_fma(B[3].im, c[2][0].im, B[2].re, B[3].im, RMODE);

//     mpfr_fma(B[3].im, c[2][1], B[1].im, B[3].im, RMODE);
     mpfr_fma(B[3].im, c[2][1].re, B[1].im, B[3].im, RMODE);
     mpfr_fma(B[3].im, c[2][1].im, B[1].re, B[3].im, RMODE);

//     mpfr_fma(B[3].im, c[2][2], B[0].im, B[3].im, RMODE);
     mpfr_fma(B[3].im, c[2][2].re, B[0].im, B[3].im, RMODE);
     mpfr_fma(B[3].im, c[2][2].im, B[0].re, B[3].im, RMODE);

     mpfr_fms(B[3].im, in[jj].im, B[1].re, B[3].im, RMODE);
     mpfr_fma(B[3].im, in[jj].re, B[1].im, B[3].im, RMODE);
     //------------------------------------------------
//     mpfr_mul(C[3].re, c[2][0], C[2].re, RMODE);
     mpfr_mul(C[3].re, c[2][0].im, C[2].im, RMODE);
     mpfr_fms(C[3].re, c[2][0].re, C[2].re, C[3].re, RMODE);

//     mpfr_fma(C[3].re, c[2][1], C[1].re, C[3].re, RMODE);
     mpfr_fms(C[3].re, c[2][1].im, C[1].im, C[3].re, RMODE);
     mpfr_fms(C[3].re, c[2][1].re, C[1].re, C[3].re, RMODE);

//     mpfr_fma(C[3].re, c[2][2], C[0].re, C[3].re, RMODE);
     mpfr_fms(C[3].re, c[2][2].im, C[0].im, C[3].re, RMODE);
     mpfr_fms(C[3].re, c[2][2].re, C[0].re, C[3].re, RMODE);

     mpfr_fma(C[3].re, in[jj].im, C[1].im, C[3].re, RMODE);
     mpfr_fms(C[3].re, in[jj].re, C[1].re, C[3].re, RMODE);                    // c3
     mpfr_add(C[3].re, A[1].re, C[3].re, RMODE);
//     mpfr_mul(C[3].im, c[2][0], C[2].im, RMODE);
     mpfr_mul(C[3].im, c[2][0].re, C[2].im, RMODE);
     mpfr_fma(C[3].im, c[2][0].im, C[2].re, C[3].im, RMODE);

//     mpfr_fma(C[3].im, c[2][1], C[1].im, C[3].im, RMODE);
     mpfr_fma(C[3].im, c[2][1].re, C[1].im, C[3].im, RMODE);
     mpfr_fma(C[3].im, c[2][1].im, C[1].re, C[3].im, RMODE);

//     mpfr_fma(C[3].im, c[2][2], C[0].im, C[3].im, RMODE);
     mpfr_fma(C[3].im, c[2][2].re, C[0].im, C[3].im, RMODE);
     mpfr_fma(C[3].im, c[2][2].im, C[0].re, C[3].im, RMODE);

     mpfr_fms(C[3].im, in[jj].im, C[1].re, C[3].im, RMODE);
     mpfr_fma(C[3].im, in[jj].re, C[1].im, C[3].im, RMODE);
     mpfr_add(C[3].im, A[1].im, C[3].im, RMODE);
     //------------------------------------------------
//     mpfr_mul(D[3].re, c[2][0], D[2].re, RMODE);
     mpfr_mul(D[3].re, c[2][0].im, D[2].im, RMODE);
     mpfr_fms(D[3].re, c[2][0].re, D[2].re, D[3].re, RMODE);

//     mpfr_fma(D[3].re, c[2][1], D[1].re, D[3].re, RMODE);
     mpfr_fms(D[3].re, c[2][1].im, D[1].im, D[3].re, RMODE);
     mpfr_fms(D[3].re, c[2][1].re, D[1].re, D[3].re, RMODE);

//     mpfr_fma(D[3].re, c[2][2], D[0].re, D[3].re, RMODE);
     mpfr_fms(D[3].re, c[2][2].im, D[0].im, D[3].re, RMODE);
     mpfr_fms(D[3].re, c[2][2].re, D[0].re, D[3].re, RMODE);

     mpfr_fma(D[3].re, in[jj].im, D[1].im, D[3].re, RMODE);
     mpfr_fms(D[3].re, in[jj].re, D[1].re, D[3].re, RMODE);                    // d3
     mpfr_add(D[3].re, B[1].re, D[3].re, RMODE);
//     mpfr_mul(D[3].im, c[2][0], D[2].im, RMODE);
     mpfr_mul(D[3].im, c[2][0].re, D[2].im, RMODE);
     mpfr_fma(D[3].im, c[2][0].im, D[2].re, D[3].im, RMODE);

//     mpfr_fma(D[3].im, c[2][1], D[1].im, D[3].im, RMODE);
     mpfr_fma(D[3].im, c[2][1].re, D[1].im, D[3].im, RMODE);
     mpfr_fma(D[3].im, c[2][1].im, D[1].re, D[3].im, RMODE);

//     mpfr_fma(D[3].im, c[2][2], D[0].im, D[3].im, RMODE);
     mpfr_fma(D[3].im, c[2][2].re, D[0].im, D[3].im, RMODE);
     mpfr_fma(D[3].im, c[2][2].im, D[0].re, D[3].im, RMODE);

     mpfr_fms(D[3].im, in[jj].im, D[1].re, D[3].im, RMODE);
     mpfr_fma(D[3].im, in[jj].re, D[1].im, D[3].im, RMODE);
     mpfr_add(D[3].im, B[1].im, D[3].im, RMODE);
     //------------------------------------------------
     //------------------------------------------------
//     mpfr_mul(Pout[jj].re, c[3][0], A[3].re, RMODE);
     mpfr_mul(Pout[jj].re, c[3][0].im, A[3].im, RMODE);
     mpfr_fms(Pout[jj].re, c[3][0].re, A[3].re, Pout[jj].re, RMODE);

//     mpfr_fma(Pout[jj].re, c[3][1], A[2].re, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, c[3][1].im, A[2].im, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, c[3][1].re, A[2].re, Pout[jj].re, RMODE);

//     mpfr_fma(Pout[jj].re, c[3][2], A[1].re, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, c[3][2].im, A[1].im, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, c[3][2].re, A[1].re, Pout[jj].re, RMODE);

//     mpfr_fma(Pout[jj].re, c[3][3], A[0].re, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, c[3][3].im, A[0].im, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, c[3][3].re, A[0].re, Pout[jj].re, RMODE);

     mpfr_fma(Pout[jj].re, in[jj].im, A[2].im, Pout[jj].re, RMODE);
     mpfr_fms(Pout[jj].re, in[jj].re, A[2].re, Pout[jj].re, RMODE);                    // P(x+iy)
//     mpfr_mul(Pout[jj].im, c[3][0], A[3].im, RMODE);
     mpfr_mul(Pout[jj].im, c[3][0].re, A[3].im, RMODE);
     mpfr_fma(Pout[jj].im, c[3][0].im, A[3].re, Pout[jj].im, RMODE);

//     mpfr_fma(Pout[jj].im, c[3][1], A[2].im, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, c[3][1].re, A[2].im, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, c[3][1].im, A[2].re, Pout[jj].im, RMODE);

//     mpfr_fma(Pout[jj].im, c[3][2], A[1].im, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, c[3][2].re, A[1].im, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, c[3][2].im, A[1].re, Pout[jj].im, RMODE);

//     mpfr_fma(Pout[jj].im, c[3][3], A[0].im, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, c[3][3].re, A[0].im, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, c[3][3].im, A[0].re, Pout[jj].im, RMODE);

     mpfr_fms(Pout[jj].im, in[jj].im, A[2].re, Pout[jj].im, RMODE);
     mpfr_fma(Pout[jj].im, in[jj].re, A[2].im, Pout[jj].im, RMODE);
     //-------------------------------------------------
//     mpfr_mul(Qout[jj].re, c[3][0], B[3].re, RMODE);
     mpfr_mul(Qout[jj].re, c[3][0].im, B[3].im, RMODE);
     mpfr_fms(Qout[jj].re, c[3][0].re, B[3].re, Qout[jj].re, RMODE);

//     mpfr_fma(Qout[jj].re, c[3][1], B[2].re, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, c[3][1].im, B[2].im, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, c[3][1].re, B[2].re, Qout[jj].re, RMODE);

//     mpfr_fma(Qout[jj].re, c[3][2], B[1].re, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, c[3][2].im, B[1].im, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, c[3][2].re, B[1].re, Qout[jj].re, RMODE);

//     mpfr_fma(Qout[jj].re, c[3][3], B[0].re, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, c[3][3].im, B[0].im, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, c[3][3].re, B[0].re, Qout[jj].re, RMODE);

     mpfr_fma(Qout[jj].re, in[jj].im, B[2].im, Qout[jj].re, RMODE);
     mpfr_fms(Qout[jj].re, in[jj].re, B[2].re, Qout[jj].re, RMODE);                    // Q(x+iy)
//     mpfr_mul(Qout[jj].im, c[3][0], B[3].im, RMODE);
     mpfr_mul(Qout[jj].im, c[3][0].re, B[3].im, RMODE);
     mpfr_fma(Qout[jj].im, c[3][0].im, B[3].re, Qout[jj].im, RMODE);

//     mpfr_fma(Qout[jj].im, c[3][1], B[2].im, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, c[3][1].re, B[2].im, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, c[3][1].im, B[2].re, Qout[jj].im, RMODE);

//     mpfr_fma(Qout[jj].im, c[3][2], B[1].im, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, c[3][2].re, B[1].im, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, c[3][2].im, B[1].re, Qout[jj].im, RMODE);

//     mpfr_fma(Qout[jj].im, c[3][3], B[0].im, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, c[3][3].re, B[0].im, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, c[3][3].im, B[0].re, Qout[jj].im, RMODE);

     mpfr_fms(Qout[jj].im, in[jj].im, B[2].re, Qout[jj].im, RMODE);
     mpfr_fma(Qout[jj].im, in[jj].re, B[2].im, Qout[jj].im, RMODE);
     //--------------------------------------------------
//     mpfr_mul(Ppout[jj].re, c[3][0], C[3].re, RMODE);
     mpfr_mul(Ppout[jj].re, c[3][0].im, C[3].im, RMODE);
     mpfr_fms(Ppout[jj].re, c[3][0].re, C[3].re, Ppout[jj].re, RMODE);

//     mpfr_fma(Ppout[jj].re, c[3][1], C[2].re, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, c[3][1].im, C[2].im, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, c[3][1].re, C[2].re, Ppout[jj].re, RMODE);

//     mpfr_fma(Ppout[jj].re, c[3][2], C[1].re, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, c[3][2].im, C[1].im, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, c[3][2].re, C[1].re, Ppout[jj].re, RMODE);

//     mpfr_fma(Ppout[jj].re, c[3][3], C[0].re, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, c[3][3].im, C[0].im, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, c[3][3].re, C[0].re, Ppout[jj].re, RMODE);

     mpfr_fma(Ppout[jj].re, in[jj].im, C[2].im, Ppout[jj].re, RMODE);
     mpfr_fms(Ppout[jj].re, in[jj].re, C[2].re, Ppout[jj].re, RMODE);                    // P'(x+iy)
     mpfr_add(Ppout[jj].re, A[2].re, Ppout[jj].re, RMODE);
//     mpfr_mul(Ppout[jj].im, c[3][0], C[3].im, RMODE);
     mpfr_mul(Ppout[jj].im, c[3][0].re, C[3].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[3][0].im, C[3].re, Ppout[jj].im, RMODE);

//     mpfr_fma(Ppout[jj].im, c[3][1], C[2].im, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[3][1].re, C[2].im, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[3][1].im, C[2].re, Ppout[jj].im, RMODE);

//     mpfr_fma(Ppout[jj].im, c[3][2], C[1].im, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[3][2].re, C[1].im, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[3][2].im, C[1].re, Ppout[jj].im, RMODE);

//     mpfr_fma(Ppout[jj].im, c[3][3], C[0].im, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[3][3].re, C[0].im, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, c[3][3].im, C[0].re, Ppout[jj].im, RMODE);

     mpfr_fms(Ppout[jj].im, in[jj].im, C[2].re, Ppout[jj].im, RMODE);
     mpfr_fma(Ppout[jj].im, in[jj].re, C[2].im, Ppout[jj].im, RMODE);
     mpfr_add(Ppout[jj].im, A[2].im, Ppout[jj].im, RMODE);
     //--------------------------------------------------
//     mpfr_mul(Qpout[jj].re, c[3][0], D[3].re, RMODE);
     mpfr_mul(Qpout[jj].re, c[3][0].im, D[3].im, RMODE);
     mpfr_fms(Qpout[jj].re, c[3][0].re, D[3].re, Qpout[jj].re, RMODE);

//     mpfr_fma(Qpout[jj].re, c[3][1], D[2].re, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, c[3][1].im, D[2].im, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, c[3][1].re, D[2].re, Qpout[jj].re, RMODE);

//     mpfr_fma(Qpout[jj].re, c[3][2], D[1].re, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, c[3][2].im, D[1].im, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, c[3][2].re, D[1].re, Qpout[jj].re, RMODE);

//     mpfr_fma(Qpout[jj].re, c[3][3], D[0].re, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, c[3][3].im, D[0].im, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, c[3][3].re, D[0].re, Qpout[jj].re, RMODE);

     mpfr_fma(Qpout[jj].re, in[jj].im, D[2].im, Qpout[jj].re, RMODE);
     mpfr_fms(Qpout[jj].re, in[jj].re, D[2].re, Qpout[jj].re, RMODE);                    // Q'(x+iy)
     mpfr_add(Qpout[jj].re, B[2].re, Qpout[jj].re, RMODE);
//     mpfr_mul(Qpout[jj].im, c[3][0], D[3].im, RMODE);
     mpfr_mul(Qpout[jj].im, c[3][0].re, D[3].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[3][0].im, D[3].re, Qpout[jj].im, RMODE);

//     mpfr_fma(Qpout[jj].im, c[3][1], D[2].im, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[3][1].re, D[2].im, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[3][1].im, D[2].re, Qpout[jj].im, RMODE);

//     mpfr_fma(Qpout[jj].im, c[3][2], D[1].im, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[3][2].re, D[1].im, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[3][2].im, D[1].re, Qpout[jj].im, RMODE);

//     mpfr_fma(Qpout[jj].im, c[3][3], D[0].im, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[3][3].re, D[0].im, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, c[3][3].im, D[0].re, Qpout[jj].im, RMODE);

     mpfr_fms(Qpout[jj].im, in[jj].im, D[2].re, Qpout[jj].im, RMODE);
     mpfr_fma(Qpout[jj].im, in[jj].re, D[2].im, Qpout[jj].im, RMODE);
     mpfr_add(Qpout[jj].im, B[2].im, Qpout[jj].im, RMODE);
    }
  } else {
   for (int jj = 0; jj < NN; jj++) {
    mpfr_set_ui(A[0].re, 0, RMODE);	  mpfr_set_zero(C[0].re, 1);
    mpfr_set_ui(A[0].im, 0, RMODE);	  mpfr_set_zero(C[0].im, 1);
    mpfr_set_si(A[1].re, -1, RMODE);	  mpfr_set_zero(C[1].re, 1);
    mpfr_set_ui(A[1].im, 0, RMODE);	  mpfr_set_zero(C[1].im, 1);
    mpfr_set_ui(B[0].re, 1, RMODE);	  mpfr_set_zero(D[0].re, 1);
    mpfr_set_ui(B[0].im, 0, RMODE);	  mpfr_set_zero(D[0].im, 1);
    /*mpfr_neg(B[1].re, c[0][0], RMODE);*/mpfr_set_zero(D[1].re, 1);
    /*mpfr_set_ui(B[1].im, 0, RMODE); */  mpfr_set_zero(D[1].im, 1);
    mpfr_neg(B[1].re, c[0][0].re, RMODE);
    mpfr_neg(B[1].im, c[0][0].im, RMODE);
    //-----------------------------------------------
//    mpfr_mul(A[2].re, c[1][0], A[1].re, RMODE);
    mpfr_mul(A[2].re, c[1][0].im, A[1].im, RMODE);
    mpfr_fms(A[2].re, c[1][0].re, A[1].re, A[2].re, RMODE);

//    mpfr_fma(A[2].re, c[1][1], A[0].re, A[2].re, RMODE);
    mpfr_fms(A[2].re, c[1][1].im, A[0].im, A[2].re, RMODE);
    mpfr_fms(A[2].re, c[1][1].re, A[0].re, A[2].re, RMODE);

    mpfr_fma(A[2].re, in[jj].im, A[0].im, A[2].re, RMODE);
    mpfr_fms(A[2].re, in[jj].re, A[0].re, A[2].re, RMODE);
//    mpfr_mul(A[2].im, c[1][0], A[1].im, RMODE);
    mpfr_mul(A[2].im, c[1][0].re, A[1].im, RMODE);
    mpfr_fma(A[2].im, c[1][0].im, A[1].re, A[2].im, RMODE);

//    mpfr_fma(A[2].im, c[1][1], A[0].im, A[2].im, RMODE);  // a2
    mpfr_fma(A[2].im, c[1][1].re, A[0].im, A[2].im, RMODE); 
    mpfr_fma(A[2].im, c[1][1].im, A[0].re, A[2].im, RMODE); 

    mpfr_fms(A[2].im, in[jj].im, A[0].re, A[2].im, RMODE);
    mpfr_fma(A[2].im, in[jj].re, A[0].im, A[2].im, RMODE);
    //------------------------------------------------
//    mpfr_mul(B[2].re, c[1][0], B[1].re, RMODE);
    mpfr_mul(B[2].re, c[1][0].im, B[1].im, RMODE);
    mpfr_fms(B[2].re, c[1][0].re, B[1].re, B[2].re, RMODE);

//    mpfr_fma(B[2].re, c[1][1], B[0].re, B[2].re, RMODE);
    mpfr_fms(B[2].re, c[1][1].im, B[0].im, B[2].re, RMODE);
    mpfr_fms(B[2].re, c[1][1].re, B[0].re, B[2].re, RMODE);

    mpfr_fma(B[2].re, in[jj].im, B[0].im, B[2].re, RMODE);
    mpfr_fms(B[2].re, in[jj].re, B[0].re, B[2].re, RMODE);        
//    mpfr_mul(B[2].im, c[1][0], B[1].im, RMODE);
    mpfr_mul(B[2].im, c[1][0].re, B[1].im, RMODE);
    mpfr_fma(B[2].im, c[1][0].im, B[1].re, B[2].im, RMODE);

//    mpfr_fma(B[2].im, c[1][1], B[0].im, B[2].im, RMODE);  // b2
    mpfr_fma(B[2].im, c[1][1].re, B[0].im, B[2].im, RMODE);  // b2
    mpfr_fma(B[2].im, c[1][1].im, B[0].re, B[2].im, RMODE);  // b2

    mpfr_fms(B[2].im, in[jj].im, B[0].re, B[2].im, RMODE);
    mpfr_fma(B[2].im, in[jj].re, B[0].im, B[2].im, RMODE);
    //------------------------------------------------
//    mpfr_mul(C[2].re, c[1][0], C[1].re, RMODE);
    mpfr_mul(C[2].re, c[1][0].im, C[1].im, RMODE);
    mpfr_fms(C[2].re, c[1][0].re, C[1].re, C[2].re, RMODE);

//    mpfr_fma(C[2].re, c[1][1], C[0].re, C[2].re, RMODE);
    mpfr_fms(C[2].re, c[1][1].im, C[0].im, C[2].re, RMODE);
    mpfr_fms(C[2].re, c[1][1].re, C[0].re, C[2].re, RMODE);

    mpfr_fma(C[2].re, in[jj].im, C[0].im, C[2].re, RMODE);
    mpfr_fms(C[2].re, in[jj].re, C[0].re, C[2].re, RMODE);  // C2
    mpfr_add(C[2].re, A[0].re, C[2].re, RMODE);
//    mpfr_mul(C[2].im, c[1][0], C[1].im, RMODE);
    mpfr_mul(C[2].im, c[1][0].re, C[1].im, RMODE);
    mpfr_fma(C[2].im, c[1][0].im, C[1].re, C[2].im, RMODE);

//    mpfr_fma(C[2].im, c[1][1], C[0].im, C[2].im, RMODE);
    mpfr_fma(C[2].im, c[1][1].re, C[0].im, C[2].im, RMODE);
    mpfr_fma(C[2].im, c[1][1].im, C[0].re, C[2].im, RMODE);

    mpfr_fms(C[2].im, in[jj].im, C[0].re, C[2].im, RMODE);
    mpfr_fma(C[2].im, in[jj].re, C[0].im, C[2].im, RMODE);
    mpfr_add(C[2].im, A[0].im, C[2].im, RMODE);
    //-----------------------------------------------------
//    mpfr_mul(D[2].re, c[1][0], D[1].re, RMODE);
    mpfr_mul(D[2].re, c[1][0].im, D[1].im, RMODE);
    mpfr_fms(D[2].re, c[1][0].re, D[1].re, D[2].re, RMODE);

//    mpfr_fma(D[2].re, c[1][1], D[0].re, D[2].re, RMODE);
    mpfr_fms(D[2].re, c[1][1].im, D[0].im, D[2].re, RMODE);
    mpfr_fms(D[2].re, c[1][1].re, D[0].re, D[2].re, RMODE);

    mpfr_fma(D[2].re, in[jj].im, D[0].im, D[2].re, RMODE);
    mpfr_fms(D[2].re, in[jj].re, D[0].re, D[2].re, RMODE);   // D2
    mpfr_add(D[2].re, B[0].re, D[2].re, RMODE);
//    mpfr_mul(D[2].im, c[1][0], D[1].im, RMODE);
    mpfr_mul(D[2].im, c[1][0].re, D[1].im, RMODE);
    mpfr_fma(D[2].im, c[1][0].im, D[1].re, D[2].im, RMODE);

//    mpfr_fma(D[2].im, c[1][1], D[0].im, D[2].im, RMODE);
    mpfr_fma(D[2].im, c[1][1].re, D[0].im, D[2].im, RMODE);
    mpfr_fma(D[2].im, c[1][1].im, D[0].re, D[2].im, RMODE);

    mpfr_fms(D[2].im, in[jj].im, D[0].re, D[2].im, RMODE);
    mpfr_fma(D[2].im, in[jj].re, D[0].im, D[2].im, RMODE);
    mpfr_add(D[2].im, B[0].im, D[2].im, RMODE);
    //------------------------------------------------
    //------------------------------------------------
//    mpfr_mul(A[3].re, c[2][0], A[2].re, RMODE);
    mpfr_mul(A[3].re, c[2][0].im, A[2].im, RMODE);
    mpfr_fms(A[3].re, c[2][0].re, A[2].re, A[3].re, RMODE);

//    mpfr_fma(A[3].re, c[2][1], A[1].re, A[3].re, RMODE);
    mpfr_fms(A[3].re, c[2][1].im, A[1].im, A[3].re, RMODE);
    mpfr_fms(A[3].re, c[2][1].re, A[1].re, A[3].re, RMODE);

//    mpfr_fma(A[3].re, c[2][2], A[0].re, A[3].re, RMODE);
    mpfr_fms(A[3].re, c[2][2].im, A[0].im, A[3].re, RMODE);
    mpfr_fms(A[3].re, c[2][2].re, A[0].re, A[3].re, RMODE);

    mpfr_fma(A[3].re, in[jj].im, A[1].im, A[3].re, RMODE);
    mpfr_fms(A[3].re, in[jj].re, A[1].re, A[3].re, RMODE);                    // a3
//    mpfr_mul(A[3].im, c[2][0], A[2].im, RMODE);
    mpfr_mul(A[3].im, c[2][0].re, A[2].im, RMODE);
    mpfr_fma(A[3].im, c[2][0].im, A[2].re, A[3].im, RMODE);

//    mpfr_fma(A[3].im, c[2][1], A[1].im, A[3].im, RMODE);
    mpfr_fma(A[3].im, c[2][1].re, A[1].im, A[3].im, RMODE);
    mpfr_fma(A[3].im, c[2][1].im, A[1].re, A[3].im, RMODE);

//    mpfr_fma(A[3].im, c[2][2], A[0].im, A[3].im, RMODE);
    mpfr_fma(A[3].im, c[2][2].re, A[0].im, A[3].im, RMODE);
    mpfr_fma(A[3].im, c[2][2].im, A[0].re, A[3].im, RMODE);

    mpfr_fms(A[3].im, in[jj].im, A[1].re, A[3].im, RMODE);
    mpfr_fma(A[3].im, in[jj].re, A[1].im, A[3].im, RMODE);
    //------------------------------------------------
//    mpfr_mul(B[3].re, c[2][0], B[2].re, RMODE);
    mpfr_mul(B[3].re, c[2][0].im, B[2].im, RMODE);
    mpfr_fms(B[3].re, c[2][0].re, B[2].re, B[3].re, RMODE);

//    mpfr_fma(B[3].re, c[2][1], B[1].re, B[3].re, RMODE);
    mpfr_fms(B[3].re, c[2][1].im, B[1].im, B[3].re, RMODE);
    mpfr_fms(B[3].re, c[2][1].re, B[1].re, B[3].re, RMODE);

//    mpfr_fma(B[3].re, c[2][2], B[0].re, B[3].re, RMODE);
    mpfr_fms(B[3].re, c[2][2].im, B[0].im, B[3].re, RMODE);
    mpfr_fms(B[3].re, c[2][2].re, B[0].re, B[3].re, RMODE);

    mpfr_fma(B[3].re, in[jj].im, B[1].im, B[3].re, RMODE);
    mpfr_fms(B[3].re, in[jj].re, B[1].re, B[3].re, RMODE);                    // b3

//    mpfr_mul(B[3].im, c[2][0], B[2].im, RMODE);
    mpfr_mul(B[3].im, c[2][0].re, B[2].im, RMODE);
    mpfr_fma(B[3].im, c[2][0].im, B[2].re, B[3].im, RMODE);

//    mpfr_fma(B[3].im, c[2][1], B[1].im, B[3].im, RMODE);
    mpfr_fma(B[3].im, c[2][1].re, B[1].im, B[3].im, RMODE);
    mpfr_fma(B[3].im, c[2][1].im, B[1].re, B[3].im, RMODE);

//    mpfr_fma(B[3].im, c[2][2], B[0].im, B[3].im, RMODE);
    mpfr_fma(B[3].im, c[2][2].re, B[0].im, B[3].im, RMODE);
    mpfr_fma(B[3].im, c[2][2].im, B[0].re, B[3].im, RMODE);

    mpfr_fms(B[3].im, in[jj].im, B[1].re, B[3].im, RMODE);
    mpfr_fma(B[3].im, in[jj].re, B[1].im, B[3].im, RMODE);
    //------------------------------------------------
//    mpfr_mul(C[3].re, c[2][0], C[2].re, RMODE);
    mpfr_mul(C[3].re, c[2][0].im, C[2].im, RMODE);
    mpfr_fms(C[3].re, c[2][0].re, C[2].re, C[3].re, RMODE);

//    mpfr_fma(C[3].re, c[2][1], C[1].re, C[3].re, RMODE);
    mpfr_fms(C[3].re, c[2][1].im, C[1].im, C[3].re, RMODE);
    mpfr_fms(C[3].re, c[2][1].re, C[1].re, C[3].re, RMODE);

//    mpfr_fma(C[3].re, c[2][2], C[0].re, C[3].re, RMODE);
    mpfr_fms(C[3].re, c[2][2].im, C[0].im, C[3].re, RMODE);
    mpfr_fms(C[3].re, c[2][2].re, C[0].re, C[3].re, RMODE);

    mpfr_fma(C[3].re, in[jj].im, C[1].im, C[3].re, RMODE);
    mpfr_fms(C[3].re, in[jj].re, C[1].re, C[3].re, RMODE);                    // c3
    mpfr_add(C[3].re, A[1].re, C[3].re, RMODE);
//    mpfr_mul(C[3].im, c[2][0], C[2].im, RMODE);
    mpfr_mul(C[3].im, c[2][0].re, C[2].im, RMODE);
    mpfr_fma(C[3].im, c[2][0].im, C[2].re, C[3].im, RMODE);

//    mpfr_fma(C[3].im, c[2][1], C[1].im, C[3].im, RMODE);
    mpfr_fma(C[3].im, c[2][1].re, C[1].im, C[3].im, RMODE);
    mpfr_fma(C[3].im, c[2][1].im, C[1].re, C[3].im, RMODE);

//    mpfr_fma(C[3].im, c[2][2], C[0].im, C[3].im, RMODE);
    mpfr_fma(C[3].im, c[2][2].re, C[0].im, C[3].im, RMODE);
    mpfr_fma(C[3].im, c[2][2].im, C[0].re, C[3].im, RMODE);

    mpfr_fms(C[3].im, in[jj].im, C[1].re, C[3].im, RMODE);
    mpfr_fma(C[3].im, in[jj].re, C[1].im, C[3].im, RMODE);
    mpfr_add(C[3].im, A[1].im, C[3].im, RMODE);
    //------------------------------------------------
//    mpfr_mul(D[3].re, c[2][0], D[2].re, RMODE);
    mpfr_mul(D[3].re, c[2][0].im, D[2].im, RMODE);
    mpfr_fms(D[3].re, c[2][0].re, D[2].re, D[3].re, RMODE);

//    mpfr_fma(D[3].re, c[2][1], D[1].re, D[3].re, RMODE);
    mpfr_fms(D[3].re, c[2][1].im, D[1].im, D[3].re, RMODE);
    mpfr_fms(D[3].re, c[2][1].re, D[1].re, D[3].re, RMODE);

//    mpfr_fma(D[3].re, c[2][2], D[0].re, D[3].re, RMODE);
    mpfr_fms(D[3].re, c[2][2].im, D[0].im, D[3].re, RMODE);
    mpfr_fms(D[3].re, c[2][2].re, D[0].re, D[3].re, RMODE);

    mpfr_fma(D[3].re, in[jj].im, D[1].im, D[3].re, RMODE);
    mpfr_fms(D[3].re, in[jj].re, D[1].re, D[3].re, RMODE);                    // d3
    mpfr_add(D[3].re, B[1].re, D[3].re, RMODE);
//    mpfr_mul(D[3].im, c[2][0], D[2].im, RMODE);
    mpfr_mul(D[3].im, c[2][0].re, D[2].im, RMODE);
    mpfr_fma(D[3].im, c[2][0].im, D[2].re, D[3].im, RMODE);

//    mpfr_fma(D[3].im, c[2][1], D[1].im, D[3].im, RMODE);
    mpfr_fma(D[3].im, c[2][1].re, D[1].im, D[3].im, RMODE);
    mpfr_fma(D[3].im, c[2][1].im, D[1].re, D[3].im, RMODE);

//    mpfr_fma(D[3].im, c[2][2], D[0].im, D[3].im, RMODE);
    mpfr_fma(D[3].im, c[2][2].re, D[0].im, D[3].im, RMODE);
    mpfr_fma(D[3].im, c[2][2].im, D[0].re, D[3].im, RMODE);

    mpfr_fms(D[3].im, in[jj].im, D[1].re, D[3].im, RMODE);
    mpfr_fma(D[3].im, in[jj].re, D[1].im, D[3].im, RMODE);
    mpfr_add(D[3].im, B[1].im, D[3].im, RMODE);
    //-----------------------------------------------
    //-----------------------------------------------
//    mpfr_mul(A[4].re, c[3][0], A[3].re, RMODE);
    mpfr_mul(A[4].re, c[3][0].im, A[3].im, RMODE);
    mpfr_fms(A[4].re, c[3][0].re, A[3].re, A[4].re, RMODE);

//    mpfr_fma(A[4].re, c[3][1], A[2].re, A[4].re, RMODE);
    mpfr_fms(A[4].re, c[3][1].im, A[2].im, A[4].re, RMODE);
    mpfr_fms(A[4].re, c[3][1].re, A[2].re, A[4].re, RMODE);

//    mpfr_fma(A[4].re, c[3][2], A[1].re, A[4].re, RMODE);
    mpfr_fms(A[4].re, c[3][2].im, A[1].im, A[4].re, RMODE);
    mpfr_fms(A[4].re, c[3][2].re, A[1].re, A[4].re, RMODE);

//    mpfr_fma(A[4].re, c[3][3], A[0].re, A[4].re, RMODE);
    mpfr_fms(A[4].re, c[3][3].im, A[0].im, A[4].re, RMODE);
    mpfr_fms(A[4].re, c[3][3].re, A[0].re, A[4].re, RMODE);

    mpfr_fma(A[4].re, in[jj].im, A[2].im, A[4].re, RMODE);
    mpfr_fms(A[4].re, in[jj].re, A[2].re, A[4].re, RMODE);                    // a4
//    mpfr_mul(A[4].im, c[3][0], A[3].im, RMODE);
    mpfr_mul(A[4].im, c[3][0].re, A[3].im, RMODE);
    mpfr_fma(A[4].im, c[3][0].im, A[3].re, A[4].im, RMODE);

//    mpfr_fma(A[4].im, c[3][1], A[2].im, A[4].im, RMODE);
    mpfr_fma(A[4].im, c[3][1].re, A[2].im, A[4].im, RMODE);
    mpfr_fma(A[4].im, c[3][1].im, A[2].re, A[4].im, RMODE);

//    mpfr_fma(A[4].im, c[3][2], A[1].im, A[4].im, RMODE);
    mpfr_fma(A[4].im, c[3][2].re, A[1].im, A[4].im, RMODE);
    mpfr_fma(A[4].im, c[3][2].im, A[1].re, A[4].im, RMODE);

//    mpfr_fma(A[4].im, c[3][3], A[0].im, A[4].im, RMODE);
    mpfr_fma(A[4].im, c[3][3].re, A[0].im, A[4].im, RMODE);
    mpfr_fma(A[4].im, c[3][3].im, A[0].re, A[4].im, RMODE);

    mpfr_fms(A[4].im, in[jj].im, A[2].re, A[4].im, RMODE);
    mpfr_fma(A[4].im, in[jj].re, A[2].im, A[4].im, RMODE);
    //--------------------------------------------------
//    mpfr_mul(B[4].re, c[3][0], B[3].re, RMODE);
    mpfr_mul(B[4].re, c[3][0].im, B[3].im, RMODE);
    mpfr_fms(B[4].re, c[3][0].re, B[3].re, B[4].re, RMODE);

//    mpfr_fma(B[4].re, c[3][1], B[2].re, B[4].re, RMODE);
    mpfr_fms(B[4].re, c[3][1].im, B[2].im, B[4].re, RMODE);
    mpfr_fms(B[4].re, c[3][1].re, B[2].re, B[4].re, RMODE);

//    mpfr_fma(B[4].re, c[3][2], B[1].re, B[4].re, RMODE);
    mpfr_fms(B[4].re, c[3][2].im, B[1].im, B[4].re, RMODE);
    mpfr_fms(B[4].re, c[3][2].re, B[1].re, B[4].re, RMODE);

//    mpfr_fma(B[4].re, c[3][3], B[0].re, B[4].re, RMODE);
    mpfr_fms(B[4].re, c[3][3].im, B[0].im, B[4].re, RMODE);
    mpfr_fms(B[4].re, c[3][3].re, B[0].re, B[4].re, RMODE);

    mpfr_fma(B[4].re, in[jj].im, B[2].im, B[4].re, RMODE);
    mpfr_fms(B[4].re, in[jj].re, B[2].re, B[4].re, RMODE);                    // b4
//    mpfr_mul(B[4].im, c[3][0], B[3].im, RMODE);
    mpfr_mul(B[4].im, c[3][0].re, B[3].im, RMODE);
    mpfr_fma(B[4].im, c[3][0].im, B[3].re, B[4].im, RMODE);

//    mpfr_fma(B[4].im, c[3][1], B[2].im, B[4].im, RMODE);
    mpfr_fma(B[4].im, c[3][1].re, B[2].im, B[4].im, RMODE);
    mpfr_fma(B[4].im, c[3][1].im, B[2].re, B[4].im, RMODE);

//    mpfr_fma(B[4].im, c[3][2], B[1].im, B[4].im, RMODE);
    mpfr_fma(B[4].im, c[3][2].re, B[1].im, B[4].im, RMODE);
    mpfr_fma(B[4].im, c[3][2].im, B[1].re, B[4].im, RMODE);

//    mpfr_fma(B[4].im, c[3][3], B[0].im, B[4].im, RMODE);
    mpfr_fma(B[4].im, c[3][3].re, B[0].im, B[4].im, RMODE);
    mpfr_fma(B[4].im, c[3][3].im, B[0].re, B[4].im, RMODE);

    mpfr_fms(B[4].im, in[jj].im, B[2].re, B[4].im, RMODE);
    mpfr_fma(B[4].im, in[jj].re, B[2].im, B[4].im, RMODE);
    //--------------------------------------------------
//    mpfr_mul(C[4].re, c[3][0], C[3].re, RMODE);
    mpfr_mul(C[4].re, c[3][0].im, C[3].im, RMODE);
    mpfr_fms(C[4].re, c[3][0].re, C[3].re, C[4].re, RMODE);

//    mpfr_fma(C[4].re, c[3][1], C[2].re, C[4].re, RMODE);
    mpfr_fms(C[4].re, c[3][1].im, C[2].im, C[4].re, RMODE);
    mpfr_fms(C[4].re, c[3][1].re, C[2].re, C[4].re, RMODE);

//    mpfr_fma(C[4].re, c[3][2], C[1].re, C[4].re, RMODE);
    mpfr_fms(C[4].re, c[3][2].im, C[1].im, C[4].re, RMODE);
    mpfr_fms(C[4].re, c[3][2].re, C[1].re, C[4].re, RMODE);

//    mpfr_fma(C[4].re, c[3][3], C[0].re, C[4].re, RMODE);
    mpfr_fms(C[4].re, c[3][3].im, C[0].im, C[4].re, RMODE);
    mpfr_fms(C[4].re, c[3][3].re, C[0].re, C[4].re, RMODE);

    mpfr_fma(C[4].re, in[jj].im, C[2].im, C[4].re, RMODE);
    mpfr_fms(C[4].re, in[jj].re, C[2].re, C[4].re, RMODE);                    // C4
    mpfr_add(C[4].re, A[2].re, C[4].re, RMODE);
//    mpfr_mul(C[4].im, c[3][0], C[3].im, RMODE);
    mpfr_mul(C[4].im, c[3][0].re, C[3].im, RMODE);
    mpfr_fma(C[4].im, c[3][0].im, C[3].re, C[4].im, RMODE);

//    mpfr_fma(C[4].im, c[3][1], C[2].im, C[4].im, RMODE);
    mpfr_fma(C[4].im, c[3][1].re, C[2].im, C[4].im, RMODE);
    mpfr_fma(C[4].im, c[3][1].im, C[2].re, C[4].im, RMODE);

//    mpfr_fma(C[4].im, c[3][2], C[1].im, C[4].im, RMODE);
    mpfr_fma(C[4].im, c[3][2].re, C[1].im, C[4].im, RMODE);
    mpfr_fma(C[4].im, c[3][2].im, C[1].re, C[4].im, RMODE);

//    mpfr_fma(C[4].im, c[3][3], C[0].im, C[4].im, RMODE);
    mpfr_fma(C[4].im, c[3][3].re, C[0].im, C[4].im, RMODE);
    mpfr_fma(C[4].im, c[3][3].im, C[0].re, C[4].im, RMODE);

    mpfr_fms(C[4].im, in[jj].im, C[2].re, C[4].im, RMODE);
    mpfr_fma(C[4].im, in[jj].re, C[2].im, C[4].im, RMODE);
    mpfr_add(C[4].im, A[2].im, C[4].im, RMODE);
    //--------------------------------------------------
//    mpfr_mul(D[4].re, c[3][0], D[3].re, RMODE);
    mpfr_mul(D[4].re, c[3][0].im, D[3].im, RMODE);
    mpfr_fms(D[4].re, c[3][0].re, D[3].re, D[4].re, RMODE);

//    mpfr_fma(D[4].re, c[3][1], D[2].re, D[4].re, RMODE);
    mpfr_fms(D[4].re, c[3][1].im, D[2].im, D[4].re, RMODE);
    mpfr_fms(D[4].re, c[3][1].re, D[2].re, D[4].re, RMODE);

//    mpfr_fma(D[4].re, c[3][2], D[1].re, D[4].re, RMODE);
    mpfr_fms(D[4].re, c[3][2].im, D[1].im, D[4].re, RMODE);
    mpfr_fms(D[4].re, c[3][2].re, D[1].re, D[4].re, RMODE);

//    mpfr_fma(D[4].re, c[3][3], D[0].re, D[4].re, RMODE);
    mpfr_fms(D[4].re, c[3][3].im, D[0].im, D[4].re, RMODE);
    mpfr_fms(D[4].re, c[3][3].re, D[0].re, D[4].re, RMODE);

    mpfr_fma(D[4].re, in[jj].im, D[2].im, D[4].re, RMODE);
    mpfr_fms(D[4].re, in[jj].re, D[2].re, D[4].re, RMODE);                    // D4
    mpfr_add(D[4].re, B[2].re, D[4].re, RMODE);
//    mpfr_mul(D[4].im, c[3][0], D[3].im, RMODE);
    mpfr_mul(D[4].im, c[3][0].re, D[3].im, RMODE);
    mpfr_fma(D[4].im, c[3][0].im, D[3].re, D[4].im, RMODE);

//    mpfr_fma(D[4].im, c[3][1], D[2].im, D[4].im, RMODE);
    mpfr_fma(D[4].im, c[3][1].re, D[2].im, D[4].im, RMODE);
    mpfr_fma(D[4].im, c[3][1].im, D[2].re, D[4].im, RMODE);

//    mpfr_fma(D[4].im, c[3][2], D[1].im, D[4].im, RMODE);
    mpfr_fma(D[4].im, c[3][2].re, D[1].im, D[4].im, RMODE);
    mpfr_fma(D[4].im, c[3][2].im, D[1].re, D[4].im, RMODE);

//    mpfr_fma(D[4].im, c[3][3], D[0].im, D[4].im, RMODE);
    mpfr_fma(D[4].im, c[3][3].re, D[0].im, D[4].im, RMODE);
    mpfr_fma(D[4].im, c[3][3].im, D[0].re, D[4].im, RMODE);

    mpfr_fms(D[4].im, in[jj].im, D[2].re, D[4].im, RMODE);
    mpfr_fma(D[4].im, in[jj].re, D[2].im, D[4].im, RMODE);
    mpfr_add(D[4].im, B[2].im, D[4].im, RMODE);
    //------------------------------------------------
    //------------------------------------------------
    for (int kk = 4; kk < 2*d; kk++) {
      mpfr_set(A[0].re, A[1].re, RMODE); mpfr_set(A[0].im, A[1].im, RMODE);
      mpfr_set(A[1].re, A[2].re, RMODE); mpfr_set(A[1].im, A[2].im, RMODE);
      mpfr_set(A[2].re, A[3].re, RMODE); mpfr_set(A[2].im, A[3].im, RMODE);
      mpfr_set(A[3].re, A[4].re, RMODE); mpfr_set(A[3].im, A[4].im, RMODE);
      //-------------------------------------------------------------------
      mpfr_set(B[0].re, B[1].re, RMODE); mpfr_set(B[0].im, B[1].im, RMODE);
      mpfr_set(B[1].re, B[2].re, RMODE); mpfr_set(B[1].im, B[2].im, RMODE);
      mpfr_set(B[2].re, B[3].re, RMODE); mpfr_set(B[2].im, B[3].im, RMODE);
      mpfr_set(B[3].re, B[4].re, RMODE); mpfr_set(B[3].im, B[4].im, RMODE);
      //-------------------------------------------------------------------
      mpfr_set(C[0].re, C[1].re, RMODE); mpfr_set(C[0].im, C[1].im, RMODE);
      mpfr_set(C[1].re, C[2].re, RMODE); mpfr_set(C[1].im, C[2].im, RMODE);
      mpfr_set(C[2].re, C[3].re, RMODE); mpfr_set(C[2].im, C[3].im, RMODE);
      mpfr_set(C[3].re, C[4].re, RMODE); mpfr_set(C[3].im, C[4].im, RMODE);
      //-------------------------------------------------------------------
      mpfr_set(D[0].re, D[1].re, RMODE); mpfr_set(D[0].im, D[1].im, RMODE);
      mpfr_set(D[1].re, D[2].re, RMODE); mpfr_set(D[1].im, D[2].im, RMODE);
      mpfr_set(D[2].re, D[3].re, RMODE); mpfr_set(D[2].im, D[3].im, RMODE);
      mpfr_set(D[3].re, D[4].re, RMODE); mpfr_set(D[3].im, D[4].im, RMODE);
      //-------------------------------------------------------------------
      //-------------------------------------------------------------------
//      mpfr_mul(A[4].re, c[kk][0], A[3].re, RMODE);
      mpfr_mul(A[4].re, c[kk][0].im, A[3].im, RMODE);
      mpfr_fms(A[4].re, c[kk][0].re, A[3].re, A[4].re, RMODE);

//      mpfr_fma(A[4].re, c[kk][1], A[2].re, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[kk][1].im, A[2].im, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[kk][1].re, A[2].re, A[4].re, RMODE);

//      mpfr_fma(A[4].re, c[kk][2], A[1].re, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[kk][2].im, A[1].im, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[kk][2].re, A[1].re, A[4].re, RMODE);

//      mpfr_fma(A[4].re, c[kk][3], A[0].re, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[kk][3].im, A[0].im, A[4].re, RMODE);
      mpfr_fms(A[4].re, c[kk][3].re, A[0].re, A[4].re, RMODE);

      mpfr_fma(A[4].re, in[jj].im, A[2].im, A[4].re, RMODE);
      mpfr_fms(A[4].re, in[jj].re, A[2].re, A[4].re, RMODE);                    // ak
//      mpfr_mul(A[4].im, c[kk][0], A[3].im, RMODE);
      mpfr_mul(A[4].im, c[kk][0].re, A[3].im, RMODE);
      mpfr_fma(A[4].im, c[kk][0].im, A[3].re, A[4].im, RMODE);

//      mpfr_fma(A[4].im, c[kk][1], A[2].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[kk][1].re, A[2].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[kk][1].im, A[2].re, A[4].im, RMODE);

//      mpfr_fma(A[4].im, c[kk][2], A[1].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[kk][2].re, A[1].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[kk][2].im, A[1].re, A[4].im, RMODE);

//      mpfr_fma(A[4].im, c[kk][3], A[0].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[kk][3].re, A[0].im, A[4].im, RMODE);
      mpfr_fma(A[4].im, c[kk][3].im, A[0].re, A[4].im, RMODE);

      mpfr_fms(A[4].im, in[jj].im, A[2].re, A[4].im, RMODE);	
      mpfr_fma(A[4].im, in[jj].re, A[2].im, A[4].im, RMODE);
      //--------------------------------------------------
//      mpfr_mul(B[4].re, c[kk][0], B[3].re, RMODE);
      mpfr_mul(B[4].re, c[kk][0].im, B[3].im, RMODE);
      mpfr_fms(B[4].re, c[kk][0].re, B[3].re, B[4].re, RMODE);

//      mpfr_fma(B[4].re, c[kk][1], B[2].re, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[kk][1].im, B[2].im, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[kk][1].re, B[2].re, B[4].re, RMODE);

//      mpfr_fma(B[4].re, c[kk][2], B[1].re, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[kk][2].im, B[1].im, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[kk][2].re, B[1].re, B[4].re, RMODE);

//      mpfr_fma(B[4].re, c[kk][3], B[0].re, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[kk][3].im, B[0].im, B[4].re, RMODE);
      mpfr_fms(B[4].re, c[kk][3].re, B[0].re, B[4].re, RMODE);

      mpfr_fma(B[4].re, in[jj].im, B[2].im, B[4].re, RMODE);
      mpfr_fms(B[4].re, in[jj].re, B[2].re, B[4].re, RMODE);                    // bk
//      mpfr_mul(B[4].im, c[kk][0], B[3].im, RMODE);
      mpfr_mul(B[4].im, c[kk][0].re, B[3].im, RMODE);
      mpfr_fma(B[4].im, c[kk][0].im, B[3].re, B[4].im, RMODE);

//      mpfr_fma(B[4].im, c[kk][1], B[2].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[kk][1].re, B[2].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[kk][1].im, B[2].re, B[4].im, RMODE);

//      mpfr_fma(B[4].im, c[kk][2], B[1].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[kk][2].re, B[1].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[kk][2].im, B[1].re, B[4].im, RMODE);

//      mpfr_fma(B[4].im, c[kk][3], B[0].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[kk][3].re, B[0].im, B[4].im, RMODE);
      mpfr_fma(B[4].im, c[kk][3].im, B[0].re, B[4].im, RMODE);

      mpfr_fms(B[4].im, in[jj].im, B[2].re, B[4].im, RMODE);
      mpfr_fma(B[4].im, in[jj].re, B[2].im, B[4].im, RMODE);
      //--------------------------------------------------
//      mpfr_mul(C[4].re, c[kk][0], C[3].re, RMODE);
      mpfr_mul(C[4].re, c[kk][0].im, C[3].im, RMODE);
      mpfr_fms(C[4].re, c[kk][0].re, C[3].re, C[4].re, RMODE);

//      mpfr_fma(C[4].re, c[kk][1], C[2].re, C[4].re, RMODE);
      mpfr_fms(C[4].re, c[kk][1].im, C[2].im, C[4].re, RMODE);
      mpfr_fms(C[4].re, c[kk][1].re, C[2].re, C[4].re, RMODE);

//      mpfr_fma(C[4].re, c[kk][2], C[1].re, C[4].re, RMODE);
      mpfr_fms(C[4].re, c[kk][2].im, C[1].im, C[4].re, RMODE);
      mpfr_fms(C[4].re, c[kk][2].re, C[1].re, C[4].re, RMODE);

//      mpfr_fma(C[4].re, c[kk][3], C[0].re, C[4].re, RMODE);
      mpfr_fms(C[4].re, c[kk][3].im, C[0].im, C[4].re, RMODE);
      mpfr_fms(C[4].re, c[kk][3].re, C[0].re, C[4].re, RMODE);

      mpfr_fma(C[4].re, in[jj].im, C[2].im, C[4].re, RMODE);
      mpfr_fms(C[4].re, in[jj].re, C[2].re, C[4].re, RMODE);                    // Ck
      mpfr_add(C[4].re, A[2].re, C[4].re, RMODE);
//      mpfr_mul(C[4].im, c[kk][0], C[3].im, RMODE);
      mpfr_mul(C[4].im, c[kk][0].re, C[3].im, RMODE);
      mpfr_fma(C[4].im, c[kk][0].im, C[3].re, C[4].im, RMODE);

//      mpfr_fma(C[4].im, c[kk][1], C[2].im, C[4].im, RMODE);
      mpfr_fma(C[4].im, c[kk][1].re, C[2].im, C[4].im, RMODE);
      mpfr_fma(C[4].im, c[kk][1].im, C[2].re, C[4].im, RMODE);

//      mpfr_fma(C[4].im, c[kk][2], C[1].im, C[4].im, RMODE);
      mpfr_fma(C[4].im, c[kk][2].re, C[1].im, C[4].im, RMODE);
      mpfr_fma(C[4].im, c[kk][2].im, C[1].re, C[4].im, RMODE);

//      mpfr_fma(C[4].im, c[kk][3], C[0].im, C[4].im, RMODE);
      mpfr_fma(C[4].im, c[kk][3].re, C[0].im, C[4].im, RMODE);
      mpfr_fma(C[4].im, c[kk][3].im, C[0].re, C[4].im, RMODE);

      mpfr_fms(C[4].im, in[jj].im, C[2].re, C[4].im, RMODE);
      mpfr_fma(C[4].im, in[jj].re, C[2].im, C[4].im, RMODE);
      mpfr_add(C[4].im, A[2].im, C[4].im, RMODE);
      //--------------------------------------------------
//      mpfr_mul(D[4].re, c[kk][0], D[3].re, RMODE);
      mpfr_mul(D[4].re, c[kk][0].im, D[3].im, RMODE);
      mpfr_fms(D[4].re, c[kk][0].re, D[3].re, D[4].re, RMODE);

//      mpfr_fma(D[4].re, c[kk][1], D[2].re, D[4].re, RMODE);
      mpfr_fms(D[4].re, c[kk][1].im, D[2].im, D[4].re, RMODE);
      mpfr_fms(D[4].re, c[kk][1].re, D[2].re, D[4].re, RMODE);

//      mpfr_fma(D[4].re, c[kk][2], D[1].re, D[4].re, RMODE);
      mpfr_fms(D[4].re, c[kk][2].im, D[1].im, D[4].re, RMODE);
      mpfr_fms(D[4].re, c[kk][2].re, D[1].re, D[4].re, RMODE);

//      mpfr_fma(D[4].re, c[kk][3], D[0].re, D[4].re, RMODE);
      mpfr_fms(D[4].re, c[kk][3].im, D[0].im, D[4].re, RMODE);
      mpfr_fms(D[4].re, c[kk][3].re, D[0].re, D[4].re, RMODE);

      mpfr_fma(D[4].re, in[jj].im, D[2].im, D[4].re, RMODE);
      mpfr_fms(D[4].re, in[jj].re, D[2].re, D[4].re, RMODE);                    // Dk
      mpfr_add(D[4].re, B[2].re, D[4].re, RMODE);
//      mpfr_mul(D[4].im, c[kk][0], D[3].im, RMODE);
      mpfr_mul(D[4].im, c[kk][0].re, D[3].im, RMODE);
      mpfr_fma(D[4].im, c[kk][0].im, D[3].re, D[4].im, RMODE);

//      mpfr_fma(D[4].im, c[kk][1], D[2].im, D[4].im, RMODE);
      mpfr_fma(D[4].im, c[kk][1].re, D[2].im, D[4].im, RMODE);
      mpfr_fma(D[4].im, c[kk][1].im, D[2].re, D[4].im, RMODE);

//      mpfr_fma(D[4].im, c[kk][2], D[1].im, D[4].im, RMODE);
      mpfr_fma(D[4].im, c[kk][2].re, D[1].im, D[4].im, RMODE);
      mpfr_fma(D[4].im, c[kk][2].im, D[1].re, D[4].im, RMODE);

//      mpfr_fma(D[4].im, c[kk][3], D[0].im, D[4].im, RMODE);
      mpfr_fma(D[4].im, c[kk][3].re, D[0].im, D[4].im, RMODE);
      mpfr_fma(D[4].im, c[kk][3].im, D[0].re, D[4].im, RMODE);

      mpfr_fms(D[4].im, in[jj].im, D[2].re, D[4].im, RMODE);
      mpfr_fma(D[4].im, in[jj].re, D[2].im, D[4].im, RMODE);
      mpfr_add(D[4].im, B[2].im, D[4].im, RMODE);
      //--------------------------------------------------------------------	
      //--------------------------------------------------------------------
    }
    mpfr_set(Pout[jj].re, A[4].re, RMODE);    mpfr_set(Pout[jj].im, A[4].im, RMODE);
    mpfr_set(Qout[jj].re, B[4].re, RMODE);    mpfr_set(Qout[jj].im, B[4].im, RMODE);
    mpfr_set(Ppout[jj].re, C[4].re, RMODE);   mpfr_set(Ppout[jj].im, C[4].im, RMODE);
    mpfr_set(Qpout[jj].re, D[4].re, RMODE);   mpfr_set(Qpout[jj].im, D[4].im, RMODE);
   }
  }
}



void finish_program(){ //work on it last
 for (int j = 0; j < N; j++) {
   mpfr_clear(u[j]);
   mpfr_clear(y[j]);
 }
}

void write_real(mpfr_t *out, char* str){
  FILE *fout = fopen(str,"w");
  for (int jj = 0; jj < N; jj++) mpfr_fprintf(fout, "%.32Re\t%.32Re\n", u[jj], out[jj]);
  fclose(fout);
}

void write_cmplx(mpfc_t *out, char* str){
  FILE *fout = fopen(str,"w");
  for (int jj = 0; jj < N; jj++) {
    mpfr_fprintf(fout, "%.300Re\t%.300Re\t%.300Re\n", u[jj], out[jj].re, out[jj].im);
  }
  fclose(fout);
}

void dotpc(mpfc_t *out, mpfc_t *m1, mpfc_t *m2) {
  mpfr_set_ui(out->re, 0, RMODE);
  mpfr_set_ui(out->im, 0, RMODE);
  for (int jj = 1; jj < N; jj++ ) {
    mpfr_mul(stk1, m1[jj].re, m2[jj].re, RMODE);
    mpfr_fma(stk1, m1[jj].im, m2[jj].im, stk1, RMODE);
    mpfr_fma(out->re, M[jj], stk1, out->re, RMODE);
    mpfr_mul(stk1, m1[jj].re, m2[jj].im, RMODE);
    mpfr_fms(stk1, m1[jj].im, m2[jj].re, stk1, RMODE);
    mpfr_fma(out->im, M[jj], stk1, out->im, RMODE);
  }
}

void dotpr(mpfr_t *out, mpfc_t *m1, mpfc_t *m2) {
  mpfr_set_ui(*out, 0, RMODE);
  for (int jj = 1; jj < N; jj++ ) {
    mpfr_mul(stk1, m1[jj].re, m2[jj].re, RMODE);
    mpfr_fma(stk1, m1[jj].im, m2[jj].im, stk1, RMODE);
    mpfr_fma(*out, M[jj], stk1, *out, RMODE);
  }
}

void dotpr_eye(mpfr_t *out, mpfc_t *m1, mpfc_t *m2) {
  mpfr_set_ui(*out, 0, RMODE);
  for (int jj = 1; jj < N; jj++ ) {
    mpfr_fma(*out, m1[jj].re, m2[jj].re, *out, RMODE);
    mpfr_fma(*out, m1[jj].im, m2[jj].im, *out, RMODE);
  }
}

void absolute2(mpfr_t *out, mpfc_t* in) {
  for (int jj = 0; jj < N; jj++) {
    mpfr_hypot(out[jj], in[jj].re, in[jj].im, RMODE);
    mpfr_mul(out[jj], out[jj], out[jj], RMODE);
  }
}

void assign_real(mpfr_t *out, mpfr_t *in) {
  for (int jj = 0; jj < N; jj++ ) mpfr_set(out[jj], in[jj], RMODE);
}

void assign_cmplx(mpfc_t *out, mpfc_t *in) {
  for (int jj = 0; jj < N; jj++ ) {
    mpfr_set(out[jj].re, in[jj].re, RMODE);
    mpfr_set(out[jj].im, in[jj].im, RMODE);
  }
}

void invert(mpfc_t *out, mpfc_t *in){
  for (int jj = 0; jj < N; jj++) {
    mpfr_hypot(stk1, in[jj].re, in[jj].im, RMODE);
    mpfr_mul(stk1, stk1, stk1, RMODE);

    mpfr_set(out[jj].re, in[jj].re, RMODE);
    mpfr_neg(out[jj].im, in[jj].im, RMODE);
    mpfr_div(out[jj].re, out[jj].re, stk1, RMODE);
    mpfr_div(out[jj].im, out[jj].im, stk1, RMODE);
  }
}

void err_msg(char* str)
{
  printf("%s\n",str);
  exit(1);
}

void compute_error(mpfr_t *err) {
  compute_approximation();
  write_cmplx(X, "PovQ.txt");
  for (int jj = 0; jj < N; jj++) {
    mpfr_sub(Y[jj].re, W[jj].re, X[jj].re, RMODE);
    mpfr_sub(Y[jj].im, W[jj].im, X[jj].im, RMODE);
  }
  dotpr_eye(&stk2, Y, Y);
  dotpr_eye(&stk3, W, W);
  mpfr_div(*err, stk2, stk3, RMODE);
  mpfr_sqrt(*err, *err, RMODE);
  mpfr_sqrt(stk3, stk3, RMODE);
  mpfr_printf("Approximation error is %.10Re\n", *err); 
}

void aberth_iterations(){
  mpfc_t *zP, *zPp, *zQ, *zQp, *logD, *Acor;
  mpfr_t arg, nrm, Isum; 
  mpfc_t w, w1;
  FILE *fhroots, *fhroots2;
  int sumF = 0;
  char rname[120];

  mpfr_init(arg);    mpfr_init(nrm);    mpfr_init(Isum);
  mpfr_init(w.re);   mpfr_init(w.im);
  mpfr_init(w1.re);  mpfr_init(w1.im);
  zP   =  malloc(d*sizeof(mpfc_t));
  zQ   =  malloc(d*sizeof(mpfc_t));
  zPp  =  malloc(d*sizeof(mpfc_t));
  zQp  =  malloc(d*sizeof(mpfc_t));
  logD =  malloc(d*sizeof(mpfc_t));
  Acor =  malloc(d*sizeof(mpfc_t));
  mpfr_printf("Circle radius set to %.5Re\n", R);
  mpfr_printf("Iterating until |Q(z_k)| < %.5Re for every k\n", tolEA);
  for (int ll = 0; ll < d; ll++) {
   mpfr_init(zP[ll].re);	   mpfr_init(zP[ll].im);
   mpfr_init(logD[ll].re);	   mpfr_init(logD[ll].im);
   mpfr_init(Acor[ll].re);	   mpfr_init(Acor[ll].im);
   mpfr_init(zQ[ll].re);	   mpfr_init(zQ[ll].im);
   mpfr_init(zPp[ll].re);	   mpfr_init(zPp[ll].im);
   mpfr_init(zQp[ll].re);	   mpfr_init(zQp[ll].im);
   mpfr_set_ui(arg, 2*ll + 1, RMODE);
   mpfr_div_ui(arg, arg, d, RMODE);
   mpfr_mul(arg, arg, pie, RMODE);
   mpfr_sin_cos(z[ll].im, z[ll].re, arg, RMODE);
   mpfr_mul(z[ll].re, R, z[ll].re, RMODE);
   mpfr_sub(z[ll].re, z[ll].re, R, RMODE);
   mpfr_mul(z[ll].im, z[ll].im, R, RMODE);
  }
  int m = 0; mpfr_set_zero(Isum, 1);
  while (sumF != d) {
    if ( (m % 5) == 0 ) {
      sprintf(rname, "roots_%03d.txt", m);
      fhroots = fopen(rname, "w");
      fprintf(fhroots, "# 1. Re z_k 2. Im z_k 3. |Q(z_k)|\n# Iteration = %d  ", m);
    }
    polyevalN(zP, zPp, zQ, zQp, z, d);
    sumF = 0; mpfr_set_zero(Isum, 1);
    for (int ll = 0; ll < d; ll++) {  
      mpfr_hypot(nrm, zQ[ll].re, zQ[ll].im, RMODE);
      mpfr_add(Isum, Isum, nrm, RMODE);
      sumF = sumF + abFLAG[ll];
      mpfr_printf("|Q(z_%d)| = %.8Re\n", ll, nrm); 
      if ( (ll==0)&&((m%5)==0) ) mpfr_fprintf(fhroots, "Sum|Q(z_k)| = %.8Re\n\n", Isum);
      if ((m%5)==0) mpfr_fprintf(fhroots, "%.12Re\t%.12Re\n", z[ll].re, z[ll].im); 
      if (mpfr_cmp(nrm, tolEA) <= 0) {
        abFLAG[ll] = 1;
      } else {
        //-------------------------------------------------------------------------- 
        mpfr_sqr(arg, zQp[ll].im, RMODE);
        mpfr_fma(arg, zQp[ll].re, zQp[ll].re, arg, RMODE);
        mpfr_mul(logD[ll].re, zQ[ll].re, zQp[ll].re, RMODE);
        mpfr_fma(logD[ll].re, zQ[ll].im, zQp[ll].im, logD[ll].re, RMODE);   
        mpfr_div(logD[ll].re, logD[ll].re, arg, RMODE);
        mpfr_mul(logD[ll].im, zQ[ll].re, zQp[ll].im, RMODE);
        mpfr_fms(logD[ll].im, zQ[ll].im, zQp[ll].re, logD[ll].im, RMODE);
        mpfr_div(logD[ll].im, logD[ll].im, arg, RMODE);
        //--------------------------------------------------------------------------
        mpfr_set_zero(Acor[ll].re, 1);
        mpfr_set_zero(Acor[ll].im, 1);
        for (int kk = 0; kk < d; kk++) {
          if (kk != ll) {
            mpfr_sub(w.re, z[ll].re, z[kk].re, RMODE);
	    mpfr_sub(w.im, z[ll].im, z[kk].im, RMODE);
            mpfr_sqr(arg, w.re, RMODE);
            mpfr_fma(arg, w.im, w.im, arg, RMODE);
            mpfr_div(w.re, w.re, arg, RMODE);
            mpfr_div(w.im, w.im, arg, RMODE);
            mpfr_add(Acor[ll].re, Acor[ll].re, w.re, RMODE);
            mpfr_sub(Acor[ll].im, Acor[ll].im, w.im, RMODE);
          }
        }   
        mpfr_sqr(arg, logD[ll].im, RMODE); 
        mpfr_fma(arg, logD[ll].re, logD[ll].re, arg, RMODE);
        mpfr_fms(w.re, arg, Acor[ll].re, logD[ll].re, RMODE);
        mpfr_fma(w.im, arg, Acor[ll].im, logD[ll].im, RMODE);
        mpfr_mul(w1.re, logD[ll].re, Acor[ll].re, RMODE);
        mpfr_fms(w1.re, logD[ll].im, Acor[ll].im, w1.re, RMODE);
        mpfr_add_ui(w1.re, w1.re, 1, RMODE);
        mpfr_sqr(w1.re, w1.re, RMODE);      
        mpfr_mul(w1.im, logD[ll].re, Acor[ll].im, RMODE);
        mpfr_fma(w1.im, logD[ll].im, Acor[ll].re, w1.im, RMODE);
        mpfr_sqr(w1.im, w1.im, RMODE);
        mpfr_add(arg, w1.re, w1.im, RMODE);
        mpfr_div(w.re, w.re, arg, RMODE);      
        mpfr_div(w.im, w.im, arg, RMODE);
        mpfr_add(z[ll].re, z[ll].re, w.re, RMODE);
        mpfr_sub(z[ll].im, z[ll].im, w.im, RMODE);
      }
    }
    if ( (m%5) == 0) {
        mpfr_printf("At E-A iteration %3d :", m);
        mpfr_printf("Sum|Q(z_k)| = %.8Re\n", Isum);
      }
    if ((m%5)==0) fclose(fhroots); 
    m++;
    if ( m > N_MAXIT) break;
  }
  fhroots = fopen("all_roots.txt", "w");
  fhroots2 = fopen("pade.txt", "w");
  fprintf(fhroots, "# 1. Re z 2. Im z 3. Re Res 4. Im Res\n");
  fprintf(fhroots2, "# 1. Im z 3. Res \n");
  mpfr_fprintf(fhroots, "# Sum|Q(z_k)| = %.8Re\t", Isum);
  mpfr_fprintf(fhroots2, "# Sum|Q(z_k)| = %.8Re\t", Isum);
  if (phaseFLAG == 1) {
    sort_by_phase();
    fprintf(fhroots, "Sorted by phase\n\n");
    fprintf(fhroots2, "Sorted by phase\n\n");
  } else if (phaseFLAG == 0){
    sort_by_abs();
    fprintf(fhroots, "Sorted by absolute value\n\n");
    fprintf(fhroots2, "Sorted by absolute value\n\n");
  } else {
    sort_by_real();
    fprintf(fhroots, "Sorted by real part value\n\n");
    fprintf(fhroots2, "Sorted by real part value\n\n");
  }
  polyevalN(zP, zPp, zQ, zQp, z, d); 
  for (int jj = 0; jj < d; jj++) {
    /*mpfr_printf("P = %.4Re  +  %.4Re\nQ = %.4Re  +  %.4Re\nP' = %.4Re  +  %.4Re\nQ' = %.4Re  +  %.4Re\n\n", 
          zP[jj].re, zP[jj].im, zQ[jj].re, zQ[jj].im, zPp[jj].re, zPp[jj].im, zQp[jj].re, zQp[jj].im);*/
    mpfr_mul(arg, zQp[jj].re, zQp[jj].re, RMODE);
    mpfr_fma(arg, zQp[jj].im, zQp[jj].im, arg, RMODE);

    mpfr_mul(w.re, zP[jj].re, zQp[jj].re, RMODE);
    mpfr_fma(w.re, zP[jj].im, zQp[jj].im, w.re, RMODE);
    mpfr_div(res[jj].re, w.re, arg, RMODE);

    mpfr_mul(w.im, zP[jj].re, zQp[jj].im, RMODE);
    mpfr_fms(w.im, zP[jj].im, zQp[jj].re, w.im, RMODE);
    mpfr_div(res[jj].im, w.im, arg, RMODE);
   
    mpfr_div(zsc[jj].re, z[jj].re, L, RMODE);
    mpfr_div(zsc[jj].im, z[jj].im, L, RMODE);

    mpfr_fprintf(fhroots, "%.120Re\t%.120Re\t%.120Re\t%.120Re\n", z[jj].re, z[jj].im, res[jj].re, res[jj].im);
    mpfr_fprintf(fhroots2, "%.200Re\t%.200Re\n", zsc[jj].re, res[jj].re);
  }
  fclose(fhroots);
  fclose(fhroots2);
  printf("On final iteration:\n");
  mpfr_printf("At E-A iteration %3d:", m);
  mpfr_printf("Sum|Q(z_k)| = %.8Re\n", Isum);
}

void newton_suppression() {
  mpfc_t *zP, *zPp, *zQ, *zQp, *logD, *Acor; //, *root;
  mpfr_t arg, nrm, Isum, aux; 
  mpfc_t w, w1;
  FILE *fhroots = fopen("roots_00.txt","w");
  FILE *fhroots2;
  //int sumF = 0;
  //char rname[120];
  
  fclose(fhroots);
  mpfr_init(arg);    mpfr_init(nrm);    mpfr_init(Isum);
  mpfr_init(w.re);   mpfr_init(w.im);   mpfr_init(aux);
  mpfr_init(w1.re);  mpfr_init(w1.im);
  zP   =  malloc(d*sizeof(mpfc_t));
  zQ   =  malloc(d*sizeof(mpfc_t));
  zPp  =  malloc(d*sizeof(mpfc_t));
  zQp  =  malloc(d*sizeof(mpfc_t));
  logD =  malloc(d*sizeof(mpfc_t));
  Acor =  malloc(d*sizeof(mpfc_t));

  mpfr_printf("Circle radius set to %.5Re\n", R);
  mpfr_printf("Iterating until |Q(z_k)| < %.5Re for every k\n", tolEA);

  for (int ll = 0; ll < d; ll++) {
    mpfr_init(zP[ll].re);	   	   mpfr_init(zP[ll].im);
    mpfr_init(logD[ll].re);		   mpfr_init(logD[ll].im);
    mpfr_init(Acor[ll].re);		   mpfr_init(Acor[ll].im);
    mpfr_init(zQ[ll].re);	   	   mpfr_init(zQ[ll].im);
    mpfr_init(zPp[ll].re);		   mpfr_init(zPp[ll].im);
    mpfr_init(zQp[ll].re);		   mpfr_init(zQp[ll].im);
    mpfr_set_ui(arg, 2*ll + 1, RMODE);
    mpfr_div_ui(arg, arg, d, RMODE);
    mpfr_mul(arg, arg, pie, RMODE);
    //mpfr_sin_cos(z[ll].im, z[ll].re, arg, RMODE);
    //mpfr_mul(z[ll].re, R, z[ll].re, RMODE);
    //mpfr_sub(z[ll].re, z[ll].re, R, RMODE);
    //mpfr_mul(z[ll].im, z[ll].im, R, RMODE);
    mpfr_set_d(z[ll].im,  0.25, RMODE);
    mpfr_set_d(z[ll].re, -1.00, RMODE);
  }

  int iter; 
  mpfr_set_zero(Isum, 1);
  for (int m = 0; m < d; m++) {
    //polyevalN(zP, zPp, zQ, zQp, &z[m], 1);
    //mpfr_set_zero(Isum, 1);
    mpfr_set_ui(nrm, 1, RMODE);
    iter = 0;
    while (mpfr_cmp(nrm, tolEA) > 0 ) {  
      polyevalN(&zP[m], &zPp[m], &zQ[m], &zQp[m], &z[m], 1);
      mpfr_hypot(nrm, zQ[m].re, zQ[m].im, RMODE);
      //-------------------------------------------------------------------------- 
      mpfr_sqr(arg, zQp[m].im, RMODE);
      mpfr_fma(arg, zQp[m].re, zQp[m].re, arg, RMODE);
      mpfr_mul(logD[m].re, zQ[m].re, zQp[m].re, RMODE);
      mpfr_fma(logD[m].re, zQ[m].im, zQp[m].im, logD[m].re, RMODE);         // logD = Q/Q'
      mpfr_div(logD[m].re, logD[m].re, arg, RMODE);
      mpfr_mul(logD[m].im, zQ[m].re, zQp[m].im, RMODE);
      mpfr_fms(logD[m].im, zQ[m].im, zQp[m].re, logD[m].im, RMODE);
      mpfr_div(logD[m].im, logD[m].im, arg, RMODE);
      //--------------------------------------------------------------------------
      mpfr_set_zero(Acor[m].re, 1);
      mpfr_set_zero(Acor[m].im, 1);
      for (int kk = 0; kk < m; kk++) {
          mpfr_sub(w.re, z[m].re, z[kk].re, RMODE);
	  mpfr_sub(w.im, z[m].im, z[kk].im, RMODE);
          mpfr_sqr(arg, w.re, RMODE);
          mpfr_fma(arg, w.im, w.im, arg, RMODE);
          mpfr_div(w.re, w.re, arg, RMODE);
          mpfr_div(w.im, w.im, arg, RMODE);
          mpfr_add(Acor[m].re, Acor[m].re, w.re, RMODE);
          mpfr_sub(Acor[m].im, Acor[m].im, w.im, RMODE);
      }   
      mpfr_sqr(arg, logD[m].im, RMODE); 
      mpfr_fma(arg, logD[m].re, logD[m].re, arg, RMODE);
      mpfr_fms(w.re, arg, Acor[m].re, logD[m].re, RMODE);  // w.re = N.re - |N|^2 A.re
      mpfr_fma(w.im, arg, Acor[m].im, logD[m].im, RMODE);  // w.im = N.im + |N|^2 A.im
      //--------------------------------------------------------------------------
      mpfr_mul(w1.re, logD[m].re, Acor[m].re, RMODE);
      mpfr_fms(w1.re, logD[m].im, Acor[m].im, w1.re, RMODE);
      mpfr_add_ui(w1.re, w1.re, 1, RMODE);
      mpfr_sqr(w1.re, w1.re, RMODE);      
      mpfr_mul(w1.im, logD[m].re, Acor[m].im, RMODE);
      mpfr_fma(w1.im, logD[m].im, Acor[m].re, w1.im, RMODE);
      mpfr_sqr(w1.im, w1.im, RMODE);
      mpfr_add(arg, w1.re, w1.im, RMODE);  // arg = |1 - NA|^2
      //--------------------------------------------------------------------------
      mpfr_div(w.re, w.re, arg, RMODE);      
      mpfr_div(w.im, w.im, arg, RMODE);
      mpfr_add(z[m].re, z[m].re, w.re, RMODE); // z.re = z.re - w.re/arg
      mpfr_sub(z[m].im, z[m].im, w.im, RMODE); // z.im = z.im - w.im/arg
      //--------------------------------------------------------------------------
      iter++;
      //mpfr_printf("|Q(z_%d)| = %.6Re at Newton step %d\n", m, nrm, iter); 
      if (iter > N_MAXIT) {
        mpfr_printf("Warning! Max Iterations Reached for root %d\n", m); 
	break;
      } 
    }
    for (int l = 0; l < EXTRA_NEWTON; l++) {
      polyevalN(&zP[m], &zPp[m], &zQ[m], &zQp[m], &z[m], 1);
      mpfr_hypot(nrm, zQ[m].re, zQ[m].im, RMODE);
      //-------------------------------------------------------------------------- 
      mpfr_sqr(arg, zQp[m].im, RMODE);
      mpfr_fma(arg, zQp[m].re, zQp[m].re, arg, RMODE);
      mpfr_mul(logD[m].re, zQ[m].re, zQp[m].re, RMODE);
      mpfr_fma(logD[m].re, zQ[m].im, zQp[m].im, logD[m].re, RMODE);         // logD = Q/Q'
      mpfr_div(logD[m].re, logD[m].re, arg, RMODE);
      mpfr_mul(logD[m].im, zQ[m].re, zQp[m].im, RMODE);
      mpfr_fms(logD[m].im, zQ[m].im, zQp[m].re, logD[m].im, RMODE);
      mpfr_div(logD[m].im, logD[m].im, arg, RMODE);
      //--------------------------------------------------------------------------
      mpfr_set_zero(Acor[m].re, 1);
      mpfr_set_zero(Acor[m].im, 1);
      for (int kk = 0; kk < m; kk++) {
          mpfr_sub(w.re, z[m].re, z[kk].re, RMODE);
	  mpfr_sub(w.im, z[m].im, z[kk].im, RMODE);
          mpfr_sqr(arg, w.re, RMODE);
          mpfr_fma(arg, w.im, w.im, arg, RMODE);
          mpfr_div(w.re, w.re, arg, RMODE);
          mpfr_div(w.im, w.im, arg, RMODE);
          mpfr_add(Acor[m].re, Acor[m].re, w.re, RMODE);
          mpfr_sub(Acor[m].im, Acor[m].im, w.im, RMODE);
      }   
      mpfr_sqr(arg, logD[m].im, RMODE); 
      mpfr_fma(arg, logD[m].re, logD[m].re, arg, RMODE);
      mpfr_fms(w.re, arg, Acor[m].re, logD[m].re, RMODE);  // w.re = N.re - |N|^2 A.re
      mpfr_fma(w.im, arg, Acor[m].im, logD[m].im, RMODE);  // w.im = N.im + |N|^2 A.im
      //--------------------------------------------------------------------------
      mpfr_mul(w1.re, logD[m].re, Acor[m].re, RMODE);
      mpfr_fms(w1.re, logD[m].im, Acor[m].im, w1.re, RMODE);
      mpfr_add_ui(w1.re, w1.re, 1, RMODE);
      mpfr_sqr(w1.re, w1.re, RMODE);      
      mpfr_mul(w1.im, logD[m].re, Acor[m].im, RMODE);
      mpfr_fma(w1.im, logD[m].im, Acor[m].re, w1.im, RMODE);
      mpfr_sqr(w1.im, w1.im, RMODE);
      mpfr_add(arg, w1.re, w1.im, RMODE);  // arg = |1 - NA|^2
      //--------------------------------------------------------------------------
      mpfr_div(w.re, w.re, arg, RMODE);      
      mpfr_div(w.im, w.im, arg, RMODE);
      mpfr_add(z[m].re, z[m].re, w.re, RMODE); // z.re = z.re - w.re/arg
      mpfr_sub(z[m].im, z[m].im, w.im, RMODE); // z.im = z.im - w.im/arg
      //--------------------------------------------------------------------------
    }


    mpfr_add(Isum, Isum, nrm, RMODE);
    mpfr_printf("|Q(z_%d)| = %.6Re at Newton step %d. Doing extra %d iterations to avoid loss of precision.\n", m, nrm, iter, EXTRA_NEWTON);
  }


  fhroots = fopen("all_roots.txt", "w");
  fhroots2 = fopen("pade.txt", "w");
  fprintf(fhroots, "# 1. Re z 2. Im z 3. Re Res 4. Im Res\n");
  fprintf(fhroots2, "# 1. Im z 2. Re Res\n");
  mpfr_fprintf(fhroots, "# Sum|Q(z_k)| = %.8Re\t", Isum);
  mpfr_fprintf(fhroots2, "# Sum|Q(z_k)| = %.8Re\t", Isum);
  if (phaseFLAG == 1) {
    sort_by_phase();
    fprintf(fhroots, "Sorted by phase\n\n");
    fprintf(fhroots2, "Sorted by phase\n\n");
  } else if (phaseFLAG == 0) {
    sort_by_abs();
    fprintf(fhroots, "Sorted by absolute value\n\n");
    fprintf(fhroots2, "Sorted by absolute value\n\n");
  } else {
    sort_by_real();
    fprintf(fhroots, "Sorted by real part value\n\n");
    fprintf(fhroots2, "Sorted by real part value\n\n");
  }
  polyevalN(zP, zPp, zQ, zQp, z, d); 
  for (int jj = 0; jj < d; jj++) {

    mpfr_mul(arg, zQp[jj].re, zQp[jj].re, RMODE);
    mpfr_fma(arg, zQp[jj].im, zQp[jj].im, arg, RMODE);

    mpfr_mul(w.re, zP[jj].re, zQp[jj].re, RMODE);
    mpfr_fma(w.re, zP[jj].im, zQp[jj].im, w.re, RMODE);
    mpfr_div(res[jj].re, w.re, arg, RMODE);

    mpfr_mul(w.im, zP[jj].re, zQp[jj].im, RMODE);
    mpfr_fms(w.im, zP[jj].im, zQp[jj].re, w.im, RMODE);
    mpfr_div(res[jj].im, w.im, arg, RMODE);

    mpfr_div(zsc[jj].re, z[jj].re, L, RMODE);
    mpfr_div(zsc[jj].im, z[jj].im, L, RMODE);

    mpfr_fprintf(fhroots, "%.120Re\t%.120Re\t%.120Re\t%.120Re\n", z[jj].re, z[jj].im, res[jj].re, res[jj].im);

    mpfr_div(ressc[jj].re, res[jj].re, L, RMODE);
    mpfr_neg(ressc[jj].re, ressc[jj].re, RMODE);
    mpfr_fprintf(fhroots2, "%.200Re\t%.200Re\n", zsc[jj].re, ressc[jj].re);
  }
  fclose(fhroots);
  fclose(fhroots2);

  FILE *fhslau = fopen("FqCqdebug.txt","w");
  fprintf(fhslau, " rhoBd:\n\tN/A\n ell:\n\t\tN/A\n Nsing:\n\t\t%ld\n Npair:\n\t\tN/A\n ===============\n\n",d);
  fprintf(fhslau, " Beta locations:\n"); 
  for (int jj = 0; jj < d; jj++) {
    mpfr_fprintf(fhslau, "%2d\t%.32Re\t%.32Re\n", jj+1, z[d-jj-1].re, z[d-jj-1].im);
  }
  fprintf(fhslau, "\n Gamma strengths:\n");
  for (int jj = 0; jj < d; jj++) {
    mpfr_neg(res[d-jj-1].re, res[d-jj-1].re, RMODE);
    mpfr_neg(res[d-jj-1].im, res[d-jj-1].im, RMODE);
    mpfr_fprintf(fhslau, "%2d\t%.32Re\t%.32Re\n", jj+1, res[d-jj-1].re, res[d-jj-1].im);
    mpfr_neg(res[d-jj-1].re, res[d-jj-1].re, RMODE);
    mpfr_neg(res[d-jj-1].im, res[d-jj-1].im, RMODE);
  }
  fclose(fhslau);

  FILE *fhpavel = fopen("formatted.txt","w");
  fprintf(fhpavel, "# 1. chi 2. gamma_re 3. d_chi 4. rho\n");
  fprintf(fhpavel, "# N = %ld\n\n", N);
  for (int jj = 0; jj < d; jj++) {
    mpfr_neg(z[d-jj-1].re, z[d-jj-1].re, RMODE);
    mpfr_fprintf(fhpavel, "%.120Re\t%.120Re\t 0.\t0.\n", z[d-jj-1].re, res[d-jj-1].re);
    mpfr_neg(z[d-jj-1].re, z[d-jj-1].re, RMODE);
  }
  fclose(fhpavel);
}

void sort_by_phase() {
  mpfr_t swap, mx;
  mpfr_t *list;
  
  list = malloc(d*sizeof(mpfr_t));
  mpfr_init(swap); mpfr_init(mx); 
  for (int jj = 0; jj < d; jj++) {
    mpfr_init(list[jj]);
    mpfr_neg(mx, z[jj].re, RMODE);
    mpfr_atan2(list[jj], z[jj].im, mx, RMODE);
  }
  for (int jj = 0; jj < d-1; jj++) {
    for (int kk = 0; kk < d-jj-1; kk++) {
      mpfr_neg(mx, z[kk].re, RMODE);
      mpfr_atan2(list[kk], z[kk].im, mx, RMODE);
      mpfr_neg(mx, z[kk+1].re, RMODE);
      mpfr_atan2(list[kk+1], z[kk+1].im, mx, RMODE);
      if (mpfr_cmp(list[kk], list[kk+1]) > 0) {
        mpfr_set(swap, z[kk].re, RMODE);
        mpfr_set(z[kk].re, z[kk+1].re, RMODE);
        mpfr_set(z[kk+1].re, swap, RMODE);
	mpfr_set(swap, z[kk].im, RMODE);
        mpfr_set(z[kk].im, z[kk+1].im, RMODE);
        mpfr_set(z[kk+1].im, swap, RMODE);
      }
    }  
  }
}

void sort_by_abs() {
  mpfr_t swap;
  mpfr_t *list;  
  list = malloc(d*sizeof(mpfr_t));
  mpfr_init(swap);
  for (int jj = 0; jj < d; jj++) {
    mpfr_init(list[jj]);
    mpfr_hypot(list[jj], z[jj].im, z[jj].re, RMODE);
  }
  for (int jj = 0; jj < d-1; jj++) {
    for (int kk = 0; kk < d-jj-1; kk++) {
      mpfr_hypot(list[kk], z[kk].im, z[kk].re, RMODE);
      mpfr_hypot(list[kk+1], z[kk+1].im, z[kk+1].re, RMODE);
      if (mpfr_cmp(list[kk], list[kk+1]) > 0) {
        mpfr_set(swap, z[kk].re, RMODE);
        mpfr_set(z[kk].re, z[kk+1].re, RMODE);
        mpfr_set(z[kk+1].re, swap, RMODE);
	mpfr_set(swap, z[kk].im, RMODE);
        mpfr_set(z[kk].im, z[kk+1].im, RMODE);
        mpfr_set(z[kk+1].im, swap, RMODE);
      }
    }  
  }
}

void sort_by_real() {
  mpfr_t swap;
  mpfr_t *list;  
  list = malloc(d*sizeof(mpfr_t));
  mpfr_init(swap);
  for (int jj = 0; jj < d; jj++) {
    mpfr_init_set(list[jj], z[jj].im, RMODE);
    //mpfr_hypot(list[jj], z[jj].im, z[jj].re, RMODE);
  }
  for (int jj = 0; jj < d-1; jj++) {
    for (int kk = 0; kk < d-jj-1; kk++) {
      //mpfr_hypot(list[kk], z[kk].im, z[kk].re, RMODE);
      //mpfr_hypot(list[kk+1], z[kk+1].im, z[kk+1].re, RMODE);
      if (mpfr_cmp(list[kk], list[kk+1]) > 0) {
        mpfr_set(swap, z[kk].re, RMODE);
        mpfr_set(z[kk].re, z[kk+1].re, RMODE);
        mpfr_set(z[kk+1].re, swap, RMODE);
	mpfr_set(swap, z[kk].im, RMODE);
        mpfr_set(z[kk].im, z[kk+1].im, RMODE);
        mpfr_set(z[kk+1].im, swap, RMODE);
      }
    }  
  }
}

void compare_pade () {
  //FILE *fhpade = fopen("err_pade.dat", "w");
  mpfc_t w, w1;
  mpfr_t arg;
 
  mpfr_init(w.re);	mpfr_init(w.im);
  mpfr_init(w1.re);	mpfr_init(w1.im);
  mpfr_init(arg);
  //fprintf(fhpade, "# 1. u 2. re Error 3. im Error\n\n");
  for (int jj = 0; jj < N; jj ++) {
    mpfr_set_zero(X[jj].re, 1);
    mpfr_set_zero(X[jj].im, 1);
    for (int kk = 0; kk < d; kk++) {
      mpfr_neg(w.re, z[kk].re, RMODE);
      mpfr_sub(w.im, y[jj], z[kk].im, RMODE);
      mpfr_mul(arg, w.re, w.re, RMODE);
      mpfr_fma(arg, w.im, w.im, arg, RMODE);
      mpfr_mul(w1.re, res[kk].im, w.im, RMODE);
      mpfr_fma(w1.re, res[kk].re, w.re, w1.re, RMODE);
      mpfr_mul(w1.im, res[kk].re, w.im, RMODE);
      mpfr_fms(w1.im, res[kk].im, w.re, w1.im, RMODE); 
      mpfr_div(w.re, w1.re, arg, RMODE);
      mpfr_div(w.im, w1.im, arg, RMODE);
      mpfr_add(X[jj].re, X[jj].re, w.re, RMODE);
      mpfr_add(X[jj].im, X[jj].im, w.im, RMODE);
    }
    mpfr_sqr(arg, y[jj], RMODE);
    mpfr_add_ui(arg, arg, 1, RMODE);
    mpfr_div_ui(arg, arg, 2, RMODE);
    mpfr_sub(Y[jj].re, W[jj].re, X[jj].re, RMODE);
    mpfr_sub(Y[jj].im, W[jj].im, X[jj].im, RMODE);
  }
  //write_cmplx(Y, "pade_err.dat");
  //write_cmplx(X, "pade_approx.dat");
  dotpr_eye(&stk2, Y, Y);
  dotpr_eye(&stk3, W, W);
  mpfr_div(arg, stk2, stk3, RMODE);
  mpfr_sqrt(arg, arg, RMODE);
  mpfr_printf("Relative error of Pade is %.8Re\n", arg);
  FILE *fherr = fopen("compress.log","w");
  mpfr_fprintf(fherr, "compress.f relative error: %.8Re\n", arg);
  fclose(fherr);
  
}

void test_funcs() {
  mpfr_t err;
  mpfc_t result;
  mpfr_t aux, arg, xx, yy, zz;
  mpfr_t Abs1, Phi1, Phi2;
  mpfr_t PieOverTwo;

  mpfr_init(PieOverTwo);
  mpfr_const_pi(PieOverTwo, RMODE);
  mpfr_div_ui(PieOverTwo, PieOverTwo, 2, RMODE);
  mpfr_init(zz);
  mpfr_init(yy);
  mpfr_init(xx);
  mpfr_init(aux);
  mpfr_init(arg);
  mpfr_init(err);
  mpfr_init(result.re);
  mpfr_init(result.im);
  if (DEMO_SQRT_SYMM) {
    mpfr_t 	A, B;
    mpfr_t	d1, d2;
    mpfr_init_set_ui(A, 1, RMODE);
    mpfr_init_set_ui(B, 2, RMODE);
    mpfr_init(d1);
    mpfr_init(d2);
    mpfr_init(Abs1);
    mpfr_init(Phi1);
    mpfr_init(Phi2);
    for (int jj = 0; jj < N; jj++) {
      // set phase
      mpfr_add(aux, A, B, RMODE);
      mpfr_mul(Phi1, y[jj], aux, RMODE); 

      mpfr_sqr(Phi2, y[jj], RMODE);
      mpfr_mul(aux, A, B, RMODE);
      mpfr_sub(Phi2, Phi2, aux, RMODE);
      mpfr_atan2(arg, Phi1, Phi2, RMODE);
      mpfr_div_ui(arg, arg, 2, RMODE);
      mpfr_add(arg, arg, PieOverTwo, RMODE);
      mpfr_cos(xx, arg, RMODE); 
      mpfr_sin(yy, arg, RMODE); 
      // set abs-value
      mpfr_sqr(aux, y[jj], RMODE);
      mpfr_fms(Abs1, A, B, aux, RMODE);
      mpfr_add(aux, A, B, RMODE);
      mpfr_mul(aux, aux, y[jj], RMODE);
      mpfr_sqr(aux, aux, RMODE);
      mpfr_add(Abs1, Abs1, aux, RMODE);
      mpfr_sqrt(Abs1, Abs1, RMODE);
      mpfr_sqrt(Abs1, Abs1, RMODE);
      // set real/imaginary parts
      mpfr_mul(xx, xx, Abs1, RMODE);
      mpfr_mul(yy, yy, Abs1, RMODE);
      mpfr_hypot(d1, y[jj], A, RMODE);
      mpfr_hypot(d2, y[jj], B, RMODE);
      mpfr_mul(aux, d1, d2, RMODE);
      //mpfr_div(W[jj].re, xx, aux, RMODE); 
      //mpfr_div(W[jj].im, yy, aux, RMODE); 
      mpfr_neg(W[jj].im, W[jj].im, RMODE); 
    }
    printf("Initialized to 1/sqrt((z-ia)(ib-z))\n");
  }
/*
  for (int jj = 0; jj < N; jj++) {
    mpfr_sqr(aux, y[jj], RMODE);
    mpfr_sqr(aux, aux, RMODE);  // y^4
    mpfr_add_ui(aux, aux, 4, RMODE);
    mpfr_sqrt(aux, aux, RMODE);
    mpfr_sqrt(aux, aux, RMODE);  // (4 + y^4)^(1/4)
    mpfr_ui_div(W[jj].re, 1, aux, RMODE);
    mpfr_ui_div(W[jj].im, 1, aux, RMODE);
    mpfr_mul_ui(yy, y[jj], 2, RMODE);
    mpfr_sqr(xx, y[jj], RMODE);  
    mpfr_ui_sub(xx, 2, xx, RMODE);
    mpfr_atan2(arg, yy, xx, RMODE);
    mpfr_div_ui(arg, arg, 2, RMODE); // arg = 0.5atan(2y/(4-y^2))
    mpfr_cos(xx, arg, RMODE);
    mpfr_sin(yy, arg, RMODE); 
    mpfr_mul(W[jj].re, W[jj].re, xx, RMODE);
    mpfr_mul(W[jj].im, W[jj].im, yy, RMODE);
    mpfr_neg(W[jj].im, W[jj].im, RMODE);
  }*/
  
  /*for (int jj = 0; jj < N; jj++) {
    mpfr_add_ui(xx, y[jj], 1, RMODE);
    mpfr_sub_ui(yy, y[jj], 1, RMODE);
    mpfr_set(result.re, xx, RMODE);
    mpfr_set(result.im, yy, RMODE);

    mpfr_sqr(xx, xx, RMODE);
    mpfr_sqr(yy, yy, RMODE);

    mpfr_add_ui(xx, xx, 1, RMODE);
    mpfr_add_ui(yy, yy, 1, RMODE);
    mpfr_ui_div(xx, 1, xx, RMODE);
    mpfr_ui_div(yy, 1, yy, RMODE);
    mpfr_add(W[jj].re, xx, yy, RMODE);
    mpfr_mul(W[jj].im, result.re, xx, RMODE);
    mpfr_fma(W[jj].im, result.im, yy, W[jj].im, RMODE);
    mpfr_neg(W[jj].im, W[jj].im, RMODE);

    mpfr_sqr(xx, y[jj], RMODE);
    mpfr_add_ui(xx, xx, 1, RMODE);
    mpfr_ui_div(xx, 2, xx, RMODE);
    mpfr_add(W[jj].re, W[jj].re, xx, RMODE);

    mpfr_mul(xx, xx, y[jj], RMODE);
    mpfr_sub(W[jj].im, W[jj].im, xx, RMODE);
  }*/

  /*for (int jj = 0; jj < N; jj++) {
    mpfr_ui_sub(xx, 1, y[jj], RMODE);
    mpfr_add_ui(yy, y[jj], 1, RMODE);
    mpfr_set(result.re, xx, RMODE); // result.re = 1-y
    mpfr_set(result.im, yy, RMODE); // result.im = 1+y
    mpfr_sqr(xx, xx, RMODE);
    mpfr_sqr(yy, yy, RMODE);
    mpfr_sqr(zz, y[jj], RMODE);
    mpfr_add_ui(xx, xx, 1, RMODE);
    mpfr_add_ui(yy, yy, 1, RMODE);
    mpfr_add_ui(zz, zz, 1, RMODE);
    mpfr_ui_div(xx, 1, xx, RMODE);  // xx = 1/(1 + (1-y)^2)
    mpfr_ui_div(yy, 1, yy, RMODE);  // yy = 1/(1 + (1+y)^2)
    mpfr_ui_div(zz, 1, zz, RMODE);  // zz = 1/(1 + y^2)
    mpfr_add_ui(aux, result.re, 1, RMODE);
    mpfr_mul(W[jj].re, aux, xx, RMODE);
    mpfr_add_ui(aux, result.im, 1, RMODE);
    mpfr_fma(W[jj].re, aux, yy, W[jj].re, RMODE);
    mpfr_add(W[jj].re, W[jj].re, zz, RMODE);
    mpfr_mul(W[jj].im, y[jj], xx, RMODE);
    mpfr_fma(W[jj].im, y[jj], yy, W[jj].im, RMODE);
    mpfr_fma(W[jj].im, y[jj], zz, W[jj].im, RMODE);
    mpfr_neg(W[jj].im, W[jj].im, RMODE);
  }*/
  //write_cmplx(W, "3poles.txt");
  dotpr_eye(&aux, W, W);
  dotpr_eye(&arg, eye, eye);
  mpfr_div(aux, aux, arg, RMODE);
  //printf("Initialized to 1/sqrt(z^2+1)\n");
  mpfr_printf("Input data norm is %.10Re\n", aux);
  generate_Q0();
  int l = 0;
  while (l < N_AGH) {
    set_weight();    
    gram_schmidt();  
    update_Q();
    compute_error(&err);
    l++;
  }
  write_cmplx(P, "P.txt");
  write_cmplx(Q, "Q.txt");
  newton_suppression();
  compare_pade();
}



void my_mpfr_version() {
  printf("MPFR library: %-12s\nMPFR header: %s (based on %d.%d.%d)\n",
         mpfr_get_version (), MPFR_VERSION_STRING, MPFR_VERSION_MAJOR,
         MPFR_VERSION_MINOR, MPFR_VERSION_PATCHLEVEL);
}
