#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <math.h>
//#include <complex.h>
#include <mpfr.h>

#define RMODE		MPFR_RNDN
#define N_MAXIT		1000
#define EXTRA_NEWTON	32
#define AD_ABERTH	10
#define ROOT_OUT	1
#define N_AGH		14
#define DEMO_SQRT_SYMM	0

typedef struct profile {
  mpfr_t *raw_u, *raw_x, *raw_y;
  mpfr_t *m;
  long int n;
} vect;

typedef struct multiprecision_complex {
  mpfr_t re, im;
} mpfc_t;

extern void my_mpfr_version();
extern void init_program(int *argc, char **argv);
extern void finish_program();
extern void write_real(mpfr_t *out, char* str);
extern void write_cmplx(mpfc_t *out, char* str);
extern void err_msg(char* str);
extern void dotpr(mpfr_t *out, mpfc_t *m1, mpfc_t *m2);
extern void dotpc(mpfc_t *out, mpfc_t *m1, mpfc_t *m2);
extern void test_funcs();
extern void assign_cmplx(mpfc_t *out, mpfc_t *in);
extern void invert(mpfc_t *out, mpfc_t *in);
extern void absolute2(mpfr_t *out, mpfc_t* in);
extern void sort_by_phase();
extern void sort_by_abs();
extern void sort_by_real();
