/*
  DISCLAIMER:
  This software was produced by the National Institute of Standards
  and Technology (NIST), an agency of the U.S. government, and by
  statute is not subject to copyright in the United States.
  Recipients of this software assume all responsibility associated
  with its operation, modification, maintenance, and subsequent
  redistribution.
*/

#ifndef VARIATES_H
#define VARIATES_H

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

typedef struct {
  long int seed;
} unit_random_struct;

/*
  If you want separate chunks of random number generators, here are
  the seeds of 12 evenly-spread ranges with 178,956,970 in each:

  1101211447
  2021127233
  1925176231
  1304948567
  1375081611
  774234184
  676806766
  934251302
  1589551955
  1316071563
  1713378112
  573050001
*/

extern void unit_random_init(unit_random_struct *r);
extern void unit_random_seed(unit_random_struct *r, long int s);
extern long int unit_random_integer_min(unit_random_struct *r);
extern long int unit_random_integer_max(unit_random_struct *r);
extern long int unit_random_integer(unit_random_struct *r);
extern double unit_random_real(unit_random_struct *r);

typedef struct {
  unit_random_struct u;
  double min;
  double diff;
} uniform_random_struct;

extern void uniform_random_init(uniform_random_struct *r, double a, double b);
extern void uniform_random_set(uniform_random_struct *r, double a, double b);
extern void uniform_random_seed(uniform_random_struct *r, long int s);
extern double uniform_random_real(uniform_random_struct *r);

typedef struct {
  unit_random_struct u1, u2;
  double x1, x2;
  double mean;
  double sd;
  char return_x2;
} normal_random_struct;

extern void normal_random_init(normal_random_struct *r, double mean, double sd);
extern void normal_random_set(normal_random_struct *r, double mean, double sd);
extern void normal_random_seed(normal_random_struct *r, long int s1, long int s2);
extern double normal_random_real(normal_random_struct *r);

typedef struct {
  unit_random_struct u;
  double sd;
} exponential_random_struct;

extern void exponential_random_init(exponential_random_struct *r, double sd);
extern void exponential_random_set(exponential_random_struct *r, double sd);
extern void exponential_random_seed(exponential_random_struct *r, long int s);
extern double exponential_random_real(exponential_random_struct *r);

typedef struct {
  unit_random_struct u;
  double alpha_inv;
  double beta;
  char degen;
} weibull_random_struct;

extern void weibull_random_init(weibull_random_struct *r, double alpha, double beta);
extern void weibull_random_set(weibull_random_struct *r, double alpha, double beta);
extern void weibull_random_seed(weibull_random_struct *r, long int s);
extern double weibull_random_real(weibull_random_struct *r);

typedef struct {
  unit_random_struct u1, u2;
  double alpha, alpha_inv;
  double beta;
  double a, b, q, theta, d;
  char range;
} gamma_random_struct;

extern void gamma_random_init(gamma_random_struct *r, double alpha, double beta);
extern void gamma_random_set(gamma_random_struct *r, double alpha, double beta);
extern void gamma_random_seed(gamma_random_struct *r, long int s1, long int s2);
extern double gamma_random_real(gamma_random_struct *r);

typedef struct {
  gamma_random_struct g;
} pearson_v_random_struct;

extern void pearson_v_random_init(pearson_v_random_struct *r, double alpha, double beta);
extern void pearson_v_random_set(pearson_v_random_struct *r, double alpha, double beta);
extern void pearson_v_random_seed(pearson_v_random_struct *r, long int s1, long int s2);
extern double pearson_v_random_real(pearson_v_random_struct *r);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif	/* VARIATES_H */

