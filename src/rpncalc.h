#ifndef RPNCALC_H
#define RPNCALC_H

#include "variates.h"		/* xxx_random_struct */

enum {RPN_OK, RPN_ERROR, RPN_HELP, RPN_QUIT};

/*
  User-sized stack of doubles
 */

typedef struct {
  double *stack;
  double mem;			/* 1-value memory */
  double sumx, sumy, sumxx, sumyy, sumxy, n; /* statistics vars */
  int size;			/* stack[size-1] = last one */
  int next;			/* index of next to push, also num in stack */
  int base;			/* base used for numbers */
  int askprec;			/* asked-for precision */
  int prec;			/* actual precision, <= sig figs */
  int sigfig;			/* max significant figures, from base */
  int angle_unit;		/* 0 = rad, 1 = deg */
  uniform_random_struct urand;
  normal_random_struct nrand;
  exponential_random_struct erand;
} DS;

extern int ds_init(DS *ds, double *stack, int size);
extern int ds_clear(DS *ds);
extern int ds_allclear(DS *ds);
extern int ds_push(DS *ds, double val);
extern int ds_pop(DS *ds, double *val);
extern int ds_swap(DS *ds);
extern int ds_replace(DS *ds, int howmany, double val);
extern int ds_fromtop(DS *ds, int down, double *val);
extern int ds_setbase(DS *ds, int base);
extern int ds_base(DS *ds);
extern int ds_prec(DS *ds);

extern int convert_d_to_s(char *buf, double x, int base, int prec, int len);

/*
  Using a calc you've inited and plan to continue using,
  evaluate an RPN string. The calc's stack will be left intact
  so you can call this again and preserve intermediate results.
  You will need to call rpncalc_pop() to get the result.
 */
extern int rpncalc_eval(DS *ds, char *ptr);

/*
  This is useful if you just want the result once, and don't want to
  create and reuse a calculator. No need to call rpncalc_pop(); this 
  will be done for you and the result stored in val.
*/
extern int rpncalc_eval_full(char *ptr, double *val);

#endif /* RPNCALC_H */

