#include <math.h>		/* all math functions */
#include <float.h>		/* DBL_MIN */
#include <errno.h>		/* errno */
#include "rpncalc.h"		/* our decls */
#include "ptime.h"		/* ptime */
#include "variates.h"		/* uniform_random, ... */

/*
  Reverse Polish Notation calculator.

  To do:

  Exponential notation, e.g., 1.23e4. See FIXME note.

  rot, ndup, etc. stack operations

  Expressions, e.g., 
  expr +
  Then "1 2" will implicity tack on expr, yielding 3. expr <blank> clears.
  This lets you pipe in columns of numbers and compute a value.
*/

/* useful constants */
static const double CONST_E =  2.7182818284590452354;
static const double CONST_PI = 3.1415926535897932385;
/* 1/ln(10), for computing log10(x) = ln(x)/ln(10) */
static const double CONST_LN10_INV = 0.43429448190325182765;
static const double CONST_SPEED_OF_LIGHT = 299792458;

/* useful conversion factors */
static const double CONV_MI_TO_M = 5280.0 * 12.0 * 0.0254;
static const double CONV_FT_TO_M = 12.0 * 0.0254;
static const double CONV_IN_TO_MM = 25.4;

#define TODEG(x) ((x) * 57.295779513082320875)
#define TORAD(x) ((x) * 0.017453292519943295770)
#define TOFARENHEIT(c) (9.0/5.0*(c)+32.0)
#define TOCELSIUS(f) (5.0/9.0*((f)-32.0))

#define DOUBLE_BITS 53		/* bits in the fraction part of a double */

static int sigfig(int base)
{
  return (int) (DOUBLE_BITS * log(2)/log(base));
}

int ds_init(DS *ds, double *stack, int size)
{
  if (size <= 0) return RPN_ERROR;

  ds->stack = stack;
  ds->mem = 0.0;
  ds->sumx = ds->sumy = ds->sumxx = ds->sumyy = ds->sumxy = ds->n = 0.0;
  ds->size = size;
  ds->next = 0;
  ds->base = 10;
  ds->sigfig = sigfig(ds->base);
  ds->askprec = ds->sigfig;
  ds->prec = ds->sigfig;
  ds->angle_unit = 0;		/* radians */
  uniform_random_init(&ds->urand, 0, 1);
  normal_random_init(&ds->nrand, 0, 1);
  exponential_random_init(&ds->erand, 1);

  return RPN_OK;
}

int ds_clear(DS *ds)
{
  ds->next = 0;

  return RPN_OK;
}

/* clear memory, but leave base, degrees, etc. alone */
int ds_allclear(DS *ds)
{
  ds->mem = 0.0;
  ds->sumx = ds->sumy = ds->sumxx = ds->sumyy = ds->sumxy = ds->n = 0.0;

  return ds_clear(ds);
}

int ds_push(DS *ds, double val)
{
  if (ds->next == ds->size) {
    /* full */
    return RPN_ERROR;
  }

  ds->stack[ds->next++] = val;

  return RPN_OK;
}

int ds_pop(DS *ds, double *val)
{
  if (ds->next == 0) {
    /* empty */
    return RPN_ERROR;
  }

  *val = ds->stack[--ds->next];

  return RPN_OK;
}

int ds_dup(DS *ds)
{
  if (ds->next == 0 ||
      ds->next == ds->size) {
    /* empty or full */
    return RPN_ERROR;
  }

  ds->stack[ds->next] = ds->stack[ds->next - 1];
  ds->next++;

  return RPN_OK;
}

int ds_swap(DS *ds)
{
  double temp;

  if (ds->next < 2) {
    /* don't have 2 to swap */
    return RPN_ERROR;
  }

  temp = ds->stack[ds->next - 1];
  ds->stack[ds->next - 1] = ds->stack[ds->next - 2];
  ds->stack[ds->next - 2] = temp;

  return RPN_OK;
}

int ds_rot(DS *ds)
{
  double temp;
  int t;

  if (ds->next < 2) {
    /* the rot of empty or one element is just the same stack */
    return RPN_OK;
  }

  temp = ds->stack[0];
  for (t = 0; t < ds->next - 1; t++) {
    ds->stack[t] = ds->stack[t + 1];
  }
  ds->stack[ds->next -1] = temp;

  return RPN_OK;
}

int ds_drop(DS *ds)
{
  if (ds->next == 0) {
    /* empty */
    return RPN_ERROR;
  }

  ds->next--;

  return RPN_OK;
}

int ds_replace(DS *ds, int howmany, double val)
{
  if (ds->next < howmany) {
    /* too few to replace */
    return RPN_ERROR;
  }

  ds->next -= (howmany - 1);
  ds->stack[ds->next - 1] = val;

  return RPN_OK;
}

/*
  ds_fromtop() puts stack element 'down' from top of stack into 'val'.
  If down is 0, val gets the top of the stack. Stack is unchanged.
 */
int ds_fromtop(DS *ds, int down, double *val)
{
  if (down < 0 ||
      down >= ds->next) {
    /* not enough on stack */
    return RPN_ERROR;
  }

  *val = ds->stack[ds->next - 1 - down];

  return RPN_OK;
}

int ds_setbase(DS *ds, int base)
{
  if (base < 2) return RPN_ERROR;
  if (base > 36) return RPN_ERROR;

  ds->base = base;
  ds->sigfig = sigfig(base);
  ds->prec = ds->askprec;
  if (ds->prec > ds->sigfig) ds->prec = ds->sigfig;

  return RPN_OK;
}

int ds_base(DS *ds)
{
  return ds->base;
}

int ds_setprec(DS *ds, int prec)
{
  ds->askprec = prec;
  ds->prec = prec < 0 ? 0 : prec > ds->sigfig ? ds->sigfig : prec;

  return RPN_OK;
}

int ds_prec(DS *ds)
{
  return ds->prec;
}

/*
  Return sample standard deviation,

  .     1
  sqrt(--- sum((xi - mean)^2))
  .    N-1
*/

double ds_stddev(DS *ds)
{
  double sumx = 0.0, sumxx = 0.0, mean;
  int t;

  if (ds->next < 2) return 0.0;

  for (t = 0; t < ds->next; t++) {
    sumx += ds->stack[t];
    sumxx += ds->stack[t] * ds->stack[t];
  }
  mean = sumx / ds->next;

  return sqrt((sumxx - 2.0*mean*sumx + ds->next*mean*mean) / (ds->next-1));
}

double ds_stddev_x(DS *ds)
{
  double mean;

  if (ds->n < 2) return 0.0;

  mean = ds->sumx / ds->n;

  return sqrt((ds->sumxx - 2.0*mean*ds->sumx + ds->n*mean*mean) / (ds->n-1));
}

double ds_stddev_y(DS *ds)
{
  double mean;

  if (ds->n < 2) return 0.0;

  mean = ds->sumy / ds->n;

  return sqrt((ds->sumyy - 2.0*mean*ds->sumy + ds->n*mean*mean) / (ds->n-1));
}

/*
  Given xy points in stat regs, computes linear regression coefficients

  y = ax + b, correlation coeff r
 */

double ds_leastsq_a(DS *ds)
{
  double denom;

  denom = ds->n*ds->sumxx - ds->sumx*ds->sumx;

  if (denom == 0.0) return 0.0;

  return (ds->n*ds->sumxy - ds->sumx*ds->sumy) / denom;
}

double ds_leastsq_b(DS *ds)
{
  double denom;

  denom = ds->n*ds->sumxx - ds->sumx*ds->sumx;

  if (denom == 0.0) return 0.0;

  return (ds->sumxx*ds->sumy - ds->sumx*ds->sumxy) / denom;
}

double ds_leastsq_r(DS *ds)
{
  double denom;

  denom = (ds->n*ds->sumxx - ds->sumx*ds->sumx) * (ds->n*ds->sumyy - ds->sumy*ds->sumy);

  if (denom <= 0.0) return 0.0;

  return (ds->n*ds->sumxy - ds->sumx*ds->sumy) / sqrt(denom);
}

#define isspace(c) ((c) == ' ' || (c) == '\t' || (c) == '\n' || (c) == '\r')
#define isnullspace(c) ((c) == 0 || (c) == ' ' || (c) == '\t' || (c) == '\n' || (c) == '\r')

/*
  Returns 0 if chars of "s" match chars of "c", up to null
  or whitespace, e.g., these match:

  "atan2 hello" "atan2"
  "atan2 hello" "atan2 "
  "atan2" "atan2 hello"

  and these don't:

  "atan2" "atan"
  "atan" "atan2"
  "atan2" "atan hello"

  Used in place of hashes when a hashed values collide. Not yet used.
*/

#if 0
static int strspcmp(char *s, char *c)
{
  while (*s == *c) {
    if (*s == 0) return 0;	/* they match completely */
    s++, c++;
  }

  if (*s == 0 && isspace(*c)) return 0;
  if (*c == 0 && isspace(*s)) return 0;
  return 1;
}
#endif

int compute_hash(char *buffer)
{
  int hash;
  int len;

  if (isnullspace(buffer[0])) {
    return 0;
  }
  hash = ((int) buffer[0]) << 24;

  if (isnullspace(buffer[1])) {
    return hash + 1;
  }
  hash += ((int) buffer[1]) << 16;

  if (isnullspace(buffer[2])) {
    return hash + 2;
  }
  hash += ((int) buffer[2]) << 8;

  if (isnullspace(buffer[3])) {
    return hash + 3;
  }

  len = 4;
  while (! isnullspace(buffer[len])) len++;
  hash += len;

  return hash;
}

#define compute_hash_1(a) ((((int) (a)) << 24) + 1)
#define compute_hash_2(a,b) ((((int) (a)) << 24) + (((int) (b)) << 16) + 2)
#define compute_hash_3(a,b,c) ((((int) (a)) << 24) + (((int) (b)) << 16) + (((int) (c)) << 8) + 3)
#define compute_hash_n(a,b,c,n) ((n) == 3 ? compute_hash_3(a,b,c) : (((int) (a)) << 24) + (((int) (b)) << 16) +(((int) (c)) << 8) + n)

int decompute_hash(int hash, char *buffer)
{
  int len;
  int t;

  len = hash & 0xFF;
  buffer[2] = (hash >> 8) & 0xFF;
  buffer[1] = (hash >> 16) & 0xFF;
  buffer[0] = (hash >> 24) & 0xFF;

  for (t = 3; t < len; t++) {
    buffer[t] = 'X';
  }
  buffer[len] = 0;

  return 0;
}

#define round(x) (x) < 0 ? (int) ((x) - 0.5) : (int) ((x) + 0.5)

static int factorial(double x, double *f)
{
  double cum = x;
  int r = round(x);

  /* check for integer */
  if (fabs(x - r) > DBL_MIN) return RPN_ERROR;

  if (r < 0) return RPN_ERROR;
  if (r < 2) {
    *f = 1;
    return RPN_OK;
  }

  while (r-- > 2) cum *= r;
  *f = cum;
  return RPN_OK;
}

static int rpncalc_op(DS *ds, char *op)
{
  double val;
  double top, next;
  int hash;
  int t;
  int i, j;

  hash = compute_hash(op);

  switch (hash) {
  case compute_hash_1('c'):	/* c */
    return ds_clear(ds);
  case compute_hash_2('a','c'):	/* ac */
    return ds_allclear(ds);
  case compute_hash_3('d','e','c'): /* dec */
    return ds_setbase(ds, 10);
  case compute_hash_3('h','e','x'): /* hex */
    return ds_setbase(ds, 16);
  case compute_hash_3('b','i','n'): /* bin */
    return ds_setbase(ds, 2);
  case compute_hash_3('d','u','p'): /* dup */
    return ds_dup(ds);
  case compute_hash_n('s','w','a',4): /* swap */
    return ds_swap(ds);
  case compute_hash_3('r','o','t'): /* rot */
    return ds_rot(ds);
  case compute_hash_n('d','r','o',4): /* drop */
  case compute_hash_1('.'):	/* . short for drop */
    return ds_drop(ds);
  case compute_hash_n('d','e','p',5): /* depth */
    return ds_push(ds, ds->next);
  case compute_hash_3('a','v','g'): /* avg */
    if (0 == ds->next) return RPN_ERROR;
    val = 0;
    for (t = 0; t < ds->next; t++) {
      val += ds->stack[t];
    }
    return ds_push(ds, val / ds->next);
  case compute_hash_3('s','t','d'): /* std */
    if (0 == ds->next) return RPN_ERROR;
    return ds_push(ds, ds_stddev(ds));
  case compute_hash_n('s','t','a',4): /* stat */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    if (ds->next % 2) return RPN_ERROR;
    for (t = 0; t < ds->next; t += 2) {
      ds->sumx += ds->stack[t];
      ds->sumy += ds->stack[t + 1];
      ds->sumxx += ds->stack[t] * ds->stack[t];
      ds->sumyy += ds->stack[t + 1] * ds->stack[t + 1];
      ds->sumxy += ds->stack[t] * ds->stack[t + 1];
      ds->n++;
    }
    ds->next = 0;
    return RPN_OK;
  case compute_hash_1('n'):	/* n, number of stat points */
    return ds_push(ds, ds->n);
  case compute_hash_2('s','x'): /* sx, sum of x */
    return ds_push(ds, ds->sumx);
  case compute_hash_2('s','y'): /* sy, sum of y */
    return ds_push(ds, ds->sumy);
  case compute_hash_3('s','x','x'): /* sxx, sum of x^2 */
    return ds_push(ds, ds->sumxx);
  case compute_hash_3('s','y','y'): /* syy, sum of y^2 */
    return ds_push(ds, ds->sumyy);
  case compute_hash_3('s','x','y'): /* sxy, sum of x*y */
    return ds_push(ds, ds->sumxy);
  case compute_hash_2('m','x'):	/* mx, mean of x */
    return ds->n == 0 || ds_push(ds, ds->sumx / ds->n);
  case compute_hash_2('m','y'):	/* my, mean of y */
    return ds->n == 0 || ds_push(ds, ds->sumy / ds->n);
  case compute_hash_3('s','d','x'): /* sdx, stddev of x */
    return ds->n < 2 || ds_push(ds, ds_stddev_x(ds));
  case compute_hash_3('s','d','y'): /* sdy, stddev of y */
    return ds->n < 2 || ds_push(ds, ds_stddev_y(ds));
  case compute_hash_1('a'):	/* a in linear regression ax+b */
    return ds->n < 2 || ds_push(ds, ds_leastsq_a(ds));
  case compute_hash_1('b'):	/* b in linear regression ax+b */
    return ds->n < 2 || ds_push(ds, ds_leastsq_b(ds));
  case compute_hash_1('r'):	/* r, correlation coefficient */
    return ds->n < 2 || ds_push(ds, ds_leastsq_r(ds));
  case compute_hash_n('=','b','a',5): /* =base */
    return ds_fromtop(ds, 0, &top) == RPN_OK ? ds_drop(ds), ds_setbase(ds, top) : RPN_ERROR;
  case compute_hash_n('=','p','r',5): /* =prec */
    return ds_fromtop(ds, 0, &top) == RPN_OK ? ds_drop(ds), ds_setprec(ds, top) : RPN_ERROR;
  case compute_hash_n('?','b','a',5): /* ?base */
    return ds_push(ds, ds->base);
  case compute_hash_n('?','p','r',5): /* ?prec */
    return ds_push(ds, ds->prec);
  case compute_hash_3('?','s','f'): /* ?sf, significant figures */
    return ds_push(ds, ds->sigfig);
  case compute_hash_3('s','t','o'): /* sto, X->M */
    return ds_fromtop(ds, 0, &top) == RPN_OK ? ds_drop(ds), ds->mem = top, 0 : RPN_ERROR;
  case compute_hash_3('r','c','l'): /* rcl, RCL */
    return ds_push(ds, ds->mem);
  case compute_hash_3('s','u','m'): /* sum, M+ */
    return ds_fromtop(ds, 0, &top) == RPN_OK ? ds_drop(ds), ds->mem += top, 0 : RPN_ERROR;
  case compute_hash_3('e','x','c'): /* exc, EXC */
    return ds_fromtop(ds, 0, &top) == RPN_OK ? ds_replace(ds, 1, ds->mem), ds->mem = top, 0 : RPN_ERROR;
  case compute_hash_2('-','+'):	/* -+ */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, -top);
  case compute_hash_2('+','-'):	/* +- */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, -top);
  case compute_hash_3('i','n','v'): /* inv */
    return ds_fromtop(ds, 0, &top) || ! (fabs(top) > DBL_MIN) || ds_replace(ds, 1, 1.0 / top);
  case compute_hash_2('s','q'):	/* sq */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, top * top);
  case compute_hash_n('s','q','r',4): /* sqr */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, sqrt(top));
  case compute_hash_1('!'):	/* ! */
    return ds_fromtop(ds, 0, &top) || factorial(top, &val) || ds_replace(ds, 1, val);
  case compute_hash_3('s','i','n'): /* sin */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, sin(ds->angle_unit == 0 ? top : TORAD(top)));
  case compute_hash_3('c','o','s'): /* cos */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, cos(ds->angle_unit == 0 ? top : TORAD(top)));
  case compute_hash_3('t','a','n'): /* tan */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, tan(ds->angle_unit == 0 ? top : TORAD(top)));
  case compute_hash_n('s','i','n',4): /* sinh */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    errno = 0;
    val = sinh(top);
    return errno != 0 || ds_replace(ds, 1, val);
  case compute_hash_n('c','o','s',4): /* cosh */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    errno = 0;
    val = cosh(top);
    return errno != 0 || ds_replace(ds, 1, val);
  case compute_hash_n('t','a','n',4): /* tanh */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    errno = 0;
    val = tanh(top);
    return errno != 0 || ds_replace(ds, 1, val);
  case compute_hash_n('a','s','i',4): /* asin */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    errno = 0;
    val = asin(top);
    return errno != 0 || ds_replace(ds, 1, ds->angle_unit == 0 ? val : TODEG(val));
  case compute_hash_n('a','c','o',4): /* acos */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    errno = 0;
    val = acos(top);
    return errno != 0 || ds_replace(ds, 1, ds->angle_unit == 0 ? val : TODEG(val));
  case compute_hash_n('a','t','a',4): /* atan */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    val = atan(top);
    return ds_replace(ds, 1, ds->angle_unit == 0 ? val : TODEG(val));
  case compute_hash_n('a','t','a',5): /* atan2 */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    val = atan2(next, top);
    return ds_replace(ds, 2, ds->angle_unit == 0 ? val : TODEG(val));
  case compute_hash_3('e','x','p'): /* exp */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, exp(top));
  case compute_hash_2('l','n'):	/* ln */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    errno = 0;
    val = log(top);
    return errno != 0 || ds_replace(ds, 1, val);
  case compute_hash_3('l','o','g'): /* log */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    errno = 0;
    val = log(top);
    if (errno) return RPN_ERROR;
    val *= CONST_LN10_INV;
    return ds_replace(ds, 1, val);
  case compute_hash_n('l','o','g',4): /* logn */
    return ds_fromtop(ds, 1, &next) || ds_fromtop(ds, 0, &top) || next <= 0.0 || top <= 0.0 || ds_replace(ds, 2, log(next)/log(top));
  case compute_hash_3('a','b','s'): /* abs */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, fabs(top));
  case compute_hash_n('r','o','u',5): /* round */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, round(top));
  case compute_hash_3('r', 'a', 'd'): /* rad */
    ds->angle_unit = 0;
    return RPN_OK;
  case compute_hash_3('d','e','g'): /* deg */
    ds->angle_unit = 1;
    return RPN_OK;
  case compute_hash_n('t','o','d',5): /* todeg */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, TODEG(top));
  case compute_hash_n('t','o','r',5): /* torad */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, TORAD(top));
  case compute_hash_3('t','o','f'): /* tof, to farenheit */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, TOFARENHEIT(top));
  case compute_hash_3('t','o','c'): /* toc, to celsius */
    return ds_fromtop(ds, 0, &top) || ds_replace(ds, 1, TOCELSIUS(top));
  case compute_hash_n('x','s','t',5): /* xstat */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    for (t = 0; t < ds->next; t++) {
      ds->sumx += ds->n;
      ds->sumy += ds->stack[t];
      ds->sumxx += ds->n * ds->n;
      ds->sumyy += ds->stack[t] * ds->stack[t];
      ds->sumxy += ds->n * ds->stack[t];
      ds->n++;
    }
    ds->next = 0;
    return RPN_OK;
  case compute_hash_1('+'):	/* + */
    return ds_fromtop(ds, 1, &next) || ds_fromtop(ds, 0, &top) || ds_replace(ds, 2, next + top);
  case compute_hash_1('-'):	/* - */
    return ds_fromtop(ds, 1, &next) || ds_fromtop(ds, 0, &top) || ds_replace(ds, 2, next - top);
  case compute_hash_1('*'):	/* * */
  case compute_hash_1('x'):	/* x */
    return ds_fromtop(ds, 1, &next) || ds_fromtop(ds, 0, &top) || ds_replace(ds, 2, next * top);
  case compute_hash_1('/'):	/* / */
    return ds_fromtop(ds, 1, &next) || ds_fromtop(ds, 0, &top) || ! (fabs(top) > DBL_MIN) || ds_replace(ds, 2, next / top);
  case compute_hash_3('d','i','v'): /* div */
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    i = round(next), j = round(top);
    if (j == 0) return RPN_ERROR;
    return ds_replace(ds, 2, i / j);
  case compute_hash_3('m','o','d'): /* mod */
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    i = round(next), j = round(top);
    if (j == 0) return RPN_ERROR;
    return ds_replace(ds, 2, i % j);
  case compute_hash_n('f','m','o',4): /* fmod */
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    return ds_replace(ds, 2, fmod(next, top));
  case compute_hash_n('f','l','o',5): /* floor */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    return ds_replace(ds, 1, floor(top));
  case compute_hash_n('c','e','i',4): /* ceil */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    return ds_replace(ds, 1, ceil(top));
  case compute_hash_3('p','o','w'): /* pow */
  case compute_hash_1('^'):	/* ^ short for pow */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    errno = 0;
    val = pow(next, top);
    return errno || ds_replace(ds, 2, val);
  case compute_hash_2('>','>'):	      /* >> */
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    i = round(next), j = round(top);
    return ds_replace(ds, 2, i >> j);
  case compute_hash_2('<','<'):	      /* >> */
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    i = round(next), j = round(top);
    return ds_replace(ds, 2, i << j);
  case compute_hash_1('|'):	      /* | */
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    i = round(next), j = round(top);
    return ds_replace(ds, 2, i | j);
  case compute_hash_1('&'):	      /* | */
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    i = round(next), j = round(top);
    return ds_replace(ds, 2, i & j);
  case compute_hash_1('~'):	      /* | */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    i = round(top);
    return ds_replace(ds, 1, ~i);
  case compute_hash_n('t','o','x',4): /* toxy, r-theta to x-y conversion */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    ds->next -= 2;
    return ds_push(ds, next * cos(ds->angle_unit == 0 ? top : TORAD(top))) || ds_push(ds, next * sin(ds->angle_unit == 0 ? top : TORAD(top)));
  case compute_hash_n('t','o','r',4): /* tort, x-y to r-theta conversion */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    ds->next -= 2;
    val = atan2(top, next);
    return ds_push(ds, sqrt(next*next + top*top)) || ds_push(ds, ds->angle_unit == 0 ? val : TODEG(val));

  case compute_hash_n('m','i','2',4): /* mi2m */
    return ds_push(ds, CONV_MI_TO_M);
  case compute_hash_n('f','t','2',4): /* ft2m */
    return ds_push(ds, CONV_FT_TO_M);
  case compute_hash_n('i','n','2',5): /* in2mm */
    return ds_push(ds, CONV_IN_TO_MM);

    /* time */
  case compute_hash_n('t','i','m',4): /* time */
    return ds_push(ds, ptime());

    /* random variates */
  case compute_hash_n('=','u','r',6): /* =urand, set a and b */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    ds->next -= 2;
    uniform_random_set(&ds->urand, next, top);
    return RPN_OK;
  case compute_hash_n('=','n','r',6): /* =nrand, set mean and sd */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    if (ds_fromtop(ds, 1, &next)) return RPN_ERROR;
    ds->next -= 2;
    normal_random_set(&ds->nrand, next, top);
    return RPN_OK;
  case compute_hash_n('=','e','r',6): /* =erand, set sd */
    if (ds_fromtop(ds, 0, &top)) return RPN_ERROR;
    ds->next -= 1;
    exponential_random_set(&ds->erand, top);
    return RPN_OK;
  case compute_hash_n('u','r','a',5): /* urand, uniform random number */
    return ds_push(ds, uniform_random_real(&ds->urand));
  case compute_hash_n('n','r','a',5): /* nrand, normal random number */
    return ds_push(ds, normal_random_real(&ds->nrand));
  case compute_hash_n('e','r','a',5): /* erand, exponential random number */
    return ds_push(ds, exponential_random_real(&ds->erand));

    /* useful constants */
  case compute_hash_2('p','i'):	/* pi */
    return ds_push(ds, CONST_PI);
  case compute_hash_1('e'):	/* e */
    return ds_push(ds, CONST_E);
  case compute_hash_2('v','c'):	/* vc, speed of light */
    return ds_push(ds, CONST_SPEED_OF_LIGHT);
  }

  return RPN_ERROR;
}

int isdigitbase(char digit, int base)
{
  if (base <= 10) {
    if (digit >= '0' && digit - '0' < base) return 1;
    return 0;
  }

  if (base <= 36) {
    if (digit >= '0' && digit <= '9') return 1;
    if (digit >= 'A' && digit - 'A' + 10 < base) return 1;
    return 0;
  }

  return 0;			/* can't handle > base 36 */
}

double todoublebase(char digit, int base)
{
  if (base <= 10) {
    if (digit >= '0' && digit - '0' < base) return (double) (digit - '0');
    return 0.0;
  }

  if (base <= 36) {
    if (digit >= '0' && digit <= '9') return (double) (digit - '0');
    if (digit >= 'A' && digit - 'A' + 10 < base) return (double) (digit - 'A' + 10);
    return 0.0;
  }

  return 0.0;			/* can't handle > base 36 */
}

char tocharbase(int digit, int base)
{
  if (digit < 0) return '0';

  if (base <= 10) {
    if (digit < base) return digit + '0';
    return '0';
  }

  if (base <= 36) {
    if (digit < 10) return digit + '0';
    if (digit < base) return digit + 'A' - 10;
  }

  return '0';
}

static char *skipwhite(char *buffer)
{
  while (isspace(*buffer)) buffer++;

  return buffer;
}

static char *skipnonwhite(char *buffer)
{
  while (! isspace(*buffer) && *buffer != 0) buffer++;

  return buffer;
}

/*
  a replacement for strtod, sort of
 */
static int convert_s_to_d(char *ptr, double *x, int base)
{
  double num = 0.0, fracnum;
  int started = 0;
  int gotnum = 0;
  int infrac = 0;
  int minus = 0;
  int inexp = 0;
  int t;
  char c;

  if (0 == *(ptr = skipwhite(ptr))) return RPN_ERROR;

  while (0 != (c = *ptr++) &&
	 !isspace(c)) {
    if (c == '-') {
      if (started) return RPN_ERROR;	/* already in number */
      minus = 1, started = 1;
      continue;
    }

    if (c == '+') {
      if (started) return RPN_ERROR;	/* already in number */
      started = 1;
      continue;
    }

    if (c == 'e') {
      if (! started || inexp) return RPN_ERROR;
      inexp = 1;
    }

    /*
     FIXME-- you added 'inexp' flag, now build cumulative power of base
     as multiplier below, and then use it
    */

    if (c == '.') {
      if (infrac) return RPN_ERROR;	/* already in fraction */
      infrac = 1, started = 1;
      continue;
    }

    if (isdigitbase(c, base)) {
      started = 1;
      if (infrac) {
	/*
	  FIXME-- a better way to do this is to work from the end toward
	  the radix point, and divide cumulative number once each digit.
	  This saves a bunch of divisions. 
	 */
	t = infrac;
	fracnum = todoublebase(c, base);
	while (t-- > 0) {
	  fracnum /= (double) base;
	}
	num += fracnum;
	infrac++;
      } else {
	num *= (double) base;
	num += todoublebase(c, base);
      }
      gotnum = 1;
      continue;
    }

    /* else it's a bad character */
    return RPN_ERROR;
  }

  *x = minus ? -num : num;

  return (gotnum == 1 ? RPN_OK : RPN_ERROR);
}

int convert_d_to_s(char *buf, double x, int base, int prec, int n)
{
  double base_to_prec;
  double roundinc;
  double frac;
  double whole;
  int digit;
  int t;
  char c;
  char *ptr;
  char *wholestart;
  char *point;
  char *lastzero;
  char temp;

  ptr = buf;

  /* set the negative sign and make x nonnegative for subsequent work */
  if (x < 0.0) {
    *ptr++ = '-', n--;
    if (n <= 0) {*ptr = 0; return RPN_ERROR;}
    x = -x;
  } 
  wholestart = ptr;

  /* round to the last digit by adding half of the value of the last digit */
  if (prec < 0) prec = 0;
  base_to_prec = 1.0;
  for (t = 0; t < prec; t++) {
    base_to_prec *= base;
  }
  roundinc = 0.5 / base_to_prec;
  x += roundinc;

  whole = floor(x);
  if (whole <= 0.0) {
    /* print a leading zero for numbers < 1 */
    *ptr++ = '0', n--;
    if (n <= 0) {*ptr = 0; return RPN_ERROR;}
  }
  frac = x - whole;
  while (whole >= 1.0) {
    digit = (int) (whole - (floor(whole / base) * base));
    *ptr++ = tocharbase(digit, base), n--;
    if (n <= 0) {*ptr = 0; return RPN_ERROR;}
    whole /= base;
    prec--;			/* TEST */
  }

  point = ptr;			/* set spot for decimal point */
  /* swap start with end of integer part, moving inward */
  ptr--;
  while (wholestart < ptr) {
    temp = *wholestart;
    *wholestart = *ptr;
    *ptr = temp;
    wholestart++, ptr--;
  }
  ptr = point;

  if (prec <= 0) {		/* TEST-- was == 0 */
    *ptr = 0;
    return RPN_OK;
  }
  lastzero = ptr;
  *ptr++ = '.', n--;
  if (n <= 0) {*ptr = 0; return RPN_ERROR;}

  while (prec-- > 0) {
    frac *= base;
    digit = (int) frac;
    frac = frac - digit;
    c = tocharbase(digit, base);
    *ptr++ = c, n--;
    if (c != '0') lastzero = ptr;
    if (n <= 0) {*lastzero = 0; return RPN_ERROR;}
  }
  *lastzero = 0;

  return RPN_OK;
}

/*
  this is used to evaluate an incremental calc, where the stack
  is preserved between calls. To get the value out, do an
  rpncalc_pop() when you want the final value
 */
int rpncalc_eval(DS *ds, char *ptr)
{
  double x;
  int base;

  while (0 != *(ptr = skipwhite(ptr))) {
    /*
      Handle operators first, then numbers. Since the operator names
      are all lower case and the interpreter is case-sensitive, to
      input numbers that may be confused with operators (e.g., dec),
      numbers should be uppercase, e.g., DEC for 0xDEC.
    */
    if ('?' == ptr[0] && (0 == ptr[1] || isspace(ptr[1]))) return RPN_HELP;
    if ('q' == ptr[0] && (0 == ptr[1] || isspace(ptr[1]))) return RPN_QUIT;
    base = ds_base(ds);
    if (0 != rpncalc_op(ds, ptr)) {
      if (0 == convert_s_to_d(ptr, &x, base)) {
	/* it's a number, so push it */
	ds_push(ds, x);
      } else {
	/* it's not an operator or number */
	return RPN_ERROR;
      }
    }
    /* else it's an operator, we just handled it */

    /* go on to the next one */
    ptr = skipnonwhite(ptr);
  }

  return RPN_OK;
};

/*
  this is useful if you just want the result once, and don't want to
  create and reuse a calculator. No rpncalc_pop() is necessary since
  it's done for you and the result stored in val
*/
int rpncalc_eval_full(char *ptr, double *val)
{
  DS ds;
  double stack[10];

  ds_init(&ds, stack, sizeof(stack));

  return rpncalc_eval(&ds, ptr) || ds_pop(&ds, val);
}
