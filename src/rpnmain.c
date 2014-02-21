/*
  rpnmain.c

  Front end to rpncalc.c, for Reverse Polish Notation (postfix)
  expression evaluation. 

  To Do: enable infix notation.

  Here's is how to convert from an infix notation to postfix notation:

  1. Initialize an empty stack (string stack), prepare input infix
  expression and clear RPN string.

  2. Repeat until we reach end of infix expression:

  I. Get token (operand or operator); skip white spaces.

  II. If token is:

  a. Left parenthesis: Push it into stack.

  b. Right parenthesis: Keep popping from the stack and appending to
  RPN string until we reach the left parenthesis.  If stack becomes
  empty and we didn't reach the left parenthesis then break out with
  error "Unbalanced parenthesis".

  c. Operator: If stack is empty or operator has a higher precedence
  than the top of the stack then push operator into stack. Else if
  operator has lower precedence then we keep popping and appending to
  RPN string, this is repeated until operator in stack has lower
  precedence than the current operator.

  d. An operand: we simply append it to RPN string.

  III. When the infix expression is finished, we start popping off the
  stack and appending to RPN string till stack becomes empty.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_READLINE_READLINE_H && HAVE_LIBREADLINE
#define USE_READLINE 1
#endif

#if HAVE_READLINE_HISTORY_H && HAVE_LIBHISTORY
#define USE_HISTORY 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_READLINE
#include <readline/readline.h>
#endif
#ifdef USE_HISTORY
#include <readline/history.h>
#endif
#include "rpncalc.h"

static void print_help(void)
{
  printf("Use Reverse Polish Notation (RPN), 1 2 + instead of 1 + 2.\n");
  printf("Numbers are pushed onto the stack for use by operators.\n");
  printf("Operators are lower case, numbers are uppercase for bases > 10.\n");
  printf("The stack is shown after each line, left-to-right is bottom-to-top\n");
  printf("Operators (X means top of stack, X Y mean next and top, respectively):\n");
  printf("c            clear (except memory)\n");
  printf("ac           all clear (including memory)\n");
  printf("dec          use decimal base 10\n");
  printf("sto          copy X into memory and drop it\n");
  printf("rcl          push contents of memory onto stack\n");
  printf("sum          add X to memory and drop it\n");
  printf("exc          exchange X with memory\n");
  printf("hex          use hexadecimal base 16\n");
  printf("bin          use binary base 2\n");
  printf("dup          duplicate X\n");
  printf("swap         swap X and Y\n");
  printf("drop         drop X\n");
  printf("depth        push the depth of the stack onto stack\n");

  printf("statistics functions:\n");
  printf("avg          push average of numbers on stack\n");
  printf("std          push std dev of numbers on stack\n");
  printf("stat         X Y ... points go into cumulative statistics\n");
  printf("xstat        X ... singles go into cumulative statistics\n");
  printf("n            push number of stat points\n");
  printf("sx           push sum of x values of the stat points\n");
  printf("sy           push sum of y values of the stat points\n");
  printf("sxx          push sum of squares of x values of the stat points\n");
  printf("syy          push sum of squares of y values of the stat points\n");
  printf("sxy          push sum of products of x and y values of the stat points\n");
  printf("mx           push mean of x values of the stat points\n");
  printf("my           push mean of y values of the stat points\n");
  printf("sdx          push std dev of x values of the stat points\n");
  printf("sdy          push std dev of y values of the stat points\n");
  printf("a            push linear regression 'a' value of ax+b\n");
  printf("b            push linear regression 'b' value of ax+b\n");
  printf("r            push correlation coefficient of linear regression\n");

  printf("sqrt         replace X with its square root\n");
  printf("sq           replace X with its square\n");
  printf("inv          replace X with its inverse, 1/X\n");

  printf("=base        set the base to X\n");
  printf("=prec        set the precision to X\n");
  printf("?base        push the base\n");
  printf("?prec        push the precision\n");
  printf("?sf          push the number of significant figures\n");

  printf(">>           replace X Y with X shifted right by Y\n");
  printf("<<           replace X Y with X shifted left by Y\n");
  printf("&            replace X Y with X bitwise-and Y\n");
  printf("|            replace X Y with X bitwise-or Y\n");
  printf("~            replace X with its bitwise negation\n");

  printf("sin          replace X (in radians) with its sine\n");
  printf("cos          replace X (in radians) with its cosine\n");
  printf("tan          replace X (in radians) with its tangent\n");
  printf("atan2        replace X Y with arctangent(x/y)\n");

  printf("=urand       set uniform random generator (a,b) to X Y\n");
  printf("=nrand       set normal random generator mean, sd to X Y\n");
  printf("=erand       set exponential random generator sd to X\n");
  printf("urand        generate uniform random number using set (a,b)\n");
  printf("nrand        generate normal randomd number using set mean, sd\n");
  printf("erand        generate exponential random number using set sd\n");

  printf("pi           push pi\n");
  printf("e            push e, the base of the natural log\n");
  printf("vc           push speed of light\n");
}

/*
  RPN calculator test example

  Syntax: rpn {<expression>}

  If expression is provided, evaluate this, otherwise read from stdin.
*/

int main(int argc, char *argv[])
{
  enum {BUFFERSIZE = 256};
  char buffer[BUFFERSIZE];
  int bufferleft = BUFFERSIZE;
  char *line;
  DS ds;
  enum {STACKSIZE = 10};
  double stack[STACKSIZE];
  int base;
  int prec;
  int t;
  int retval;

  ds_init(&ds, stack, STACKSIZE);

#ifdef USE_HISTORY
  using_history();
#endif

  *buffer = 0;
  if (argc > 1) {
    for (t = 1; t < argc; t++) {
      if (strlen(argv[t]) < bufferleft - 1) {
	strcat(buffer, argv[t]);
	strcat(buffer, " ");
	bufferleft -= (strlen(argv[t]) + 1);
      }
    }
  }

  do {
    if (argc > 1) {
      line = buffer;
    } else {
#ifdef USE_READLINE
      line = readline("");
      if (NULL == line) {
	break;
      }
#else
      if (NULL == fgets(buffer, BUFFERSIZE, stdin)) {
	return 0;		/* end of file */
      }
      line = buffer;
#endif
#ifdef USE_HISTORY
      if (*line) {
	add_history(line);
      }
#endif
    }

    retval = rpncalc_eval(&ds, line);
    if (RPN_OK == retval) {
      if (ds.next == 0) {
	printf("(empty)\n");
      } else {
	prec = ds_prec(&ds);
	base = ds_base(&ds);
	for (t = 0; t < ds.next; t++) {
	  retval = convert_d_to_s(buffer, ds.stack[t], base, prec, BUFFERSIZE);
	  if (RPN_OK != retval) {
	    printf("error\n");
	  } else {
	    buffer[BUFFERSIZE-1] = 0;
	    printf("%s ", buffer);
	  }
	}
	printf("\n");
      }
    } else if (RPN_ERROR == retval) {
      printf("error\n");
    } else if (RPN_HELP == retval) {
      print_help();
    } else if (RPN_QUIT == retval) {
      break;
    }

    if (argc > 1) {
      break;
    } else {
#ifdef USE_READLINE
      free(line);
#endif
    }
  } while (! feof(stdin));

  return RPN_ERROR == retval ? 1 : 0;
}
