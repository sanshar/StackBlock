#include <stdio.h>
#include <iostream>
#include "pario.h"
#include <stdio.h>
#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>


void print_trace(int nSig)
{
  printf("print_trace: got signal %d\n", nSig);

  void           *array[32];    /* Array to store backtrace symbols */
  size_t          size;     /* To store the exact no of values stored */
  char          **strings;    /* To store functions from the backtrace list in ARRAY */
  size_t          nCnt;

  size = backtrace(array, 32);

  strings = backtrace_symbols(array, size);

  /* prints each string of function names of trace*/
  for (nCnt = 0; nCnt < size; nCnt++)
    fprintf(stderr, "%s\n", strings[nCnt]);


  abort();
}

#ifdef MOLPRO
#include "global/CxOutputStream.h"
blockout Bout(&xout);
blockerr Berr(&xerr);
#else
blockout Bout;
blockerr Berr;
#endif

std::ostream &bout = *(Bout.outstream);
std::ostream &berr = *(Berr.errstream);
