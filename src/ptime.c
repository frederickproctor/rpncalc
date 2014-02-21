#include "ptime.h"

#ifdef WIN32

#include <windows.h>
double ptime(void)
{
  FILETIME f;
  ULARGE_INTEGER t;
  GetSystemTimeAsFileTime(&f);
  t.QuadPart = f.dwHighDateTime;
  t.QuadPart <<= 32;
  t.QuadPart |= f.dwLowDateTime;
  /* now time is in 100 nanosecond units */
  return ((double) t.QuadPart) * 1.0e-7;
}

#else

#include <time.h>
double ptime(void)
{
  struct timespec tv;

  if (0 == clock_gettime(CLOCK_REALTIME, &tv)) {
    return ((double) tv.tv_sec) + ((double) tv.tv_nsec) * 1.0e-9;
  }
  return 0.0;
}
#endif

