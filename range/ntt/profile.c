/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


double get_time()
{
  struct timeval time;
  gettimeofday(&time, NULL);
  return time.tv_sec + time.tv_usec / 1000000.0;
}


#if PROFILE

void prof_init(prof_t* prof)
{
  prof->count = 0;
}


void prof_clear(prof_t* prof)
{
}


void prof_event(prof_t* prof, char* label)
{
  assert(prof->count < MAX_PROF);
  double time = get_time();
  prof->labels[prof->count] = label;
  prof->times[prof->count] = time;
  prof->count++;
}


void prof_report(prof_t* prof)
{
  for (size_t i = 1; i < prof->count; i++)
    {
      double diff = prof->times[i] - prof->times[i - 1];
      printf("%7.2lfs: %s\n", diff, prof->labels[i - 1]);
    }
  printf("%7.2lfs: TOTAL\n", prof->times[prof->count - 1] - prof->times[0]);
}

#endif
