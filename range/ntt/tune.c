/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/

/*
  usage:

  tune max_threads > tunetab.c

  diagnostics written to stderr
*/


#include "ntt-internal.h"
#include <math.h>


int compare_double(const void* a, const void* b)
{
  double aa = * (double*) a;
  double bb = * (double*) b;
  if (aa < bb)
    return -1;
  if (aa > bb)
    return 1;
  return 0;
}


// wall time in seconds for multiplication of size n limbs using GMP
// (which == 0) or NTT with given number of threads (which == 1)
// (not accurate for very small limb count)
double time_func(size_t n, int which, int num_threads)
{
  mp_limb_t* rop = safe_malloc(2 * n * sizeof(mp_limb_t));
  mp_limb_t* op1 = safe_malloc(n * sizeof(mp_limb_t));
  mp_limb_t* op2 = safe_malloc(n * sizeof(mp_limb_t));

  mpn_random(rop, 2 * n);
  mpn_random(op1, n);
  mpn_random(op2, n);

#define REPEAT 3

  double t[REPEAT];

  for (int i = 0; i < REPEAT; i++)
    {
      // increase repeat count until it takes at least 0.1s
      size_t count = 1;
      double len = 0.0;

      while (1)
	{
	  double start = get_time();
	  for (size_t j = 0; j < count; j++)
	    {
	      if (which)
		ntt_mpn_mul_lowmem(rop, op1, n, 1, op2, n, 1, num_threads);
	      else
		mpn_mul_n(rop, op1, op2, n);
	    }
	  double end = get_time();
	  len = end - start;
	  if (len >= 0.1)
	    break;
	  count *= 2;
	}

      t[i] = len / count;
    }

  safe_free(rop);
  safe_free(op1);
  safe_free(op2);

  qsort(t, REPEAT, sizeof(double), compare_double);
  return t[REPEAT / 2];
}


double time_gmp(size_t n)
{
  return time_func(n, 0, 0);
}


double time_ntt(size_t n, int num_threads)
{
  return time_func(n, 1, num_threads);
}


void usage()
{
  fprintf(stderr, "syntax: tune max_threads\n");
  exit(0);
}


int main(int argc, char* argv[])
{
  ntt_init();

  if (argc != 2)
    usage();

  int max_threads = atoi(argv[1]);

  double step = 1.618;
  double min_n = 1e2;
  double max_n = 1e7;

  size_t* threshold = safe_malloc((max_threads + 1) * sizeof(size_t));
  for (int num_threads = 1; num_threads <= max_threads; num_threads++)
    threshold[num_threads] = SIZE_MAX;

  int done = 0;

  for (size_t n = min_n; n < max_n && done < max_threads;
       n = (size_t) ceil(n * step))
    {
      fprintf(stderr, "\nn = %8zu\n", n);
      fflush(stderr);

      double t_gmp = time_gmp(n);
      fprintf(stderr, "  gmp              = %.3le sec\n", t_gmp);
      fflush(stderr);

      for (int num_threads = 1; num_threads <= max_threads; num_threads++)
	if (threshold[num_threads] == SIZE_MAX)
	  {
	    double t_ntt = time_ntt(n, num_threads);
	    fprintf(stderr, "  ntt (%2d threads) = %.3le sec\n",
		    num_threads, t_ntt);
	    if (t_ntt < t_gmp)
	      {
		threshold[num_threads] = n;
		done++;
	      }
	}
    }
  fprintf(stderr, "\n");
  fflush(stderr);

  printf("#include \"ntt-internal.h\"\n");
  printf("\n");
  printf("int tune_tab_size = %d;\n", max_threads + 1);
  printf("\n");
  printf("size_t tune_tab[] = {0");
  for (int num_threads = 1; num_threads <= max_threads; num_threads++)
    {
      if (threshold[num_threads] == SIZE_MAX)
	printf(", SIZE_MAX");
      else
	printf(", %zu", threshold[num_threads]);
    }
  printf("};\n");

  return 0;
}
