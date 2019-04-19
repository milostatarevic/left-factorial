/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


void fatal_error(char* fmt, ...)
{
  va_list list;
  va_start(list, fmt);
  fprintf(stderr, "NTT: ");
  vfprintf(stderr, fmt, list);
  fprintf(stderr, "\n");
  va_end(list);
  exit(1);
}


uint64_t global_p[MAX_NUM_PRIMES] = {
  PRIME0,
  PRIME1,
  PRIME2,
  PRIME3,
  PRIME4,
  PRIME5,
  //  PRIME6,
  //  PRIME7
};

uint64_t global_pinv[MAX_NUM_PRIMES];
uint64_t global_pinvb[MAX_NUM_PRIMES];
uint64_t global_w[MAX_NUM_PRIMES];



void ntt_init()
{
  global_bit_bound[0] = 0.0;

  // compute global precomputed data
  for (unsigned i = 0; i < MAX_NUM_PRIMES; i++)
    {
      uint64_t p = global_p[i];
      global_pinvb[i] = calc_pinvb(p);
      global_pinv[i] = calc_pinv(p);

      // compute global_w
      uint64_t w;
      for (uint64_t x = 2; ; x++)
	{
	  w = pow_mod(x, (p-1) / MAX_K, p);
	  if (pow_mod(w, MAX_K / 2, p) != 1 && pow_mod(w, MAX_K / 3, p) != 1)
	    break;
	}
      global_w[i] = w;

      // compute global_crt_* stuff
      if (i >= 1)
	{
	  if (i == 1)
	    global_crt_u[1][0] = global_p[0];
	  else
	    global_crt_u[i][i-1] = mpn_mul_1(global_crt_u[i], global_crt_u[i-1],
					     i-1, global_p[i-1]);
	  uint64_t s = mpn_mod_1(global_crt_u[i], i, p);
	  s = inv_mod(s, p);
	  s = sub_mod(0, s, p);
	  s = mul_mod(s, pow_mod(2, 64*(i-1), p), p);
	  global_crt_s[i] = s;
	  global_crt_spinv[i] = calc_ypinv(s, p, calc_pinv(p));
	}

      // compute global_bit_bound
      global_bit_bound[i + 1] = global_bit_bound[i] +
	0.999 * log((double) p) / log(2.0);
    }
}


uint64_t calc_w(unsigned i, size_t K)
{
  return pow_mod_pinv(global_w[i], MAX_K / K, global_p[i], global_pinv[i]);
}
