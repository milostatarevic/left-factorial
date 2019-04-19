/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


uint64_t pow_mod(uint64_t x, uint64_t n, uint64_t p)
{
  // binary powering
  uint64_t y = 1;
  while (n != 0)
    {
      if (n & 1)
	y = mul_mod(y, x, p);
      x = mul_mod(x, x, p);
      n >>= 1;
    }
  return y;
}


uint64_t pow_mod_pinv(uint64_t x, uint64_t n, uint64_t p, uint64_t pinv)
{
  // binary powering
  uint64_t y = 1;
  while (n != 0)
    {
      if (n & 1)
	y = mul_mod_pinv(y, x, p, pinv);
      x = mul_mod_pinv(x, x, p, pinv);
      n >>= 1;
    }
  return y;
}


uint64_t inv_mod(uint64_t x, uint64_t p)
{
  // use GMP
  // todo: perhaps use mpn?
  mpz_t xx, pp;
  mpz_init_set_ui(xx, x);
  mpz_init_set_ui(pp, p);
  mpz_invert(xx, xx, pp);
  uint64_t r = mpz_get_ui(xx);
  mpz_clear(xx);
  mpz_clear(pp);
  return r;
}


uint64_t calc_pinvb(uint64_t p)
{
  uint64_t pinvb = p;          // correct mod 2^3
  pinvb *= (2 - p * pinvb);    // correct mod 2^6
  pinvb *= (2 - p * pinvb);    // correct mod 2^12
  pinvb *= (2 - p * pinvb);    // correct mod 2^24
  pinvb *= (2 - p * pinvb);    // correct mod 2^48
  pinvb *= (2 - p * pinvb);    // correct mod 2^64
  return pinvb;
}


uint64_t calc_pinv(uint64_t p)
{
  // use GMP
  // todo: perhaps use mpn?
  mpz_t x;
  mpz_init_set_ui(x, 1);
  mpz_mul_2exp(x, x, 126);
  mpz_tdiv_q_ui(x, x, p);
  uint64_t r = mpz_get_ui(x);   // get low 64 bits
  mpz_clear(x);
  return r;
}
