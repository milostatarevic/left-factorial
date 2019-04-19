/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


/*
  these random-number generation routines are NOT THREADSAFE.
*/

// uniformly random 32 bits
INLINE uint32_t random32()
{
  static uint64_t state = 0;
  state = state * 6364136223846793005UL + 1442695040888963407UL;
  return state >> 32;
}

// uniformly random 64 bits
uint64_t random64();

// uniform random integer in 0 <= x < p
// assumes 2^61 < p < 2^62
uint64_t random_mod(uint64_t p);
