/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include <gmp.h>


/*
  Initialises the library; must be called before using any of the
  functions below.
*/
extern void ntt_init();


/*
  Drop-in replacement for mpn_mul.
*/
void ntt_mpn_mul(mp_limb_t* rop, mp_limb_t* op1, size_t n1,
		 mp_limb_t* op2, size_t n2);


/*
  Same as ntt_mpn_mul, with extra functionality.

  Unlike mpn_mul, the output region may overlap either or both input regions.

  Must have n1 >= 1 and n2 >= 1, but unlike mpn_mul, it is not required that
  n1 <= n2.

  If preserve_op1 == 1, then the contents of op1 are not modified (unless of
  course op1 overlaps the output region). If preserve_op1 == 0, op1 may be
  destroyed, and this may reduce total memory consumption. Similarly for
  preserve_op2.

  However, if either preserve_op1 == 0 or preserve_op2 == 0, then {op1,n1} must
  be either disjoint from or identical to {op,n2}.

  If optimise_memory == 1, it will attempt to use less memory, possibly with
  some increase in time.

  num_threads specifies the desired number of threads to use. (This might
  require OpenMP's nested parallelism enabled, if the caller is itself a thread
  in an OpenMP team?)
*/
void ntt_mpn_mul_bonus(mp_limb_t* rop,
		       mp_limb_t* op1, size_t n1, int preserve_op1,
		       mp_limb_t* op2, size_t n2, int preserve_op2,
		       int optimise_memory, int num_threads);
