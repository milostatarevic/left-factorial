/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


void wtab_init(wtab_t* wtab, size_t K, unsigned k2, unsigned k3,
	       uint64_t w, uint64_t p)
{
  uint64_t pinv = calc_pinv(p);

  wtab->w = safe_malloc(2 * (k2 + k3 + 1) * sizeof(uint64_t*));
  wtab->wpinv = wtab->w + (k2 + k3 + 1);

  size_t space = 2 * K * sizeof(uint64_t) + 2 * (k2 + k3 + 1) * CACHE_LINE;
  wtab->buf = (uint64_t*) safe_malloc(space);
  uint64_t* ptr = wtab->buf;

  // wtab->w[j] is always SHIFT1 words modulo the cache line size.
  // wtab->wpinv[j] is always SHIFT2 words modulo the cache line size.
  // The data buffers in fft() are always aligned with cache lines.
  // This scheme hopefully reduces e.g. bank conflicts.
#define SHIFT1 2
#define SHIFT2 4

  // for storing the cube root of unity:
  wtab->w[0] = ptr;
  wtab->wpinv[0] = ptr + 1;
  ptr += 2;

  size_t L = K;

  // radix-3 layers
  for (unsigned i = 0; i < k3; i++)
    {
      // round up to SHIFT1
      ptr = (uint64_t*) next_cache_line(ptr - SHIFT1) + SHIFT1;
      uint64_t* ptr_wtab;
      ptr_wtab = wtab->w[k2 + k3 - i] = ptr;

      uint64_t wpinv = calc_ypinv(w, p, pinv);
      uint64_t wsqr = mul_mod_ypinv(w, w, wpinv, p);
      uint64_t wsqrpinv = calc_ypinv(wsqr, p, pinv);

      L = L / 3;

      // compute values for wtab
      if (i == 0)
	{
	  ptr[0] = 1;
	  for (size_t j = 1; j < L; j++)    // powers of w
	    ptr[j] = mul_mod_ypinv(ptr[j - 1], w, wpinv, p);
	  ptr[L] = 1;
	  for (size_t j = 1; j < L; j++)    // powers of w^2
	    ptr[j + L] = mul_mod_ypinv(ptr[j + L - 1], wsqr, wsqrpinv, p);

	  // get cube root of unity while we're here
	  uint64_t u = mul_mod_ypinv(ptr[L - 1], w, wpinv, p);
	  wtab->w[0][0] = u;
	  wtab->wpinv[0][0] = calc_ypinv(u, p, pinv);
	}
      else
	{
	  // copy from previous iteration
	  uint64_t* prev = wtab->w[k2 + k3 - i + 1];
	  for (size_t j = 0; j < 2 * L; j++)
	    ptr[j] = prev[3 * j];
	}
      ptr += 2 * L;

      // round up to SHIFT2
      ptr = (uint64_t*) next_cache_line(ptr - SHIFT2) + SHIFT2;
      wtab->wpinv[k2 + k3 - i] = ptr;

      // compute values for wtabinv
      if (i == 0)
	{
	  for (size_t j = 0; j < 2 * L; j++)
	    ptr[j] = calc_ypinv(ptr_wtab[j], p, pinv);
	}
      else
	{
	  // copy from previous iteration
	  uint64_t* prev = wtab->wpinv[k2 + k3 - i + 1];
	  for (size_t j = 0; j < 2 * L; j++)
	    ptr[j] = prev[3 * j];
	}
      ptr += 2 * L;
      
      // w := w^3
      w = mul_mod_ypinv(wsqr, w, wpinv, p);
    }

  // radix-2 layers
  for (unsigned i = 0; i < k2; i++)
    {
      // round up to SHIFT1
      ptr = (uint64_t*) next_cache_line(ptr - SHIFT1) + SHIFT1;
      uint64_t* ptr_wtab;
      ptr_wtab = wtab->w[k2 - i] = ptr;

      uint64_t wpinv = calc_ypinv(w, p, pinv);

      L = L / 2;

      // compute values for wtab
      if (i == 0)
	{
	  ptr[0] = 1;
	  for (size_t j = 1; j < L; j++)    // powers of w
	    ptr[j] = mul_mod_ypinv(ptr[j - 1], w, wpinv, p);
	}
      else
	{
	  // copy from previous iteration
	  uint64_t* prev = wtab->w[k2 - i + 1];
	  for (size_t j = 0; j < L; j++)
	    ptr[j] = prev[2 * j];
	}
      ptr += L;

      // round up to SHIFT2
      ptr = (uint64_t*) next_cache_line(ptr - SHIFT2) + SHIFT2;
      wtab->wpinv[k2 - i] = ptr;

      // compute values for wtabinv
      if (i == 0)
	{
	  for (size_t j = 0; j < L; j++)
	    ptr[j] = calc_ypinv(ptr_wtab[j], p, pinv);
	}
      else
	{
	  // copy from previous iteration
	  uint64_t* prev = wtab->wpinv[k2 - i + 1];
	  for (size_t j = 0; j < L; j++)
	    ptr[j] = prev[2 * j];
	}
      ptr += L;

      // w := w^3
      w = mul_mod_ypinv(w, w, wpinv, p);
    }

  assert(w == 1);
  assert((ptr - wtab->buf) * sizeof(uint64_t) <= space);
}


void wtab_clear(wtab_t* wtab)
{
  safe_free(wtab->buf);
  safe_free(wtab->w);
}



/*
  One radix-3 layer of fft_base.
  Inputs is K/(3*L) blocks of size 3*L.
  u = cube root of unity
*/
void fft_base_radix3_layer(uint64_t* op, size_t K, size_t L,
                           uint64_t* w, uint64_t* wpinv,
                           uint64_t u, uint64_t upinv, uint64_t p)
{
  for (size_t i = K; i > 0; i -= 3*L, op += 3*L)
    {
      for (size_t j = L - 1; j > 0; j--)
        {
          // radix-3 butterfly
          uint64_t x1 = op[j + L];              // x1 in [0, 2p)
          uint64_t x2 = op[j + 2*L];            // x2 in [0, 2p)
          uint64_t t0 = add_mod(x1, x2, 2*p);   // x1 + x2 in [0, 2p)
          uint64_t t1 = x1 - x2 + 2*p;          // x1 - x2 in [0, 4p)
          uint64_t t2 = mul_mod_ypinv_lazy(t1, u, upinv, p);
                                                // u*x1 - u*x2 in [0, 2p)
          uint64_t t3 = sub_mod(x2, t2, 2*p);   // -u*x1 - u^2*x2 in [0, 2p)
          uint64_t t4 = sub_mod(t0, t3, 2*p);   // -u^2*x1 - u*x2 in [0, 2p)
          uint64_t x0 = op[j];                  // x0 in [0, 2p)
          uint64_t y0 = add_mod(x0, t0, 2*p);   // x0 + x1 + x2 in [0, 2p)
          op[j] = y0;
          uint64_t y1 = x0 - t3 + 2*p;          // x0 + u*x1 + u^2*x2 in [0, 4p)
          op[j + L] = mul_mod_ypinv_lazy(y1, w[j], wpinv[j], p);
          uint64_t y2 = x0 - t4 + 2*p;          // x0 + u^2*x1 + u*x2 in [0, 4p)
	  op[j + 2*L] = mul_mod_ypinv_lazy(y2, w[j + L], wpinv[j + L], p);
        }
      // last butterfly
      uint64_t x1 = op[L];                  // x1 in [0, 2p)
      uint64_t x2 = op[2*L];                // x2 in [0, 2p)
      uint64_t t0 = add_mod(x1, x2, 2*p);   // x1 + x2 in [0, 2p)
      uint64_t t1 = x1 - x2 + 2*p;          // x1 - x2 in [0, 4p)
      uint64_t t2 = mul_mod_ypinv_lazy(t1, u, upinv, p);
                                            // u*x1 - u*x2 in [0, 2p)
      uint64_t t3 = sub_mod(x2, t2, 2*p);   // -u*x1 - u^2*x2 in [0, 2p)
      uint64_t t4 = sub_mod(t0, t3, 2*p);   // -u^2*x1 - u*x2 in [0, 2p)
      uint64_t x0 = op[0];                  // x0 in [0, 2p)
      op[0] = add_mod(x0, t0, 2*p);         // x0 + x1 + x2 in [0, 2p)
      op[L] = sub_mod(x0, t3, 2*p);         // x0 + u*x1 + u^2*x2 in [0, 2p)
      op[2*L] = sub_mod(x0, t4, 2*p);       // x0 + u^2*x1 + u*x2 in [0, 2p)
    }
}



/*
  One radix-2 layer of fft_base.
  Inputs is K/(2*L) blocks of size 2*L.
  Assumes L is divisible by 2.
*/
void fft_base_radix2_layer(uint64_t* op, size_t K, size_t L,
                           uint64_t* w, uint64_t* wpinv, uint64_t p)
{
  for (size_t i = K; i > 0; i -= 2*L, op += 2*L)
    {
      for (size_t j = L - 2; j > 0; j -= 2)
        {
          {
            // one butterfly
            uint64_t x0 = op[j + 1];              // x0 in [0, 2p)
            uint64_t x1 = op[j + L + 1];          // x1 in [0, 2p)
            op[j + 1] = add_mod(x0, x1, 2*p);     // x0 + x1 in [0, 2p)
            uint64_t t = x0 - x1 + 2*p;           // x0 - x1 in [0, 4p)
            op[j + L + 1] =
              mul_mod_ypinv_lazy(t, w[j + 1], wpinv[j + 1], p);
          }
          {
            // another butterfly
            uint64_t x0 = op[j];                  // x0 in [0, 2p)
            uint64_t x1 = op[j + L];              // x1 in [0, 2p)
            op[j] = add_mod(x0, x1, 2*p);         // x0 + x1 in [0, 2p)
            uint64_t t = x0 - x1 + 2*p;           // x0 - x1 in [0, 4p)
            op[j + L] = mul_mod_ypinv_lazy(t, w[j], wpinv[j], p);
          }
        }
      // last two butterflies
      {
        uint64_t x0 = op[1];                  // x0 in [0, 2p)
        uint64_t x1 = op[L + 1];              // x1 in [0, 2p)
        op[1] = add_mod(x0, x1, 2*p);         // x0 + x1 in [0, 2p)
        uint64_t t = x0 - x1 + 2*p;           // x0 - x1 in [0, 4p)
        op[L + 1] = mul_mod_ypinv_lazy(t, w[1], wpinv[1], p);
      }
      {
        uint64_t x0 = op[0];                  // x0 in [0, 2p)
        uint64_t x1 = op[L];                  // x1 in [0, 2p)
        op[0] = add_mod(x0, x1, 2*p);         // x0 + x1 in [0, 2p)
        op[L] = sub_mod(x0, x1, 2*p);         // x0 - x1 in [0, 2p)
      }
    }
}



/*
  Last two radix-2 layers of fft_base.
  Inputs is K/4 blocks of size 4.
  u = 4th root of unity
*/
void fft_base_radix2_last_2_layers(uint64_t* op, size_t K,
                                   uint64_t u, uint64_t upinv, uint64_t p)
{
  for (size_t i = K; i > 0; i -= 4, op += 4)
    {
      uint64_t x0 = op[0];                       // x0 in [0, 2p)
      uint64_t x1 = op[1];                       // x1 in [0, 2p)
      uint64_t x2 = op[2];                       // x2 in [0, 2p)
      uint64_t x3 = op[3];                       // x3 in [0, 2p)

      uint64_t y0 = add_mod(x0, x2, 2*p);        // y0 = x0 + x2 in [0, 2p)
      uint64_t y2 = sub_mod(x0, x2, 2*p);        // y2 = x0 - x2 in [0, 2p)
      uint64_t y1 = add_mod(x1, x3, 2*p);        // y1 = x1 + x3 in [0, 2p)
      uint64_t t = x1 - x3 + 2*p;                // x1 - x3 in [0, 4p)
      uint64_t y3 = mul_mod_ypinv_lazy(t, u, upinv, p);
                                                 // y3 = u*(x1 - x3) in [0, 2p)
      op[0] = add_mod(y0, y1, 2*p);
      op[1] = sub_mod(y0, y1, 2*p);
      op[2] = add_mod(y2, y3, 2*p);
      op[3] = sub_mod(y2, y3, 2*p);
    }
}



/*
  Last radix-2 layer of fft_base.
  Inputs is K/2 blocks of size 2.
*/
void fft_base_radix2_last_layer(uint64_t* op, size_t K, uint64_t p)
{
  if (K % 4 != 0)
    {
      uint64_t x0 = op[0];
      uint64_t x1 = op[1];
      op[0] = add_mod(x0, x1, 2*p);
      op[1] = sub_mod(x0, x1, 2*p);
      op += 2;
      K -= 2;
    }

  // now assume K divisible by 4
  for (size_t i = K; i > 0; i -= 4, op += 4)
    {
      {
        uint64_t x0 = op[0];
        uint64_t x1 = op[1];
        op[0] = add_mod(x0, x1, 2*p);
        op[1] = sub_mod(x0, x1, 2*p);
      }
      {
        uint64_t x0 = op[2];
        uint64_t x1 = op[3];
        op[2] = add_mod(x0, x1, 2*p);
        op[3] = sub_mod(x0, x1, 2*p);
      }
    }
}



void fft_base(uint64_t* op, size_t K, unsigned k2, unsigned k3, size_t n,
	      wtab_t* wtab, uint64_t p)
{
  for (size_t i = n; i < K; i++)
    op[i] = 0;

  if (K >= FFT_SPLIT_THRESHOLD)
    {
      // perform one FFT layer and then recurse into each half/third
      if (k3 >= 1)
	{
	  size_t L = K / 3;
	  uint64_t u = wtab->w[0][0];
	  uint64_t upinv = wtab->wpinv[0][0];
          fft_base_radix3_layer(op, K, L, wtab->w[k2 + k3],
				wtab->wpinv[k2 + k3], u, upinv, p);
	  fft_base(op, L, k2, k3 - 1, L, wtab, p);
	  fft_base(op + L, L, k2, k3 - 1, L, wtab, p);
	  fft_base(op + 2*L, L, k2, k3 - 1, L, wtab, p);
	}
      else
	{
	  size_t L = K / 2;
          fft_base_radix2_layer(op, K, L, wtab->w[k2], wtab->wpinv[k2], p);
	  fft_base(op, L, k2 - 1, k3, L, wtab, p);
	  fft_base(op + L, L, k2 - 1, k3, L, wtab, p);
	}

      return;
    }

  // base case iterative FFT

  size_t L = K;

  // radix-3 layers
  if (k3 >= 1)
    {
      // u = cube root of unity
      uint64_t u = wtab->w[0][0];
      uint64_t upinv = wtab->wpinv[0][0];

      for (unsigned i = 0; i < k3; i++)
        {
          L /= 3;
          fft_base_radix3_layer(op, K, L, wtab->w[k2 + k3 - i],
				wtab->wpinv[k2 + k3 - i], u, upinv, p);
        }
    }

  // radix-2 layers
  if (k2 >= 2)
    {
      for (unsigned i = 0; i < k2 - 2; i++)
        {
          L /= 2;
          fft_base_radix2_layer(op, K, L, wtab->w[k2 - i],
				wtab->wpinv[k2 - i], p);
        }

      uint64_t u = wtab->w[2][1];    // 4th root of unity
      uint64_t upinv = wtab->wpinv[2][1];
      fft_base_radix2_last_2_layers(op, K, u, upinv, p);
    }
  else if (k2 == 1)
    fft_base_radix2_last_layer(op, K, p);
}


size_t fft_base_pattern(size_t K, unsigned k2, unsigned k3)
{
  // all radix-3 layers first
  return ((size_t) 1 << k3) - 1;
}



/*
  One radix-3 layer of ifft_base.
  Inputs is K/(3*L) blocks of size 3*L.
  u = cube root of unity
*/
void ifft_base_radix3_layer(uint64_t* op, size_t K, size_t L,
			    uint64_t* w, uint64_t* wpinv,
			    uint64_t u, uint64_t upinv, uint64_t p)
{
  for (size_t i = K; i > 0; i -= 3*L, op += 3*L)
    {
      for (size_t j = L - 1; j > 0; j--)
        {
          // radix-3 butterfly
          uint64_t x0 = op[j];                  // x0 in [0, 4p)
	  x0 -= (x0 >= 2*p) ? (2*p) : 0;        // x0 in [0, 2p)
          uint64_t x1 = mul_mod_ypinv_lazy(op[j + L], w[j], wpinv[j], p);
	                                        // x1 in [0, 2p) 
	  uint64_t x2 =                         // x2 in [0, 2p)
	    mul_mod_ypinv_lazy(op[j + 2*L], w[j + L], wpinv[j + L], p);
          uint64_t t0 = add_mod(x1, x2, 2*p);   // x1 + x2 in [0, 2p)
          uint64_t t1 = x1 - x2 + 2*p;          // x1 - x2 in [0, 4p)
          uint64_t t2 = mul_mod_ypinv_lazy(t1, u, upinv, p);
                                                // u*x1 - u*x2 in [0, 2p)
          uint64_t t3 = sub_mod(x2, t2, 2*p);   // -u*x1 - u^2*x2 in [0, 2p)
          uint64_t t4 = sub_mod(t0, t3, 2*p);   // -u^2*x1 - u*x2 in [0, 2p)
	  op[j] = x0 + t0;                      // x0 + x1 + x2 in [0, 4p)
	  op[j + L] = x0 - t3 + 2*p;            // x0 + u*x1 + u^2*x2 in [0, 4p)
	  op[j + 2*L] = x0 - t4 + 2*p;          // x0 + u^2*x1 + u*x2 in [0, 4p)
        }
      // last butterfly
      uint64_t x0 = op[0];                  // x0 in [0, 4p)
      x0 -= (x0 >= 2*p) ? (2*p) : 0;        // x0 in [0, 2p)
      uint64_t x1 = op[L];                  // x1 in [0, 4p)
      x1 -= (x1 >= 2*p) ? (2*p) : 0;        // x1 in [0, 2p)
      uint64_t x2 = op[2*L];                // x2 in [0, 4p)
      x2 -= (x2 >= 2*p) ? (2*p) : 0;        // x2 in [0, 2p)
      uint64_t t0 = add_mod(x1, x2, 2*p);   // x1 + x2 in [0, 2p)
      uint64_t t1 = x1 - x2 + 2*p;          // x1 - x2 in [0, 4p)
      uint64_t t2 = mul_mod_ypinv_lazy(t1, u, upinv, p); 
                                            // u*x1 - u*x2 in [0, 2p)
      uint64_t t3 = sub_mod(x2, t2, 2*p);   // -u*x1 - u^2*x2 in [0, 2p)
      uint64_t t4 = sub_mod(t0, t3, 2*p);   // -u^2*x1 - u*x2 in [0, 2p)
      op[0] = x0 + t0;                      // x0 + x1 + x2 in [0, 4p)
      op[L] = x0 - t3 + 2*p;                // x0 + u*x1 + u^2*x2 in [0, 4p)
      op[2*L] = x0 - t4 + 2*p;              // x0 + u^2*x1 + u*x2 in [0, 4p)
    }
}



/*
  One radix-2 layer of ifft_base.
  Inputs is K/(2*L) blocks of size 2*L.
  Assumes L is divisible by 2.
*/
void ifft_base_radix2_layer(uint64_t* op, size_t K, size_t L,
			    uint64_t* w, uint64_t* wpinv, uint64_t p)
{
  for (size_t i = K; i > 0; i -= 2*L, op += 2*L)
    {
      for (size_t j = L - 2; j > 0; j -= 2)
        {
          {
            // one butterfly
            uint64_t x0 = op[j + 1];              // x0 in [0, 4p)
	    x0 -= (x0 >= 2*p) ? (2*p) : 0;        // x0 in [0, 2p)
            uint64_t x1 =
	      mul_mod_ypinv_lazy(op[j + L + 1], w[j + 1], wpinv[j + 1], p);
                                                  // x1 in [0, 2p)
	    op[j + 1] = x0 + x1;                  // x0 + x1 in [0, 4p)
	    op[j + L + 1] = x0 - x1 + 2*p;        // x0 - x1 in [0, 4p)
          }
          {
            // another butterfly
            uint64_t x0 = op[j];                  // x0 in [0, 4p)
	    x0 -= (x0 >= 2*p) ? (2*p) : 0;        // x0 in [0, 2p)
            uint64_t x1 = mul_mod_ypinv_lazy(op[j + L], w[j], wpinv[j], p);
                                                  // x1 in [0, 2p)
	    op[j] = x0 + x1;                      // x0 + x1 in [0, 4p)
	    op[j + L] = x0 - x1 + 2*p;            // x0 - x1 in [0, 4p)
          }
        }
      // last two butterflies
      {
	uint64_t x0 = op[1];                  // x0 in [0, 4p)
	x0 -= (x0 >= 2*p) ? (2*p) : 0;        // x0 in [0, 2p)
	uint64_t x1 = mul_mod_ypinv_lazy(op[L + 1], w[1], wpinv[1], p);
                                              // x1 in [0, 2p)
	op[1] = x0 + x1;                      // x0 + x1 in [0, 4p)
	op[L + 1] = x0 - x1 + 2*p;            // x0 - x1 in [0, 4p)
      }
      {
	uint64_t x0 = op[0];                  // x0 in [0, 4p)
	x0 -= (x0 >= 2*p) ? (2*p) : 0;        // x0 in [0, 2p)
	uint64_t x1 = op[L];                  // x1 in [0, 4p)
	x1 -= (x1 >= 2*p) ? (2*p) : 0;        // x1 in [0, 2p)
	op[0] = x0 + x1;                      // x0 + x1 in [0, 4p)
	op[L] = x0 - x1 + 2*p;                // x0 - x1 in [0, 4p)
      }
    }
}



/*
  Last two radix-2 layers of ifft_base.
  Inputs is K/4 blocks of size 4.
  u = 4th root of unity
*/
void ifft_base_radix2_last_2_layers(uint64_t* op, size_t K,
				    uint64_t u, uint64_t upinv, uint64_t p)
{
  for (size_t i = K; i > 0; i -= 4, op += 4)
    {
      uint64_t x0 = op[0];
      x0 -= (x0 >= 2*p) ? (2*p) : 0;             // x0 in [0, 2p)
      uint64_t x1 = op[1];
      x1 -= (x1 >= 2*p) ? (2*p) : 0;             // x1 in [0, 2p)
      uint64_t x2 = op[2];
      x2 -= (x2 >= 2*p) ? (2*p) : 0;             // x2 in [0, 2p)
      uint64_t x3 = op[3];
      x3 -= (x3 >= 2*p) ? (2*p) : 0;             // x3 in [0, 2p)

      uint64_t y0 = add_mod(x0, x1, 2*p);        // y0 = x0 + x1 in [0, 2p)
      uint64_t y1 = sub_mod(x0, x1, 2*p);        // y1 = x0 - x1 in [0, 2p)
      uint64_t y2 = add_mod(x2, x3, 2*p);        // y2 = x2 + x3 in [0, 2p)
      uint64_t t = x2 - x3 + 2*p;                // t = x2 - x3 in [0, 4p)
      uint64_t y3 = mul_mod_ypinv_lazy(t, u, upinv, p);
                                                 // y3 = u*(x2 - x3) in [0, 2p)

      op[0] = y0 + y2;                           // y0 + y2 in [0, 4p)
      op[2] = y0 - y2 + 2*p;                     // y0 - y2 in [0, 4p)
      op[1] = y1 + y3;                           // y1 + y3 in [0, 4p)
      op[3] = y1 - y3 + 2*p;                     // y1 - y3 in [0, 4p)
    }
}



/*
  Last radix-2 layer of ifft_base.
  Inputs is K/2 blocks of size 2.
*/
void ifft_base_radix2_last_layer(uint64_t* op, size_t K, uint64_t p)
{
  if (K % 4 != 0)
    {
      uint64_t x0 = op[0];
      x0 -= (x0 >= 2*p) ? (2*p) : 0;
      uint64_t x1 = op[1];
      x1 -= (x1 >= 2*p) ? (2*p) : 0;
      op[0] = x0 + x1;
      op[1] = x0 - x1 + 2*p;
      op += 2;
      K -= 2;
    }

  // now assume K divisible by 4
  for (size_t i = K; i > 0; i -= 4, op += 4)
    {
      {
        uint64_t x0 = op[0];
	x0 -= (x0 >= 2*p) ? (2*p) : 0;
        uint64_t x1 = op[1];
	x1 -= (x1 >= 2*p) ? (2*p) : 0;
        op[0] = x0 + x1;
        op[1] = x0 - x1 + 2*p;
      }
      {
        uint64_t x0 = op[2];
	x0 -= (x0 >= 2*p) ? (2*p) : 0;
        uint64_t x1 = op[3];
	x1 -= (x1 >= 2*p) ? (2*p) : 0;
        op[2] = x0 + x1;
        op[3] = x0 - x1 + 2*p;
      }
    }
}



void ifft_base(uint64_t* op, size_t K, unsigned k2, unsigned k3,
	       wtab_t* wtab, uint64_t p)
{
  if (K >= FFT_SPLIT_THRESHOLD)
    {
      // perform one FFT layer and then recurse into each half/third
      if (k3 >= 1)
	{
	  size_t L = K / 3;
	  ifft_base(op, L, k2, k3 - 1, wtab, p);
	  ifft_base(op + L, L, k2, k3 - 1, wtab, p);
	  ifft_base(op + 2*L, L, k2, k3 - 1, wtab, p);
	  uint64_t u = wtab->w[0][0];
	  uint64_t upinv = wtab->wpinv[0][0];
          ifft_base_radix3_layer(op, K, L, wtab->w[k2 + k3],
				 wtab->wpinv[k2 + k3], u, upinv, p);
	}
      else
	{
	  size_t L = K / 2;
	  ifft_base(op, L, k2 - 1, k3, wtab, p);
	  ifft_base(op + L, L, k2 - 1, k3, wtab, p);
          ifft_base_radix2_layer(op, K, L, wtab->w[k2], wtab->wpinv[k2], p);
	}

      return;
    }

  // base case iterative FFT

  size_t L = 1;

  // radix-2 layers
  if (k2 >= 2)
    {
      uint64_t u = wtab->w[2][1];    // 4th root of unity
      uint64_t upinv = wtab->wpinv[2][1];

      ifft_base_radix2_last_2_layers(op, K, u, upinv, p);
      L *= 4;

      for (int i = k2 - 3; i >= 0; i--)
        {
          ifft_base_radix2_layer(op, K, L, wtab->w[k2 - i],
				 wtab->wpinv[k2 - i], p);
          L *= 2;
        }
    }
  else if (k2 == 1)
    {
      ifft_base_radix2_last_layer(op, K, p);
      L *= 2;
    }

  // radix-3 layers
  if (k3 >= 1)
    {
      // u = cube root of unity
      uint64_t u = wtab->w[0][0];
      uint64_t upinv = wtab->wpinv[0][0];

      for (int i = k3 - 1; i >= 0; i--)
        {
          ifft_base_radix3_layer(op, K, L, wtab->w[k2 + k3 - i],
				 wtab->wpinv[k2 + k3 - i], u, upinv, p);
	  L *= 3;
        }
    }
}



void conv_base(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
	       uint64_t* op1, size_t n1, uint64_t* temp1,
	       uint64_t* op2, size_t n2, uint64_t* temp2,
	       wtab_t* wtab, wtab_t* winvtab, uint64_t p, uint64_t pinvb)
{
  int squaring = (op1 == op2) && (n1 == n2);

  if (temp1 != op1)
    memcpy(temp1, op1, n1 * sizeof(uint64_t));
  if (!squaring && temp2 != op2)
    memcpy(temp2, op2, n2 * sizeof(uint64_t));
  if (squaring)
    temp2 = temp1;

  if (K >= FFT_SPLIT_THRESHOLD)
    {
      // recursively split into two/three convolutions

      for (size_t i = n1; i < K; i++)
	temp1[i] = 0;
      if (!squaring)
	for (size_t i = n2; i < K; i++)
	  temp2[i] = 0;

      if (k3 >= 1)
	{
	  size_t L = K / 3;
	  uint64_t u = wtab->w[0][0];
	  uint64_t upinv = wtab->wpinv[0][0];

          fft_base_radix3_layer(temp1, K, L, wtab->w[k2 + k3],
				wtab->wpinv[k2 + k3], u, upinv, p);
	  if (!squaring)
	    fft_base_radix3_layer(temp2, K, L, wtab->w[k2 + k3],
				  wtab->wpinv[k2 + k3], u, upinv, p);

	  conv_base(rop, L, k2, k3 - 1, temp1, L, temp1,
		    temp2, L, temp2, wtab, winvtab, p, pinvb);
	  conv_base(rop + L, L, k2, k3 - 1, temp1 + L, L, temp1 + L,
		    temp2 + L, L, temp2 + L, wtab, winvtab, p, pinvb);
	  conv_base(rop + 2*L, L, k2, k3 - 1, temp1 + 2*L, L, temp1 + 2*L,
		    temp2 + 2*L, L, temp2 + 2*L, wtab, winvtab, p, pinvb);

	  u = winvtab->w[0][0];
	  upinv = winvtab->wpinv[0][0];
	  ifft_base_radix3_layer(rop, K, L, winvtab->w[k2 + k3],
				 winvtab->wpinv[k2 + k3], u, upinv, p);
	}
      else
	{
	  size_t L = K / 2;

          fft_base_radix2_layer(temp1, K, L, wtab->w[k2], wtab->wpinv[k2], p);
	  if (!squaring)
	    fft_base_radix2_layer(temp2, K, L, wtab->w[k2], wtab->wpinv[k2], p);

	  conv_base(rop, L, k2 - 1, k3, temp1, L, temp1,
		    temp2, L, temp2, wtab, winvtab, p, pinvb);
	  conv_base(rop + L, L, k2 - 1, k3, temp1 + L, L, temp1 + L,
		    temp2 + L, L, temp2 + L, wtab, winvtab, p, pinvb);

	  ifft_base_radix2_layer(rop, K, L, winvtab->w[k2],
				 winvtab->wpinv[k2], p);
	}

      return;
    }

  // base case: FFT each input, multiply, and IFFT

  fft_base(temp1, K, k2, k3, n1, wtab, p);
  if (!squaring)
    fft_base(temp2, K, k2, k3, n2, wtab, p);

  pointwise_multiply(rop, temp1, temp2, K, p, pinvb, 1);

  ifft_base(rop, K, k2, k3, winvtab, p);
}
