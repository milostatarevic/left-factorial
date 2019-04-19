/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


void fft(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
	 uint64_t* op, size_t n, uint64_t w, uint64_t p,
	 int num_threads)
{
  if (K < FFT_ARRAY_THRESHOLD)
    {
      wtab_t wtab;
      wtab_init(&wtab, K, k2, k3, w, p);
      if (rop != op)
	memcpy(rop, op, n * sizeof(uint64_t));
      fft_base(rop, K, k2, k3, n, &wtab, p);
      wtab_clear(&wtab);
    }
  else
    {
      fft_array(rop, K, k2, k3, op, n, w, p, num_threads);
    }
}


size_t fft_pattern(size_t K, unsigned k2, unsigned k3)
{
  if (K < FFT_ARRAY_THRESHOLD)
    return fft_base_pattern(K, k2, k3);
  else
    return fft_array_pattern(K, k2, k3);
}


void fft_vbuf(vbuf_t* rop, size_t K, unsigned k2, unsigned k3,
	      vbuf_t* op, size_t n, uint64_t w, uint64_t p,
	      int num_threads)
{
  if (K <= rop->block_size)
    {
      fft(vbuf_get_block(rop, 0), K, k2, k3,
	  vbuf_get_block(op, 0), n, w, p, num_threads);
    }
  else
    {
      fft_array_vbuf(rop, K, k2, k3, op, n, w, p, num_threads);
    }
}


size_t fft_vbuf_pattern(size_t K, unsigned k2, unsigned k3, size_t block_size)
{
  if (K <= block_size)
    return fft_pattern(K, k2, k3);
  else
    return fft_array_vbuf_pattern(K, k2, k3, block_size);
}


void ifft(uint64_t* op, size_t K, unsigned k2, unsigned k3,
	  uint64_t w, uint64_t p, int num_threads)
{
  if (K < FFT_ARRAY_THRESHOLD)
    {
      wtab_t wtab;
      wtab_init(&wtab, K, k2, k3, w, p);
      ifft_base(op, K, k2, k3, &wtab, p);
      wtab_clear(&wtab);
    }
  else
    {
      ifft_array(op, K, k2, k3, w, p, num_threads);
    }
}


void ifft_vbuf(vbuf_t* op, size_t K, unsigned k2, unsigned k3,
	       uint64_t w, uint64_t p, int num_threads)
{
  if (K <= op->block_size)
    {
      ifft(vbuf_get_block(op, 0), K, k2, k3, w, p, num_threads);
    }
  else
    {
      ifft_array_vbuf(op, K, k2, k3, w, p, num_threads);
    }
}



void pointwise_multiply(uint64_t* rop, uint64_t* op1, uint64_t* op2, size_t n,
			uint64_t p, uint64_t pinvb, int num_threads)
{
#pragma omp parallel for num_threads(num_threads)
  for (ptrdiff_t j = n - 1; j >= 0; j--)
    rop[j] = mul_mod_pinvb_lazy(op1[j], op2[j], p, pinvb);
}



void conv(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
	  uint64_t* op1, size_t n1, uint64_t* temp1,
	  uint64_t* op2, size_t n2, uint64_t* temp2,
	  uint64_t w, uint64_t p, int num_threads)
{
  if (K < FFT_ARRAY_THRESHOLD)
    {
      wtab_t wtab, winvtab;
      wtab_init(&wtab, K, k2, k3, w, p);
      wtab_init(&winvtab, K, k2, k3, inv_mod(w, p), p);
      conv_base(rop, K, k2, k3, op1, n1, temp1, op2, n2, temp2,
		&wtab, &winvtab, p, calc_pinvb(p));
      wtab_clear(&wtab);
      wtab_clear(&winvtab);
    }
  else
    {
      conv_array(rop, K, k2, k3, op1, n1, temp1, op2, n2, temp2,
		 w, p, num_threads);
    }
}
