/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


/*
  Interpret src and dst as rectangular arrays, with dst_stride and src_stride
  entries per row respectively. Copy r rows * c columns from src to dst.
*/
void copy_block(uint64_t* dst, size_t dst_stride, 
		uint64_t* src, size_t src_stride,
		size_t r, size_t c)
{
  for (size_t j = 0; j < r; j++)
    for (size_t i = 0; i < c; i++)
      dst[j*dst_stride + i] = src[j*src_stride + i];
}


// similar to copy_block, but transposes the array.
void transpose_block(uint64_t* dst, size_t dst_stride, 
		     uint64_t* src, size_t src_stride,
		     size_t r, size_t c)
{
  for (size_t j = 0; j < r; j++)
    for (size_t i = 0; i < c; i++)
      dst[j + i*dst_stride] = src[j*src_stride + i];
}


void fft_cols_vbuf(vbuf_t* rop,
		   size_t R, unsigned r2, unsigned r3,
		   size_t C, unsigned c2, unsigned c3,
		   vbuf_t* op, size_t n, uint64_t w, uint64_t p,
		   int num_threads)
{
  uint64_t pinv = calc_pinv(p);
  uint64_t wC = pow_mod_pinv(w, C, p, pinv);

  // use fft_base for column transforms if they are small enough
  int use_fft_base = (R < FFT_ARRAY_THRESHOLD);

  // compute root tables if necessary
  wtab_t wtab;
  if (use_fft_base)
    wtab_init(&wtab, R, r2, r3, wC, p);

  // number of columns to extract on each iteration
  const size_t clump = CACHE_LINE / sizeof(uint64_t);

  // number of rows in each memory pool block
  size_t block_rows = rop->block_size / C;
  assert(rop->block_size == C * block_rows);

#pragma omp parallel num_threads(num_threads)
  {
    uint64_t temp[clump*clump];
    uint64_t* cols = (uint64_t*) safe_malloc(R * clump * sizeof(uint64_t));

    // nR = number of complete rows before n entries
    size_t nR = n / C;
    size_t nR2 = MIN(R, nR + 1);
    // nC = number of columns in row #nR
    size_t nC = n - nR * C;

#pragma omp for schedule(dynamic,1)
    // handle columns in blocks
    for (ptrdiff_t i = 0; i < C; i += clump)
      {
	// c = number of columns to handle on this iteration
	size_t c = MIN(clump, C - i);

	// transfer c columns from op to cols
	{
	  for (size_t j = 0; j < nR; j += block_rows)
	    {
	      // r = number of rows to transfer on this iteration
	      size_t r = MIN(block_rows, nR - j);

	      for (size_t jj = 0; jj < r; jj += clump)
		{
		  size_t rr = MIN(clump, r - jj);

		  // transpose from op to cols via temp
		  copy_block(temp, clump, vbuf_get(op, i + (j + jj) * C),
			     C, rr, c);
		  transpose_block(cols + j + jj, R, temp, clump, rr, c);
		}
	    }

	  // handle last partial row
	  if (nR < R)
	    {
	      size_t cc = nC - i;
	      cc = MAX((ptrdiff_t) cc, 0);
	      cc = MIN(cc, c);
	      uint64_t* src = vbuf_get(op, i + nR * C);
	      for (size_t ii = 0; ii < cc; ii++)
		cols[nR + ii*R] = src[ii];
	      for (size_t ii = cc; ii < c; ii++)
		cols[nR + ii*R] = 0;
	    }
	}
	
	// transform columns
	for (size_t ii = 0; ii < c; ii++)
	  {
	    if (use_fft_base)
	      fft_base(cols + ii*R, R, r2, r3, nR2, &wtab, p);
	    else
	      fft_array(cols + ii*R, R, r2, r3, cols + ii*R, nR2, wC, p, 1);
	  }

	// transfer columns from cols to rop
	{
	  // handle rows in blocks
	  for (size_t j = 0; j < R; j += block_rows)
	    {
	      // r = number of rows to transfer on this iteration
	      size_t r = MIN(block_rows, R - j);

	      for (size_t jj = 0; jj < r; jj += clump)
		{
		  size_t rr = MIN(clump, r - jj);

		  // transpose from cols to rop via temp
		  transpose_block(temp, clump, cols + j + jj, R, c, rr);
		  copy_block(vbuf_get(rop, i + (j + jj) * C),
			     C, temp, clump, rr, c);
		}
	    }
	}
      }

    free(cols);
  }

  if (use_fft_base)
    wtab_clear(&wtab);
}


void fft_cols(uint64_t* rop,
	      size_t R, unsigned r2, unsigned r3,
	      size_t C, unsigned c2, unsigned c3,
	      uint64_t* op, size_t n, uint64_t w, uint64_t p,
	      int num_threads)
{
  // create wrapper vbufs and use vbuf version
  vbuf_t vbuf_rop, vbuf_op;
  vbuf_init_wrap(&vbuf_op, R * C, op);
  vbuf_init_wrap(&vbuf_rop, R * C, rop);
  fft_cols_vbuf(&vbuf_rop, R, r2, r3, C, c2, c3,
		&vbuf_op, n, w, p, num_threads);
  vbuf_clear(&vbuf_op);
  vbuf_clear(&vbuf_rop);
}



// helper function for fft_rows and fft_rows_vbuf
void compute_twiddles(uint64_t* rop, size_t R, unsigned r2, unsigned r3,
		      uint64_t w, uint64_t p, uint64_t pinv)
{
  size_t pattern = fft_pattern(R, r2, r3);

  uint64_t wpow[r2 + r3];
  uint64_t wpowpinv[r2 + r3];
  uint64_t t = w;
  for (unsigned i = 0; i < r2 + r3; i++)
    {
      uint64_t tpinv = calc_ypinv(t, p, pinv);
      wpow[i] = t;
      wpowpinv[i] = tpinv;
      uint64_t tsqr = mul_mod_ypinv(t, t, tpinv, p);
      t = ((pattern >> i) & 1) ? mul_mod_ypinv(tsqr, t, tpinv, p) : tsqr;
    }

  // twiddle factors stored with implied factor of 2^64
  rop[0] = (-((uint64_t) 1)) % p + 1;
  size_t L = 1;
  for (int i = r2 + r3 - 1; i >= 0; i--)
    {
      uint64_t u = wpow[i];
      uint64_t upinv = wpowpinv[i];
      if ((pattern >> i) & 1)
	{
	  for (size_t j = 0; j < 2*L; j++)
	    rop[j + L] = mul_mod_ypinv_lazy(rop[j], u, upinv, p);
	  L *= 3;
	}
      else
	{
	  for (size_t j = 0; j < L; j++)
	    rop[j + L] = mul_mod_ypinv_lazy(rop[j], u, upinv, p);
	  L *= 2;
	}
    }
}


// multiply op[i] by (t/2^64)^i for 1 <= i < C
// inputs and outputs in [0, 2p)
void apply_twiddles(uint64_t* op, size_t C, uint64_t t,
		    uint64_t p, uint64_t pinvb)
{
  uint64_t u = t;
  for (size_t j = 1; j < C; j++)
    {
      op[j] = mul_mod_pinvb_lazy(op[j], u, p, pinvb);
      u = mul_mod_pinvb_lazy(u, t, p, pinvb);
    }
}


void fft_rows_vbuf(vbuf_t* op,
		   size_t R, unsigned r2, unsigned r3,
		   size_t C, unsigned c2, unsigned c3,
		   uint64_t w, uint64_t p,
		   int num_threads)
{
  assert(op->block_size % C == 0);
  uint64_t pinv = calc_pinv(p);
  uint64_t pinvb = calc_pinvb(p);
  uint64_t wR = pow_mod_pinv(w, R, p, pinv);

  // use fft_base for row transforms if they are small enough
  int use_fft_base = (C < FFT_ARRAY_THRESHOLD);

  // compute root tables if necessary
  wtab_t wtab;
  if (use_fft_base)
    wtab_init(&wtab, C, c2, c3, wR, p);

  // compute twiddle[i] = basic twiddle factor to use for row #i
  uint64_t* t = safe_malloc(R * sizeof(uint64_t));
  compute_twiddles(t, R, r2, r3, w, p, pinv);

  // transform rows
#pragma omp parallel for schedule(dynamic,1) num_threads(num_threads)
  for (ptrdiff_t i = 0; i < R; i++)
    {
      uint64_t* x = vbuf_get(op, i * C);

      // apply twiddle factors
      if (i != 0)
	apply_twiddles(x, C, t[i], p, pinvb);

      // transform row
      if (use_fft_base)
	fft_base(x, C, c2, c3, C, &wtab, p);
      else
	fft_array(x, C, c2, c3, x, C, wR, p, 1);
    }

  safe_free(t);
  if (use_fft_base)
    wtab_clear(&wtab);
}



void fft_rows(uint64_t* op,
	      size_t R, unsigned r2, unsigned r3,
	      size_t C, unsigned c2, unsigned c3,
	      uint64_t w, uint64_t p,
	      int num_threads)
{
  // create wrapper vbufs and use vbuf version
  vbuf_t vbuf;
  vbuf_init_wrap(&vbuf, R * C, op);
  fft_rows_vbuf(&vbuf, R, r2, r3, C, c2, c3, w, p, num_threads);
  vbuf_clear(&vbuf);
}


size_t fft_cols_rows_pattern(size_t R, unsigned r2, unsigned r3,
			     size_t C, unsigned c2, unsigned c3)
{
  return (fft_pattern(C, c2, c3) << (r2 + r3)) + fft_pattern(R, r2, r3);
}


size_t fft_cols_rows_vbuf_pattern(size_t R, unsigned r2, unsigned r3,
				  size_t C, unsigned c2, unsigned c3,
				  size_t block_size)
{
  return fft_cols_rows_pattern(R, r2, r3, C, c2, c3);
}



// select array decomposition for fft_array
void fft_array_params(size_t* R, unsigned* r2, unsigned* r3,
		      size_t* C, unsigned* c2, unsigned* c3,
		      unsigned k2, unsigned k3)
{
  *r2 = k2 / 2;
  *r3 = k3 / 2;
  *c2 = k2 - *r2;
  *c3 = k3 - *r3;
  *R = pow23(*r2, *r3);
  *C = pow23(*c2, *c3);
}


void fft_array(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
	       uint64_t* op, size_t n, uint64_t w, uint64_t p,
	       int num_threads)
{
  size_t R, C;
  unsigned r2, r3, c2, c3;
  fft_array_params(&R, &r2, &r3, &C, &c2, &c3, k2, k3);
  fft_cols(rop, R, r2, r3, C, c2, c3, op, n, w, p, num_threads);
  fft_rows(rop, R, r2, r3, C, c2, c3, w, p, num_threads);
}


size_t fft_array_pattern(size_t K, unsigned k2, unsigned k3)
{
  size_t R, C;
  unsigned r2, r3, c2, c3;
  fft_array_params(&R, &r2, &r3, &C, &c2, &c3, k2, k3);
  return fft_cols_rows_pattern(R, r2, r3, C, c2, c3);
}



// select array decomposition for fft_array_vbuf
void fft_array_vbuf_params(size_t* R, unsigned* r2, unsigned* r3,
			   size_t* C, unsigned* c2, unsigned* c3,
			   unsigned k2, unsigned k3, size_t block_size)
{
  unsigned _r2, _r3, _c2, _c3;
  size_t _R, _C;
  fft_array_params(&_R, &_r2, &_r3, &_C, &_c2, &_c3, k2, k3);

  // now ensure block_size divisible by C, i.e. divisible by 2^c2 and 3^c3
  while (block_size % pow23(_c2, 0))
    _c2--, _r2++, _C /= 2, _R *= 2;

  while (block_size % pow23(0, _c3))
    _c3--, _r3++, _C /= 3, _R *= 3;

  *R = _R;
  *r2 = _r2;
  *r3 = _r3;
  *C = _C;
  *c2 = _c2;
  *c3 = _c3;
}



void fft_array_vbuf(vbuf_t* rop, size_t K, unsigned k2, unsigned k3,
		    vbuf_t* op, size_t n, uint64_t w, uint64_t p,
		    int num_threads)
{
  size_t R, C;
  unsigned r2, r3, c2, c3;
  fft_array_vbuf_params(&R, &r2, &r3, &C, &c2, &c3, k2, k3, rop->block_size);
  fft_cols_vbuf(rop, R, r2, r3, C, c2, c3, op, n, w, p, num_threads);
  fft_rows_vbuf(rop, R, r2, r3, C, c2, c3, w, p, num_threads);
}


size_t fft_array_vbuf_pattern(size_t K, unsigned k2, unsigned k3,
			      size_t block_size)
{
  size_t R, C;
  unsigned r2, r3, c2, c3;
  fft_array_vbuf_params(&R, &r2, &r3, &C, &c2, &c3, k2, k3, block_size);
  return fft_cols_rows_vbuf_pattern(R, r2, r3, C, c2, c3, block_size);
}




void ifft_cols_vbuf(vbuf_t* op,
		    size_t R, unsigned r2, unsigned r3,
		    size_t C, unsigned c2, unsigned c3,
		    uint64_t w, uint64_t p, int num_threads)
{
  uint64_t pinv = calc_pinv(p);
  uint64_t wC = pow_mod_pinv(w, C, p, pinv);

  // use ifft_base for column transforms if they are small enough
  int use_ifft_base = (R < FFT_ARRAY_THRESHOLD);

  // compute root tables if necessary
  wtab_t wtab;
  if (use_ifft_base)
    wtab_init(&wtab, R, r2, r3, wC, p);

  // number of columns to extract on each iteration
  const size_t clump = CACHE_LINE / sizeof(uint64_t);

  // number of rows in each memory pool block
  size_t block_rows = op->block_size / C;
  assert(op->block_size == C * block_rows);

#pragma omp parallel num_threads(num_threads)
  {
    uint64_t temp[clump*clump];
    uint64_t* cols = (uint64_t*) safe_malloc(R * clump * sizeof(uint64_t));

#pragma omp for schedule(dynamic,1)
    // handle columns in blocks
    for (ptrdiff_t i = 0; i < C; i += clump)
      {
	// c = number of columns to handle on this iteration
	size_t c = MIN(clump, C - i);

	// transfer c columns from op to cols
	{
	  for (size_t j = 0; j < R; j += block_rows)
	    {
	      // r = number of rows to transfer on this iteration
	      size_t r = MIN(block_rows, R - j);

	      for (size_t jj = 0; jj < r; jj += clump)
		{
		  size_t rr = MIN(clump, r - jj);

		  // transpose from op to cols via temp
		  copy_block(temp, clump, vbuf_get(op, i + (j + jj) * C),
			     C, rr, c);
		  transpose_block(cols + j + jj, R, temp, clump, rr, c);
		}
	    }
	}
	
	// transform columns
	for (size_t ii = 0; ii < c; ii++)
	  {
	    if (use_ifft_base)
	      ifft_base(cols + ii*R, R, r2, r3, &wtab, p);
	    else
	      ifft_array(cols + ii*R, R, r2, r3, wC, p, 1);
	  }

	// transfer columns from cols to rop
	{
	  // handle rows in blocks
	  for (size_t j = 0; j < R; j += block_rows)
	    {
	      // r = number of rows to transfer on this iteration
	      size_t r = MIN(block_rows, R - j);

	      for (size_t jj = 0; jj < r; jj += clump)
		{
		  size_t rr = MIN(clump, r - jj);

		  // transpose from cols to op via temp
		  transpose_block(temp, clump, cols + j + jj, R, c, rr);
		  copy_block(vbuf_get(op, i + (j + jj) * C),
			     C, temp, clump, rr, c);
		}
	    }
	}
      }

    free(cols);
  }

  if (use_ifft_base)
    wtab_clear(&wtab);
}


void ifft_cols(uint64_t* op,
	       size_t R, unsigned r2, unsigned r3,
	       size_t C, unsigned c2, unsigned c3,
	       uint64_t w, uint64_t p, int num_threads)
{
  // create wrapper vbufs and use vbuf version
  vbuf_t vbuf;
  vbuf_init_wrap(&vbuf, R * C, op);
  ifft_cols_vbuf(&vbuf, R, r2, r3, C, c2, c3, w, p, num_threads);
  vbuf_clear(&vbuf);
}



void ifft_rows_vbuf(vbuf_t* op,
		    size_t R, unsigned r2, unsigned r3,
		    size_t C, unsigned c2, unsigned c3,
		    uint64_t w, uint64_t p, int num_threads)
{
  assert(op->block_size % C == 0);
  uint64_t pinv = calc_pinv(p);
  uint64_t pinvb = calc_pinvb(p);
  uint64_t wR = pow_mod_pinv(w, R, p, pinv);

  // use ifft_base for row transforms if they are small enough
  int use_ifft_base = (C < FFT_ARRAY_THRESHOLD);

  // compute root tables if necessary
  wtab_t wtab;
  if (use_ifft_base)
    wtab_init(&wtab, C, c2, c3, wR, p);

  // compute t[i] = basic twiddle factor to use for row #i
  uint64_t* t = safe_malloc(R * sizeof(uint64_t));
  compute_twiddles(t, R, r2, r3, w, p, pinv);

  // transform rows
#pragma omp parallel for schedule(dynamic,1) num_threads(num_threads)
  for (ptrdiff_t i = 0; i < R; i++)
    {
      uint64_t* x = vbuf_get(op, i * C);

      // transform row
      if (use_ifft_base)
	ifft_base(x, C, c2, c3, &wtab, p);
      else
	ifft_array(x, C, c2, c3, wR, p, 1);

      // apply twiddle factors
      if (i != 0)
	apply_twiddles(x, C, t[i], p, pinvb);
    }

  safe_free(t);
  if (use_ifft_base)
    wtab_clear(&wtab);
}



void ifft_rows(uint64_t* op,
	       size_t R, unsigned r2, unsigned r3,
	       size_t C, unsigned c2, unsigned c3,
	       uint64_t w, uint64_t p, int num_threads)
{
  // create wrapper vbuf and use vbuf version
  vbuf_t vbuf;
  vbuf_init_wrap(&vbuf, R * C, op);
  ifft_rows_vbuf(&vbuf, R, r2, r3, C, c2, c3, w, p, num_threads);
  vbuf_clear(&vbuf);
}


void ifft_array(uint64_t* op, size_t K, unsigned k2, unsigned k3,
		uint64_t w, uint64_t p, int num_threads)
{
  size_t R, C;
  unsigned r2, r3, c2, c3;
  fft_array_params(&R, &r2, &r3, &C, &c2, &c3, k2, k3);
  ifft_rows(op, R, r2, r3, C, c2, c3, w, p, num_threads);
  ifft_cols(op, R, r2, r3, C, c2, c3, w, p, num_threads);
}


void ifft_array_vbuf(vbuf_t* op, size_t K, unsigned k2, unsigned k3,
		     uint64_t w, uint64_t p, int num_threads)
{
  size_t R, C;
  unsigned r2, r3, c2, c3;
  fft_array_vbuf_params(&R, &r2, &r3, &C, &c2, &c3, k2, k3, op->block_size);
  ifft_rows_vbuf(op, R, r2, r3, C, c2, c3, w, p, num_threads);
  ifft_cols_vbuf(op, R, r2, r3, C, c2, c3, w, p, num_threads);
}



void conv_cols_rows_vbuf(vbuf_t* rop,
			 size_t R, unsigned r2, unsigned r3,
			 size_t C, unsigned c2, unsigned c3,
			 vbuf_t* op1, size_t n1, vbuf_t* temp1,
			 vbuf_t* op2, size_t n2, vbuf_t* temp2,
			 uint64_t w, uint64_t p, int num_threads)
{
  assert(rop->block_size % C == 0);

  int squaring = (op1 == op2) && (n1 == n2);

  uint64_t pinv = calc_pinv(p);
  uint64_t pinvb = calc_pinvb(p);
  uint64_t wR = pow_mod_pinv(w, R, p, pinv);
  uint64_t winv = inv_mod(w, p);
  uint64_t winvR = pow_mod_pinv(winv, R, p, pinv);

  // use fft_base for row transforms if they are small enough
  int use_fft_base = (C < FFT_ARRAY_THRESHOLD);

  // transform columns
  fft_cols_vbuf(temp1, R, r2, r3, C, c2, c3, op1, n1, w, p, num_threads);
  if (!squaring)
    fft_cols_vbuf(temp2, R, r2, r3, C, c2, c3, op2, n2, w, p, num_threads);

  // todo: could use omp parallel sections for following tables...
  
  // compute root tables if necessary
  wtab_t wtab, winvtab;
  if (use_fft_base)
    {
      wtab_init(&wtab, C, c2, c3, wR, p);
      wtab_init(&winvtab, C, c2, c3, winvR, p);
    }

  // compute t[i] = basic twiddle factor to use for row #i
  uint64_t* t = safe_malloc(R * sizeof(uint64_t));
  uint64_t* tinv = safe_malloc(R * sizeof(uint64_t));
  compute_twiddles(t, R, r2, r3, w, p, pinv);
  compute_twiddles(tinv, R, r2, r3, winv, p, pinv);

  // convolve rows
#pragma omp parallel for schedule(dynamic,1) num_threads(num_threads)
  for (ptrdiff_t i = 0; i < R; i++)
    {
      uint64_t* x1 = vbuf_get(temp1, i * C);
      uint64_t* x2 = squaring ? x1 : vbuf_get(temp2, i * C);
      uint64_t* y = vbuf_get(rop, i * C);

      // apply twiddle factors
      if (i != 0)
	{
	  apply_twiddles(x1, C, t[i], p, pinvb);
	  if (!squaring)
	    apply_twiddles(x2, C, t[i], p, pinvb);
	}

      // convolve this row
      if (use_fft_base)
	conv_base(y, C, c2, c3, x1, C, x1, x2, C, x2,
		  &wtab, &winvtab, p, pinvb);
      else
	conv_array(y, C, c2, c3, x1, C, x1, x2, C, x2, wR, p, 1);

      // apply outgoing twiddle factors
      if (i != 0)
	apply_twiddles(y, C, tinv[i], p, pinvb);
    }

  safe_free(t);
  safe_free(tinv);
  if (use_fft_base)
    {
      wtab_clear(&wtab);
      wtab_clear(&winvtab);
    }

  // inverse transform columns
  ifft_cols_vbuf(rop, R, r2, r3, C, c2, c3, winv, p, num_threads);
}



void conv_cols_rows(uint64_t* rop,
		    size_t R, unsigned r2, unsigned r3,
		    size_t C, unsigned c2, unsigned c3,
		    uint64_t* op1, size_t n1, uint64_t* temp1,
		    uint64_t* op2, size_t n2, uint64_t* temp2,
		    uint64_t w, uint64_t p, int num_threads)
{
  // create wrapper vbufs and use vbuf version
  int squaring = (op1 == op2) && (n1 == n2);

  vbuf_t vbuf_rop, vbuf_op1, vbuf_op2, vbuf_temp1, vbuf_temp2;
  vbuf_init_wrap(&vbuf_rop, R * C, rop);
  vbuf_init_wrap(&vbuf_op1, R * C, op1);
  vbuf_init_wrap(&vbuf_temp1, R * C, temp1);
  if (!squaring)
    {
      vbuf_init_wrap(&vbuf_op2, R * C, op2);
      vbuf_init_wrap(&vbuf_temp2, R * C, temp2);
    }

  conv_cols_rows_vbuf(&vbuf_rop, R, r2, r3, C, c2, c3,
		      &vbuf_op1, n1, &vbuf_temp1,
		      squaring ? &vbuf_op1 : &vbuf_op2, n2,
		      squaring ? &vbuf_temp1 : &vbuf_temp2,
		      w, p, num_threads);

  vbuf_clear(&vbuf_rop);
  vbuf_clear(&vbuf_op1);
  vbuf_clear(&vbuf_temp1);
  if (!squaring)
    {
      vbuf_clear(&vbuf_op2);
      vbuf_clear(&vbuf_temp2);
    }
}



void conv_array(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
		uint64_t* op1, size_t n1, uint64_t* temp1,
		uint64_t* op2, size_t n2, uint64_t* temp2,
		uint64_t w, uint64_t p, int num_threads)
{
  size_t R, C;
  unsigned r2, r3, c2, c3;
  fft_array_params(&R, &r2, &r3, &C, &c2, &c3, k2, k3);
  conv_cols_rows(rop, R, r2, r3, C, c2, c3,
		 op1, n1, temp1, op2, n2, temp2, w, p, num_threads);
}


void conv_array_vbuf(vbuf_t* rop, size_t K, unsigned k2, unsigned k3,
		     vbuf_t* op1, size_t n1, vbuf_t* temp1,
		     vbuf_t* op2, size_t n2, vbuf_t* temp2,
		     uint64_t w, uint64_t p, int num_threads)
{
  size_t R, C;
  unsigned r2, r3, c2, c3;
  fft_array_vbuf_params(&R, &r2, &r3, &C, &c2, &c3, k2, k3, rop->block_size);
  conv_cols_rows_vbuf(rop, R, r2, r3, C, c2, c3,
		      op1, n1, temp1, op2, n2, temp2, w, p, num_threads);
}
