/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "ntt-internal.h"


// special case of split for 0 < r < 64
void split_1(mp_limb_t* rop, size_t m, mp_limb_t* op, size_t n, unsigned r)
{
  mp_limb_t mask = ((mp_limb_t) 1 << r) - 1;
  size_t dst = 0;

  // the loop below reads from bit index 64*src + src_bit
  size_t src = 0;
  unsigned src_bit = 0;

  while (dst < m && src < n - 1)
    {
      mp_limb_t x0 = op[src];
      mp_limb_t x1 = op[src + 1];
      x0 = (x0 >> src_bit) | (x1 << 1) << (63 - src_bit);
      rop[dst++] = x0 & mask;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m && src < n)
    {
      mp_limb_t x0 = op[src] >> src_bit;
      rop[dst++] = x0 & mask;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m)
    rop[dst++] = 0;
}


// special case of split for 64 < r < 128
void split_2(mp_limb_t* rop, size_t m, mp_limb_t* op, size_t n, unsigned r)
{
  mp_limb_t mask = ((mp_limb_t) 1 << (r - 64)) - 1;
  size_t dst = 0;

  // the loop below reads from bit index 64*src + src_bit
  size_t src = 0;
  unsigned src_bit = 0;

  m *= 2;

  while (dst < m && (ptrdiff_t) src < (ptrdiff_t) n - 2)
    {
      mp_limb_t x0 = op[src];
      mp_limb_t x1 = op[src + 1];
      mp_limb_t x2 = op[src + 2];
      x0 = (x0 >> src_bit) | (x1 << 1) << (63 - src_bit);
      x1 = (x1 >> src_bit) | (x2 << 1) << (63 - src_bit);
      rop[dst++] = x0;
      rop[dst++] = x1 & mask;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m && src < n - 1)
    {
      mp_limb_t x0 = op[src];
      mp_limb_t x1 = op[src + 1];
      x0 = (x0 >> src_bit) | (x1 << 1) << (63 - src_bit);
      x1 >>= src_bit;
      rop[dst++] = x0;
      rop[dst++] = x1 & mask;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m && src < n)
    {
      mp_limb_t x0 = op[src] >> src_bit;
      rop[dst++] = x0;
      rop[dst++] = 0;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m)
    {
      rop[dst++] = 0;
      rop[dst++] = 0;
    }
}


// special case of split for 128 < r < 192
void split_3(mp_limb_t* rop, size_t m, mp_limb_t* op, size_t n, unsigned r)
{
  mp_limb_t mask = ((mp_limb_t) 1 << (r - 128)) - 1;
  size_t dst = 0;

  // the loop below reads from bit index 64*src + src_bit
  size_t src = 0;
  unsigned src_bit = 0;

  m *= 3;

  while (dst < m && (ptrdiff_t) src < (ptrdiff_t) n - 3)
    {
      mp_limb_t x0 = op[src];
      mp_limb_t x1 = op[src + 1];
      mp_limb_t x2 = op[src + 2];
      mp_limb_t x3 = op[src + 3];
      x0 = (x0 >> src_bit) | (x1 << 1) << (63 - src_bit);
      x1 = (x1 >> src_bit) | (x2 << 1) << (63 - src_bit);
      x2 = (x2 >> src_bit) | (x3 << 1) << (63 - src_bit);
      rop[dst++] = x0;
      rop[dst++] = x1;
      rop[dst++] = x2 & mask;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m && (ptrdiff_t) src < (ptrdiff_t) n - 2)
    {
      mp_limb_t x0 = op[src];
      mp_limb_t x1 = op[src + 1];
      mp_limb_t x2 = op[src + 2];
      x0 = (x0 >> src_bit) | (x1 << 1) << (63 - src_bit);
      x1 = (x1 >> src_bit) | (x2 << 1) << (63 - src_bit);
      x2 >>= src_bit;
      rop[dst++] = x0;
      rop[dst++] = x1;
      rop[dst++] = x2 & mask;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m && src < n - 1)
    {
      mp_limb_t x0 = op[src];
      mp_limb_t x1 = op[src + 1];
      x0 = (x0 >> src_bit) | (x1 << 1) << (63 - src_bit);
      x1 >>= src_bit;
      rop[dst++] = x0;
      rop[dst++] = x1;
      rop[dst++] = 0;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m && src < n)
    {
      mp_limb_t x0 = op[src] >> src_bit;
      rop[dst++] = x0;
      rop[dst++] = 0;
      rop[dst++] = 0;
      src_bit += r;
      src += src_bit / 64;
      src_bit &= 63;
    }

  while (dst < m)
    {
      rop[dst++] = 0;
      rop[dst++] = 0;
      rop[dst++] = 0;
    }
}


void split(mp_limb_t* rop, size_t m, mp_limb_t* op, size_t n, unsigned r)
{
  if (n == 0)
    {
      unsigned s = (r + 63) / 64;
      zero(rop, m * s);
    }
  else if (r % 64 == 0)
    {
      // just split along limb boundaries
      unsigned t = r / 64;
      size_t d = MIN(t * m, n);
      memcpy(rop, op, d * sizeof(mp_limb_t));
      zero(rop + d, t * m - d);
    }
  else if (r < 64)
    split_1(rop, m, op, n, r);
  else if (r < 128)
    split_2(rop, m, op, n, r);
  else if (r < 192)
    split_3(rop, m, op, n, r);
  else
    fatal_error("split: unexpected r = %u", r);
}


// special case of reduce() for t == 1
void reduce_1(uint64_t* rop, mp_limb_t* op, size_t m,
	      uint64_t u, uint64_t p, uint64_t pinv, uint64_t pinvb)
{
  uint64_t upinv = calc_ypinv(u, p, pinv);

  for (ptrdiff_t i = m - 1; i >= 0; i--)
    rop[i] = mul_mod_ypinv_lazy(op[i], u, upinv, p);
}


// special case of reduce() for t == 2
void reduce_2(uint64_t* rop, mp_limb_t* op, size_t m,
	      uint64_t u, uint64_t p, uint64_t pinv, uint64_t pinvb)
{
  // incorporate extra factor of 2^64
  u = mul_mod_pinv(u, -p, p, pinv);
  uint64_t upinv = calc_ypinv(u, p, pinv);

  for (ptrdiff_t i = m - 1; i >= 0; i--)
    {
      uint64_t x0 = op[2*i];
      uint64_t x1 = op[2*i + 1];
      x1 -= (x1 >= 2*p) ? (2*p) : 0;
      uint64_t h1 = mulhi(x0 * pinvb, p);
      uint64_t y1 = x1 - h1 + p;
      rop[i] = mul_mod_ypinv_lazy(y1, u, upinv, p);
    }
}


// special case of reduce() for t == 3
void reduce_3(uint64_t* rop, mp_limb_t* op, size_t m,
	      uint64_t u, uint64_t p, uint64_t pinv, uint64_t pinvb)
{
  // incorporate extra factor of 2^128
  u = mul_mod_pinv(u, -p, p, pinv);
  u = mul_mod_pinv(u, -p, p, pinv);
  uint64_t upinv = calc_ypinv(u, p, pinv);

  for (ptrdiff_t i = m - 1; i >= 0; i--)
    {
      uint64_t x0 = op[3*i];
      uint64_t x1 = op[3*i + 1];
      uint64_t x2 = op[3*i + 2];
      x2 -= (x2 >= 2*p) ? (2*p) : 0;
      uint64_t h1 = mulhi(x0 * pinvb, p);
      uint64_t y1 = x1 - h1;
      y1 += (x1 < h1) ? p : 0;
      uint64_t h2 = mulhi(y1 * pinvb, p);
      uint64_t y2 = x2 - h2 + p;
      rop[i] = mul_mod_ypinv_lazy(y2, u, upinv, p);
    }
}


void reduce(uint64_t* rop, mp_limb_t* op, size_t m, unsigned t,
	    uint64_t u, uint64_t p, uint64_t pinv, uint64_t pinvb)
{
  if (t == 1)
    reduce_1(rop, op, m, u, p, pinv, pinvb);
  else if (t == 2)
    reduce_2(rop, op, m, u, p, pinv, pinvb);
  else if (t == 3)
    reduce_3(rop, op, m, u, p, pinv, pinvb);
  else
    fatal_error("reduce: unexpected t = %u", t);
}



void split_reduce_vbuf(vbuf_t* rop, size_t m, mp_limb_t* op, size_t n,
		       unsigned r, uint64_t* u, unsigned num_primes,
		       int destroy_op, int num_threads)
{
  size_t block_size = rop->block_size;
  assert(block_size % 64 == 0);

  size_t block_size2 = CACHE_SIZE / (2 * sizeof(mp_limb_t) * num_primes);
  block_size2 = MIN(block_size, block_size2);
  block_size2 -= block_size2 % 64;   // ensure multiple of 64

  size_t t = (r + 63) / 64;

  // set up a vbuf to manage memory in op, so that we can use vbuf_free_region()
  vbuf_t op_vbuf;
  size_t op_blocks = n / block_size;
  if (op_blocks == 0)
    destroy_op = 0;
  if (destroy_op)
    {
      vbuf_init(&op_vbuf, op_blocks, rop->pool);
      for (size_t i = 0; i < op_blocks; i++)
	op_vbuf.blocks[i] = (uint64_t*) (op + i * block_size);
    }

#pragma omp parallel num_threads(num_threads)
  {
    mp_limb_t temp[block_size2 * t];
    uint64_t* ptr[num_primes];

#pragma omp for schedule(dynamic,1)
    for (ptrdiff_t j = 0; j < m; j += block_size)
      {
	// s = number of output coefficients to write in this block
	size_t s = MIN(block_size, m - j);

	for (unsigned i = 0; i < num_primes; i++)
	  ptr[i] = vbuf_get(&rop[i], j);

	for (ptrdiff_t jj = 0; jj < s; jj += block_size2)
	  {
	    // ss = number of output coefficients to write in this block
	    size_t ss = MIN(block_size2, s - jj);

	    // k = starting limb in op for this block
	    size_t k = (j + jj) / 64 * r;

	    split(temp, ss, op + k, (k < n) ? (n - k) : 0, r);

	    for (unsigned i = 0; i < num_primes; i++)
	      reduce(ptr[i] + jj, temp, ss, t, u[i], global_p[i],
		     global_pinv[i], global_pinvb[i]);
	  }
	
	// give memory from op to pool if possible
	if (destroy_op)
	  {
	    size_t k1 = j / 64 * r;
	    size_t k2 = MIN((j + s) / 64 * r, n);
	    if (k1 < k2)
	      vbuf_free_region(&op_vbuf, k1, k2);
	  }
      }
  }

  if (destroy_op)
    {
      for (size_t i = 0; i < op_blocks; i++)
	op_vbuf.blocks[i] = NULL;
      vbuf_clear(&op_vbuf);
    }
}



void split_reduce(uint64_t** rop, size_t m, mp_limb_t* op, size_t n,
		  unsigned r, uint64_t* u, unsigned num_primes, int num_threads)
{
  // create wrapper vbufs and use vbuf version
  vbuf_t vbuf_rop[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    vbuf_init_wrap(&vbuf_rop[i], (m / 64 + 1) * 64, rop[i]);

  split_reduce_vbuf(vbuf_rop, m, op, n, r, u, num_primes, 0, num_threads);

  for (unsigned i = 0; i < num_primes; i++)
    vbuf_clear(&vbuf_rop[i]);
}



uint64_t global_crt_s[MAX_NUM_PRIMES];
uint64_t global_crt_spinv[MAX_NUM_PRIMES];
mp_limb_t global_crt_u[MAX_NUM_PRIMES][MAX_NUM_PRIMES];


void crt(mp_limb_t* rop, uint64_t** op, size_t m, unsigned num_primes)
{
  if (num_primes == 1)
    {
      // only one prime, simply reduce data mod p
      uint64_t p = global_p[0];
      uint64_t* src = op[0];
      for (ptrdiff_t j = m - 1; j >= 0; j--)
	{
	  uint64_t x = src[j];
	  x -= (x >= 2*p) ? (2*p) : 0;
	  x -= (x >= p) ? p : 0;
	  rop[j] = x;
	}
      return;
    }

  // first pass: CRT data from first two primes
  {
    uint64_t* op0 = op[0];
    uint64_t* op1 = op[1];
    mp_limb_t* dest = rop + (m-1) * num_primes;
    uint64_t p0 = global_p[0];
    uint64_t p1 = global_p[1];
    uint64_t s = global_crt_s[1];
    uint64_t spinv = global_crt_spinv[1];

    for (ptrdiff_t j = m - 1; j >= 0; j--, dest -= num_primes)
      {
	uint64_t x = op0[j];                 // [0, 4*p0)
	x -= (x >= 2*p0) ? (2*p0) : 0;
	x -= (x >= p0) ? p0 : 0;             // [0, p0)
	uint64_t y = op1[j];                 // [0, 4*p1)
	uint64_t r = x - y;
	r += (x < y) ? (4*p1) : 0;           // [0, 2^64)
	r = mul_mod_ypinv(r, s, spinv, p1);  // r = (y - x)/p0 mod p1
	uint64_t s1, s0;
#if AVOID_128_BIT
	MUL128(s1, s0, r, p0);
	s0 += x;
	s1 += (s0 < x);
#else
	uint128_t s = (uint128_t) r * p0 + x;
              // s = x0 mod p0, x1 mod p1, s in [0, p0*p1)
	s0 = (uint64_t) s;
	s1 = s >> 64;
#endif
	dest[0] = s0;
	dest[1] = s1;
      }
  }

  // remaining primes
  for (unsigned i = 2; i < num_primes; i++)
    {
      uint64_t* src = op[i];
      mp_limb_t* dest = rop + (m - 1) * num_primes;
      uint64_t p = global_p[i];
      uint64_t pinvb = global_pinvb[i];
      mp_limb_t* u = global_crt_u[i];
      uint64_t s = global_crt_s[i];
      uint64_t spinv = global_crt_spinv[i];

      for (ptrdiff_t j = m - 1; j >= 0; j--, dest -= num_primes)
	{
	  // Let X be the previous residue, stored in the first i limbs of dest,
	  // i.e. X = op[k][j] mod p[k] for 0 <= k < i, and
	  // 0 <= X < p[0]*...*p[i-1].

	  // Compute r = (X - op[i][j]) / (2^64)^(i-1) mod p[i], in [0, 4*p[i]).
	  uint64_t x = dest[0];     // [0, 2^64)
	  uint64_t y = src[j];      // [0, 4p)
	  uint64_t r = x - y;
	  r += (x < y) ? (4*p) : 0;
	  for (unsigned k = 1; k < i; k++)
	    {
	      uint64_t q = r * pinvb;
	      uint64_t h = mulhi(q, p);                  // [0, p)
	      uint64_t x = dest[k];                      // [0, 2^64)
	      r = x - h;
	      r += (x < h) ? p : 0;
	    }

	  // r := (op[i][j] - X) / (p[0]*...*p[i-1]) mod p[i], in [0, p[i])
	  r = mul_mod_ypinv(r, s, spinv, p);

	  // add r*p[0]*...*p[i-1] to X
	  dest[i] = mpn_addmul_1(dest, u, i, r);
	}
    }
}



// same as crt(), but reads from a vbuf, starting at position src
void crt_from_vbuf(mp_limb_t* rop, vbuf_t* op, size_t src,
		   size_t m, unsigned num_primes)
{
  size_t block_size = op->block_size;
  uint64_t* ptr[num_primes];

  while (m > 0)
    {
      // k = # input coefficients to handle on this iteration
      size_t k = MIN(m, (src / block_size + 1) * block_size - src);

      for (unsigned i = 0; i < num_primes; i++)
	ptr[i] = vbuf_get(&op[i], src);

      crt(rop, ptr, k, num_primes);

      rop += num_primes * k;
      src += k;
      m -= k;
    }
}


void recompose(mp_limb_t* rop, size_t n, mp_limb_t* op, size_t m, unsigned s,
	       unsigned r, unsigned r0)
{
  // reduce to case r0 < 64
  while (r0 >= 64)
    {
      r0 -= 64;
      *rop++ = 0;
      if (--n == 0)
	return;
    }

  mp_limb_t temp[s + 1];

  // on each iteration of main loop, dst increases by either u or u + 1
  size_t u = r / 64;
  assert(s >= u);

  ptrdiff_t src = 0;
  ptrdiff_t dst = 0;
  unsigned dst_bit = r0;

  zero(rop + dst, MIN(s - u + 1, n - dst));

  // at the beginning of each iteration, the first dst + s - u + 1 limbs of
  // rop contain valid data
  while (src < m && dst < (ptrdiff_t) (n - s - 1))
    {
      if (dst_bit != 0)
	temp[s] = mpn_lshift(temp, op, s, dst_bit);
      else
	{
	  memcpy(temp, op, s * sizeof(mp_limb_t));
	  temp[s] = 0;
	}

      mp_limb_t cy = mpn_add_n(rop + dst, rop + dst, temp, s - u + 1);
      if (u > 0)
	mpn_add_1(rop + dst + s - u + 1, temp + s - u + 1, u, cy);
      rop[dst + s + 1] = 0;

      dst_bit += r;
      dst += dst_bit / 64;
      dst_bit &= 63;
      src++;
      op += s;
    }

  // same loop as above, but don't write beyond end of rop
  while (src < m && dst < n)
    {
      if (dst_bit != 0)
	temp[s] = mpn_lshift(temp, op, s, dst_bit);
      else
	{
	  memcpy(temp, op, s * sizeof(mp_limb_t));
	  temp[s] = 0;
	}

      mp_limb_t cy = mpn_add_n(rop + dst, rop + dst, temp,
			       MIN(s - u + 1, n - dst));
      if (u > 0 && dst + s - u + 1 < n)
	mpn_add_1(rop + dst + s - u + 1, temp + s - u + 1,
		  MIN(u, n - (dst + s - u + 1)), cy);
      if (dst + s + 1 < n)
	rop[dst + s + 1] = 0;

      dst_bit += r;
      dst += dst_bit / 64;
      dst_bit &= 63;
      src++;
      op += s;
    }

  // zero rest of output
  if (dst + s - u + 1 < n)
    zero(rop + dst + s - u + 1, n - (dst + s - u + 1));
}


// zeroes out n words of op, starting at index dst
void vbuf_zero(vbuf_t* op, size_t dst, size_t n)
{
  size_t block_size = op->block_size;
  size_t i = dst / block_size;
  size_t j = dst - i * block_size;

  while (n > 0)
    {
      size_t k = MIN(n, block_size - j);
      zero(vbuf_get_block(op, i) + j, k);
      n -= k;
      j = 0;
      i++;
    }
}


// adds {op,m} into {rop+dst,n}.
// must have 0 <= m <= n
void vbuf_mpn_add(vbuf_t* rop, size_t dst, size_t n, mp_limb_t* op, size_t m)
{
  size_t block_size = rop->block_size;
  size_t i = dst / block_size;
  size_t j = dst - i * block_size;
  mp_limb_t cy = 0;

  while (m > 0)
    {
      size_t k = MIN(m, block_size - j);
      mp_limb_t* ptr = (mp_limb_t*) vbuf_get_block(rop, i) + j;
      cy = mpn_add_1(ptr, ptr, k, cy);
      cy += mpn_add_n(ptr, ptr, op, k);
      op += k;
      m -= k;
      n -= k;
      j += k;
      if (j == block_size)
	{
	  j = 0;
	  i++;
	}
    }

  // propagate carry
  while (n > 0 && cy)
    {
      size_t k = MIN(n, block_size - j);
      mp_limb_t* ptr = (mp_limb_t*) vbuf_get_block(rop, i) + j;
      cy = mpn_add_1(ptr, ptr, k, cy);
      n -= k;
      j = 0;
      i++;
    }
}


void crt_recompose_vbuf(vbuf_t* rop, size_t n, unsigned r,
			vbuf_t* op, size_t m, int destroy_op,
			unsigned num_primes, int num_threads)
{
  size_t block_size = rop->block_size;
  size_t block_size2 = CACHE_SIZE / (2 * sizeof(mp_limb_t) * num_primes);

  size_t u = 2 * num_primes + 2;    // size of each saved chunk
  size_t num_blocks = n / block_size + 1;
  size_t num_blocks2 = block_size / block_size2 + 1;
  mp_limb_t* saved = safe_malloc(u * num_blocks * num_blocks2
				 * sizeof(mp_limb_t));

  // split up output into memory pool blocks
#pragma omp parallel num_threads(num_threads)
  {
    mp_limb_t* temp = safe_malloc((64 * block_size2 / r + 10)
				  * num_primes * sizeof(mp_limb_t));

#pragma omp for schedule(dynamic,1)     
    for (ptrdiff_t j = 0; j < n; j += block_size)
      {
	mp_limb_t* saved_ptr = saved + u * (j / block_size) * num_blocks2;

	// k = number of output limbs to write on this iteration
	size_t k = MIN(n - j, block_size);

	// split up further into blocks for better locality
	for (size_t jj = 0; jj < k; jj += block_size2, saved_ptr += u)
	  {
	    // h = number of output limbs to write on this iteration
	    size_t h = MIN(k - jj, block_size2);

	    // coefficients between t1 and t2 are contained entirely in the
	    // output region of length h starting at limb j + jj
	    size_t t1 = (64 * (j + jj) + r - 1) / r;
	    size_t t2 = 64 * (j + jj + h) - 62 * num_primes - 1 + r;
	    if ((ptrdiff_t) t2 < 0)
	      t2 = 0;
	    t2 /= r;
	    if (t2 < t1)
	      t2 = t1;

	    // crt/recompose those coefficients, writing exactly h limbs
	    if (t1 < m)
	      {
		crt_from_vbuf(temp, op, t1, MIN(t2, m) - t1, num_primes);
		recompose((mp_limb_t*) vbuf_get(rop, j + jj), h, temp,
			  MIN(t2, m) - t1, num_primes, r,
			  t1 * r - 64 * (j + jj));
	      }
	    else
	      vbuf_zero(rop, j + jj, h);

	    // t3 = value of t1 for next block
	    size_t t3 = (64 * (j + jj + h) + r - 1) / r;

	    // crt/recompose coefficients between t2 and t3, writing to "saved",
	    // to be added back later
	    if (t2 < m)
	      {
		crt_from_vbuf(temp, op, t2, MIN(t3, m) - t2, num_primes);
		recompose(saved_ptr, u, temp, MIN(t3, m) - t2, num_primes,
			  r, (t2 * r) % 64);
	      }
	    else
	      zero(saved_ptr, u);

	    // return memory to pool if possible
	    if (destroy_op && t1 < MIN(t3, m))
	      {
		for (unsigned i = 0; i < num_primes; i++)
		  vbuf_free_region(&op[i], t1, MIN(t3, m));
	      }
	  }
      }

    safe_free(temp);
  }

  // add back saved blocks
  for (size_t j = 0; j < n; j += block_size)
    {
      mp_limb_t* saved_ptr = saved + u * (j / block_size) * num_blocks2;
      size_t k = MIN(n - j, block_size);
      for (size_t jj = 0; jj < k; jj += block_size2, saved_ptr += u)
	{
	  size_t h = MIN(k - jj, block_size2);
	  size_t t1 = (64 * (j + jj) + r - 1) / r;
	  size_t t2 = 64 * (j + jj + h) - 62 * num_primes - 1 + r;
	  if ((ptrdiff_t) t2 < 0)
	    t2 = 0;
	  t2 /= r;
	  if (t2 < t1)
	    t2 = t1;
	  size_t c = (r * t2) / 64;
	  if (c < n)
	    vbuf_mpn_add(rop, c, n - c, saved_ptr, MIN(u, n - c));
	}
    }

  safe_free(saved);
}



void crt_recompose(mp_limb_t* rop, size_t n, unsigned r,
		   uint64_t** op, size_t m, unsigned num_primes,
		   int num_threads)
{
  // create wrapper vbufs and use vbuf version
  vbuf_t vbuf_rop, vbuf_op[num_primes];
  vbuf_init_wrap(&vbuf_rop, n, (uint64_t*) rop);
  for (unsigned i = 0; i < num_primes; i++)
    vbuf_init_wrap(&vbuf_op[i], m, op[i]);

  crt_recompose_vbuf(&vbuf_rop, n, r, vbuf_op, m, 0, num_primes, num_threads);

  vbuf_clear(&vbuf_rop);
  for (unsigned i = 0; i < num_primes; i++)
    vbuf_clear(&vbuf_op[i]);
}



void ntt_mpn_mul(mp_limb_t* rop, mp_limb_t* op1, size_t n1,
		 mp_limb_t* op2, size_t n2)
{
  ntt_mpn_mul_bonus(rop, op1, n1, 1, op2, n2, 1, 0, 1);
}



// assumes n1 >= n2
void ntt_mpn_mul_gmp(mp_limb_t* rop, mp_limb_t* op1, size_t n1,
		     mp_limb_t* op2, size_t n2)
{
  size_t n3 = n1 + n2;
  if (overlap(rop, n3, op1, n1) || overlap(rop, n3, op2, n2))
    {
      // input overlaps output; put output into a temporary
      mp_limb_t* temp = safe_malloc(n3 * sizeof(mp_limb_t));
      mpn_mul(temp, op1, n1, op2, n2);
      memcpy(rop, temp, n3 * sizeof(mp_limb_t));
      safe_free(temp);
    }
  else
    // no overlap; call GMP directly
    mpn_mul(rop, op1, n1, op2, n2);
}


double global_bit_bound[MAX_NUM_PRIMES + 1];


/*
  Chooses parameters for NTT multiplication:
  k2, k3 = describe transform length
  num_primes = number of FFT primes
  r = we evaluate at 2^r
  lowmem = 1 if we are called by ntt_mpn_mul_lowmem
  assumes n1 >= n2
*/
void ntt_mpn_mul_params(unsigned* k2, unsigned* k3, unsigned* num_primes,
			unsigned* r, size_t n1, size_t n2, int lowmem)
{
  unsigned _k2, _k3, _num_primes, _r;
  size_t m1, m2, m, K;
  double cost, best = -1.0;

  for (_num_primes = 3; _num_primes <= MAX_NUM_PRIMES; _num_primes++)
    {
      // find largest possible r for given num_primes
      _r = (62 * _num_primes) / 2;

      while (1)
	{
	  // length of polynomials
	  m2 = n2 * 64 / _r + 1;
	  // u = maximum bitsize of coefficients of product
	  double u = 2 * _r + log((double) m2) / log(2.0);
	  if (u < global_bit_bound[_num_primes])
	    break;
	  _r--;
	}

      m1 = n1 * 64 / _r + 1;
      m = m1 + m2;

      // try various allowable transform sizes >= m
      for (_k3 = 0; _k3 <= MAX_K3; _k3++)
	{
	  _k2 = lowmem ? 6 : 0;
	  K = pow23(_k2, _k3);
	  while (K < m)
	    K *= 2, _k2++;

	  cost = 1.0 * K * _num_primes;

	  if (best == -1.0 || cost < best)
	    {
	      // best combination found so far
	      best = cost;
	      *num_primes = _num_primes;
	      *k2 = _k2;
	      *k3 = _k3;
	      *r = _r;
	    }
	}
    }
}



/*
  Computes 2^32/sqrt(K) mod p[i].
  This is the scale factor that must be applied at the split/reduce stage.
  The 1/sqrt(K) covers the FFT * IFFT = K scaling, and the 2^32 covers the
  factor of 1/2^64 introduced by pointwise_multiply.
*/
uint64_t scale_factor(unsigned i, unsigned k2, unsigned k3)
{
  uint64_t p = global_p[i];
  uint64_t pinv = global_pinv[i];
  uint64_t a = pow23(k2 / 2, k3 / 2);
  uint64_t z24 = calc_w(i, 24);
  uint64_t z12 = mul_mod_pinv(z24, z24, p, pinv);
  uint64_t z8 = mul_mod_pinv(z24, z12, p, pinv);
  uint64_t z4 = mul_mod_pinv(z8, z8, p, pinv);
  uint64_t z6 = mul_mod_pinv(z12, z12, p, pinv);
  uint64_t z3 = mul_mod_pinv(z6, z6, p, pinv);
  if (k2 & 1)
    {
      // sqrt(2) = zeta_8 * (1 - zeta_4)
      uint64_t c = mul_mod_pinv(z8, p + 1 - z4, p, pinv);
      a = mul_mod_pinv(a, c, p, pinv);
    }
  if (k3 & 1)
    {
      // sqrt(3) = zeta_4 * (2 * zeta_3 + 1)
      uint64_t c = mul_mod_pinv(2 * z3 + 1, z4, p, pinv);
      a = mul_mod_pinv(a, c, p, pinv);
    }
  // at this stage a = sqrt(K) mod p
  a = inv_mod(a, p);
  a = mul_mod_pinv(a, (uint64_t) 1 << 32, p, pinv);
  return a;
}



// ntt_mpn_mul_bonus without memory-saving techniques
// assumes n1 >= n2
void ntt_mpn_mul_plain(mp_limb_t* rop,
		       mp_limb_t* op1, size_t n1,
		       mp_limb_t* op2, size_t n2,
		       int num_threads)
{
  int squaring = (op1 == op2) && (n1 == n2);

  // choose transform parameters
  unsigned k2, k3, num_primes, r;
  ntt_mpn_mul_params(&k2, &k3, &num_primes, &r, n1, n2, 0);

  size_t K = pow23(k2, k3);
  size_t m1 = n1 * 64 / r + 1;    // polynomial lengths
  size_t m2 = n2 * 64 / r + 1;
  size_t m = m1 + m2;
  assert(m <= K);

  // compute scale factor for each prime
  uint64_t u[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    u[i] = scale_factor(i, k2, k3);

  // split/reduce op1 and op2
  uint64_t* F1[num_primes];
  uint64_t* F2[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    {
      F1[i] = safe_malloc(K * sizeof(uint64_t));
      F2[i] = F1[i] + m1;
    }
  split_reduce(F1, m1, op1, n1, r, u, num_primes, num_threads);
  if (!squaring)
    split_reduce(F2, m2, op2, n2, r, u, num_primes, num_threads);

  // convolution mod p for each prime
  uint64_t* temp;
  if (!squaring)
    temp = safe_malloc(K * sizeof(uint64_t));
  for (unsigned i = 0; i < num_primes; i++)
    {
      conv(F1[i], K, k2, k3,
	   squaring ? F1[i] : F2[i], m2, squaring ? F1[i] : temp,
	   F1[i], m1, F1[i], calc_w(i, K), global_p[i], num_threads);
    }
  if (!squaring)
    safe_free(temp);

  // CRT and recompose
  crt_recompose(rop, n1 + n2, r, F1, m, num_primes, num_threads);

  // clean up
  for (unsigned i = 0; i < num_primes; i++)
    safe_free(F1[i]);
}



// low-memory version of ntt_mpn_mul_bonus
// assumes n1 >= n2
void ntt_mpn_mul_lowmem(mp_limb_t* rop,
			mp_limb_t* op1, size_t n1, int preserve_op1,
			mp_limb_t* op2, size_t n2, int preserve_op2,
			int num_threads)
{
  prof_t prof;
  prof_init(&prof);
  prof_event(&prof, "begin");

  int squaring = (op1 == op2) && (n1 == n2);

  // choose transform parameters
  unsigned k2, k3, num_primes, r;
  ntt_mpn_mul_params(&k2, &k3, &num_primes, &r, n1, n2, 1);

  size_t K = pow23(k2, k3);
  size_t m1 = n1 * 64 / r + 1;    // polynomial lengths
  size_t m2 = n2 * 64 / r + 1;
  size_t m = m1 + m2;
  assert(m <= K);

  // compute scale factor for each prime
  uint64_t u[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    u[i] = scale_factor(i, k2, k3);

  // choose memory pool block size
  size_t block_size;
  if (K >= FFT_ARRAY_THRESHOLD)
    {
      // if we are going to use an array decomposition for the FFTs, ensure
      // the block size is a large enough multiple of the number of columns
      unsigned r2 = k2 / 2;
      unsigned r3 = k3 / 2;
      unsigned c2 = k2 - r2;
      unsigned c3 = k3 - r3;
      size_t R = pow23(r2, r3);
      size_t C = pow23(c2, c3);
      block_size = 8 * C;
    }
  else
    {
      block_size = K;
    }

  size_t num_blocks = K / block_size;

  pool_t pool;
  pool_init(&pool, block_size);

  // contribute memory from rop to the pool (except those bits contained in
  // op1 or op2)
  for (size_t i = 0; i + block_size <= n1 + n2; i += block_size)
    {
      if (!overlap(rop + i, block_size, op1, n1) &&
	  !overlap(rop + i, block_size, op2, n2))
	pool_insert(&pool, (uint64_t*) (rop + i));
    }

  prof_event(&prof, "split/reduce");

  // split/reduce op1 and op2
  vbuf_t F1[num_primes];
  vbuf_t F2[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    {
      vbuf_init(&F1[i], num_blocks, &pool);
      if (!squaring)
	vbuf_init(&F2[i], num_blocks, &pool);
    }
  split_reduce_vbuf(F1, m1, op1, n1, r, u, num_primes,
		    !preserve_op1, num_threads);
  if (!squaring)
    split_reduce_vbuf(F2, m2, op2, n2, r, u, num_primes,
		      !preserve_op2, num_threads);

  // convolution mod p for each prime
  for (unsigned i = 0; i < num_primes; i++)
    {
      prof_event(&prof, "convolution");

      conv_array_vbuf(&F1[i], K, k2, k3,
		      squaring ? &F1[i] : &F2[i], m2,
		      squaring ? &F1[i] : &F2[i],
		      &F1[i], m1, &F1[i],
		      calc_w(i, K), global_p[i], num_threads);
      if (!squaring)
	vbuf_clear(&F2[i]);
    }

  prof_event(&prof, "CRT/recompose");

  // CRT and recompose
  vbuf_t out;
  vbuf_init(&out, (n1 + n2) / block_size + 1, &pool);
  crt_recompose_vbuf(&out, n1 + n2, r, F1, m, 1, num_primes, num_threads);

  prof_event(&prof, "reassemble");

  for (unsigned i = 0; i < num_primes; i++)
    vbuf_clear(&F1[i]);

  // we're basically done, but some of the blocks of "out" might be contained
  // in rop, we need to move them out of the way
#pragma omp parallel for num_threads(num_threads) schedule(dynamic,1)
  for (ptrdiff_t i = 0; i < n1 + n2; i += block_size)
    {
      mp_limb_t* src = (mp_limb_t*) vbuf_get_block(&out, i / block_size);
      mp_limb_t* dst = src;
      while (overlap(dst, block_size, rop, n1 + n2))
	dst = (mp_limb_t*) pool_request(&pool);
      if (dst != src)
	memcpy(dst, src, block_size * sizeof(mp_limb_t));
      out.blocks[i / block_size] = (uint64_t*) dst;
    }

  // and finally copy the output into rop
#pragma omp parallel for num_threads(num_threads) schedule(dynamic,1)
  for (ptrdiff_t i = 0; i < n1 + n2; i += block_size)
    {
      mp_limb_t* src = (mp_limb_t*) vbuf_get_block(&out, i / block_size);
      size_t s = MIN(n1 + n2 - i, block_size);
      memcpy(rop + i, src, s * sizeof(mp_limb_t));
    }

  vbuf_clear(&out);
  pool_clear(&pool);

  prof_event(&prof, "finish");
  prof_report(&prof);
  prof_clear(&prof);
}



void ntt_mpn_mul_bonus(mp_limb_t* rop,
		       mp_limb_t* op1, size_t n1, int preserve_op1,
		       mp_limb_t* op2, size_t n2, int preserve_op2,
		       int optimise_memory, int num_threads)
{
  // swap inputs to make n1 >= n2
  if (n1 < n2)
    {
      { mp_limb_t* temp = op1; op1 = op2; op2 = temp; }
      { size_t temp = n1; n1 = n2; n2 = temp; }
      { int temp = preserve_op1; preserve_op1 = preserve_op2;
	preserve_op2 = temp; }
    }

  size_t threshold = tune_tab[MIN(num_threads, tune_tab_size - 1)];
  if ((n1 + n2) / 2 < threshold)
    {
      // multiplication is small; fall back on GMP
      ntt_mpn_mul_gmp(rop, op1, n1, op2, n2);
      return;
    }

  if (optimise_memory)
    ntt_mpn_mul_lowmem(rop, op1, n1, preserve_op1, op2, n2, preserve_op2,
		       num_threads);
  else
    ntt_mpn_mul_plain(rop, op1, n1, op2, n2, num_threads);
}
