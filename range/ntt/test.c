/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/


#include "test.h"


#define TEST_START printf("%s... ", __func__); fflush(stdout);
#define TEST_OK printf("ok\n"); fflush(stdout);


void print_array(uint64_t* x, size_t K)
{
  printf("[");
  if (K >= 1)
    printf("0x%016llx", x[0]);
  for (size_t i = 1; i < K; i++)
    printf(", 0x%016llx", x[i]);
  printf("]");
}


// global random state
gmp_randstate_t state;


#define NUM_SENTRIES 4
#define SENTRY 0x0123456789abcdefULL


// allocate block of size uint64_t's and write sentries before and after block
void* sentry_malloc(size_t size)
{
  uint64_t* ptr = safe_malloc((size + 2 * NUM_SENTRIES) * sizeof(uint64_t));
  for (size_t i = 0; i < NUM_SENTRIES; i++)
    ptr[i] = ptr[size + NUM_SENTRIES + i] = SENTRY;
  return ptr + NUM_SENTRIES;
}

// free block allocated via sentry_malloc
void sentry_free(void* ptr)
{
  safe_free((uint64_t*) ptr - NUM_SENTRIES);
}

// check that the sentries weren't modified
void sentry_check(void* ptr, size_t size)
{
  uint64_t* x = (uint64_t*) ptr - NUM_SENTRIES;
  for (size_t i = 0; i < NUM_SENTRIES; i++)
    assert(x[i] == SENTRY);
  for (size_t i = 0; i < NUM_SENTRIES; i++)
    assert(x[size + NUM_SENTRIES + i] == SENTRY);
}


uint64_t random64()
{
  uint64_t x0 = random32();
  uint64_t x1 = random32();
  return x0 | (x1 << 32);
}


// assumes 2^61 < p < 2^62
uint64_t random_mod(uint64_t p)
{
  uint64_t x;
  do x = random64() >> 2; while (x >= p);
  return x;
}


int probably_prime(uint64_t p)
{
  mpz_t y;
  mpz_init_set_ui(y, p);
  int result = mpz_probab_prime_p(y, 20);
  mpz_clear(y);
  return result;
}


void test_calc_pinvb()
{
  TEST_START;

  uint64_t mask = (int64_t) (-1);
  for (unsigned k = 1; k <= 64; k++, mask >>= 1)
    for (int trial = 0; trial < 100; trial++)
      {
	uint64_t p = (random64() & mask) | 1;
	uint64_t pinvb = calc_pinvb(p);
	assert(p * pinvb == 1);
      }

  TEST_OK;
}


void test_mul_mod_pinv()
{
  TEST_START;

  mpz_t z;
  mpz_init(z);

  for (int i = 0; i < 1000; i++)
    {
      uint64_t p = (random64() >> 3) | 0x2000000000000001ULL;
      uint64_t pinv = calc_pinv(p);
      for (int j = 0; j < 10000; j++)
	{
	  uint64_t x = random64();
	  uint64_t y = random_mod(p);
	  uint64_t r = mul_mod_pinv(x, y, p, pinv);
	  mpz_set_ui(z, x);
	  mpz_mul_ui(z, z, y);
	  mpz_tdiv_r_ui(z, z, p);
	  assert(mpz_get_ui(z) == r);
	}
    }

  mpz_clear(z);

  TEST_OK;
}


void test_calc_ypinv()
{
  TEST_START;

  mpz_t z;
  mpz_init(z);

  for (int i = 0; i < 1000; i++)
    {
      uint64_t p = (random64() >> 3) | 0x2000000000000001ULL;
      uint64_t pinv = calc_pinv(p);
      for (int j = 0; j < 10000; j++)
	{
	  uint64_t y = random_mod(p);
	  uint64_t ypinv = calc_ypinv(y, p, pinv);
	  mpz_set_ui(z, y);
	  mpz_mul_2exp(z, z, 64);
	  mpz_tdiv_q_ui(z, z, p);
	  assert(mpz_get_ui(z) == ypinv);
	}
    }

  mpz_clear(z);

  TEST_OK;
}


void test_mul_mod_ypinv()
{
  TEST_START;

  mpz_t z;
  mpz_init(z);

  for (int i = 0; i < 1000; i++)
    {
      uint64_t p = (random64() >> 3) | 0x2000000000000001ULL;
      uint64_t pinv = calc_pinv(p);
      for (int j = 0; j < 10000; j++)
	{
	  uint64_t x = random64();
	  uint64_t y = random_mod(p);
	  uint64_t ypinv = calc_ypinv(y, p, pinv);
	  uint64_t r = mul_mod_ypinv(x, y, ypinv, p);
	  mpz_set_ui(z, x);
	  mpz_mul_ui(z, z, y);
	  mpz_tdiv_r_ui(z, z, p);
	  assert(mpz_get_ui(z) == r);
	}
    }

  mpz_clear(z);

  TEST_OK;
}



/*
  reference implementation of fft() for testing
  inputs and outputs in [0, p)
  inplace operation only
*/
void fft_ref(uint64_t* op, size_t K, unsigned k2, unsigned k3, size_t n,
	     uint64_t w, uint64_t p)
{
  for (size_t i = n; i < K; i++)
    op[i] = 0;

  if (k3 >= 1)
    {
      // decompose into 3 rows, K/3 columns

      uint64_t u = pow_mod(w, K/3, p);    // cube root of 1
      uint64_t t = 1;                     // w^i

      // transform columns and apply twiddle factors
      for (size_t i = 0; i < K/3; i++)
	{
	  // radix-3 butterfly
	  uint64_t x0 = op[i];
	  uint64_t x1 = op[i + K/3];
	  uint64_t x2 = op[i + 2*K/3];
	  uint64_t x1u = mul_mod(x1, u, p);
	  uint64_t x1uu = mul_mod(x1u, u, p);
	  uint64_t x2u = mul_mod(x2, u, p);
	  uint64_t x2uu = mul_mod(x2u, u, p);
	  uint64_t y0 = add_mod(x0, add_mod(x1, x2, p), p);
	  uint64_t y1 = add_mod(x0, add_mod(x1u, x2uu, p), p);
	  uint64_t y2 = add_mod(x0, add_mod(x1uu, x2u, p), p);
	  op[i] = y0;
	  op[i + K/3] = mul_mod(y1, t, p);
	  op[i + 2*K/3] = mul_mod(y2, mul_mod(t, t, p), p);
	  t = mul_mod(t, w, p);
	}

      // transform rows
      uint64_t wcube = mul_mod(mul_mod(w, w, p), w, p);
      fft_ref(op, K/3, k2, k3 - 1, K/3, wcube, p);
      fft_ref(op + K/3, K/3, k2, k3 - 1, K/3, wcube, p);
      fft_ref(op + 2*K/3, K/3, k2, k3 - 1, K/3, wcube, p);
    }
  else if (k2 >= 1)
    {
      // decompose into 2 rows, K/2 columns

      // transform columns and apply twiddle factors
      uint64_t t = 1;                     // w^i
      for (size_t i = 0; i < K/2; i++)
	{
	  // radix-2 butterfly
	  uint64_t x0 = op[i];
	  uint64_t x1 = op[i + K/2];
	  op[i] = add_mod(x0, x1, p);
	  op[i + K/2] = mul_mod(sub_mod(x0, x1, p), t, p);
	  t = mul_mod(t, w, p);
	}

      // transform rows
      uint64_t wsqr = mul_mod(w, w, p);
      fft_ref(op, K/2, k2 - 1, k3, K/2, wsqr, p);
      fft_ref(op + K/2, K/2, k2 - 1, k3, K/2, wsqr, p);
    }
}

size_t fft_ref_pattern(size_t K, unsigned k2, unsigned k3)
{
  // all radix-3 layers first
  return ((size_t) 1 << k3) - 1;
}



/*
  Rearranges output of fft() from scrambled order into correct order.
  Input is op, output rop, no overlap allowed.
*/
void unscramble(uint64_t* rop, uint64_t* op, size_t K, unsigned k2, unsigned k3,
		size_t pattern)
{
  for (size_t j = 0; j < K; j++)
    {
      size_t _j = j, _pattern = pattern;
      size_t i = 0;
      for (unsigned t = 0; t < k2 + k3; t++, _pattern >>= 1)
	{
	  if (_pattern & 1)
	    i = 3*i + (_j % 3), _j /= 3;
	  else
	    i = 2*i + (_j % 2), _j /= 2;
	}
      rop[j] = op[i];
    }
}


void testcase_fft_base(size_t K, unsigned k2, unsigned k3, size_t n,
		       uint64_t w, uint64_t p)
{
  uint64_t* op = safe_malloc(K * sizeof(uint64_t));
  uint64_t* ref = safe_malloc(K * sizeof(uint64_t));
  uint64_t* temp = safe_malloc(K * sizeof(uint64_t));
  uint64_t* out = safe_malloc(K * sizeof(uint64_t));

  wtab_t wtab;
  wtab_init(&wtab, K, k2, k3, w, p);

  // random data to transform
  for (size_t i = 0; i < n; i++)
    op[i] = random_mod(p);

  // compute FFT using reference implementation
  for (size_t i = 0; i < n; i++)
    temp[i] = op[i];
  fft_ref(temp, K, k2, k3, n, w, p);
  unscramble(ref, temp, K, k2, k3, fft_ref_pattern(K, k2, k3));

  // compute FFT using fft_base
  for (size_t i = 0; i < n; i++)
    op[i] += (random32() & 1) ? p : 0;
  for (size_t i = n; i < K; i++)
    op[i] = random64();
  fft_base(op, K, k2, k3, n, &wtab, p);
  unscramble(out, op, K, k2, k3, fft_base_pattern(K, k2, k3));

  // compare results
  for (size_t i = 0; i < K; i++)
    {
      assert(out[i] < 2*p);
      assert(out[i] % p == ref[i]);
    }

  wtab_clear(&wtab);
  safe_free(op);
  safe_free(ref);
  safe_free(temp);
  safe_free(out);
}


void test_fft_base()
{
  TEST_START;

  for (int trial = 0; trial < 10; trial++)
    for (unsigned k2 = 0; k2 <= MAX_K2; k2++)
      for (unsigned k3 = 0; k3 <= MAX_K3; k3++)
	{
	  size_t K = pow23(k2, k3);
	  if (K < 5 * FFT_SPLIT_THRESHOLD)
	    {
	      unsigned i = random32() % MAX_NUM_PRIMES;
	      size_t n = random64() % K + 1;
	      uint64_t p = global_p[i];
	      size_t w = calc_w(i, K);
	      testcase_fft_base(K, k2, k3, n, w, p);
	    }
	}

  TEST_OK;
}


void testcase_ifft_base(size_t K, unsigned k2, unsigned k3,
			uint64_t w, uint64_t p)
{
  uint64_t* op = safe_malloc(K * sizeof(uint64_t));
  uint64_t* temp = safe_malloc(K * sizeof(uint64_t));

  uint64_t winv = inv_mod(w, p);

  wtab_t wtab1, wtab2;
  wtab_init(&wtab1, K, k2, k3, w, p);
  wtab_init(&wtab2, K, k2, k3, winv, p);

  // random data to transform
  for (size_t i = 0; i < K; i++)
    op[i] = random_mod(p);
  // expected result of fft_base followed by ifft_base
  for (size_t i = 0; i < K; i++)
    temp[i] = mul_mod(op[i], K, p);

  // compute FFT using fft_base
  fft_base(op, K, k2, k3, K, &wtab1, p);

  // then perform inverse FFT using ifft_base
  for (size_t i = 0; i < K; i++)
    op[i] = (op[i] % p) + (random32() & 3) * p;
  ifft_base(op, K, k2, k3, &wtab2, p);

  // check it really is the inverse
  for (size_t i = 0; i < K; i++)
    {
      assert(op[i] < 4*p);
      assert(op[i] % p == temp[i]);
    }

  wtab_clear(&wtab1);
  wtab_clear(&wtab2);
  safe_free(op);
  safe_free(temp);
}


void test_ifft_base()
{
  TEST_START;

  for (int trial = 0; trial < 10; trial++)
    for (unsigned k2 = 0; k2 <= MAX_K2; k2++)
      for (unsigned k3 = 0; k3 <= MAX_K3; k3++)
	{
	  size_t K = pow23(k2, k3);
	  if (K < 5 * FFT_SPLIT_THRESHOLD)
	    {
	      unsigned i = random32() % MAX_NUM_PRIMES;
	      size_t n = random64() % K + 1;
	      uint64_t p = global_p[i];
	      size_t w = calc_w(i, K);
	      testcase_ifft_base(K, k2, k3, w, p);
	    }
	}

  TEST_OK;
}


void testcase_fft_cols_rows(size_t R, unsigned r2, unsigned r3,
			    size_t C, unsigned c2, unsigned c3,
			    size_t n, uint64_t w, uint64_t p,
			    int inplace, int num_threads)
{
  unsigned k2 = r2 + c2;
  unsigned k3 = r3 + c3;
  size_t K = pow23(k2, k3);

  uint64_t* op = safe_malloc(K * sizeof(uint64_t));
  uint64_t* temp = safe_malloc(K * sizeof(uint64_t));
  uint64_t* ref = safe_malloc(K * sizeof(uint64_t));
  uint64_t* out = safe_malloc(K * sizeof(uint64_t));
  uint64_t* dest = inplace ? op : temp;

  // random data to transform
  for (size_t i = 0; i < n; i++)
    op[i] = random_mod(p);

  // compute FFT using reference implementation
  for (size_t i = 0; i < n; i++)
    temp[i] = op[i];
  fft_ref(temp, K, k2, k3, n, w, p);
  unscramble(ref, temp, K, k2, k3, fft_ref_pattern(K, k2, k3));

  // compute FFT using fft_cols and fft_rows
  for (size_t i = 0; i < n; i++)
    op[i] += (random32() & 1) ? p : 0;
  for (size_t i = n; i < K; i++)
    op[i] = random64();
  fft_cols(dest, R, r2, r3, C, c2, c3, op, n, w, p, num_threads);
  fft_rows(dest, R, r2, r3, C, c2, c3, w, p, num_threads);
  unscramble(out, dest, K, k2, k3, fft_cols_rows_pattern(R, r2, r3, C, c2, c3));

  // compare results
  for (size_t i = 0; i < K; i++)
    {
      assert(out[i] < 2*p);
      assert(out[i] % p == ref[i]);
    }

  safe_free(op);
  safe_free(ref);
  safe_free(temp);
  safe_free(out);
}


void test_fft_cols_rows()
{
  TEST_START;

  for (int trial = 0; trial < 10; trial++)
    for (unsigned c2 = 0; c2 < 5; c2++)
      for (unsigned c3 = 0; c3 < 3; c3++)
	for (unsigned r2 = 0; r2 < 5; r2++)
	  for (unsigned r3 = 0; r3 < 3; r3++)
	    {
	      unsigned i = random32() % MAX_NUM_PRIMES;
	      unsigned k2 = c2 + r2;
	      unsigned k3 = c3 + r3;
	      size_t C = pow23(c2, c3);
	      size_t R = pow23(r2, r3);
	      size_t K = pow23(k2, k3);
	      size_t n = random64() % K + 1;
	      uint64_t p = global_p[i];
	      size_t w = calc_w(i, K);
	      int inplace = random32() & 1;
	      int num_threads = random32() % 4 + 1;
	      testcase_fft_cols_rows(R, r2, r3, C, c2, c3, n, w, p,
				     inplace, num_threads);
	    }

  TEST_OK;
}


void testcase_fft_cols_rows_vbuf(size_t R, unsigned r2, unsigned r3,
				 size_t C, unsigned c2, unsigned c3,
				 size_t block_size, size_t n,
				 uint64_t w, uint64_t p,
				 int inplace, int num_threads)
{
  assert(block_size % C == 0);
  unsigned k2 = r2 + c2;
  unsigned k3 = r3 + c3;
  size_t K = pow23(k2, k3);

  pool_t pool;
  pool_init(&pool, block_size);

  uint64_t* op = safe_malloc(K * sizeof(uint64_t));
  uint64_t* ref = safe_malloc(K * sizeof(uint64_t));
  uint64_t* temp = safe_malloc(K * sizeof(uint64_t));
  uint64_t* out = safe_malloc(K * sizeof(uint64_t));
  vbuf_t op_vbuf, rop_vbuf;
  vbuf_init(&op_vbuf, K / block_size + 1, &pool);
  vbuf_init(&rop_vbuf, K / block_size + 1, &pool);
  vbuf_t* dest_vbuf = inplace ? &op_vbuf : &rop_vbuf;

  // random data to transform
  for (size_t i = 0; i < n; i++)
    op[i] = random_mod(p);

  // compute FFT using reference implementation
  memcpy(temp, op, n * sizeof(uint64_t));
  fft_ref(temp, K, k2, k3, n, w, p);
  unscramble(ref, temp, K, k2, k3, fft_ref_pattern(K, k2, k3));

  // compute FFT using fft_cols_vbuf and fft_rows_vbuf
  for (size_t i = 0; i < n; i++)
    op[i] += (random32() & 1) ? p : 0;
  for (size_t i = n; i < K; i++)
    op[i] = random64();
  for (size_t i = 0; i < K; i += block_size)
    memcpy(vbuf_get(&op_vbuf, i), op + i,
	   MIN(block_size, K - i) * sizeof(uint64_t));
  fft_cols_vbuf(dest_vbuf, R, r2, r3, C, c2, c3,
		&op_vbuf, n, w, p, num_threads);
  fft_rows_vbuf(dest_vbuf, R, r2, r3, C, c2, c3, w, p, num_threads);
  for (size_t i = 0; i < K; i += block_size)
    memcpy(temp + i, vbuf_get(dest_vbuf, i),
	   MIN(block_size, K - i) * sizeof(uint64_t));
  unscramble(out, temp, K, k2, k3,
	     fft_cols_rows_vbuf_pattern(R, r2, r3, C, c2, c3, block_size));

  // compare results
  for (size_t i = 0; i < K; i++)
    {
      assert(out[i] < 2*p);
      assert(out[i] % p == ref[i]);
    }

  vbuf_clear(&op_vbuf);
  vbuf_clear(&rop_vbuf);
  pool_clear(&pool);
  safe_free(op);
  safe_free(ref);
  safe_free(temp);
  safe_free(out);
}


void test_fft_cols_rows_vbuf()
{
  TEST_START;

  for (unsigned c2 = 0; c2 < 5; c2++)
    for (unsigned c3 = 0; c3 < 3; c3++)
	{
	  size_t C = pow23(c2, c3);
	  for (size_t block_size = C; block_size <= 10*C; block_size += C)
	    for (unsigned r2 = 0; r2 < 5; r2++)
	      for (unsigned r3 = 0; r3 < 3; r3++)
		{
		  unsigned i = random32() % MAX_NUM_PRIMES;
		  unsigned k2 = c2 + r2;
		  unsigned k3 = c3 + r3;
		  size_t R = pow23(r2, r3);
		  size_t K = pow23(k2, k3);
		  size_t n = random64() % K + 1;
		  uint64_t p = global_p[i];
		  size_t w = calc_w(i, K);
		  int inplace = random32() & 1;
		  int num_threads = random32() % 4 + 1;
		  testcase_fft_cols_rows_vbuf(R, r2, r3, C, c2, c3, block_size,
					      n, w, p, inplace, num_threads);
		}
	}

  TEST_OK;
}


void testcase_ifft_cols_rows_vbuf(size_t R, unsigned r2, unsigned r3,
				  size_t C, unsigned c2, unsigned c3,
				  size_t block_size, uint64_t w, uint64_t p,
				  int num_threads)
{
  assert(block_size % C == 0);
  unsigned k2 = r2 + c2;
  unsigned k3 = r3 + c3;
  size_t K = pow23(k2, k3);
  uint64_t winv = inv_mod(w, p);

  pool_t pool;
  pool_init(&pool, block_size);

  uint64_t* op = safe_malloc(K * sizeof(uint64_t));
  uint64_t* temp = safe_malloc(K * sizeof(uint64_t));
  vbuf_t vbuf;
  vbuf_init(&vbuf, K / block_size + 1, &pool);

  // random data to transform
  for (size_t i = 0; i < K; i++)
    op[i] = random_mod(p);

  // compute FFT using fft_cols_vbuf and fft_rows_vbuf
  for (size_t i = 0; i < K; i += block_size)
    memcpy(vbuf_get(&vbuf, i), op + i,
	   MIN(block_size, K - i) * sizeof(uint64_t));
  fft_cols_vbuf(&vbuf, R, r2, r3, C, c2, c3, &vbuf, K, w, p, num_threads);
  fft_rows_vbuf(&vbuf, R, r2, r3, C, c2, c3, w, p, num_threads);

  // compute IFFT using ifft_cols_vbuf and ifft_rows_vbuf
  ifft_rows_vbuf(&vbuf, R, r2, r3, C, c2, c3, winv, p, num_threads);
  ifft_cols_vbuf(&vbuf, R, r2, r3, C, c2, c3, winv, p, num_threads);
  for (size_t i = 0; i < K; i += block_size)
    memcpy(temp + i, vbuf_get(&vbuf, i),
	   MIN(block_size, K - i) * sizeof(uint64_t));

  // compare results
  for (size_t i = 0; i < K; i++)
    {
      assert(temp[i] < 4*p);
      assert(temp[i] % p == mul_mod(op[i], K, p));
    }

  vbuf_clear(&vbuf);
  pool_clear(&pool);
  safe_free(op);
  safe_free(temp);
}



void test_ifft_cols_rows_vbuf()
{
  TEST_START;

  for (unsigned c2 = 0; c2 < 5; c2++)
    for (unsigned c3 = 0; c3 < 3; c3++)
	{
	  size_t C = pow23(c2, c3);
	  for (size_t block_size = C; block_size <= 10*C; block_size += C)
	    for (unsigned r2 = 0; r2 < 5; r2++)
	      for (unsigned r3 = 0; r3 < 3; r3++)
		{
		  unsigned i = random32() % MAX_NUM_PRIMES;
		  unsigned k2 = c2 + r2;
		  unsigned k3 = c3 + r3;
		  size_t R = pow23(r2, r3);
		  size_t K = pow23(k2, k3);
		  uint64_t p = global_p[i];
		  size_t w = calc_w(i, K);
		  int num_threads = random32() % 4 + 1;
		  testcase_ifft_cols_rows_vbuf(R, r2, r3, C, c2, c3, block_size,
					       w, p, num_threads);
		}
	}

  TEST_OK;
}



void testcase_split(size_t n, size_t m, unsigned r)
{
  unsigned t = (r + 63) / 64;    // number of limbs per output coefficient

  mp_limb_t* x = safe_malloc((n + 1) * sizeof(mp_limb_t));
  mp_limb_t* y = sentry_malloc(t * m);

  mpz_t xx, yy, z1, z2;
  mpz_init(xx);
  mpz_init(yy);
  mpz_init(z1);
  mpz_init(z2);

  for (size_t i = 0; i < n; i++)
    x[i] = random64();

  if (n > 0)
    mpz_import(xx, n, -1, sizeof(mp_limb_t), 0, 0, x);

  split(y, m, x, n, r);
  sentry_check(y, t * m);

  mpz_import(yy, t * m, -1, sizeof(mp_limb_t), 0, 0, y);

  for (size_t j = 0; j < m; j++)
    {
      mpz_tdiv_r_2exp(z1, xx, r);
      mpz_tdiv_q_2exp(xx, xx, r);
      mpz_tdiv_r_2exp(z2, yy, 64 * t);
      mpz_tdiv_q_2exp(yy, yy, 64 * t);
      assert(mpz_cmp(z1, z2) == 0);
    }

  mpz_clear(z1);
  mpz_clear(z2);
  mpz_clear(xx);
  mpz_clear(yy);
  safe_free(x);
  sentry_free(y);
}


void test_split()
{
  TEST_START;

  for (int trial = 0; trial < 100000; trial++)
    {
      size_t n = random64() % 200;
      size_t m = random64() % 200 + 1;
      size_t r = random64() % 192 + 1;
      testcase_split(n, m, r);
    }

  TEST_OK;
}


void testcase_reduce(size_t m, unsigned t, uint64_t p)
{
  mp_limb_t* x = safe_malloc(t * m * sizeof(mp_limb_t));
  uint64_t* y = safe_malloc(m * sizeof(mp_limb_t));

  for (size_t i = 0; i < t * m; i++)
    x[i] = random64();

  uint64_t u = random_mod(p);
  reduce(y, x, m, t, u, p, calc_pinv(p), calc_pinvb(p));

  for (size_t i = 0; i < m; i++)
    {
      uint64_t r = mul_mod(mpn_mod_1(x + i*t, t, p), u, p);
      assert(y[i] < 2*p);
      assert(y[i] % p == r);
    }

  safe_free(x);
  safe_free(y);
}


void test_reduce()
{
  TEST_START;

  for (int trial = 0; trial < 100000; trial++)
    for (unsigned m = 1; m < 5; m++)
      for (unsigned t = 1; t <= 3; t++)
	{
	  unsigned i = random32() % MAX_NUM_PRIMES;
	  uint64_t p = global_p[i];
	  testcase_reduce(m, t, p);
	}

  TEST_OK;
}



void testcase_split_reduce_vbuf(size_t m, size_t n, unsigned r,
				unsigned num_primes, size_t block_size,
				int num_threads)
{
  // random input
  mp_limb_t* op = safe_malloc((n + 1) * sizeof(mp_limb_t));
  for (size_t i = 0; i < n; i++)
    op[i] = random64();

  // split and reduce using basic functions
  unsigned t = (r + 63) / 64;
  mp_limb_t* temp = safe_malloc(m * t * sizeof(mp_limb_t));
  split(temp, m, op, n, r);
  uint64_t* ref[num_primes];
  uint64_t u[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    {
      u[i] = random_mod(global_p[i]);
      ref[i] = safe_malloc(m * sizeof(uint64_t));
      reduce(ref[i], temp, m, t, u[i],
	     global_p[i], global_pinv[i], global_pinvb[i]);
    }
  safe_free(temp);

  // now try it using split_reduce_vbuf
  pool_t pool;
  vbuf_t out[num_primes];

  pool_init(&pool, block_size);
  for (unsigned i = 0; i < num_primes; i++)
    vbuf_init(&out[i], m / block_size + 1, &pool);

  split_reduce_vbuf(out, m, op, n, r, u, num_primes, 1, num_threads);

  // compare results
  for (unsigned i = 0; i < num_primes; i++)
    for (size_t j = 0; j < m; j += block_size)
      {
	uint64_t* ptr = vbuf_get_block(&out[i], j / block_size);
	for (size_t h = 0; h < MIN(m - j, block_size); h++)
	  assert(ptr[h] == ref[i][j + h]);
      }

  // clean up
  for (unsigned i = 0; i < num_primes; i++)
    {
      vbuf_clear(&out[i]);
      safe_free(ref[i]);
    }
  pool_clear(&pool);
  safe_free(op);
}



void test_split_reduce_vbuf()
{
  TEST_START;

  for (int trial = 0; trial < 100000; trial++)
    {
      size_t block_size = 64 * (random32() % 10 + 1);
      unsigned num_primes = random32() % MAX_NUM_PRIMES + 1;
      unsigned r = random32() % 192 + 1;
      size_t m = random64() % (10 * block_size) + 1;
      size_t n = random64() % (10 * block_size);
      int num_threads = random32() % 4 + 1;
      testcase_split_reduce_vbuf(m, n, r, num_primes, block_size, num_threads);
    }

  TEST_OK;
}


void testcase_crt(size_t m, unsigned num_primes)
{
  uint64_t* op[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    op[i] = safe_malloc(m * sizeof(uint64_t));

  for (unsigned i = 0; i < num_primes; i++)
    for (size_t j = 0; j < m; j++)
      op[i][j] = random_mod(global_p[i]) + (random32() % 3) * global_p[i];

  mp_limb_t* rop = safe_malloc(m * num_primes * sizeof(mp_limb_t));

  crt(rop, op, m, num_primes);

  // q = product of p's
  mp_limb_t q[num_primes];
  q[0] = 1;
  for (unsigned i = 1; i < num_primes; i++)
    q[i] = 0;
  for (unsigned i = 0; i < num_primes; i++)
    mpn_mul_1(q, q, num_primes, global_p[i]);

  // check output is correct
  for (size_t j = 0; j < m; j++)
    {
      for (unsigned i = 0; i < num_primes; i++)
	{
	  assert(mpn_cmp(rop + j * num_primes, q, num_primes) < 0);
	  uint64_t x = mpn_mod_1(rop + j * num_primes, num_primes, global_p[i]);
	  uint64_t y = op[i][j] % global_p[i];
	  assert(x == y);
	}
    }

  safe_free(rop);
  for (unsigned i = 0; i < num_primes; i++)
    safe_free(op[i]);
}


void test_crt()
{
  TEST_START;

  for (int trial = 0; trial < 1000000; trial++)
    {
      size_t m = random32() % 10 + 1;
      unsigned num_primes = random32() % MAX_NUM_PRIMES + 1;
      testcase_crt(m, num_primes);
    }

  TEST_OK;
}


void testcase_recompose(size_t n, size_t m, unsigned s, unsigned r, unsigned r0)
{
  mp_limb_t* op = safe_malloc((m + 1) * s * sizeof(mp_limb_t));
  mpn_random2(op, (m + 1) * s);

  mpz_t x, y;
  mpz_init(x);
  mpz_init(y);
  for (size_t i = 0; i < m; i++)
    {
      mpz_import(y, s, -1, sizeof(mp_limb_t), 0, 0, op + i * s);
      mpz_mul_2exp(y, y, r0 + i * r);
      mpz_add(x, x, y);
    }

  mp_limb_t* ref = safe_malloc(n * sizeof(mp_limb_t));
  for (size_t i = 0; i < n; i++)
    {
      ref[i] = mpz_get_ui(x);
      mpz_tdiv_q_2exp(x, x, 64);
    }

  mp_limb_t* rop = sentry_malloc(n);

  recompose(rop, n, op, m, s, r, r0);
  sentry_check(rop, n);

  assert(mpn_cmp(rop, ref, n) == 0);

  mpz_clear(x);
  mpz_clear(y);
  sentry_free(rop);
  safe_free(op);
  safe_free(ref);
}


void test_recompose()
{
  TEST_START;

  for (int trial = 0; trial < 100000; trial++)
    {
      size_t n = random32() % 200 + 1;
      size_t m = random32() % 200;
      size_t s = random32() % 6 + 1;
      unsigned r = random32() % (64 * s) + 1;
      unsigned r0 = random32() % 300;
      testcase_recompose(n, m, s, r, r0);
    }

  TEST_OK;
}



void testcase_crt_recompose_vbuf(size_t n, size_t m, unsigned num_primes,
				 unsigned r, size_t block_size, int num_threads)
{
  uint64_t* op[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    op[i] = safe_malloc((m + 1) * sizeof(uint64_t));
  mp_limb_t* ref = safe_malloc(n * sizeof(mp_limb_t));
  mp_limb_t* temp = safe_malloc((m + 1) * num_primes * sizeof(mp_limb_t));

  mpz_t x;
  mpz_init(x);
  
  for (size_t j = 0; j < m; j++)
    {
      // (next line assumes FFT primes are "close" to 2^62)
      mpz_rrandomb(x, state, 62 * num_primes - 1);
      for (unsigned i = 0; i < num_primes; i++)
	op[i][j] = mpz_fdiv_ui(x, global_p[i]) + (random32() & 3) * global_p[i];
      for (unsigned i = 0; i < num_primes; i++)
	{
	  temp[j * num_primes + i] = mpz_get_ui(x);
	  mpz_tdiv_q_2exp(x, x, 64);
	}
    }

  mpz_clear(x);

  // compute recomposition using simple function
  recompose(ref, n, temp, m, num_primes, r, 0);
  safe_free(temp);

  // now try using crt_recompose_vbuf
  pool_t pool;
  pool_init(&pool, block_size);

  vbuf_t op_vbuf[num_primes];
  for (unsigned i = 0; i < num_primes; i++)
    {
      vbuf_init(&op_vbuf[i], m / block_size + 1, &pool);
      for (size_t j = 0; j < m; j++)
	vbuf_get(&op_vbuf[i], j)[0] = op[i][j];
    }

  vbuf_t rop_vbuf;
  vbuf_init(&rop_vbuf, n / block_size + 1, &pool);

  crt_recompose_vbuf(&rop_vbuf, n, r, op_vbuf, m, 1, num_primes, num_threads);
  
  // check results
  for (size_t j = 0; j < n; j++)
    assert(vbuf_get(&rop_vbuf, j)[0] == ref[j]);

  vbuf_clear(&rop_vbuf);
  for (unsigned i = 0; i < num_primes; i++)
    vbuf_clear(&op_vbuf[i]);
  pool_clear(&pool);
  for (unsigned i = 0; i < num_primes; i++)
    safe_free(op[i]);
  safe_free(ref);
}


void test_crt_recompose_vbuf()
{
  TEST_START;

  for (int trial = 0; trial < 10000; trial++)
    {
      size_t m = random64() % (1 << (random32() % 15));
      size_t n = random64() % (1 << (random32() % 15)) + 1;
      unsigned r = random32() % 192 + 1;
      unsigned num_primes;
      do num_primes = random32() % MAX_NUM_PRIMES + 1;
      while (r / 64 > num_primes);
      int num_threads = random32() % 4 + 1;
      size_t block_size = 64 * (random32() % 10 + 1);
      testcase_crt_recompose_vbuf(n, m, num_primes, r, block_size, num_threads);
    }

  TEST_OK;
}



void testcase_ntt_mpn_mul_bonus(size_t n1, int preserve_op1,
				size_t n2, int preserve_op2,
				int squaring,
				int optimise_memory, int num_threads)
{
  if (squaring)
    n2 = n1;

  // put op1, op2, rop in random locations (possibly overlapping) within
  // a buffer of length k
  size_t k = 3 * (n1 + n2);
  mp_limb_t* buf = safe_malloc(k * sizeof(mp_limb_t));
  for (size_t i = 0; i < k; i++)
    buf[i] = SENTRY;

  mp_limb_t* op1 = buf + random64() % (k - n1 + 1);
  mp_limb_t* op2 = squaring ? op1 : (buf + random64() % (k - n2 + 1));
  mp_limb_t* rop = buf + random64() % (k - (n1 + n2) + 1);

  if (random32() & 1)
    mpn_random2(op1, n1);
  else
    mpn_random(op1, n1);

  if (!squaring)
    {
      if (random32() & 1)
	mpn_random2(op2, n2);
      else
	mpn_random(op2, n2);
    }

  mp_limb_t* buf_copy = safe_malloc(k * sizeof(mp_limb_t));
  memcpy(buf_copy, buf, k * sizeof(mp_limb_t));

  mp_limb_t* ref = safe_malloc((n1 + n2) * sizeof(mp_limb_t));

  if (n1 >= n2)
    mpn_mul(ref, op1, n1, op2, n2);
  else
    mpn_mul(ref, op2, n2, op1, n1);

  ntt_mpn_mul_bonus(rop, op1, n1, preserve_op1, op2, n2, preserve_op2,
		    optimise_memory, num_threads);

  // check memory hasn't been modified where it shouldn't have been
  for (size_t i = 0; i < k; i++)
    {
      // ignore changes in output region
      if (buf + i >= rop && buf + i < rop + n1 + n2)
	continue;

      // check operands haven't been modified if appropriate preserve
      // flags are set, and check regions outside operands haven't been modified
      int in_op1 = buf + i >= op1 && buf + i < op1 + n1;
      int in_op2 = buf + i >= op2 && buf + i < op2 + n2;
      if ((in_op1 && preserve_op1) || (in_op2 && preserve_op2) ||
	  (!in_op1 && !in_op2))
	assert(buf[i] == buf_copy[i]);
    }

  // check result is correct
  assert(mpn_cmp(rop, ref, n1 + n2) == 0);

  safe_free(buf_copy);
  safe_free(buf);
  safe_free(ref);
}



void test_ntt_mpn_mul_bonus()
{
  TEST_START;

  for (int trial = 0; trial < 10000; trial++)
    {
      size_t n1 = random64() % (1 << (random32() % 15)) + 1;
      size_t n2 = random64() % (1 << (random32() % 15)) + 1;
      int preserve_op1 = random32() & 1;
      int preserve_op2 = random32() & 1;
      int optimise_memory = random32() & 1;
      int squaring = random32() & 1;
      int num_threads = random32() % 4 + 1;
      testcase_ntt_mpn_mul_bonus(n1, preserve_op1, n2, preserve_op2,
				 squaring, optimise_memory, num_threads);
    }

  TEST_OK;
}



#if 1

int main()
{
  gmp_randinit_default(state);
  ntt_init();

#if 1
  test_calc_pinvb();
  test_mul_mod_pinv();
  test_calc_ypinv();
  test_mul_mod_ypinv();
  test_fft_base();
  test_ifft_base();
  test_fft_cols_rows();
  test_fft_cols_rows_vbuf();
  test_ifft_cols_rows_vbuf();
  test_split();
  test_reduce();
  test_split_reduce_vbuf();
  test_crt();
  test_recompose();
  test_crt_recompose_vbuf();
  test_ntt_mpn_mul_bonus();
#endif

  return 0;
}


#else


void usage()
{
  printf("syntax: test [ gmp | ntt ] size num_threads\n");
  exit(0);
}


int main(int argc, char* argv[])
{
  gmp_randinit_default(state);
  ntt_init();

  if (argc != 4)
    usage();

  size_t n = atol(argv[2]);
  int num_threads = atoi(argv[3]);

  int which;
  if (strcmp(argv[1], "gmp") == 0)
    which = 0;
  else if (strcmp(argv[1], "ntt") == 0)
    which = 1;
  else
    usage();

  mp_limb_t* op1, * op2, * rop;

  if (which == 0)
    {
      op1 = safe_malloc(n * sizeof(mp_limb_t));
      op2 = safe_malloc(n * sizeof(mp_limb_t));
      rop = safe_malloc(2 * n * sizeof(mp_limb_t));
    }
  else
    {
      rop = safe_malloc(2 * n * sizeof(mp_limb_t));
      op1 = rop;
      op2 = op1 + n;
    }

  for (size_t i = 0; i < n; i++)
    op1[i] = random64();
  for (size_t i = 0; i < n; i++)
    op2[i] = random64();

  // also do multiplication mod q for consistency checking
  mp_limb_t q = 0xffffffffffffffc5ULL;
  mp_limb_t res1 = mpn_mod_1(op1, n, q);
  mp_limb_t res2 = mpn_mod_1(op2, n, q);

  printf("start...\n");
  fflush(stdout);
  double time1 = get_time();
  if (which == 0)
    mpn_mul(rop, op1, n, op2, n);
  else
    ntt_mpn_mul_bonus(rop, op1, n, 0, op2, n, 0, 1, num_threads);
  double time2 = get_time();
  printf("done\n");

  mp_limb_t res3 = mpn_mod_1(rop, 2*n, q);
#if AVOID_128_BIT
  {
    mp_limb_t z0, z1, dummy, r;
    MUL128(z1, z0, res1, res2);
    DIV128(dummy, r, z1, z0, q);
    assert(r == res3);
  }
#else
  assert((uint128_t) res1 * res2 % q == res3);
#endif

  printf("%llu * %llu = %llu (mod %llu)\n",
	 (uint64_t) res1, (uint64_t) res2, (uint64_t) res3, (uint64_t) q);
  printf("time = %.1lf\n", time2 - time1);

  if (which == 0)
    {
      free(op1);
      free(op2);
      free(rop);
    }
  else
    free(rop);
      
  return 0;
}

#endif
