/*
  Copyright (C) 2012, David Harvey

  This file is part of the ntt library (version 0.1.2).

  ntt is released under a modified BSD license. See the file COPYING in the
  source distribution for details.
*/

#include <omp.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "ntt.h"


// set this to avoid gcc 128-bit types, use inline assembly instead
#ifndef AVOID_128_BIT
#define AVOID_128_BIT 1
#endif


// enable profiling output
#define PROFILE 0


// estimated L1 cache size, in bytes
#define CACHE_SIZE 32768

// estimated cache line size, in bytes
#define CACHE_LINE 64

// transforms longer than this use the "array" algorithm
#define FFT_ARRAY_THRESHOLD 65536

// fft_base splits up transforms larger than this
#define FFT_SPLIT_THRESHOLD 1024


#if !(GMP_LIMB_BITS == 64 && GMP_NUMB_BITS == 64)
#error NTT requires GMP compiled with 64-bit limbs and no nails
#endif

#ifndef MAX
#define MAX(aaa, bbb) (((aaa) >= (bbb)) ? (aaa) : (bbb))
#endif

#ifndef MIN
#define MIN(aaa, bbb) (((aaa) <= (bbb)) ? (aaa) : (bbb))
#endif


#define INLINE static inline


#if AVOID_128_BIT

// gcc x86_64 inline assembly
#define MUL128(z1, z0, x, y) __asm__("mulq %3" : "=a" (z0), "=d" (z1) : "%0" ((uint64_t)(x)), "rm" ((uint64_t)(y)))
#define DIV128(q, r, x1, x0, y) __asm__("divq %4" : "=a" (q), "=d" (r) : "0" ((uint64_t)(x0)), "1" ((uint64_t)(x1)), "rm" ((uint64_t)(y)))

INLINE uint64_t mulhi(uint64_t x, uint64_t y)
{
  uint64_t z1, z0;
  MUL128(z1, z0, x, y);
  return z1;
}

#else

// signed and unsigned 128-bit types
typedef int int128_t __attribute__((mode(TI)));
typedef unsigned int uint128_t __attribute__((mode(TI)));

INLINE uint64_t mulhi(uint64_t x, uint64_t y)
{
  return ((uint128_t) x * y) >> 64;
}

#endif



#define fatal_error _NTT_fatal_error
void fatal_error(char* fmt, ...);


#ifdef _OPENMP
#define INIT_LOCK(x) omp_init_lock(x)
#define DESTROY_LOCK(x) omp_destroy_lock(x)
#define SET_LOCK(x) omp_set_lock(x)
#define UNSET_LOCK(x) omp_unset_lock(x)
#else
#define INIT_LOCK(x) { }
#define DESTROY_LOCK(x) { }
#define SET_LOCK(x) { }
#define UNSET_LOCK(x) { }
#endif


// *****************************************************************************
//  profiling

// current wall time in seconds
#define get_time _NTT_get_time
double get_time();


#define MAX_PROF 100

#if PROFILE

typedef struct
{
  size_t count;
  char* labels[MAX_PROF];
  double times[MAX_PROF];
} prof_t;


#define prof_init _NTT_prof_init
void prof_init(prof_t* prof);

#define prof_clear _NTT_prof_clear
void prof_clear(prof_t* prof);

#define prof_event _NTT_prof_event
void prof_event(prof_t* prof, char* label);

#define prof_report _NTT_prof_report
void prof_report(prof_t* prof);

#else

// dummy implementation

typedef int prof_t;

#define prof_init(X)
#define prof_clear(X)
#define prof_event(X, Y)
#define prof_report(X)

#endif

// *****************************************************************************
//  memory management


// wrapper for malloc
// on failure, print error and abort
#define safe_malloc _NTT_safe_malloc
void* safe_malloc(size_t n);

#define safe_realloc _NTT_safe_realloc
void* safe_realloc(void* ptr, size_t n);

#define safe_free _NTT_safe_free
void safe_free(void* ptr);


// zeroes out n words starting at ptr
// n == 0 is allowed
#define zero _NTT_zero
INLINE void zero(void* ptr, size_t n)
{
  uint64_t* x = ptr;
  for (ptrdiff_t k = n - 1; k >= 0; k--)
    x[k] = 0;
}


// rounds up ptr to the next cache line
#define next_cache_line _NTT_next_cache_line
INLINE void* next_cache_line(void* ptr)
{
  return (void*) ((((uintptr_t) ptr) + CACHE_LINE - 1)
		  & ~((uintptr_t) CACHE_LINE - 1));
}


// a node in a singly-linked list of uint64_t*.
struct list_node_struct
{
  struct list_node_struct* next;    // next node in list
  uint64_t* ptr;                    // actual data
};

typedef struct list_node_struct list_node_t;


// adds node pointing to "ptr" at beginning of list
// ptr cannot be NULL
#define list_insert _NTT_list_insert
void list_insert(list_node_t** list, uint64_t* ptr);

// extracts and removes node from beginning of list
// returns NULL if list is empty
#define list_extract _NTT_list_extract
uint64_t* list_extract(list_node_t** list);



/*
  A pool_t manages a collection of blocks of memory of a fixed size.
*/
typedef struct
{
#ifdef _OPENMP
  omp_lock_t lock;      // for access to "available" and "from_heap"
#endif

  size_t block_size;    // block size in words

  // list of available blocks
  list_node_t* available;

  // list of blocks that must be freed by pool_clear
  list_node_t* from_heap;

} pool_t;


// Initialises "pool" with given block size.
#define pool_init _NTT_pool_init
void pool_init(pool_t* pool, size_t block_size);


/*
  Destroys the pool.
  All blocks that were allocated via this pool must be returned to the pool
  before calling pool_clear.
*/
#define pool_clear _NTT_pool_clear
void pool_clear(pool_t* pool);


/*
  Adds the given block to the pool. This routine is threadsafe.

  This could be a block previously obtained via pool_request(), in which case
  ownership returns to the pool. Or it could be a block owned by the caller,
  in which case the caller must ensure it is available during the lifetime of
  the pool, and ownership returns to the caller after the pool is destroyed.
*/
#define pool_insert _NTT_pool_insert
void pool_insert(pool_t* pool, uint64_t* ptr);


// request a block from the pool. This routine is threadsafe.
#define pool_request _NTT_pool_request
uint64_t* pool_request(pool_t* pool);


/*
  A vbuf_t represents a "virtual buffer". It is a collection of blocks of memory
  interpreted as a single contiguous array. Blocks are allocated lazily, from
  a memory pool, when they are needed.
*/
typedef struct
{
  // must match block size of memory pool
  size_t block_size;

  // if wrap != NULL, then this is not really a vbuf, but rather a vbuf
  // wrapping an ordinary block of memory located at "wrap". The caller is
  // responsible for allocating and freeing "wrap". The remaining data
  // members (except block_size) are unused. Such a vbuf is created with
  // vbuf_init_wrap.
  uint64_t* wrap;

#ifdef _OPENMP
  omp_lock_t lock;   // for access to "blocks" and "regions"
#endif

  // associated memory pool
  pool_t* pool;

  // number of blocks; fixed at initialisation
  size_t n;

  // block[i] is a pointer to the i-th block, for 0 <= i < n,
  // or NULL if block is not yet allocated.
  uint64_t** blocks;

  // list of regions that have been freed by vbuf_free_region, represented as
  // pairs (start, end). They are in sorted order, and "end" for one interval
  // is always strictly less than "start" for the next interval.
  size_t* regions;
  size_t regions_size;    // number of pairs currently stored in "regions"
  size_t regions_alloc;   // number of pairs currently allocated

} vbuf_t;


// initialise a vbuf_t. Initially no blocks are allocated.
#define vbuf_init _NTT_vbuf_init
void vbuf_init(vbuf_t* vbuf, size_t n, pool_t* pool);

#define vbuf_init_wrap _NTT_vbuf_init_wrap
void vbuf_init_wrap(vbuf_t* vbuf, size_t block_size, uint64_t* wrap);

// destroys vbuf, returns all blocks to memory pool
#define vbuf_clear _NTT_vbuf_clear
void vbuf_clear(vbuf_t* vbuf);

// get pointer to i-th block. If not yet allocated, this routine performs the
// allocation. This routine is threadsafe.
#define vbuf_get_block _NTT_vbuf_get_block
uint64_t* vbuf_get_block(vbuf_t* vbuf, size_t i);

// get pointer to index #i within virtual buffer. This routine is threadsafe.
#define vbuf_get _NTT_vbuf_get
uint64_t* vbuf_get(vbuf_t* vbuf, size_t i);

// marks region (start, end) as having been consumed. When a block is fully
// consumed, it is returned to the memory pool (and should not subsequently
// be accessed)
#define vbuf_free_region _NTT_vbuf_free_region
void vbuf_free_region(vbuf_t* vbuf, size_t start, size_t end);


// returns 1 if {x1,n1} overlaps {x2,n2}
#define overlap _NTT_overlap
INLINE int overlap(mp_limb_t* x1, ptrdiff_t n1, mp_limb_t* x2, ptrdiff_t n2)
{
  ptrdiff_t d = x2 - x1;
  return (d < n1) && (d > -n2);
}


// *****************************************************************************
//  single-word modular arithmetic


// 1 if p is probably prime, according to GMP
#define probably_prime _NTT_probably_prime
int probably_prime(uint64_t p);


// x - y mod n, assuming x, y in [0, n), and n < 2^63
#define sub_mod _NTT_sub_mod
INLINE uint64_t sub_mod(uint64_t x, uint64_t y, uint64_t n)
{
  uint64_t z = x - y;
  return z + (((int64_t) z < 0) ? n : 0);
}


// x + y mod n, assuming x, y in [0, n), and n < 2^63
#define add_mod _NTT_add_mod
INLINE uint64_t add_mod(uint64_t x, uint64_t y, uint64_t n)
{
  uint64_t z = x + y;
  return z - ((z >= n) ? n : 0);
}



// !!!! the remaining functions all assume p is odd and in (2^61, 2^62).


// x*y mod p, assuming x*y < 2^64*p
// This is a VERY SLOW implementation using hardware division
#define mul_mod _NTT_mul_mod
INLINE uint64_t mul_mod(uint64_t x, uint64_t y, uint64_t p)
{
#ifdef AVOID_128_BIT
  uint64_t z0, z1, dummy, r;
  MUL128(z1, z0, x, y);
  DIV128(dummy, r, z1, z0, p);
  return r;
#else
  return ((uint128_t) x * y) % p;
#endif
}


// x^n mod p, assuming x in [0, p), n >= 0
// This is a VERY SLOW implementation using hardware division
#define pow_mod _NTT_pow_mod
uint64_t pow_mod(uint64_t x, uint64_t n, uint64_t p);


// 1/x mod p, assuming x in (0, p)
#define inv_mod _NTT_inv_mod
uint64_t inv_mod(uint64_t x, uint64_t p);


// floor(2^126 / p) - 2^64, in [0, 2^64)
#define calc_pinv _NTT_calc_pinv
uint64_t calc_pinv(uint64_t p);


// returns x*y mod p, assumimg x*y < 2^64*p, pinv = calc_pinv(p).
// The result is in [0, 4p).
#define mul_mod_pinv_lazy _NTT_mul_mod_pinv_lazy
INLINE uint64_t mul_mod_pinv_lazy(uint64_t x, uint64_t y,
				  uint64_t p, uint64_t pinv)
{
  /*
    Let z = x*y = 2^62*z1 + z0, and let q = floor(z1*pinv/2^64) + z1. Then
      0 <= z1*pinv/2^64 + z1 - q < 1,
    (Note also
      q <= z1*pinv/2^64 + z1 = z1*(pinv/2^64 + 1) <= z*(pinv/2^64 + 1)/2^62
         < 4*p*(pinv/2^64 + 1) = p*pinv/2^62 + 4*p < (2^126 - 2^64*p)/2^62 + 4p
         <= 2^64.)
    so
      0 <= z1*p*pinv/2^64 + z1*p - q*p < p.
    We also have
      0 <= 2^126/p - (pinv + 2^64) < 1,
    so
      0 <= 2^62*z1 - z1*p*pinv/2^64 - z1*p < z1*p/2^64 < p.
    Adding:
      0 <= 2^62*z1 - q*p < 2p,
    so finally
      0 <= z - q*p < 2p + z0 < 4p.
  */
  uint64_t z0, z1;
#ifdef AVOID_128_BIT
  MUL128(z1, z0, x, y);
  z1 = (z1 << 2) + (z0 >> 62);
#else
  uint128_t z = (uint128_t) x * y;
  z1 = z >> 62;
  z0 = (uint64_t) z;
#endif
  uint64_t q = mulhi(z1, pinv) + z1;
  return z0 - q * p;
}


// same as mul_mod_pinv_lazy, but with result in [0, p)
#define mul_mod_pinv _NTT_mul_mod_pinv
INLINE uint64_t mul_mod_pinv(uint64_t x, uint64_t y,
			     uint64_t p, uint64_t pinv)
{
  uint64_t r = mul_mod_pinv_lazy(x, y, p, pinv);
  r -= ((r >= (2*p)) ? (2*p) : 0);
  r -= ((r >= p) ? p : 0);
  return r;
}


// x^n mod p, assuming x in [0, p), n >= 0, pinv = calc_pinv(p)
#define pow_mod_pinv _NTT_pow_mod_pinv
uint64_t pow_mod_pinv(uint64_t x, uint64_t n, uint64_t p, uint64_t pinv);


// floor(y * 2^64 / p), assuming y in [0, p) and pinv = calc_pinv(p)
#define calc_ypinv _NTT_calc_ypinv
INLINE uint64_t calc_ypinv(uint64_t y, uint64_t p, uint64_t pinv)
{
  /*
    we have
      0 <= 2^126/p - (pinv + 2^64) < 1,
    so
      0 <= 2^64*y/p - pinv*4*y/2^64 - 4*y < 2^(-62)*y < 1.
    Let u = floor(4*y*pinv/2^64) + 4*y. Then u < 2^64, and
      0 <= 4*y*pinv/2^64 + 4*y - u < 1.
    Adding:
      0 <= 2^64*y/p - u < 2.
    Let z = floor(2^64*y/p). Then
      -1 < z - 2^64*y/p <= 0.
    Adding:
      -1 < z - u < 2.
    So the correct value for z is either u or u + 1.
  */
  uint64_t t = 4 * y;
  uint64_t u = mulhi(t, pinv) + t;
  // correction if u is one too small
  return u + ((-u * p) >= p);
}


// x*y mod p, assumimg x in [0, 2^64), y in [0, p),
// ypinv = calc_ypinv(y, p, pinv). The result is in [0, 2p).
// (This is basically Shoup's algorithm from NTL)
#define mul_mod_ypinv_lazy _NTT_mul_mod_ypinv_lazy
INLINE uint64_t mul_mod_ypinv_lazy(uint64_t x, uint64_t y,
				   uint64_t ypinv, uint64_t p)
{
  /*
    we have
      0 <= 2^64*y/p - ypinv < 1,
    so
      0 <= x*y - ypinv*x*p/2^64 < x*p/2^64 < p.
    Let q = floor(x * ypinv / 2^64). Then
      0 <= x*ypinv/2^64 - q < 1
    so
      0 <= x*ypinv*p/2^64 - q*p < p.
    Adding:
      0 <= x*y - q*p < 2p.
  */
  uint64_t q = mulhi(x, ypinv);
  return x * y - q * p;
}


// same as mul_mod_ypinv_lazy(), but with result in [0, p).
#define mul_mod_ypinv _NTT_mul_mod_ypinv
INLINE uint64_t mul_mod_ypinv(uint64_t x, uint64_t y,
			      uint64_t ypinv, uint64_t p)
{
  uint64_t r = mul_mod_ypinv_lazy(x, y, ypinv, p);
  return r - ((r >= p) ? p : 0);
}




// 1/p mod 2^64
#define calc_pinvb _NTT_calc_pinvb
uint64_t calc_pinvb(uint64_t p);


// x*y/2^64 mod p, assumimg x*y < 2^64*p, pinvb = calc_pinvb(p).
// The result is in [0, 2p).
// (This is basically Montgomery's algorithm.)
#define mul_mod_pinvb_lazy _NTT_mul_mod_pinvb_lazy
INLINE uint64_t mul_mod_pinvb_lazy(uint64_t x, uint64_t y,
				   uint64_t p, uint64_t pinvb)
{
  uint64_t z0, z1;
#if AVOID_128_BIT
  MUL128(z1, z0, x, y);
#else
  uint128_t z = (uint128_t) x * y;
  z0 = (uint64_t) z;
  z1 = z >> 64;
#endif
  uint64_t q = z0 * pinvb;
  uint64_t h = mulhi(q, p);
  return z1 - h + p;
}


// same as mul_mod_pinvb_lazy(), but with result in [0, p).
#define mul_mod_pinvb _NTT_mul_mod_pinvb
INLINE uint64_t mul_mod_pinvb(uint64_t x, uint64_t y,
			      uint64_t p, uint64_t pinvb)
{
  uint64_t z0, z1;
#if AVOID_128_BIT
  MUL128(z1, z0, x, y);
#else
  uint128_t z = (uint128_t) x * y;
  z0 = (uint64_t) z;
  z1 = z >> 64;
#endif
  uint64_t q = z0 * pinvb;
  uint64_t h = mulhi(q, p);
  int64_t r = z1 - h;
  return r + ((r < 0) ? p : 0);
}


// *****************************************************************************
//  tuning

// tune_tab[n] is the threshold for switching from GMP to NTT multiplication if
// there are n threads available. If tune_tab[n] == SIZE_MAX, then GMP is always
// used. Length of table is tune_tab_size.
#define tune_tab _NTT_tune_tab
extern size_t tune_tab[];

#define tune_tab_size _NTT_tune_tab_size
extern int tune_tab_size;


// *****************************************************************************
//  FFT utility routines


// 2^k2 * 3^k3
#define pow23 _NTT_pow23
INLINE uint64_t pow23(unsigned k2, unsigned k3)
{
  static uint64_t pow3[] = {
    1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049, 177147, 531441,
    1594323, 4782969, 14348907, 43046721, 129140163, 387420489, 1162261467};
  return pow3[k3] << k2;
}


// allowable FFT transform lengths are 2^k2 * 2^k3
// where k2 <= MAX_K2, k3 <= MAX_K3
#define MAX_K2 36
#define MAX_K3 6
#define MAX_K (pow23(MAX_K2, MAX_K3))


// a few primes p such that p = 1 mod 2^MAX_K2 * 3^MAX_K3 and 2^61 < p < 2^62
// (must be decreasing order)
#define MAX_NUM_PRIMES 6
#define PRIME0 0x3ffdad3000000001ULL
#define PRIME1 0x3ff9954000000001ULL
#define PRIME2 0x3ff57d5000000001ULL
#define PRIME3 0x3ff4f4a000000001ULL
#define PRIME4 0x3ff35a9000000001ULL
#define PRIME5 0x3ff2d1e000000001ULL
//#define PRIME6 0x3ff1ee1000000001ULL
//#define PRIME7 0x3fee5ed000000001ULL
#define global_p _NTT_global_p
extern uint64_t global_p[];

// corresponding values of pinv, pinvb
#define global_pinv _NTT_global_pinv
extern uint64_t global_pinv[];
#define global_pinvb _NTT_global_pinvb
extern uint64_t global_pinvb[];

// global_w[i] is a primitive 2^MAX_K2*3^MAX_K3-th root of unity
// modulo global_p[i].
#define global_w _NTT_global_w
extern uint64_t global_w[];

// finds standard K-th root of unity mod global_p[i], where K must divide MAX_K
#define calc_w _NTT_calc_w
uint64_t calc_w(unsigned i, size_t K);



// *****************************************************************************
//  FFT routines

/*
  Parameters for FFT routines:

    p = modulus, i.e. we are doing a number-theoretic transform mod p.
    K = transform length = 2^k2 * 3^k3, where k2 <= MAX_K2, k3 <= MAX_K3
    n = input length, 0 <= n <= K
    op = input buffer of length n. Input elements are in [0, 2p).
    rop = output buffer of length K. Output elements are in [0, 2p).
    num_threads = number of threads to use

  Some routines allow in-place operation. In this case rop == op, and op must
  have length K (the last K - n inputs are ignored).

  Some routines allow out-of-place operation. In this case {op,n} and {rop,K}
  must not overlap. The input is not modified.

  Some routines can operate on virtual buffers (see vbuf_t), which are made up
  of blocks from the memory pool. In this case K must be a multiple of the
  block size.

  Roots of unity are described in one of two ways. Either:

    w = primitive K-th root, 0 <= w < p.

  Or, for the lowest-level routines:

    wtab = tables of roots. The table wtab->w[i] is used for layer #i of the
      transform, where layer 1 is the bottom layer. Layers 1 <= i <= k2 are
      radix-2, layers k2 < i <= k2 + k3 are radix-3.

      For a radix-2 layer, let L = 2^i be the corresponding transform length,
      and let w' = w^(K/L) be the corresponding primitive L-th root.
      Then wtab->w[i][j] = w'^j for 0 <= j < L/2.

      For a radix-3 layer, let L = 2^k2 * 3^(i - k2), and let w' = w^(K/L).
      Then wtab->w[i][j] = w'^j, wtab->w[i][j+L/3] = w'^(2j) for 0 <= j < L/3.

      wtab->wpinv = corresponding tables for use with Shoup multiplication,
      i.e. wtab->wpinv[i][j] = floor(wtab->w[i][j] * 2^64 / p).

      Finally, entry 0 of both tables is a single entry describing the cube root
      of unity w^(K/3) (if K is divisible by 3, otherwise unused).

  The output is written in "scrambled" order, but the meaning of this depends
  on the decomposition into radix-2 and radix-3 transforms. The decomposition
  chosen depends on k2, k3 and on the memory pool block size. Each fft_XXX()
  routine has a corresponding fft_XXX_pattern() function that returns a
  bitstring P describing the decomposition. Bit #i of P is set if the
  corresponding layer is radix-3, otherwise cleared if it's radix-2.
*/


typedef struct
{
  uint64_t* buf;
  uint64_t** w;      // array of pointers into buf
  uint64_t** wpinv;  // array of pointers into buf
} wtab_t;


// Initialises wtab for given parameters
#define wtab_init _NTT_wtab_init
void wtab_init(wtab_t* wtab, size_t K, unsigned k2, unsigned k3,
	       uint64_t w, uint64_t p);

#define wtab_clear _NTT_wtab_clear
void wtab_clear(wtab_t* wtab);


// inplace operation only
#define fft_base _NTT_fft_base
void fft_base(uint64_t* op, size_t K, unsigned k2, unsigned k3, size_t n,
	      wtab_t* wtab, uint64_t p);

#define fft_base_pattern _NTT_fft_base_pattern
size_t fft_base_pattern(size_t K, unsigned k2, unsigned k3);


// Main FFT entry point. Both in-place and out-of-place operation allowed.
#define fft _NTT_fft
void fft(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
	 uint64_t* op, size_t n, uint64_t w, uint64_t p,
	 int num_threads);

#define fft_pattern _NTT_fft_pattern
size_t fft_pattern(size_t K, unsigned k2, unsigned k3);


// same as fft(), but operates on vbuf_t's
#define fft_vbuf _NTT_fft_vbuf
void fft_vbuf(vbuf_t* rop, size_t K, unsigned k2, unsigned k3,
	      vbuf_t* op, size_t n, uint64_t w, uint64_t p,
	      int num_threads);

#define fft_vbuf_pattern _NTT_fft_vbuf_pattern
size_t fft_vbuf_pattern(size_t K, unsigned k2, unsigned k3, size_t block_size);


/*
  fft_cols and fft_rows together perform an FFT using an array decomposition,
  where R = 2^r2 * 3^r3 = number of rows, C = 2^c2 * 3^c3 = number of columns,
  K = R * C.

  fft_cols performs the column transforms. Both in-place and out-of-place
  operation are allowed.

  fft_rows multiplies by twiddle factors and performs the row transforms.
  It only operates in-place.

  The vbuf versions require that the pool block size is divisible by C.

  Corresponding output order is obtained from fft_cols_rows_pattern().
*/
#define fft_cols _NTT_fft_cols
void fft_cols(uint64_t* rop,
	      size_t R, unsigned r2, unsigned r3,
	      size_t C, unsigned c2, unsigned c3,
	      uint64_t* op, size_t n, uint64_t w, uint64_t p,
	      int num_threads);

#define fft_rows _NTT_fft_rows
void fft_rows(uint64_t* op,
	      size_t R, unsigned r2, unsigned r3,
	      size_t C, unsigned c2, unsigned c3,
	      uint64_t w, uint64_t p,
	      int num_threads);

#define fft_cols_rows_pattern _NTT_fft_cols_rows_pattern
size_t fft_cols_rows_pattern(size_t R, unsigned r2, unsigned r3,
			     size_t C, unsigned c2, unsigned c3);


#define fft_cols_vbuf _NTT_fft_cols_vbuf
void fft_cols_vbuf(vbuf_t* rop,
		   size_t R, unsigned r2, unsigned r3,
		   size_t C, unsigned c2, unsigned c3,
		   vbuf_t* op, size_t n, uint64_t w, uint64_t p,
		   int num_threads);

#define fft_rows_vbuf _NTT_fft_rows_vbuf
void fft_rows_vbuf(vbuf_t* op,
		   size_t R, unsigned r2, unsigned r3,
		   size_t C, unsigned c2, unsigned c3,
		   uint64_t w, uint64_t p,
		   int num_threads);

#define fft_cols_rows_vbuf_pattern _NTT_fft_cols_rows_vbuf_pattern
size_t fft_cols_rows_vbuf_pattern(size_t R, unsigned r2, unsigned r3,
				  size_t C, unsigned c2, unsigned c3,
				  size_t block_size);


/*
  fft_array() chooses an appropriate array decomposition and performs the
  FFT using fft_cols() followed by fft_rows().
*/
#define fft_array _NTT_fft_array
void fft_array(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
	       uint64_t* op, size_t n, uint64_t w, uint64_t p,
	       int num_threads);

#define fft_array_pattern _NTT_fft_array_pattern
size_t fft_array_pattern(size_t K, unsigned k2, unsigned k3);


#define fft_array_vbuf _NTT_fft_array_vbuf
void fft_array_vbuf(vbuf_t* rop, size_t K, unsigned k2, unsigned k3,
		    vbuf_t* op, size_t n, uint64_t w, uint64_t p,
		    int num_threads);

#define fft_array_vbuf_pattern _NTT_fft_array_vbuf_pattern
size_t fft_array_vbuf_pattern(size_t K, unsigned k2, unsigned k3,
			      size_t block_size);



// *****************************************************************************
//  IFFT routines


/*
  The IFFT routines are basically the FFT routines in reverse.
  w = the inverse of the root used for the correpsonding FFT.
  There is no n parameter.
  All IFFTs are inplace.
  Inputs and outputs are in [0, 4p).
*/


#define ifft _NTT_ifft
void ifft(uint64_t* op, size_t K, unsigned k2, unsigned k3,
	  uint64_t w, uint64_t p, int num_threads);

#define ifft_vbuf _NTT_ifft_vbuf
void ifft_vbuf(vbuf_t* op, size_t K, unsigned k2, unsigned k3,
	       uint64_t w, uint64_t p, int num_threads);


#define ifft_base _NTT_ifft_base
void ifft_base(uint64_t* op, size_t K, unsigned k2, unsigned k3,
	       wtab_t* wtab, uint64_t p);


#define ifft_cols _NTT_ifft_cols
void ifft_cols(uint64_t* op,
	       size_t R, unsigned r2, unsigned r3,
	       size_t C, unsigned c2, unsigned c3,
	       uint64_t w, uint64_t p, int num_threads);

#define ifft_rows _NTT_ifft_rows
void ifft_rows(uint64_t* op,
	       size_t R, unsigned r2, unsigned r3,
	       size_t C, unsigned c2, unsigned c3,
	       uint64_t w, uint64_t p, int num_threads);

#define ifft_cols_vbuf _NTT_ifft_cols_vbuf
void ifft_cols_vbuf(vbuf_t* op,
		    size_t R, unsigned r2, unsigned r3,
		    size_t C, unsigned c2, unsigned c3,
		    uint64_t w, uint64_t p, int num_threads);

#define ifft_rows_vbuf _NTT_ifft_rows_vbuf
void ifft_rows_vbuf(vbuf_t* op,
		    size_t R, unsigned r2, unsigned r3,
		    size_t C, unsigned c2, unsigned c3,
		    uint64_t w, uint64_t p, int num_threads);

#define ifft_array _NTT_ifft_array
void ifft_array(uint64_t* op, size_t K, unsigned k2, unsigned k3,
		uint64_t w, uint64_t p, int num_threads);

#define ifft_array_vbuf _NTT_ifft_array_vbuf
void ifft_array_vbuf(vbuf_t* op, size_t K, unsigned k2, unsigned k3,
		     uint64_t w, uint64_t p, int num_threads);


// *****************************************************************************
//  convolution routines

/*
  rop[i] = op1[i] * op2[i] / 2^64 mod p for 0 <= i < n.
  rop must be either disjoint from or identical to op1.
  rop must be either disjoint from or identical to op2.
  Inputs in [0, 2p), outputs in [0, 2p).
*/
#define pointwise_multiply _NTT_pointwise_multiply
void pointwise_multiply(uint64_t* rop, uint64_t* op1, uint64_t* op2, size_t n,
			uint64_t p, uint64_t pinvb, int num_threads);


/*
  Computes cyclic convolution of {op1,n1} and {op2,n2},
  multiplied by 1/2^64 mod p, written to {rop,K}.
  Must have 1 <= n1 <= K and 1 <= n2 <= K, and 0 <= u < p.
  Inputs are in [0, 2p), outputs in [0, 4p).

  temp1 and temp2 are buffers of length K used to store the fourier transforms
  of op1 and op2. The routine first transforms op1, then op2 (unless {op1,n1}
  == {op,n2}, in which case the second transform is skipped), then multiplies
  them together into rop, then inverse transforms rop. So it is allowed for
  temp1 and/or temp2 to overlap op1, op2, rop, as long as the data dependencies
  in the above algorithm are not obviously violated.
*/
#define conv _NTT_conv
void conv(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
	  uint64_t* op1, size_t n1, uint64_t* temp1,
	  uint64_t* op2, size_t n2, uint64_t* temp2,
	  uint64_t w, uint64_t p, int num_threads);


#define conv_cols_rows_vbuf _NTT_conv_cols_rows_vbuf
void conv_cols_rows_vbuf(vbuf_t* rop,
			 size_t R, unsigned r2, unsigned r3,
			 size_t C, unsigned c2, unsigned c3,
			 vbuf_t* op1, size_t n1, vbuf_t* temp1,
			 vbuf_t* op2, size_t n2, vbuf_t* temp2,
			 uint64_t w, uint64_t p, int num_threads);

#define conv_cols_rows _NTT_conv_cols_rows
void conv_cols_rows(uint64_t* rop,
		    size_t R, unsigned r2, unsigned r3,
		    size_t C, unsigned c2, unsigned c3,
		    uint64_t* op1, size_t n1, uint64_t* temp1,
		    uint64_t* op2, size_t n2, uint64_t* temp2,
		    uint64_t w, uint64_t p, int num_threads);

#define conv_array _NTT_conv_array
void conv_array(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
		uint64_t* op1, size_t n1, uint64_t* temp1,
		uint64_t* op2, size_t n2, uint64_t* temp2,
		uint64_t w, uint64_t p, int num_threads);


#define conv_array_vbuf _NTT_conv_array_vbuf
void conv_array_vbuf(vbuf_t* rop, size_t K, unsigned k2, unsigned k3,
		     vbuf_t* op1, size_t n1, vbuf_t* temp1,
		     vbuf_t* op2, size_t n2, vbuf_t* temp2,
		     uint64_t w, uint64_t p, int num_threads);

#define conv_base _NTT_conv_base
void conv_base(uint64_t* rop, size_t K, unsigned k2, unsigned k3,
	       uint64_t* op1, size_t n1, uint64_t* temp1,
	       uint64_t* op2, size_t n2, uint64_t* temp2,
	       wtab_t* wtab, wtab_t* winvtab, uint64_t p, uint64_t pinvb);



// *****************************************************************************
//  integer multiplication routines


/*
  Input is integer op of n limbs.
  Let op = d_0 + d_1*2^r + ... + d_{m-1}*2^(r(m-1))  mod 2^(rm),
  where 0 <= d_i < 2^r for each i.
  Writes d_0, d_1, ..., d_{m-1} to rop, exactly m * ceil(r / 64) limbs.
  Requires 1 <= r <= 3*64, n >= 0, m >= 1.
*/
#define split _NTT_split
void split(mp_limb_t* rop, size_t m, mp_limb_t* op, size_t n, unsigned r);


/*
  Input is an array d_0, d_1, ..., d_{m-1}, each t limbs, stored
  consecutively in op (t*m limbs altogether).
  Output is rop[i] = u * d_i mod p, in the range [0, 2p), for 0 <= i < m,
  where 0 <= u < p.
  Requires 1 <= t <= 3.
*/
#define reduce _NTT_reduce
void reduce(uint64_t* rop, mp_limb_t* op, size_t m, unsigned t,
	    uint64_t u, uint64_t p, uint64_t pinv, uint64_t pinvb);


/*
  Combines operations of split_mpn() and reduce() for the first num_primes
  primes in global_p[] array. Output for prime #i is stored at rop[i].
  n == 0 is allowed.
*/
#define split_reduce _NTT_split_reduce
void split_reduce(uint64_t** rop, size_t m, mp_limb_t* op, size_t n,
		  unsigned r, uint64_t* u, unsigned num_primes,
		  int num_threads);


/*
  Same as split_reduce, but writes output array #i to rop[i].
  If destroy_op == 1, blocks of memory from "op" are contributed to the pool
  as they are consumed.
  Requires pool block size to be divisible by 64.
*/
#define split_reduce_vbuf _NTT_split_reduce_vbuf
void split_reduce_vbuf(vbuf_t* rop, size_t m, mp_limb_t* op, size_t n,
		       unsigned r, uint64_t* u, unsigned num_primes,
		       int destroy_op, int num_threads);


/*
  Global precomputed data for CRT. For 1 <= i < MAX_NUM_PRIMES, we have:

  s[i] = -(2^64)^(i-1) / (p[0]*...*p[i-1]) mod p[i]
  spinv[i] = floor(2^64 * s[i] / p[i])
  u[i] = p[0]*...*p[i-1], stored in i limbs.
*/
#define global_crt_s _NTT_global_crt_s
extern uint64_t global_crt_s[MAX_NUM_PRIMES];
#define global_crt_spinv _NTT_global_crt_spinv
extern uint64_t global_crt_spinv[MAX_NUM_PRIMES];
#define global_crt_u _NTT_global_crt_u
extern mp_limb_t global_crt_u[MAX_NUM_PRIMES][MAX_NUM_PRIMES];


/*
  Input is arrays op[0], op[1], ..., op[num_primes-1], each of length m.
  Entries of op[i] are interpreted mod global_p[i].
  Inputs are in [0, 4p).
  For each 0 <= j < m, computes CRT of op[0][j], ..., op[num_primes-1][j],
  writes result (num_primes limbs) starting at rop[j*num_primes], i.e. total
  output is num_primes * m limbs.
*/
#define crt _NTT_crt
void crt(mp_limb_t* rop, uint64_t** op, size_t m, unsigned num_primes);


// same as crt(), but reads from a vbuf, starting at position src
#define crt_from_vbuf _NTT_crt_from_vbuf
void crt_from_vbuf(mp_limb_t* rop, vbuf_t* op, size_t src,
		   size_t m, unsigned num_primes);


/*
  Input is a sequence of m integers, s limbs each, stored consecutively in op.
  Call these d[i], 0 <= i < m.
  Output is 2^r0*(d[0] + d[1]*2^r + ... + d[m-1]*2^((m-1)r)) mod (2^64)^n,
  stored at rop.
  Must have m >= 0, n >= 1, s >= r/64.
*/
#define recompose _NTT_recompose
void recompose(mp_limb_t* rop, size_t n, mp_limb_t* op, size_t m, unsigned s,
	       unsigned r, unsigned r0);



/*
  Combination of crt() and recompose().
  Input is same as crt(). Output is same as recompose().
  Processed in blocks for better locality.
*/
#define crt_recompose _NTT_crt_recompose
void crt_recompose(mp_limb_t* rop, size_t n, unsigned r,
		   uint64_t** op, size_t m, unsigned num_primes,
		   int num_threads);


/*
  Same as crt_recompose, but operates on vbuf's.
  If destroy_op == 1, then blocks of op are returned to the memory pool as they
  are consumed (in which case op must not be a "wrapped" vbuf)
*/
#define crt_recompose_vbuf _NTT_crt_recompose_vbuf
void crt_recompose_vbuf(vbuf_t* rop, size_t n, unsigned r,
			vbuf_t* op, size_t m, int destroy_op,
			unsigned num_primes, int num_threads);



// bit_bound[i] = 0.999 * log2(p[0] * ... * p[i-1])
#define global_bit_bound _NTT_global_bit_bound
extern double global_bit_bound[MAX_NUM_PRIMES + 1];


// low-memory version of ntt_mpn_mul_bonus
// assumes n1 >= n2
#define ntt_mpn_mul_lowmem _NTT_mpn_mul_lowmem
void ntt_mpn_mul_lowmem(mp_limb_t* rop,
			mp_limb_t* op1, size_t n1, int preserve_op1,
			mp_limb_t* op2, size_t n2, int preserve_op2,
			int num_threads);
