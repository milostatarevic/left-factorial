/*
  This file is part of the left-factorial repository.
  Copyright (C) 2019 Milos Tatarevic

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, write to the Free Software Foundation,
  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <stdio.h>
#include <string.h>
#include <time.h>

#include <gmp.h>
#include <flint.h>
#include <nmod_poly.h>
#include <NTL/lzz_pX.h>

#define BASE 4

using namespace NTL;

slong *t_Delta, *inv_delta, *inv_table;

void init_tables(slong n)
{
  slong table_size;

  table_size = (n + 256) * sizeof(slong);

  t_Delta = (slong *) malloc(table_size);
  inv_delta = (slong *) malloc(table_size);
  inv_table = (slong *) malloc(2 * table_size);
}

static inline slong inv(slong i, nmod_t mod)
{
  return n_invmod(n_mod2_precomp(i, mod.n, mod.ninv), mod.n);
}

static inline void clear_Delta_cache(void)
{
  t_Delta[0] = -1;
  inv_table[0] = -1;
}

static inline void clear_delta_caches(void)
{
  inv_delta[0] = -1;
  clear_Delta_cache();
}

void init_inv_delta(slong d, nmod_t mod)
{
  slong i, r;

  r = 1;
  for(i = 1; i <= d; i++)
    r = n_mulmod2_preinv(r, mod.n - i, mod.n, mod.ninv);
  inv_delta[0] = inv(r, mod);

  for(i = 1; i <= d; i ++)
  {
    r = n_mulmod2_preinv(r, inv(i - d - 1, mod), mod.n, mod.ninv);
    r = n_mulmod2_preinv(r, i, mod.n, mod.ninv);

    inv_delta[i] = inv(r, mod);
  }
}

void init_inv_tables(slong d, slong a, nmod_t mod)
{
  slong i;

  if(inv_table[0] == -1)
    for(i = 0; i <= 2 * d; i ++)
      inv_table[i] = inv(a + i - d, mod);

  if(inv_delta[0] == -1)
    init_inv_delta(d, mod);
}

slong Delta(slong a, slong k, slong d, nmod_t mod)
{
  slong j, r;

  if(t_Delta[0] == -1)
  {
    r = 1;
    for(j = 0; j <= d; j ++)
      r = n_mulmod2_preinv(r, a - j, mod.n, mod.ninv);
    t_Delta[0] = r;

    for(j = 1; j <= d; j ++)
    {
      r = n_mulmod2_preinv(r, inv_table[j - 1], mod.n, mod.ninv);
      r = n_mulmod2_preinv(r, a + j, mod.n, mod.ninv);

      t_Delta[j] = r;
    }
  }

  return t_Delta[k];
}

void _shift(mp_ptr p, slong d, slong a, nmod_t mod)
{
  slong i, k, v;
  zz_pX ntl_p, ntl_s, ntl_q;

  init_inv_tables(d, a, mod);
  for(i = 0; i <= d; i ++)
    p[i] = n_mulmod2_preinv(p[i], inv_delta[i], mod.n, mod.ninv);

  // use NTL for the fast poly mul
  ntl_q.SetLength(3 * d + 1);
  ntl_s.SetLength(2 * d + 1);
  ntl_p.SetLength(d + 1);

  for(i = 0; i <= d; i ++)
    ntl_p[i] = p[i];
  for(i = 0; i <= 2*d; i ++)
    ntl_s[i] = inv_table[i];

  mul(ntl_q, ntl_s, ntl_p);

  //
  for(k = 0, i = d; k <= d; k ++, i ++)
  {
    conv(v, ntl_q[i]);
    p[k] = n_mulmod2_preinv(Delta(a, k, d, mod), v, mod.n, mod.ninv);
  }
}

void _shift_degree(mp_ptr a0, mp_ptr a1, slong a, slong step, slong len, nmod_t mod)
{
  slong i, d, v;
  mp_ptr b0, b1;

  d = len + 1;

  b0 = _nmod_vec_init(2 * d + 1);
  b1 = _nmod_vec_init(2 * d + 1);

  for (i = 0; i < d; i ++)
  {
    a0[i + d] = b0[i] = b0[i + d] = a0[i];
    a1[i + d] = b1[i] = b1[i + d] = a1[i];
  }

  clear_delta_caches();
  _shift(a0 + d, d - 1, a, mod);
  _shift(a1 + d, d - 1, a, mod);

  v = n_mulmod2_preinv(d - 1, inv(step, mod), mod.n, mod.ninv);

  clear_Delta_cache();
  _shift(b0, d - 1, v, mod);
  _shift(b1, d - 1, v, mod);

  clear_Delta_cache();
  _shift(b0 + d, d - 1, v + d, mod);
  _shift(b1 + d, d - 1, v + d, mod);

  // merge
  for(i = 0; i <= 2 * d; i ++)
  {
    a1[i] += n_mulmod2_preinv(a0[i], b1[i], mod.n, mod.ninv);
    a1[i] = n_mod2_precomp(a1[i], mod.n, mod.ninv);

    a0[i] = n_mulmod2_preinv(a0[i], b0[i], mod.n, mod.ninv);
  }
}

void rp_sqr(slong n, nmod_t mod, slong u_start, slong b_shift, slong *r0, slong *r1)
{

  slong i, j, a_0, a_1, b_0, b_1, n_shifted;
  mp_ptr a0, a1, a0_t, a1_t;

  a0 = _nmod_vec_init(n * 4 + 1);
  a1 = _nmod_vec_init(n * 4 + 1);
  a0_t = _nmod_vec_init(n + 1);
  a1_t = _nmod_vec_init(n + 1);

  a0[0] = u_start;
  a0[1] = n + u_start;

  a1[0] = a1[1] = 1;

  for(j = 1; j < n; j *= 2)
    _shift_degree(a0, a1, j + 1, n, j, mod);

  // b_shift
  if (b_shift > 1)
  {
    _nmod_vec_set(a0_t, a0, n + 1);
    _nmod_vec_set(a1_t, a1, n + 1);

    clear_delta_caches();
    for (i = 1; i < b_shift; i ++)
    {
      _shift(a0_t, n, n + 1, mod);
      _shift(a1_t, n, n + 1, mod);

      _nmod_vec_set(a0 + i * (n + 1), a0_t, n + 1);
      _nmod_vec_set(a1 + i * (n + 1), a1_t, n + 1);
    }
  }
  n_shifted = n * b_shift;

  a_0 = a0[0];
  a_1 = a1[0];

  for (i = 1; i < n_shifted; i++)
  {
    b_0 = a0[i];
    b_1 = a1[i];

    a_1 += n_mulmod2_preinv(a_0, b_1, mod.n, mod.ninv);
    a_1 = n_mod2_preinv(a_1, mod.n, mod.ninv);

    a_0 = n_mulmod2_preinv(a_0, b_0, mod.n, mod.ninv);
  }

  *r0 = a_0;
  *r1 = a_1;
}

slong calc(slong p)
{
  slong i, t, b_len, u_start, a_0, a_1, b_0, b_1;
  mp_ptr b;
  nmod_t mod;

  nmod_init(&mod, p);
  zz_p::init(mod.n);

  b = _nmod_vec_init(50);

  t = p - 1;
  for(i = 0; t != 0; i ++)
  {
    b[i] = t % BASE;
    t = t / BASE;
  }
  b_len = i - 1;

  init_tables(1UL << b_len);

  u_start = 1;
  a_1 = a_0 = 0;
  for (i = b_len; i >= 0; i --)
  {
    if(b[i] == 0) continue;

    rp_sqr(1UL << i, mod, u_start, b[i], &b_0, &b_1);
    u_start += b[i] * (1UL << (2*i));

    if (i == b_len)
    {
      a_0 = b_0;
      a_1 = b_1;
    }
    else {
      a_1 += n_mulmod2_preinv(a_0, b_1, mod.n, mod.ninv);
      a_1 = n_mod2_precomp(a_1, mod.n, mod.ninv);

      a_0 = n_mulmod2_preinv(a_0, b_0, mod.n, mod.ninv);
    }
  }
  a_1 += mod.n - 1;
  a_1 = n_mod2_preinv(a_1, mod.n, mod.ninv);

  return a_1;
}

int main(int argc, char ** argv)
{
  slong p, r_p;

  if (argc != 2)
  {
    printf("usage: %s p\n", argv[0]);
    return 1;
  }

  p = atoll(argv[1]);
  r_p = calc(p);

  printf("%lu\t%lu\n", p, r_p);
  return 0;
}
