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

#include "defs.h"

ulong mpz_sum(mpz_t *n)
{
  ulong i, size, sum;

  size = n[0]->_mp_size;
  sum = 0;
  for (i = 0; i < size; i ++)
    sum += n[0]->_mp_d[i];
  return sum;
}

void load_value(inverse_t *inverse)
{
  FILE *f;
  ulong sum;

  f = open_file("rb", "%s/inverse.value", folders.inverse);
  mpz_init(inverse->value);
  mpz_read(f, &inverse->value);
  fclose(f);

  sum = mpz_sum(&inverse->value);
  if (inverse->value_sum == 0)
    inverse->value_sum = sum;
  else if (inverse->value_sum != sum)
    error("value sum mismatch %lu != %lu", inverse->value_sum, sum);
}

void load_inverse(inverse_t *inverse)
{
  FILE *f;
  ulong sum;

  f = open_file("rb", "%s/inverse.inverse", folders.inverse);
  mpz_init(inverse->inverse);
  mpz_read(f, &inverse->inverse);
  fclose(f);

  sum = mpz_sum(&inverse->inverse);
  if (inverse->inverse_sum == 0)
    inverse->inverse_sum = sum;
  else if (inverse->inverse_sum != sum)
    error("inverse sum mismatch %lu != %lu", inverse->inverse_sum, sum);
}

void copy_value(inverse_t *inverse)
{
   exec("cp %s/c%lu-l%lu-e%lu-t%lu %s/inverse.value",
        folders.tree, PHASE_P_TREE, inverse->h_end, 0, 0, folders.inverse);
}

void save_inverse(inverse_t *inverse)
{
  FILE *f;

  f = open_file("wb", "%s/inverse.inverse", folders.inverse);
  mpz_write(f, &inverse->inverse);
  fclose(f);
}

void init_inverse(inverse_t *inverse, ulong h_end)
{
  mpz_t n;

  inverse->value_sum = 0;
  inverse->inverse_sum = 0;
  inverse->h_end = h_end;

  copy_value(inverse);
  load_value(inverse);

  inverse->bits = 2 * inverse->value[0]._mp_size * sizeof(mp_limb_t) * 8 + 1;

  mpz_init_set_ui(n, 2);
  mpz_mul_2exp(n, n, inverse->bits - 1);

  mpz_init(inverse->inverse);
  ntt_mpz_fdiv_q(inverse->inverse, n, inverse->value);
  mpz_clear(n);

  save_inverse(inverse);
  inverse->blank = 0;
}

void reduce_mod(mpz_t *r, inverse_t *inverse)
{
  mpz_t t, d;

  if (inverse == NULL || inverse->blank || mpz_cmp(*r, inverse->value) <= 0)
    return;

  if (mpz_cmp(*r, inverse->value) == 0)
  {
    mpz_set_ui(*r, 0);
    return;
  }

  if (state.use_ntt == 0)
  {
    mpz_mod(*r, *r, inverse->value);
    return;
  }

  if (r[0]->_mp_size > inverse->bits / sizeof(mp_limb_t) / 8)
  {
    error("insufficient inverse size: r: %d", r[0]->_mp_size);
    mpz_mod(*r, *r, inverse->value);
    return;
  }

  mpz_init(t);
  mpz_init(d);

  // save extra memory
  mpz_clear(inverse->value);

  //
  load_inverse(inverse);
  m_mul(&t, r, &inverse->inverse);
  mpz_clear(inverse->inverse);

  mpz_fdiv_q_2exp(t, t, inverse->bits);
  mpz_realloc2(t, t[0]._mp_size * sizeof(mp_limb_t) * 8);

  //
  load_value(inverse);
  m_mul(&d, &t, &inverse->value);
  mpz_clear(t);

  mpz_sub(*r, *r, d);
  if (mpz_cmp(inverse->value, *r) <= 0)
    mpz_sub(*r, *r, inverse->value);
  mpz_clear(d);
}
