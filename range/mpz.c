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

mpz_t *mpz_list_init(ulong len)
{
  ulong i;
  mpz_t *list;

  list = (mpz_t *) malloc(len * sizeof(mpz_t));
  for(i = 0; i < len; i ++)
    mpz_init(list[i]);
  return list;
}

void mpz_list_free(mpz_t *list, ulong len)
{
  ulong i;

  for(i = 0; i < len; i ++)
    mpz_clear(list[i]);
  free(list);
}

void mpz_write(FILE *f, mpz_t *n)
{
  ulong size;

  size = n[0]->_mp_size;
  fwrite(&size, sizeof(ulong), 1, f);
  fwrite(n[0]->_mp_d, sizeof(mp_limb_t) * size, 1, f);
}


void mpz_read(FILE *f, mpz_t *n)
{
  ulong size;

  fread(&size, sizeof(ulong), 1, f);
  mpz_realloc2(*n, size * sizeof(mp_limb_t) * 8);

  n[0]->_mp_size = size;
  fread(n[0]->_mp_d, sizeof(mp_limb_t) * size, 1, f);
}

void fast_mod(mpz_t *r, mpz_t *p)
{
  if (p == NULL || mpz_cmp(*r, *p) < 0)
    return;

  if (state.use_ntt_mod)
    ntt_mpz_mod(*r, *r, *p);
  else
    mpz_mod(*r, *r, *p);
}

void m_mul(mpz_t *rop, mpz_t *op1, mpz_t *op2)
{
  slong i;
  ulong n_size;
  mpz_t *t;

  if (op1[0]->_mp_size < op2[0]->_mp_size)
    { t = op1; op1 = op2; op2 = t; }

  if (op1[0]->_mp_size > 0 && op2[0]->_mp_size > 0 && state.use_ntt)
  {
    if (rop[0]->_mp_d == op1[0]->_mp_d ||
        rop[0]->_mp_d == op2[0]->_mp_d)
      error("m_mul operand error!");

    n_size = op1[0]->_mp_size + op2[0]->_mp_size;
    mpz_realloc2(*rop, n_size * sizeof(mp_limb_t) * 8);
    memset(rop[0]->_mp_d, 0, n_size * sizeof(mp_limb_t));

    ntt_mpn_mul_lowmem(rop[0]->_mp_d, op1[0]->_mp_d, op1[0]->_mp_size, 1,
                                      op2[0]->_mp_d, op2[0]->_mp_size, 1,
                                      state.num_threads);

    rop[0]->_mp_size = n_size;
    for(i = n_size - 1; i >= 0 && rop[0]->_mp_d[i] == 0; i --)
      n_size --;

    if (n_size != rop[0]->_mp_size)
    {
      mpz_realloc2(*rop, n_size * sizeof(mp_limb_t) * 8);
      rop[0]->_mp_size = n_size;
    }
  }
  else
    mpz_mul(*rop, *op1, *op2);
}
