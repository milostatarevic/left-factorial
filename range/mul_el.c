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

void remove_element(ulong phase, ulong h, ulong i, ulong t)
{
  if (phase == PHASE_MUL && h == 0)
    return;

  exec("rm %s/c%lu-l%lu-e%lu-t%lu", folders.tree, phase, h, i, t);
}

void reduce_and_save_element(mpz_t *r, mpz_t *p, ulong phase, ulong h, ulong i, ulong t)
{
  if (phase != PHASE_N_TREE)
    return;

  // delete unused files
  if ((i & 1) == 1)
  {
    exec("rm %s/c%lu-l%lu-e%lu-t%lu", folders.tree, phase, h, i, t);
    return;
  }

  fast_mod(r, p);
  save_element(r, phase, h, i, t, NULL);
}

void reduce_elements(ulong phase, ulong h, ulong j, inverse_t *inverse)
{
  mpz_t a0;

  mpz_init(a0);
  load_element(&a0, phase, h, j, 0);
  save_element(&a0, phase, h + 1, j, 0, inverse);
  mpz_clear(a0);

  mpz_init(a0);
  load_element(&a0, phase, h, j, 1);
  save_element(&a0, phase, h + 1, j, 1, inverse);
  mpz_clear(a0);
}

void mul_elements_phase_mul(ulong phase, ulong h, ulong i, ulong j, inverse_t *inverse)
{
  mpz_t r, a0, a1, b0, b1;

  printf("mul :: phase: %lu - h: %lu, %lu\n", phase, h, j);

  mpz_init(r);
  mpz_init(a0);
  mpz_init(a1);
  mpz_init(b0);
  mpz_init(b1);

  if (i <= 2 * j + 1)
  {
    // 0
    load_element(&a0, phase, h, 2 * j, 0);
    save_element(&a0, phase, h + 1, j, 0, inverse);

    mpz_clear(a0);

    // 1
    load_element(&a1, phase, h, 2 * j, 1);
    save_element(&a1, phase, h + 1, j, 1, inverse);

    mpz_clear(a1);

    remove_element(phase, h, 2 * j, 0);
    remove_element(phase, h, 2 * j, 1);
  }
  else
  {
    // 0
    load_element(&a0, phase, h, 2 * j, 0);
    load_element(&b0, phase, h, 2 * j + 1, 0);

    m_mul(&r, &a0, &b0);

    mpz_clear(b0);

    if (phase == PHASE_MUL_DIFFERENCE &&
       !state.last_h[phase] && r[0]._mp_size + 1 > inverse->value[0]._mp_size)
    {
      state.last_h[phase] = h + 1;
      printf("last_h: %lu\n", state.last_h[phase]);
    }

    save_element(&r, phase, h + 1, j, 0, inverse);

    // 1
    load_element(&b1, phase, h, 2 * j + 1, 1);

    m_mul(&r, &a0, &b1);

    mpz_clear(a0);
    mpz_clear(b1);

    load_element(&a1, phase, h, 2 * j, 1);
    mpz_add(r, r, a1);

    mpz_clear(a1);

    save_element(&r, phase, h + 1, j, 1, inverse);

    remove_element(phase, h, 2 * j, 0);
    remove_element(phase, h, 2 * j, 1);
    remove_element(phase, h, 2 * j + 1, 0);
    remove_element(phase, h, 2 * j + 1, 1);
  }
  mpz_clear(r);
}

void mul_elements(ulong phase, ulong h, ulong i, ulong j, inverse_t *inverse)
{
  mpz_t r, a0, a1, b0, b1, p0;
  ulong pi;

  printf("mul :: phase: %lu - h: %lu, %lu\n", phase, h, j);

  mpz_init(r);
  mpz_init(a0);
  mpz_init(a1);
  mpz_init(b0);
  mpz_init(b1);

  if (phase == PHASE_N_TREE)
  {
    mpz_init(p0);

    pi = ((2 * j + 1 < state.level_len[h]) ? (2 * j + 1) : (2 * j));
    load_element(&p0, PHASE_P_TREE, h, pi, 0);
  }

  if (i <= 2 * j + 1)
  {
    // 0
    load_element(&a0, phase, h, 2 * j, 0);
    save_element(&a0, phase, h + 1, j, 0, inverse);

    reduce_and_save_element(&a0, &p0, phase, h, 2 * j, 0);
    mpz_clear(a0);

    // 1
    if (phase != PHASE_P_TREE)
    {
      load_element(&a1, phase, h, 2 * j, 1);
      save_element(&a1, phase, h + 1, j, 1, inverse);

      reduce_and_save_element(&a1, &p0, phase, h, 2 * j, 1);
    }
    mpz_clear(a1);
  }
  else
  {
    // 0
    load_element(&a0, phase, h, 2 * j, 0);
    load_element(&b0, phase, h, 2 * j + 1, 0);

    m_mul(&r, &a0, &b0);

    reduce_and_save_element(&b0, &p0, phase, h, 2 * j + 1, 0);
    mpz_clear(b0);

    if (phase == PHASE_N_TREE &&
       !state.last_h[phase] && r[0]._mp_size + 1 > inverse->value[0]._mp_size)
    {
      state.last_h[phase] = h + 1;
      printf("last h: %lu\n", state.last_h[phase]);
    }

    save_element(&r, phase, h + 1, j, 0, inverse);

    // 1
    if (phase == PHASE_P_TREE)
    {
      mpz_clear(a0);
      mpz_clear(a1);
      mpz_clear(b1);
    }
    else
    {
      load_element(&b1, phase, h, 2 * j + 1, 1);

      m_mul(&r, &a0, &b1);

      reduce_and_save_element(&a0, &p0, phase, h, 2 * j, 0);
      mpz_clear(a0);

      reduce_and_save_element(&b1, &p0, phase, h, 2 * j + 1, 1);
      mpz_clear(b1);

      load_element(&a1, phase, h, 2 * j, 1);
      mpz_add(r, r, a1);

      reduce_and_save_element(&a1, &p0, phase, h, 2 * j, 1);
      mpz_clear(a1);

      save_element(&r, phase, h + 1, j, 1, inverse);
    }
  }
  mpz_clear(r);

  if (phase == PHASE_N_TREE)
    mpz_clear(p0);
}

ulong mul_all_elements(ulong phase, ulong len, ulong h, inverse_t *inverse)
{
  ulong i, j, d, t, steps, vh, vj,
        from_h, from_j, j_start, counter;

  steps = 1 << MUL_STEP_DEPTH;
  counter = 0;

  if (phase == PHASE_MUL)
    len *= 2;

  for (i = len; i > 1; i = (i + 1) / 2, h ++)
  {
    printf("h: %lu, size: %lu\n", h, i);

    d = (i + 1) / 2;
    for(j = 0; j < d; j ++)
    {
      if (phase != PHASE_MUL || h)
      {
        if (phase == PHASE_MUL || phase == PHASE_MUL_DIFFERENCE)
          mul_elements_phase_mul(phase, h, i, j, inverse);
        else
          mul_elements(phase, h, i, j, inverse);
      }
      else
      {
        // save space by multiplying PHASE_MUL levels from the saved list

        reduce_elements(phase, h, j, inverse);
        printf("reduce :: phase: %lu - h: %lu, %lu\n", phase, h, j);

        from_h = 0;
        if ((j & (steps - 1)) == steps - 1)
        {
          j_start = j - steps + 1;

          for(vh = 1; vh <= MUL_STEP_DEPTH; vh ++)
            for(vj = 0; vj < 1 << (MUL_STEP_DEPTH - vh); vj ++)
              mul_elements_phase_mul(phase, vh, i, (j_start >> vh) + vj, inverse);

          from_j = counter;
          from_h = vh;
        }
        else if (j >= (d / steps) * steps)
        {
          from_j = j;
          from_h = h + 1;
        }

        if (from_h > 0)
        {
          for (t = 0; t < 2; t ++)
            if (from_h != h + 1 || from_j != counter)
              exec("mv -vn %s/c%lu-l%lu-e%lu-t%lu %s/c%lu-l%lu-e%lu-t%lu",
                   folders.tree, phase, from_h, from_j, t,
                   folders.tree, phase, h + 1, counter, t);
          counter ++;
        }
      }
    }

    if (counter > 0)
    {
      i = counter * 2;
      counter = 0;
    }

    if (phase == PHASE_MUL_DIFFERENCE && state.last_h[phase])
      break;
  }

  return h;
}
