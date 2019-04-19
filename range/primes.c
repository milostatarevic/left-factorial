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

void generate_primes()
{
  ulong i, j, p, pp, n_primes, s, segment;
  mpz_t *list_p;
  n_primes_t iter;

  list_p = mpz_list_init(state.segment_len[0]);

  n_primes_init(iter);
  p = state.p_start;
  n_primes_jump_after(iter, p + 1);

  n_primes = segment = s = 0;
  for (i = 0, j = state.p_start; j <= state.p_end; i ++, j ++)
  {
    pp = 1;
    if (p == j)
    {
      if ((p & 1) == 1) pp = j;

      p = n_primes_next(iter);
      n_primes ++;
    }

    mpz_set_ui(list_p[s], pp);
    s ++;

    if (s >= state.segment_len[0])
    {
      save_segment(list_p, PHASE_P_TREE, s, 0, segment, 0);
      s = 0;
      segment ++;
    }
  }
  n_primes_clear(iter);

  save_segment(list_p, PHASE_P_TREE, s, 0, segment, 0);
  mpz_list_free(list_p, state.segment_len[0]);
}
