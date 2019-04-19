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

#include "../ntt/ntt.h"

int gmp_ntt_threads;

void set_gmp_ntt_threads (int num_threads)
{
  gmp_ntt_threads = num_threads;
}

void ntt_mpn_mul_opt (mp_limb_t* rop, const mp_limb_t* op1, size_t n1,
                                      const mp_limb_t* op2, size_t n2)
{
  if(n1 == 0 || n2 == 0)
  {
    mpn_mul(rop, op1, n1, op2, n2);
    return;
  }

  ntt_mpn_mul_bonus(rop, (mp_limb_t*) op1, n1, 1,
                         (mp_limb_t*) op2, n2, 1, 1, gmp_ntt_threads);
}
