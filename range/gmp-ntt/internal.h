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

#include <gmp.h>
#include <gmp-impl.h>
#include <gmp/longlong.h>

void ntt_mpz_tdiv_r (mpz_ptr, mpz_srcptr, mpz_srcptr);
void ntt_mpz_tdiv_qr (mpz_ptr, mpz_ptr, mpz_srcptr, mpz_srcptr);

void ntt_mpn_mul_opt (mp_limb_t *, const mp_limb_t *, size_t, const mp_limb_t *, size_t);
void ntt_mpn_mulmod_bnm1 (mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, mp_ptr);
void ntt_mpn_tdiv_qr (mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
mp_limb_t ntt_mpn_mu_div_qr (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, mp_ptr);
