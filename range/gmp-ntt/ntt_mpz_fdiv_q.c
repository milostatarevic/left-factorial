/*
  This file is part of the left-factorial repository. The file was originally part of
  the GNU MP Library and was modified to support NTT multiplication routines.

  Copyright 1994-1996, 2000, 2001, 2005, 2012 Free Software Foundation, Inc.
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

/*
  Division rounding the quotient towards -infinity.
  The remainder gets the same sign as the denominator.
*/

#include "internal.h"

void
ntt_mpz_fdiv_q (mpz_ptr quot, mpz_srcptr dividend, mpz_srcptr divisor)
{
  mp_size_t dividend_size = SIZ (dividend);
  mp_size_t divisor_size = SIZ (divisor);
  mpz_t rem;
  TMP_DECL;

  TMP_MARK;

  MPZ_TMP_INIT (rem, ABS (divisor_size));

  ntt_mpz_tdiv_qr (quot, rem, dividend, divisor);

  if ((divisor_size ^ dividend_size) < 0 && SIZ (rem) != 0)
    mpz_sub_ui (quot, quot, 1L);

  TMP_FREE;
}
