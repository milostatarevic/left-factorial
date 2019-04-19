/*
  This file is part of the left-factorial repository. The file was originally part of
  the GNU MP Library and was modified to support NTT multiplication routines.

  Copyright 1991, 1993-1996, 2001, 2002, 2005, 2010, 2012 Free Software Foundation, Inc.
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

#include "internal.h"

void
ntt_mpz_mod (mpz_ptr rem, mpz_srcptr dividend, mpz_srcptr divisor)
{
  mp_size_t rn, bn;
  mpz_t temp_divisor;
  TMP_DECL;

  TMP_MARK;

  bn = ABSIZ(divisor);

  /* We need the original value of the divisor after the remainder has been
     preliminary calculated.  We have to copy it to temporary space if it's
     the same variable as REM.  */
  if (rem == divisor)
    {
      PTR(temp_divisor) = TMP_ALLOC_LIMBS (bn);
      MPN_COPY (PTR(temp_divisor), PTR(divisor), bn);
    }
  else
    {
      PTR(temp_divisor) = PTR(divisor);
    }
  SIZ(temp_divisor) = bn;
  divisor = temp_divisor;

  ntt_mpz_tdiv_r (rem, dividend, divisor);

  rn = SIZ (rem);
  if (rn < 0)
    mpz_add (rem, rem, divisor);

  TMP_FREE;
}
