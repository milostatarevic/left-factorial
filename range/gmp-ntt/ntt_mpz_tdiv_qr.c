/*
  This file is part of the left-factorial repository. The file was originally part of
  the GNU MP Library and was modified to support NTT multiplication routines.

  Copyright 1991, 1993, 1994, 2000, 2001, 2005, 2011, 2012 Free Software Foundation, Inc.
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
ntt_mpz_tdiv_qr (mpz_ptr quot, mpz_ptr rem, mpz_srcptr num, mpz_srcptr den)
{
  mp_size_t ql;
  mp_size_t ns, ds, nl, dl;
  mp_ptr np, dp, qp, rp;
  TMP_DECL;

  ns = SIZ (num);
  ds = SIZ (den);
  nl = ABS (ns);
  dl = ABS (ds);
  ql = nl - dl + 1;

  if (UNLIKELY (dl == 0))
    DIVIDE_BY_ZERO;

  rp = MPZ_REALLOC (rem, dl);

  if (ql <= 0)
    {
      if (num != rem)
	{
	  np = PTR (num);
	  MPN_COPY (rp, np, nl);
	  SIZ (rem) = SIZ (num);
	}
      /* This needs to follow the assignment to rem, in case the
	 numerator and quotient are the same.  */
      SIZ (quot) = 0;
      return;
    }

  qp = MPZ_REALLOC (quot, ql);

  TMP_MARK;
  np = PTR (num);
  dp = PTR (den);

  /* FIXME: We should think about how to handle the temporary allocation.
     Perhaps mpn_tdiv_qr should handle it, since it anyway often needs to
     allocate temp space.  */

  /* Copy denominator to temporary space if it overlaps with the quotient
     or remainder.  */
  if (dp == rp || dp == qp)
    {
      mp_ptr tp;
      tp = TMP_ALLOC_LIMBS (dl);
      MPN_COPY (tp, dp, dl);
      dp = tp;
    }
  /* Copy numerator to temporary space if it overlaps with the quotient or
     remainder.  */
  if (np == rp || np == qp)
    {
      mp_ptr tp;
      tp = TMP_ALLOC_LIMBS (nl);
      MPN_COPY (tp, np, nl);
      np = tp;
    }

  ntt_mpn_tdiv_qr (qp, rp, 0L, np, nl, dp, dl);

  ql -=  qp[ql - 1] == 0;
  MPN_NORMALIZE (rp, dl);

  SIZ (quot) = (ns ^ ds) >= 0 ? ql : -ql;
  SIZ (rem) = ns >= 0 ? dl : -dl;
  TMP_FREE;
}
