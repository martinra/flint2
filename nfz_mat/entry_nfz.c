/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Martin Raum

******************************************************************************/

#include "nfz_mat.h"

void
nfz_mat_entry_nfz(nfz_t e, const nfz_mat_t A, slong r, slong c,
		  const nfz_ctx_t ctx)
{
  slong length = 0;

  if (e->alloc < ctx->deg)
    fmpz_poly_realloc(e, ctx->deg);

  for (long n = 0; n < ctx->deg; ++n) {
    fmpz_set(e->coeffs + n, nfz_mat_entry(A, n, r, c));
    if (!fmpz_is_zero(e->coeffs + n))
      length = n + 1;
  }
  _fmpz_poly_set_length(e, length);
}
