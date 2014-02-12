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
_nfz_mat_interpolate_sub_row(nfz_mat_t mat, fmpz ** evl, slong r,
			     slong start_col, slong end_col,
			     const nfz_ctx_t ctx)
{
  long j, n, tmp;

  slong len = end_col - start_col;
  fmpz * intrpl = _fmpz_vec_init(ctx->deg);

  for (j = 0; j < len; ++j) {
    _nfz_interpolate(intrpl, &tmp, (const fmpz *) evl[j], ctx);
    for (n = 0; n < ctx->deg; ++n)
      fmpz_sub(nfz_mat_entry(mat, n, r, start_col + j),
	       nfz_mat_entry(mat, n, r, start_col + j),
	       intrpl + n);
  }

  _fmpz_vec_clear(intrpl, ctx->deg);
}
