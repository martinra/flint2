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

void _nfz_mat_eval_row(fmpz ** evl, nfz_mat_t mat, slong r,
		       slong start_col, slong end_col, const nfz_ctx_t ctx)
{
  long j, n;
  slong len = end_col - start_col;

  fmpz ** mat_pptr = (fmpz **) flint_malloc(ctx->deg * sizeof(fmpz *));

  for (j = 0; j < len; ++j) {
    for (n = 0; n < ctx->deg; ++n)
    mat_pptr[n] = mat->poly_coeffs[n][r] + start_col + j;

    _nfz_eval_pptr(evl + j, (const fmpz **) mat_pptr, ctx->deg, ctx);
  }

  flint_free(mat_pptr);
}
