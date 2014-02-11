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
_nfz_mat_eval(fmpz_mat_struct * evl, const nfz_mat_t A, const nfz_ctx_t ctx)
{
  long i, j, n;

  fmpz ** mat_pptr = (fmpz **) flint_malloc(ctx->deg * sizeof(fmpz *));
  fmpz ** evl_pptr = (fmpz **) flint_malloc(ctx->evl_mat->r * sizeof(fmpz *));

  for (i = 0; i < A->r; ++i)
    for (j = 0; j < A->c; ++j) {
      for (n = 0; n < ctx->deg; ++n)
	mat_pptr[n] = A->poly_coeffs[n][i] + j;
      for (n = 0; n < ctx->evl_mat->r; ++n)
	evl_pptr[n] = evl[n].rows[i] + j;

      _nfz_eval_pptr(evl_pptr, (const fmpz **) mat_pptr, ctx->deg, ctx);
    }

  flint_free(mat_pptr);
  flint_free(evl_pptr);
}
