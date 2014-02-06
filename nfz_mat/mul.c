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

#include "fmpz_mat.h"
#include "nfz.h"
#include "nfz_vec.h"

#include "nfz_mat.h"


void nfz_mat_mul(nfz_mat_t C, const nfz_mat_t A, const nfz_mat_t B,
		 const nfz_ctx_t ctx)
{
  long n;

  fmpz_mat_struct * AC_evl =
    (fmpz_mat_struct *) flint_malloc(ctx->deg * sizeof(fmpz_mat_struct));
  fmpz_mat_struct * B_evl =
    (fmpz_mat_struct *) flint_malloc(ctx->deg * sizeof(fmpz_mat_struct));
  for (n = 0; n < ctx->deg; ++n) {
    fmpz_mat_init(AC_evl + n, A->r, A->c);
    fmpz_mat_init(B_evl + n, A->r, A->c);
  }

  _nfz_mat_eval(AC_evl, A, ctx);
  _nfz_mat_eval(B_evl, B, ctx);

  for (n = 0; n < ctx->deg; ++n)
    fmpz_mat_mul(AC_evl + n, AC_evl + n, B_evl + n);

  _nfz_mat_interpolate(C, AC_evl, ctx);

  for (n = 0; n < ctx->deg; ++n) {
    fmpz_mat_clear(AC_evl + n);
    fmpz_mat_clear(B_evl + n);
  }
  flint_free(AC_evl);
  flint_free(B_evl);
}
