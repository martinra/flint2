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

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "nfz.h"
#include "nfz_vec.h"

#include "nfz_mat.h"


void
nfz_mat_scalar_divexact_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c,
			    const nfz_ctx_t ctx)
{
  long n;

  fmpz_mat_struct * mat_evl =
    (fmpz_mat_struct *) flint_malloc(ctx->deg * sizeof(fmpz_mat_struct));
  for (n = 0; n < ctx->deg; ++n)
    fmpz_mat_init(mat_evl + n, A->r, A->c);
  fmpz * c_evl = _fmpz_vec_init(ctx->deg);


  _nfz_mat_eval(mat_evl, A, ctx);
  _nfz_eval(c_evl, c->coeffs, fmpz_poly_length(c), ctx);

  for (n = 0; n < ctx->deg; ++n)
    fmpz_mat_scalar_divexact_fmpz(mat_evl + n, mat_evl + n, c_evl + n);

  _nfz_mat_interpolate(B, mat_evl, ctx);

  _fmpz_vec_clear(c_evl, ctx->deg);
  for (n = 0; n < ctx->deg; ++n)
    fmpz_mat_clear(mat_evl + n);
  flint_free(mat_evl);
}
