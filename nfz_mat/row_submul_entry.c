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

#include "nfz_vec.h"


#include "nfz_mat.h"


void _nfz_mat_row_submul_entry(nfz_mat_t mat,
			       slong r1, slong r2, slong start_col, slong end_col,
			       slong er, slong ec, const nfz_ctx_t ctx)
{
  slong i, n;
  slong len = end_col - start_col;
  fmpz ** row_evl = _nfz_vec_eval_init(len, ctx);
  fmpz * e_evl = _nfz_eval_init(ctx);

  _nfz_mat_eval_row(row_evl, mat, r2, start_col, end_col, ctx);
  _nfz_mat_eval_entry(e_evl, er, ec, ctx);

  for (i = 0; i < len; ++i)
    for (n = 0; n < ctx->evl_mat->r; ++n)
      fmpz_mul(row_evl[i] + n, row_evl[i] + n, e_evl + n);

  _nfz_mat_interpolate_sub_row(mat, row_evl, r1, start_col, end_col, ctx);

  _nfz_vec_eval_clear(row_evl, len, ctx);
  _nfz_eval_clear(e_evl, ctx);
}
