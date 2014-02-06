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
  fmpz_mat_struct * A_evl = _nfz_mat_eval_init(A->r, A->c, ctx);
  fmpz_mat_struct * B_evl = _nfz_mat_eval_init(B->r, B->c, ctx);
  fmpz_mat_struct * C_evl;
  if (C->r == A->r && C->c == A->c)
    C_evl = A_evl;
  else if (C->r == B->r && C->c == B->c)
    C_evl = B_evl;
  else
    C_evl = _nfz_mat_eval_init(C->r, C->c, ctx);


  _nfz_mat_eval(A_evl, A, ctx);
  _nfz_mat_eval(B_evl, B, ctx);

  for (long n = 0; n < ctx->deg; ++n)
    fmpz_mat_mul(C_evl + n, A_evl + n, B_evl + n);

  _nfz_mat_interpolate(C, C_evl, ctx);

  _nfz_mat_eval_clear(A_evl, ctx);
  _nfz_mat_eval_clear(B_evl, ctx);
  if ((C->r != A->r || C->c != A->c) && (C->r != B->r || C->c != B->c))
    _nfz_mat_eval_clear(C_evl, ctx);
}
