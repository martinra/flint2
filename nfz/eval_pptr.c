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

#include "fmpz_vec.h"

#include "nfz.h"

void
_nfz_eval_pptr(fmpz ** evl, const fmpz ** f, long length, const nfz_ctx_t ctx)
{
  if (length >= ctx->deg) length = ctx->deg;

  for (long i = 0; i < ctx->deg; ++i) {
    fmpz_zero(evl[i]);

    fmpz * evl_row = ctx->evl_mat->rows[i];
    for (long j = 0; j < length; ++j)
      fmpz_addmul(evl + i, f[j], evl_row + j);
  }
}
