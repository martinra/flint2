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

    Copyright (C) 2013 Martin Raum
 
******************************************************************************/

#include "fmpz_vec.h"

#include "nfz.h"

void
_nfz_eval(fmpz * evl, const nfz_t f, const nfz_ctx_t ctx)
{
  for (long i = 0; i < ctx->deg; ++i) {
    fmpz_zero(evl + i);

    fmpz * evl_row = ctx->evl_mat->rows[i];
    for (long j = 0; j < ctx->deg; ++j)
      fmpz_addmul(evl + i, f->coeffs + j, evl_row + j);
  }
}
