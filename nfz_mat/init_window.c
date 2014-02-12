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

void nfz_mat_init_window(nfz_mat_t window, const nfz_mat_t mat,
			 slong r1, slong c1, slong r2, slong c2,
			 const nfz_ctx_t ctx)
{
  long i, n;

  window->r = r2 - r1;
  window->c = c2 - c1;

  window->entries = NULL;
  window->rows = (fmpz **) flint_malloc(window->r * ctx->deg * sizeof(fmpz *));
  window->poly_coeffs = (fmpz ***) flint_malloc(ctx->deg * sizeof(fmpz **));

  for (n = 0; n < ctx->deg; ++n) {
    window->poly_coeffs[n] = window->rows + n;
    for (i = 0; i < window->r; ++i)
      window->rows[n] = mat->poly_coeffs[n][i] + c1;
  }
}
