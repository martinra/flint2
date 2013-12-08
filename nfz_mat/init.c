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

#include "nfz_mat.h"

void
nfz_mat_init(nfz_mat_t mat, slong rows, slong cols, const nfz_ctx_t ctx)
{
  if ((rows) && (cols))
    {
      slong i;

      mat->entries = (fmpz *) flint_calloc(rows * cols * (ctx->deg), sizeof(fmpz));
      mat->rows = (fmpz **) flint_calloc(rows * (ctx->deg), sizeof(fmpz *));
      mat->poly_coeffs = (fmpz ***) flint_calloc((ctx->deg), sizeof(fmpz **));

      for (i = 0; i < rows * (ctx->deg); ++i)
	mat->rows[i] = mat->entries + i * cols;
      for (i = 0; i < ctx->deg; ++i)
	mat->poly_coeffs[i] = mat->rows + i * rows;
    }
  else
    mat->entries = NULL;

  mat->r = rows;
  mat->c = cols;
}
