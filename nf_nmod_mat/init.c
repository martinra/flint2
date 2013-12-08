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

#include "nf_nmod_mat.h"

void
nf_nmod_mat_init(nf_nmod_mat_t mat, slong rows, slong cols, const nf_nmod_ctx_t ctx)
{
  if ((rows) && (cols))
    {
      slong i;

      mat->entries = (mp_limb_t *) flint_calloc(rows * cols * (ctx->deg), sizeof(mp_limb_t));
      mat->rows = (mp_limb_t **) flint_calloc(rows * (ctx->deg), sizeof(mp_limb_t *));
      mat->poly_coeffs = (mp_limb_t ***) flint_calloc((ctx->deg), sizeof(mp_limb_t **));

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
