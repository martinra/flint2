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
nfz_mat_get_nmod_mat(nf_nmod_mat_t B, const nfz_mat_t A, const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx)
{
    slong i, j;
    mp_limb_t m = ctx_nmod->mod.n;

    for (i = 0; i < A->r * ctx->deg; ++i)
      for (j = 0; j < A->c; ++j)
	B->rows[i][j] = fmpz_fdiv_ui(A->rows[i]+j, m);
}
