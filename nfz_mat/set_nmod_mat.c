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
#include "nf_nmod_mat.h"

#include "nfz_mat.h"

void nfz_mat_set_nmod_mat(nfz_mat_t A, const nf_nmod_mat_t Amod,
			  const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx)
{
  slong i, j, n;
  mp_limb_t m = nf_nmod_ctx_mod(ctx_nmod);

  for (n = 0; n < ctx->deg; ++n)
    for (i = 0; i < A->r; ++i)
      for (j = 0; j < A->c; ++j)
	fmpz_set_ui_smod(nfz_mat_entry(A, n, i, j),
			 nf_nmod_mat_entry(Amod, n, i, j), m);
}









