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
nfz_mat_CRT_nmod(nfz_mat_t out, const nfz_mat_t in1, const fmpz_t m1,
		 const nf_nmod_mat_t in2, int sign,
		 const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx)
{
  int n, r, c;
  fmpz_t fmpz_out;

  fmpz_init(fmpz_out);

  for (n = 0; n < ctx->deg; ++n)
    for (r = 0; r < out->r; ++r)
      for (c = 0; c < out->c; ++c)
	{
	  fmpz_CRT_ui(fmpz_out, nfz_mat_entry(in1, n, r, c), m1,
		      nf_nmod_mat_entry(in2, n, r, c),
		      ctx_nmod->modulus->mod.n, sign);
	  fmpz_set(nfz_mat_entry(out, n, r, c), fmpz_out);
	}

  fmpz_clear(fmpz_out);
}
