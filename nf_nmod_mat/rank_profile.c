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

slong
nf_nmod_mat_rank_profile(rank_profile_t rk_prof, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx)
{
  int i, n;

  slong c;
  mp_limb_t * e;

  e = (mp_limb_t *) flint_malloc(sizeof(mp_limb_t) * ctx->deg);

  c = 0;
  for (i = 0; i < A->r; ++i)
    {
      c = _nf_nmod_mat_first_non_zero_entry(e, A, i, c, ctx);
      if (c == -1)
	{
	  rank_profile_entry(rk_prof, i) = -1;
	  return i;
	}

      for (n = 1; n < ctx->deg; ++n)
	if (e[n] != 0)
	  return -1;

      rank_profile_entry(rk_prof, i) = c;
    }

  flint_free(e);

  return A->r;
}
