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
_nf_nmod_mat_first_non_zero_entry(mp_limb_t * entry, const nf_nmod_mat_t A, slong r, slong c, const nf_nmod_ctx_t ctx)
{
  slong i, n;

  for (i = c; i < A->r; ++i)
    for (n = 0; n < ctx->deg; ++n)
      if (nf_nmod_mat_entry(A, n, r, i) != 0)
	{
	  for (n = 0; n < ctx->deg; ++n)
	    entry[n] = nf_nmod_mat_entry(A, n, r, i);

	  return i;
	}

  return A->r;
}
