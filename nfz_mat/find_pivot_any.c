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

slong
nfz_mat_find_pivot_any(const nfz_mat_t mat,
		       slong start_row, slong end_row, slong c, const nfz_ctx_t ctx)
{
  slong i, n;

  for (i = start_row; i < end_row; ++i)
  {
    for (n = 0; n < ctx->deg; ++n)
      if (!fmpz_is_zero(nfz_mat_entry(mat, n, i, c)))
	return i;
  }

  return -1;
}
