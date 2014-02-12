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

#include "perm.h"

#include "nfz_mat.h"


#define E(n,j,k) nfz_mat_entry(B,n,j,k)

slong
nfz_mat_fflu(nfz_mat_t B, nfz_t den, slong * perm,
	     const nfz_mat_t A, int rank_check,
	     const nfz_ctx_t ctx)
{
  slong rank, pivot_row, pivot_col, pivot_row_found;
  slong Br, Bc;
  slong i;

  nfz_t pivot_entry;
  nfz_mat_t reduction_window;

  if (nfz_mat_is_empty(A, ctx))
  {
    nfz_one(den, ctx);
    return 0;
  }

  nfz_init(pivot_entry, ctx);
  nfz_init(tmp, ctx);

  nfz_mat_set(B, A, ctx);
  Br = B->r;
  Bc = B->c;
  rank = pivot_row = 0;

  for (pivot_col = 0;
       pivot_row < Br && pivot_col < Bc;
       ++pivot_col) {
    pivot_row_found = nfz_mat_find_pivot_any(B, pivot_row, Br, pivot_col, ctx);

    if (pivot_row_found == -1)
    {
      if (rank_check)
      {
	nfz_zero(den, ctx);
	rank = 0;
	break;
      }
      continue;
    }
    else if (pivot_row_found != pivot_row)
      nfz_mat_swap_rows(B, perm, pivot_row, pivot_row_found, ctx);

    ++rank;

    nfz_mat_entry_nfz(pivot_entry, B, pivot_row, pivot_col, ctx);
    nfz_mat_init_window(reduction_window, B, pivot_row + 1, pivot_col + 1, Br, Bc, ctx);


    nfz_mat_scalar_mul_nfz(reduction_window, reduction_window, pivot_entry, ctx);
    for (i = pivot_row + 1; i < Br; ++i)
      _nfz_mat_row_submul_entry(B, i, pivot_row, pivot_col + 1, Bc, i, pivot_col, ctx);
    nfz_mat_scalar_divexact_nfz(reduction_window, reduction_window, den, ctx);

    nfz_mat_clear_window(reduction_window, ctx);

    nfz_set(den, pivot_entry, ctx);
    pivot_row++;
  }

  nfz_clear(pivot_entry, ctx);

  return rank;
}
