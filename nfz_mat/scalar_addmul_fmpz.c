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

#include "fmpz_vec.h"

#include "nfz_mat.h"


void
nfz_mat_scalar_addmul_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c,
			   const nfz_ctx_t ctx)
{
  long i, n;

  for (n = 0; n < ctx->deg; ++n)
    for (i = 0; i < A->r; ++i)
      _fmpz_vec_scalar_addmul_fmpz(B->poly_coeffs[n][i], A->poly_coeffs[n][i], A->c, c);
}
