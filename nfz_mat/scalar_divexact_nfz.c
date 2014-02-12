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


void
nfz_mat_scalar_divexact_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c,
			    const nfz_ctx_t ctx)
{

  fmpz_t den;
  nfz_t c_inv;

  fmpz_init(den);
  nfz_init(c_inv, ctx);

  nfz_inv(c_inv, den, c, ctx);
  nfz_mat_scalar_mul_nfz(B, A, c_inv, ctx);
  nfz_mat_scalar_divexact_fmpz(B, B, den, ctx);

  fmpz_clear(den);
  nfz_clear(c_inv, ctx);
}
