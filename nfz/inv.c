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

#include "flint.h"

#include "nfz.h"

void
nfz_inv(nfz_t f, fmpz_t den, const nfz_t g, const nfz_ctx_t ctx)
{
  fmpz_t cont;
  fmpz_poly_t tmp;

  fmpz_poly_init(tmp);
  
  fmpz_poly_content(cont, g);
  if (!fmpz_is_one(cont)) {
    nfz_divexact_fmpz(f, g, cont, ctx);
    fmpz_poly_xgcd(den, f, tmp, f, ctx->modulus);
    fmpz_mul(den, den, cont);
  } else {
    fmpz_poly_xgcd(den, f, tmp, g, ctx->modulus);
  }
  nfz_reduce(f, ctx);

  fmpz_poly_clear(tmp);
}
