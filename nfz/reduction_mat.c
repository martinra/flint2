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

#include "fmpz_poly.h"
#include "nfz.h"

void
_nfz_reduction_mat(fmpz_mat_t mat, const fmpz_poly_t modulus, ulong deg_bd)
{
  ulong j, n;
  ulong deg;
  fmpz_t c;

  fmpz_init(c);
  deg = fmpz_poly_degree(modulus);

  /* we assume that the modulus is monic */

  for (n = 0; n < deg && n < deg_bd; ++n)
    {
      for (j = 0; j < n; ++j)
	fmpz_set_ui(fmpz_mat_entry(mat, n, j), 0);
      fmpz_set_ui(fmpz_mat_entry(mat, n, n), 1);
      for (j = n + 1; j < deg; ++j)
	fmpz_set_ui(fmpz_mat_entry(mat, n, j), 0);
    }

  for (n = deg; n < deg_bd; ++n)
    {
      for (j = 0; j < deg; ++j)
	fmpz_set_ui(fmpz_mat_entry(mat, n, j), 0);

      for (j = 0; j < deg; ++j)
	{
	  fmpz_poly_get_coeff_fmpz(c, modulus, j);
	  _fmpz_vec_scalar_submul_fmpz(mat->rows[n], mat->rows[n - deg + j], deg, c);
	}
    }

  fmpz_clear(c);
}
