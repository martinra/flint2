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

#include "nfz.h"

void
_nfz_reduction_mat(fmpz_mat_t mat, fmpz_poly_t modulus, ulong deg_bd)
{
  /* we assume that the modulus is monic */

  for (i = 0; i < deg && i < deg_bd; ++i)
    {
      for (j = 0; j < i; ++j)
	fmpz_set_ui(fmpz_mat_entry(mat, i, j), 0);
      fmpz_set_ui(fmpz_mat_entry(mat, i, i), 1);
      for (j = i + 1; j < deg; ++j)
	fmpz_set_ui(fmpz_mat_entry(mat, i, j), 0);
    }

  for (i = deg; i < deg_bd; ++i)
    {
      for (j = 0; j < deg; ++j)
	fmpz_set_ui(fmpz_mat_entry(mat, i, j), 0);

      for (j = 0; j < deg; ++j)
	_fmpz_vec_scalar_submul_fmpz(mat->r[i], mat->r[i - deg + j], deg, fmpz_poly_get_coeff(modulus, j));
    }
}
