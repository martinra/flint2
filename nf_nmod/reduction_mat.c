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

#include "nf_nmod.h"

void
_nf_nmod_reduction_mat(nmod_mat_t mat, const nmod_poly_t modulus, ulong deg_bd)
{
  ulong j, n;
  ulong deg;

  deg = nmod_poly_degree(modulus);

  /* we assume that the modulus is monic */

  for (n = 0; n < deg && n < deg_bd; ++n)
    {
      for (j = 0; j < n; ++j)
	*(mat->rows[n] + j) = 0;
      *(mat->rows[n] + n) = 1;
      for (j = n + 1; j < deg; ++j)
	*(mat->rows[n] + j) = 0;
    }

  for (n = deg; n < deg_bd; ++n)
    {
      for (j = 0; j < deg; ++j)
	*(mat->rows[n] + j) = 0;

      for (j = 0; j < deg; ++j)
	_nmod_vec_scalar_addmul_nmod(mat->rows[n], mat->rows[n - deg + j], deg, nmod_poly_get_coeff_ui(modulus, j), modulus->mod);
    }
}
